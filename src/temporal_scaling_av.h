// ============================================================
// temporal_scaling_av.h — lecture d'un fichier via libav (ffmpeg intégré)
//
// Compilé seulement avec -DTS_FFMPEG (voir compilation.bat). Remplace le
// "ffmpeg ... | temporal_scaling" : le fichier est ouvert, décodé, désentrelacé
// (filtre libavfilter au choix, ex. bwdif=1) et converti en RGB24 directement
// dans le programme, sans processus externe ni pipe. Le son est décodé dans le
// même processus (s16 48 kHz stéréo) et alimente la lecture audio.
//
//   thread demux  : lit les paquets -> files vidéo / audio
//   thread audio  : décode + rééchantillonne -> tampon circulaire PCM
//   thread lecteur (celui de gpu_stream_run) : décode la vidéo -> filtres -> RGB
//
// Les images portent leur horodatage (pts) : l'affichage se cale dessus, avec
// l'audio comme horloge maître ; les images en retard sont abandonnées dès
// la sortie du décodeur (pas de conversion ni de mise à l'échelle inutile).
// ============================================================
#ifndef TEMPORAL_SCALING_AV_H
#define TEMPORAL_SCALING_AV_H

#ifdef TS_FFMPEG

#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
#include <libavfilter/avfilter.h>
#include <libavfilter/buffersrc.h>
#include <libavfilter/buffersink.h>
#include <libavutil/opt.h>
#include <libavutil/channel_layout.h>
#include <libavutil/samplefmt.h>
#include <libavutil/mathematics.h>
#include <libswresample/swresample.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AV_AUD_RATE 48000
#define AV_AUD_CH   2
#define AV_AUD_RING (AV_AUD_RATE * AV_AUD_CH * 2 * 4)   // 4 s de PCM

// ---- file de paquets ----
typedef struct PktNode { AVPacket *pkt; struct PktNode *next; } PktNode;
typedef struct {
	PktNode *head, *tail;
	int n, cap, eof;
	pthread_mutex_t m;
	pthread_cond_t  c;
} PQueue;

static void pq_init(PQueue *q, int cap)
{
	memset(q, 0, sizeof(*q));
	q->cap = cap;
	pthread_mutex_init(&q->m, NULL);
	pthread_cond_init(&q->c, NULL);
}

static void pq_put(PQueue *q, AVPacket *pkt)   // pkt == NULL : fin de flux
{
	pthread_mutex_lock(&q->m);
	while(pkt && q->n >= q->cap) pthread_cond_wait(&q->c, &q->m);
	if(!pkt) q->eof = 1;
	else {
		PktNode *nd = (PktNode *)malloc(sizeof(*nd));
		nd->pkt = pkt; nd->next = NULL;
		if(q->tail) q->tail->next = nd; else q->head = nd;
		q->tail = nd; q->n++;
	}
	pthread_cond_broadcast(&q->c);
	pthread_mutex_unlock(&q->m);
}

// renvoie un paquet, ou NULL en fin de flux
static AVPacket *pq_get(PQueue *q)
{
	pthread_mutex_lock(&q->m);
	while(!q->head && !q->eof) pthread_cond_wait(&q->c, &q->m);
	AVPacket *p = NULL;
	if(q->head) {
		PktNode *nd = q->head;
		q->head = nd->next; if(!q->head) q->tail = NULL;
		q->n--; p = nd->pkt; free(nd);
		pthread_cond_broadcast(&q->c);
	}
	pthread_mutex_unlock(&q->m);
	return p;
}

// ---- tampon circulaire audio ----
typedef struct {
	uint8_t *buf;
	size_t rd, wr, count;
	int eof;
	pthread_mutex_t m;
	pthread_cond_t  c;
} ARing;

typedef struct {
	AVFormatContext *fmt;
	int vi, ai;
	AVCodecContext *vc, *ac;
	AVFilterGraph *vg;
	AVFilterContext *vsrc, *vsink;
	AVFrame *vframe, *fframe;
	SwrContext *swr;
	AVRational vsink_tb;
	int w, h;
	double dar;                     // rapport d'aspect d'affichage (largeur/hauteur, SAR inclus)
	int interlaced;                 // 1 si le contenu est détecté entrelacé
	int color_709, color_full;      // matrice YUV->RGB : BT.709 sinon BT.601 ; plage pleine sinon limitée
	char vf_used[512];              // chaîne de filtres effectivement appliquée
	PQueue vq, aq;
	ARing ring;
	pthread_t th_demux, th_audio;
	// état du lecteur vidéo
	int  v_draining, v_eof;
	double last_pts, frame_dur;
	long dropped;
	// audio
	volatile double audio_pts0;     // pts (s) du premier échantillon décodé, -1 = inconnu
	int has_audio;
} AvSrc;

// horloge média fournie par l'affichage (s), ou <0 si inconnue
static double (*ts_av_clock_fn)(void) = NULL;
// retard estimé (s) entre la sortie du décodeur et l'affichage (images en file x temps GPU par image)
static volatile double ts_av_lead_s = 0.0;
#ifdef _WIN32
#define TS_AV_LATE() (g_disp_late)
#else
#define TS_AV_LATE() (0.0)
#endif

static void ar_put(ARing *r, const uint8_t *src, size_t n)
{
	pthread_mutex_lock(&r->m);
	while(n > 0) {
		while(r->count >= AV_AUD_RING) pthread_cond_wait(&r->c, &r->m);
		size_t room = AV_AUD_RING - r->count;
		size_t k = n < room ? n : room;
		size_t first = AV_AUD_RING - r->wr; if(first > k) first = k;
		memcpy(r->buf + r->wr, src, first);
		memcpy(r->buf, src + first, k - first);
		r->wr = (r->wr + k) % AV_AUD_RING; r->count += k;
		src += k; n -= k;
		pthread_cond_broadcast(&r->c);
	}
	pthread_mutex_unlock(&r->m);
}

// Lecture bloquante pour la sortie audio : n octets, ou moins en fin de flux.
static size_t ts_av_audio_read(void *ctx, void *dst, size_t n)
{
	AvSrc *s = (AvSrc *)ctx;
	ARing *r = &s->ring;
	size_t got = 0;
	pthread_mutex_lock(&r->m);
	while(got < n) {
		while(r->count == 0 && !r->eof) pthread_cond_wait(&r->c, &r->m);
		if(r->count == 0) break;
		size_t k = n - got < r->count ? n - got : r->count;
		size_t first = AV_AUD_RING - r->rd; if(first > k) first = k;
		memcpy((uint8_t *)dst + got, r->buf + r->rd, first);
		memcpy((uint8_t *)dst + got + first, r->buf, k - first);
		r->rd = (r->rd + k) % AV_AUD_RING; r->count -= k; got += k;
		pthread_cond_broadcast(&r->c);
	}
	pthread_mutex_unlock(&r->m);
	return got;
}

static double ts_av_audio_pts0(void *ctx) { return ((AvSrc *)ctx)->audio_pts0; }

// ---- threads ----
static void *av_demux_thread(void *p)
{
	AvSrc *s = (AvSrc *)p;
	for(;;) {
		AVPacket *pkt = av_packet_alloc();
		if(!pkt || av_read_frame(s->fmt, pkt) < 0) { av_packet_free(&pkt); break; }
		if(pkt->stream_index == s->vi)                      pq_put(&s->vq, pkt);
		else if(pkt->stream_index == s->ai && s->has_audio) pq_put(&s->aq, pkt);
		else av_packet_free(&pkt);
	}
	pq_put(&s->vq, NULL);
	if(s->has_audio) pq_put(&s->aq, NULL);
	return NULL;
}

static void av_audio_emit(AvSrc *s, AVFrame *f)
{
	if(s->audio_pts0 < 0.0) {
		AVStream *st = s->fmt->streams[s->ai];
		int64_t pts = f->best_effort_timestamp != AV_NOPTS_VALUE ? f->best_effort_timestamp : f->pts;
		s->audio_pts0 = (pts == AV_NOPTS_VALUE) ? 0.0 : pts * av_q2d(st->time_base);
	}
	int cap = (int)av_rescale_rnd(swr_get_delay(s->swr, s->ac->sample_rate) + (f ? f->nb_samples : 0),
	                              AV_AUD_RATE, s->ac->sample_rate, AV_ROUND_UP) + 32;
	uint8_t *out = (uint8_t *)malloc((size_t)cap * AV_AUD_CH * 2);
	if(!out) return;
	int n = swr_convert(s->swr, &out, cap, f ? (const uint8_t **)f->extended_data : NULL, f ? f->nb_samples : 0);
	if(n > 0) ar_put(&s->ring, out, (size_t)n * AV_AUD_CH * 2);
	free(out);
}

static void *av_audio_thread(void *p)
{
	AvSrc *s = (AvSrc *)p;
	AVFrame *f = av_frame_alloc();
	for(;;) {
		AVPacket *pkt = pq_get(&s->aq);
		if(avcodec_send_packet(s->ac, pkt) < 0 && pkt) { av_packet_free(&pkt); continue; }
		av_packet_free(&pkt);
		while(avcodec_receive_frame(s->ac, f) >= 0) { av_audio_emit(s, f); av_frame_unref(f); }
		if(!pkt && s->aq.eof && !s->aq.head) break;
	}
	av_audio_emit(s, NULL);   // vide le rééchantillonneur
	av_frame_free(&f);
	pthread_mutex_lock(&s->ring.m);
	s->ring.eof = 1;
	pthread_cond_broadcast(&s->ring.c);
	pthread_mutex_unlock(&s->ring.m);
	return NULL;
}

static int av_open_codec(AVStream *st, AVCodecContext **out)
{
	const AVCodec *dec = avcodec_find_decoder(st->codecpar->codec_id);
	if(!dec) return -1;
	AVCodecContext *cc = avcodec_alloc_context3(dec);
	if(!cc) return -1;
	if(avcodec_parameters_to_context(cc, st->codecpar) < 0) { avcodec_free_context(&cc); return -1; }
	cc->pkt_timebase = st->time_base;
	cc->thread_count = 0;   // auto
	if(avcodec_open2(cc, dec, NULL) < 0) { avcodec_free_context(&cc); return -1; }
	*out = cc;
	return 0;
}

// Détecte si la vidéo est entrelacée : field_order du flux si connu, sinon on
// décode les premières images dans un second contexte (le principal reste intact)
// et on regarde leur drapeau "entrelacé".
static int av_probe_interlaced(const char *file)
{
	AVFormatContext *f = NULL;
	int res = 0;
	if(avformat_open_input(&f, file, NULL, NULL) < 0) return 0;
	if(avformat_find_stream_info(f, NULL) >= 0)
	{
		int vi = av_find_best_stream(f, AVMEDIA_TYPE_VIDEO, -1, -1, NULL, 0);
		if(vi >= 0)
		{
			AVStream *st = f->streams[vi];
			switch(st->codecpar->field_order) {
				case AV_FIELD_PROGRESSIVE: break;
				case AV_FIELD_TT: case AV_FIELD_BB: case AV_FIELD_TB: case AV_FIELD_BT: res = 1; break;
				default: {
					AVCodecContext *cc = NULL;
					if(av_open_codec(st, &cc) == 0) {
						AVPacket *pkt = av_packet_alloc();
						AVFrame *fr = av_frame_alloc();
						int seen = 0, inter = 0;
						for(int n = 0; n < 400 && seen < 10 && av_read_frame(f, pkt) >= 0; n++) {
							if(pkt->stream_index == vi && avcodec_send_packet(cc, pkt) >= 0) {
								while(avcodec_receive_frame(cc, fr) >= 0) {
									seen++;
									if(fr->flags & AV_FRAME_FLAG_INTERLACED) inter++;
									av_frame_unref(fr);
								}
							}
							av_packet_unref(pkt);
						}
						res = (seen > 0 && inter * 2 >= seen);
						av_packet_free(&pkt); av_frame_free(&fr);
						avcodec_free_context(&cc);
					}
				}
			}
		}
	}
	avformat_close_input(&f);
	return res;
}

// Ouvre le fichier. deint : chaîne de désentrelacement (NULL = automatique) ; vf : filtres en plus.
// La sortie du graphe est toujours en yuv420p.
static AvSrc *ts_av_open(const char *file, const char *vf, const char *deint, int want_audio, int auto_deint)
{
	av_log_set_level(AV_LOG_ERROR);
	AvSrc *s = (AvSrc *)calloc(1, sizeof(*s));
	if(!s) return NULL;
	s->last_pts = -1.0;
	s->audio_pts0 = -1.0;

	if(avformat_open_input(&s->fmt, file, NULL, NULL) < 0) { fprintf(stderr, "ffmpeg: cannot open %s\n", file); goto fail; }
	if(avformat_find_stream_info(s->fmt, NULL) < 0)        { fprintf(stderr, "ffmpeg: no stream info\n"); goto fail; }
	s->vi = av_find_best_stream(s->fmt, AVMEDIA_TYPE_VIDEO, -1, -1, NULL, 0);
	if(s->vi < 0) { fprintf(stderr, "ffmpeg: no video stream\n"); goto fail; }
	s->ai = want_audio ? av_find_best_stream(s->fmt, AVMEDIA_TYPE_AUDIO, -1, s->vi, NULL, 0) : -1;

	AVStream *vst = s->fmt->streams[s->vi];
	if(av_open_codec(vst, &s->vc) != 0) { fprintf(stderr, "ffmpeg: cannot open the video decoder\n"); goto fail; }

	// ---- graphe de filtres : [vf,]crop pair,format=yuv420p (4:2:0 conservé, le scaler fait le reste) ----
	s->vg = avfilter_graph_alloc();
	char args[512];
	AVRational sar = av_guess_sample_aspect_ratio(s->fmt, vst, NULL);
	if(sar.num <= 0 || sar.den <= 0) sar = (AVRational){1, 1};
	AVRational fr  = av_guess_frame_rate(s->fmt, vst, NULL);
	int n = snprintf(args, sizeof(args),
		"video_size=%dx%d:pix_fmt=%d:time_base=%d/%d:pixel_aspect=%d/%d",
		s->vc->width, s->vc->height, (int)s->vc->pix_fmt,
		vst->time_base.num, vst->time_base.den, sar.num, sar.den);
	if(fr.num > 0 && fr.den > 0) snprintf(args + n, sizeof(args) - n, ":frame_rate=%d/%d", fr.num, fr.den);
	if(avfilter_graph_create_filter(&s->vsrc,  avfilter_get_by_name("buffer"),     "in",  args, NULL, s->vg) < 0 ||
	   avfilter_graph_create_filter(&s->vsink, avfilter_get_by_name("buffersink"), "out", NULL, NULL, s->vg) < 0) {
		fprintf(stderr, "ffmpeg: cannot create the filter graph\n"); goto fail;
	}
	{
		char desc[1024];
		// Désentrelacement : "deint" imposé par l'utilisateur (chaîne ffmpeg, "" = aucun), sinon
		// automatique : bwdif=1 si le contenu est entrelacé (et si aucune chaîne --vf n'est donnée).
		// "vf" = filtres supplémentaires, appliqués après le désentrelacement.
		s->interlaced = av_probe_interlaced(file);
		const char *dchain = deint;
		if(!dchain) dchain = (!(vf && vf[0]) && auto_deint && s->interlaced) ? "bwdif=1" : "";
		const char *extra = (vf && vf[0]) ? vf : "";
		char used[512];
		snprintf(used, sizeof(used), "%s%s%s", dchain, (dchain[0] && extra[0]) ? "," : "", extra);
		// dimensions paires (exigé par le 4:2:0) ; yuvj420p conservé tel quel (plage pleine)
		snprintf(desc, sizeof(desc), "%s%scrop=trunc(iw/2)*2:trunc(ih/2)*2,format=pix_fmts=yuv420p|yuvj420p",
			used, used[0] ? "," : "");
		snprintf(s->vf_used, sizeof(s->vf_used), "%s", used);
		AVFilterInOut *outs = avfilter_inout_alloc(), *ins = avfilter_inout_alloc();
		outs->name = av_strdup("in");  outs->filter_ctx = s->vsrc;  outs->pad_idx = 0; outs->next = NULL;
		ins->name  = av_strdup("out"); ins->filter_ctx  = s->vsink; ins->pad_idx  = 0; ins->next  = NULL;
		int r = avfilter_graph_parse_ptr(s->vg, desc, &ins, &outs, NULL);
		avfilter_inout_free(&ins); avfilter_inout_free(&outs);
		if(r < 0 || avfilter_graph_config(s->vg, NULL) < 0) {
			fprintf(stderr, "ffmpeg: invalid filter chain \"%s\"\n", desc); goto fail;
		}
	}
	s->w = av_buffersink_get_w(s->vsink);
	s->h = av_buffersink_get_h(s->vsink);
	s->dar = (double)s->w / (double)s->h * (double)sar.num / (double)sar.den;
	{
		const AVCodecParameters *par = vst->codecpar;
		const enum AVPixelFormat pf = s->vc->pix_fmt;
		s->color_full = (par->color_range == AVCOL_RANGE_JPEG) ||
		                pf == AV_PIX_FMT_YUVJ420P || pf == AV_PIX_FMT_YUVJ422P || pf == AV_PIX_FMT_YUVJ444P;
		switch(par->color_space) {
			case AVCOL_SPC_BT709:                          s->color_709 = 1; break;
			case AVCOL_SPC_BT470BG: case AVCOL_SPC_SMPTE170M: s->color_709 = 0; break;
			default:                                       s->color_709 = (s->vc->height >= 720); break;   // HD -> 709
		}
	}
	s->vsink_tb = av_buffersink_get_time_base(s->vsink);
	AVRational ofr = av_buffersink_get_frame_rate(s->vsink);
	s->frame_dur = (ofr.num > 0 && ofr.den > 0) ? (double)ofr.den / ofr.num : 0.0;
	s->vframe = av_frame_alloc();
	s->fframe = av_frame_alloc();

	pq_init(&s->vq, 96);
	pq_init(&s->aq, 100000);   // illimitée en pratique : l'audio ne doit jamais bloquer le démultiplexage

	// ---- audio ----
	if(s->ai >= 0) {
		AVStream *ast = s->fmt->streams[s->ai];
		if(av_open_codec(ast, &s->ac) == 0) {
			AVChannelLayout in_l, out_l;
			if(av_channel_layout_copy(&in_l, &s->ac->ch_layout) < 0) in_l = (AVChannelLayout){0};
			if(in_l.order == AV_CHANNEL_ORDER_UNSPEC && in_l.nb_channels > 0)
				av_channel_layout_default(&in_l, in_l.nb_channels);
			av_channel_layout_default(&out_l, AV_AUD_CH);
			if(swr_alloc_set_opts2(&s->swr, &out_l, AV_SAMPLE_FMT_S16, AV_AUD_RATE,
			                       &in_l, s->ac->sample_fmt, s->ac->sample_rate, 0, NULL) == 0 &&
			   swr_init(s->swr) == 0) {
				s->has_audio = 1;
				s->ring.buf = (uint8_t *)malloc(AV_AUD_RING);
				pthread_mutex_init(&s->ring.m, NULL);
				pthread_cond_init(&s->ring.c, NULL);
			}
			av_channel_layout_uninit(&in_l);
			av_channel_layout_uninit(&out_l);
		}
		if(!s->has_audio) fprintf(stderr, "ffmpeg: audio disabled (no usable audio stream)\n");
	}

	pthread_create(&s->th_demux, NULL, av_demux_thread, s);
	if(s->has_audio) pthread_create(&s->th_audio, NULL, av_audio_thread, s);
	return s;

fail:
	free(s);   // fuite assumée : on quitte le programme juste après
	return NULL;
}

// Prochaine image vidéo filtrée en RGB24 dans dst (w*h*3 octets).
// pts : horodatage en secondes. Renvoie 1, ou 0 en fin de flux.
static int ts_av_read(AvSrc *s, unsigned char *dst, double *pts_out)
{
	for(;;)
	{
		int ret = av_buffersink_get_frame(s->vsink, s->fframe);
		if(ret >= 0)
		{
			double pts = (s->fframe->pts == AV_NOPTS_VALUE) ? -1.0 : s->fframe->pts * av_q2d(s->vsink_tb);
			if(pts >= 0.0) {
				if(s->last_pts >= 0.0 && pts > s->last_pts && pts - s->last_pts < 1.0)
					s->frame_dur = pts - s->last_pts;
				s->last_pts = pts;
				// image déjà en retard sur l'audio : inutile de la mettre à l'échelle
				double clk = ts_av_clock_fn ? ts_av_clock_fn() : -1.0;
				if(clk >= 0.0 && pts < clk + ts_av_lead_s + TS_AV_LATE() - 0.06) {
					s->dropped++;
					av_frame_unref(s->fframe);
					continue;
				}
			}
			// yuv420p planaire sans marge : plan Y, puis U, puis V
			{
				unsigned char *p = dst;
				for(int pl = 0; pl < 3; pl++) {
					const int pw = pl ? s->w / 2 : s->w, ph = pl ? s->h / 2 : s->h;
					for(int y = 0; y < ph; y++)
						memcpy(p + (size_t)y * pw, s->fframe->data[pl] + (size_t)y * s->fframe->linesize[pl], (size_t)pw);
					p += (size_t)pw * ph;
				}
			}
			if(pts_out) *pts_out = pts;
			av_frame_unref(s->fframe);
			return 1;
		}
		if(ret == AVERROR_EOF) return 0;
		if(ret != AVERROR(EAGAIN)) return 0;

		// le filtre veut plus d'entrée : décodeur -> filtre
		ret = avcodec_receive_frame(s->vc, s->vframe);
		if(ret >= 0) {
			s->vframe->pts = s->vframe->best_effort_timestamp;
			av_buffersrc_add_frame_flags(s->vsrc, s->vframe, AV_BUFFERSRC_FLAG_KEEP_REF);
			av_frame_unref(s->vframe);
			continue;
		}
		if(ret == AVERROR(EAGAIN)) {
			AVPacket *pkt = pq_get(&s->vq);
			if(pkt) { avcodec_send_packet(s->vc, pkt); av_packet_free(&pkt); }
			else if(!s->v_draining) { avcodec_send_packet(s->vc, NULL); s->v_draining = 1; }
			else return 0;
			continue;
		}
		if(ret == AVERROR_EOF) {
			if(!s->v_eof) { av_buffersrc_add_frame(s->vsrc, NULL); s->v_eof = 1; continue; }
			return 0;
		}
		return 0;
	}
}

#endif // TS_FFMPEG
#endif // TEMPORAL_SCALING_AV_H
