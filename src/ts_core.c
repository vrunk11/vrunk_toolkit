// ============================================================
// ts_core.c — le moteur de temporal_scaling exposé en bibliothèque (voir ts_core.h).
//
// On inclut temporal_scaling.c tel quel (son main() est renommé) : les cartes de
// mise à l'échelle, les noyaux et le backend OpenCL sont donc *exactement* ceux du
// programme en ligne de commande. TS_PLUGIN supprime les tables CPU de 64 Mo.
// ============================================================
#define TS_PLUGIN 1
#define main ts_unused_main
#include "temporal_scaling.c"
#undef main

#include "ts_core.h"

struct TsCtx {
	TsParams p;
	BresenhamMap map, map_uv;
	int          need_map_uv;
	int         *tmp_weights;
	int          uv_in_w, uv_in_h, uv_out_w, uv_out_h;
	size_t       in_len, out_frame, out_len;
	unsigned char *buf[2];       // images d'entrée contiguës Y|U|V (courante / précédente)
	int          cur;            // indice de buf[] pour l'image courante
	int          have_prev;
	int          blend;          // 1 si un fondu temporel est actif (retard d'une image)
};

static int ts_alloc_map(BresenhamMap *m, int wlen, int hlen)
{
	memset(m, 0, sizeof(*m));
	m->w_length = wlen;
	m->h_length = hlen;
	int **wa[5] = { &m->w_src, &m->w_flag, &m->w_weight, &m->w_heavy, &m->w_pos };
	int **ha[5] = { &m->h_src, &m->h_flag, &m->h_weight, &m->h_heavy, &m->h_pos };
	for(int i = 0; i < 5; i++) {
		*wa[i] = (int *)ALIGNED_MALLOC((size_t)wlen * sizeof(int), 64);
		*ha[i] = (int *)ALIGNED_MALLOC((size_t)hlen * sizeof(int), 64);
		if(!*wa[i] || !*ha[i]) return -1;
		memset(*wa[i], 0, (size_t)wlen * sizeof(int));
		memset(*ha[i], 0, (size_t)hlen * sizeof(int));
	}
	return 0;
}

static void ts_free_map(BresenhamMap *m)
{
	ALIGNED_FREE(m->w_src); ALIGNED_FREE(m->w_flag); ALIGNED_FREE(m->w_weight);
	ALIGNED_FREE(m->w_heavy); ALIGNED_FREE(m->w_pos);
	ALIGNED_FREE(m->h_src); ALIGNED_FREE(m->h_flag); ALIGNED_FREE(m->h_weight);
	ALIGNED_FREE(m->h_heavy); ALIGNED_FREE(m->h_pos);
	memset(m, 0, sizeof(*m));
}

// cartes + poids selon le mode (même logique que main() de temporal_scaling.c)
static void ts_weights_for_mode(BresenhamMap *m, int in_w, int in_h, int out_w, int out_h,
                                int ratio, int mode, double curve_base)
{
	compute_scale_map(m, in_w, in_h, out_w, out_h, ratio);
	if(mode == 3) {
		compute_curve_weights(m, out_w, out_h, ratio, CURVE_LINEAR, curve_base);
		adjust_phase_direction_weighted(m, in_w, in_h, out_w, out_h, ratio);
	} else if(mode >= 4) {
		compute_curve_weights(m, out_w, out_h, ratio, CURVE_EXPONENTIAL, curve_base);
		adjust_phase_direction_weighted(m, in_w, in_h, out_w, out_h, ratio);
	} else if(mode == 2) {
		adjust_phase_direction(m, in_w, in_h, out_w, out_h, ratio);
	}
}

TsCtx *ts_open(const TsParams *p, char *err, size_t errlen)
{
#define FAIL(...) do { if(err && errlen) snprintf(err, errlen, __VA_ARGS__); ts_close(c); return NULL; } while(0)
	TsCtx *c = (TsCtx *)calloc(1, sizeof(*c));
	if(!c) { if(err && errlen) snprintf(err, errlen, "out of memory"); return NULL; }
	c->p = *p;
	TsParams *q = &c->p;

	if(q->ratio < 1) q->ratio = 1;
	if(q->mode < 0 || q->mode > 7) FAIL("mode %d is not supported (0..7)", q->mode);
	if(q->tmp_mode < 0 || q->tmp_mode > 3) FAIL("tmp-mode %d is not supported (0..3)", q->tmp_mode);
	if(q->curve_base <= 1.0) q->curve_base = 2.0;
	if(q->in_w <= 0 || q->in_h <= 0) FAIL("invalid input size");
	if(q->chroma != TS_CHROMA_444 && (q->in_w & 1)) FAIL("input width %d must be even", q->in_w);
	if(q->chroma == TS_CHROMA_420 && (q->in_h & 1)) FAIL("input height %d must be even", q->in_h);
	if(q->out_w <= q->in_w || q->out_h <= q->in_h)
		FAIL("output %dx%d must be larger than the input %dx%d (the scaler only enlarges)", q->out_w, q->out_h, q->in_w, q->in_h);

	const int in_w = q->in_w, in_h = q->in_h, out_w = q->out_w, out_h = q->out_h, ratio = q->ratio;
	c->uv_in_w  = (q->chroma == TS_CHROMA_444) ? in_w : in_w / 2;
	c->uv_in_h  = (q->chroma == TS_CHROMA_420) ? in_h / 2 : in_h;
	c->uv_out_w = out_w;      // sortie 4:4:4 puis conversion RGB sur le GPU
	c->uv_out_h = out_h;
	const PixFmt in_fmt = (q->chroma == TS_CHROMA_420) ? PIX_420P : (q->chroma == TS_CHROMA_422) ? PIX_422P : PIX_444P;
	c->need_map_uv = !(q->chroma == TS_CHROMA_444);

	c->in_len    = (size_t)in_w * in_h + 2 * (size_t)c->uv_in_w * c->uv_in_h;
	c->out_frame = (size_t)out_w * out_h * 3;
	c->out_len   = c->out_frame * (size_t)ratio;
	c->buf[0] = (unsigned char *)ALIGNED_MALLOC(c->in_len, 64);
	c->buf[1] = (unsigned char *)ALIGNED_MALLOC(c->in_len, 64);
	c->tmp_weights = (int *)ALIGNED_MALLOC((size_t)ratio * sizeof(int), 64);
	if(!c->buf[0] || !c->buf[1] || !c->tmp_weights) FAIL("out of memory");

	if(ocl_init(q->ocl_device) != 0) FAIL("OpenCL initialisation failed (see --ocl-list of temporal_scaling.exe)");

	if(ts_alloc_map(&c->map, out_w * ratio, out_h * ratio) != 0) FAIL("out of memory");
	compute_smoothstep_lut();
	ts_weights_for_mode(&c->map, in_w, in_h, out_w, out_h, ratio, q->mode, q->curve_base);
	if(c->need_map_uv) {
		if(ts_alloc_map(&c->map_uv, c->uv_out_w * ratio, c->uv_out_h * ratio) != 0) FAIL("out of memory");
		ts_weights_for_mode(&c->map_uv, c->uv_in_w, c->uv_in_h, c->uv_out_w, c->uv_out_h, ratio, q->mode, q->curve_base);
	}

	// LUT radiales selon le mode
	if(q->mode == 3 || q->mode == 4) {
		compute_curve_rad_lut(q->mode, q->curve_base);
	} else if(q->mode >= 5) {
		compute_kaiser4_lut(q->kaiser_beta);
		compute_kaiser_rad_lut_into(kaiser_rad_lut, q->kaiser_beta, q->curve_base);
		compute_kaiser_rad_lut_into(kaiser_rad_lut_soft, q->kaiser_beta, pow(q->curve_base, 1.0 / 3.0));
		compute_curve_rad_lut(4, q->curve_base);
		if(q->mode == 7) compute_mid_rad_lut(q->kaiser_beta);
	}
	compute_temporal_weights(c->tmp_weights, ratio, (TmpMode)q->tmp_mode, q->curve_base);
	for(int r = 0; r < ratio; r++) if(c->tmp_weights[r] > 0) c->blend = 1;

	OclSetup os;
	memset(&os, 0, sizeof(os));
	os.mode        = q->mode;
	os.ratio       = ratio;
	os.bps_in      = 1;
	os.bps_out     = 1;
	os.aniso_fp    = (int)(q->aniso_strength * WEIGHT_MAX + 0.5);
	os.in_len      = c->in_len;
	os.out_len     = c->out_len;
	os.rgb_out     = 1;
	os.yuv_709     = q->yuv_709;
	os.yuv_full    = q->yuv_full;
	os.maps[0]     = &c->map;
	os.maps[1]     = c->need_map_uv ? &c->map_uv : &c->map;
	os.tmp_weights = c->tmp_weights;
	os.lut_size    = KAISER_RAD_SIZE;
	if(q->mode == 3 || q->mode == 4)  os.lutA = curve_rad_lut;
	else if(q->mode == 5)             os.lutA = kaiser_rad_lut_soft;
	else if(q->mode >= 6) {
		os.lutA = kaiser_rad_lut;
		os.krad = kaiser_rad_lut;
		os.krad_size = KAISER_RAD_SIZE;
		if(q->mode == 7) { os.lutB = curve_rad_lut; os.lutC = mid_rad_lut; }
	}

	const size_t in_y_sz   = (size_t)in_w * in_h;
	const size_t in_uv_sz  = (size_t)c->uv_in_w * c->uv_in_h;
	const size_t out_y_sz  = (size_t)out_w * out_h * ratio;
	const size_t out_uv_sz = (size_t)c->uv_out_w * c->uv_out_h * ratio;
	os.nplanes = 3;
	for(int pl = 0; pl < 3; pl++) {
		OclPlane *o = &os.planes[pl];
		o->chan    = 1;
		o->map_idx = (pl == 0) ? 0 : 1;
		o->in_w  = (pl == 0) ? in_w  : c->uv_in_w;
		o->in_h  = (pl == 0) ? in_h  : c->uv_in_h;
		o->out_w = (pl == 0) ? out_w : c->uv_out_w;
		o->out_h = (pl == 0) ? out_h : c->uv_out_h;
		o->in_off  = (pl == 0) ? 0 : (pl == 1) ? in_y_sz  : in_y_sz  + in_uv_sz;
		o->out_off = (pl == 0) ? 0 : (pl == 1) ? out_y_sz : out_y_sz + out_uv_sz;
		o->out_step = (size_t)o->out_w * o->out_h;
	}
	(void)in_fmt;
	if(ocl_setup(&os) != 0) FAIL("OpenCL setup failed (not enough GPU memory?)");
	return c;
#undef FAIL
}

const char *ts_device_name(void) { return g_ocl.dev_name; }

size_t ts_frame_bytes(const TsCtx *c) { return c->out_frame; }

int ts_push(TsCtx *c, const uint8_t *const planes[3], const int pitches[3], uint8_t *out)
{
	const TsParams *q = &c->p;
	unsigned char *dst = c->buf[c->cur];
	const int pw[3] = { q->in_w, c->uv_in_w, c->uv_in_w };
	const int ph[3] = { q->in_h, c->uv_in_h, c->uv_in_h };
	for(int pl = 0; pl < 3; pl++) {
		for(int y = 0; y < ph[pl]; y++)
			memcpy(dst + (size_t)y * pw[pl], planes[pl] + (size_t)y * pitches[pl], (size_t)pw[pl]);
		dst += (size_t)pw[pl] * ph[pl];
	}

	const unsigned char *cur = c->buf[c->cur];
	int rc;
	if(!c->blend) {
		rc = ocl_process_frame(cur, cur, out);
	} else if(!c->have_prev) {
		rc = ocl_process_frame(cur, cur, out);          // première image : sans fondu
		c->have_prev = 1;
		c->cur ^= 1;
	} else {
		const unsigned char *prev = c->buf[c->cur ^ 1];
		rc = ocl_process_frame(prev, cur, out);          // images entre précédente et courante
		c->cur ^= 1;
	}
	return rc == 0 ? 0 : -1;
}

void ts_close(TsCtx *c)
{
	if(!c) return;
	ocl_release();
	ts_free_map(&c->map);
	ts_free_map(&c->map_uv);
	ALIGNED_FREE(c->buf[0]);
	ALIGNED_FREE(c->buf[1]);
	ALIGNED_FREE(c->tmp_weights);
	free(c);
}

int ts_check(const TsParams *p, char *err, size_t errlen)
{
#define BAD(...) do { if(err && errlen) snprintf(err, errlen, __VA_ARGS__); return -1; } while(0)
	if(p->mode < 0 || p->mode > 7) BAD("mode %d is not supported (0..7)", p->mode);
	if(p->tmp_mode < 0 || p->tmp_mode > 3) BAD("tmp-mode %d is not supported (0..3)", p->tmp_mode);
	if(p->in_w <= 0 || p->in_h <= 0) BAD("invalid input size");
	if(p->chroma != TS_CHROMA_444 && (p->in_w & 1)) BAD("input width %d must be even", p->in_w);
	if(p->chroma == TS_CHROMA_420 && (p->in_h & 1)) BAD("input height %d must be even", p->in_h);
	if(p->out_w <= p->in_w || p->out_h <= p->in_h)
		BAD("output %dx%d must be larger than the input %dx%d (the scaler only enlarges)", p->out_w, p->out_h, p->in_w, p->in_h);
	if(ocl_load_library() != 0) BAD("OpenCL.dll not found");
	if(ocl_enumerate() != 0) BAD("no OpenCL device found");
	if(p->ocl_device >= g_ocl.ndev) BAD("OpenCL device %d does not exist (%d found)", p->ocl_device, g_ocl.ndev);
	if(p->ocl_device < 0) {
		int gpu = 0;
		for(int i = 0; i < g_ocl.ndev; i++) { OclDevInfo o; ocl_dev_info(i, &o); if((o.type & OCL_DEVICE_TYPE_GPU) && o.avail && o.compiler) gpu = 1; }
		if(!gpu) BAD("no usable OpenCL GPU found");
	}
	return 0;
#undef BAD
}
