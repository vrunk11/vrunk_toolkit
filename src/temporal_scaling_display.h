// ============================================================
// temporal_scaling_display.h — affichage direct (Windows, OpenGL 1.1)
//
// Affiche les frames RGB produites par temporal_scaling dans une fenêtre
// native, sans passer par un lecteur externe (VLC) : plus de pipe vidéo vers
// un autre processus, plus de seconde copie des frames.
//
// Aucune dépendance en dehors de Windows : opengl32, gdi32, user32, winmm
// (liens : -lopengl32 -lgdi32 -luser32 -lwinmm).
//
// La fenêtre doit être créée, pompée et fermée par le même thread.
// Touches : Echap = quitter, F ou F11 = plein écran, Espace = pause.
// ============================================================
#ifndef TEMPORAL_SCALING_DISPLAY_H
#define TEMPORAL_SCALING_DISPLAY_H

#ifdef _WIN32

#include <windows.h>
#include <pthread.h>
#include <GL/gl.h>
#include <stdio.h>
#include <stdint.h>

#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif

typedef struct {
	HWND   hwnd;
	HDC    dc;
	HGLRC  rc;
	GLuint tex;
	int    w, h;            // taille de la frame
	int    bps;             // octets par échantillon (1 ou 2)
	int    quit;
	int    paused;
	int    fullscreen;
	RECT   win_rect;        // position avant plein écran
	LONG   win_style;
	double qpc_freq;
	double next_t;          // date prévue de la prochaine frame
	double period;          // durée d'une frame (s)
	int    frames_shown;
	double t_start;
	long   vframe;          // rang de la prochaine frame sur la ligne de temps (mode audio)
	int    dropped;
	int    started;         // première frame présentée (mode pts)
	int    audio_was;       // l'audio jouait : à ré-ancrer quand il se termine
	double first_pts, last_clk, pause_mark;
	double t_begin;         // début de la lecture (statistiques)
} DispState;

static DispState g_disp;

// retard moyen (s) mesuré à l'affichage : le décodeur en tient compte pour choisir des images à l'heure
static volatile double g_disp_late = 0.0;

static double disp_now(void)
{
	LARGE_INTEGER c;
	QueryPerformanceCounter(&c);
	return (double)c.QuadPart / g_disp.qpc_freq;
}

// ---------- audio (waveOut, PCM s16le 48 kHz stéréo) ----------
// Le son vient d'un flux brut (stdout d'un ffmpeg). L'horloge audio (échantillons
// réellement joués) sert d'horloge maître : l'image est affichée quand l'audio
// atteint son heure, et les images en retard sont sautées pour rester synchrone.
#define AUD_RATE    48000
#define AUD_CH      2
#define AUD_BUFS    8
#define AUD_BUF_MS  40
#define AUD_BUF_BYTES (AUD_RATE * AUD_BUF_MS / 1000 * AUD_CH * 2)

typedef struct {
	FILE      *src;
	size_t   (*rd)(void *ctx, void *dst, size_t n);   // lecture personnalisée (libav) ou NULL = fread(src)
	double   (*pts0_fn)(void *ctx);                    // pts (s) du premier échantillon, <0 = inconnu
	void      *rd_ctx;
	HWAVEOUT   hwo;
	WAVEHDR    hdr[AUD_BUFS];
	char      *mem[AUD_BUFS];
	pthread_t  th;
	volatile int run;
	volatile int eof;        // le flux est épuisé ET tout a été joué
	int        active;
	int        started;
	double     offset_s;     // décalage image/son réglable (--audio-offset)
} AudState;

static AudState g_aud;

static void *aud_thread(void *p)
{
	AudState *a = (AudState *)p;
	int i = 0;
	int ended = 0;
	while(a->run)
	{
		WAVEHDR *h = &a->hdr[i];
		// attendre que le tampon soit libre (jamais écrit, ou terminé)
		while(a->run && (h->dwFlags & WHDR_PREPARED) && !(h->dwFlags & WHDR_DONE))
			Sleep(2);
		if(!a->run) break;
		size_t got = a->rd ? a->rd(a->rd_ctx, a->mem[i], AUD_BUF_BYTES)
		           : (a->src ? fread(a->mem[i], 1, AUD_BUF_BYTES, a->src) : 0);
		if(got == 0) { ended = 1; break; }
		if(got < AUD_BUF_BYTES) memset(a->mem[i] + got, 0, AUD_BUF_BYTES - got);
		if(h->dwFlags & WHDR_PREPARED) waveOutUnprepareHeader(a->hwo, h, sizeof(*h));
		h->lpData = a->mem[i];
		h->dwBufferLength = AUD_BUF_BYTES;
		h->dwFlags = 0;
		waveOutPrepareHeader(a->hwo, h, sizeof(*h));
		waveOutWrite(a->hwo, h, sizeof(*h));
		i = (i + 1) % AUD_BUFS;
	}
	if(ended) {
		// laisser jouer ce qui reste avant de déclarer la fin
		for(int k = 0; k < AUD_BUFS && a->run; k++)
			while(a->run && (a->hdr[k].dwFlags & WHDR_PREPARED) && !(a->hdr[k].dwFlags & WHDR_DONE))
				Sleep(2);
		a->eof = 1;
	}
	return NULL;
}

// src : flux PCM s16le 48000 Hz 2 canaux. Démarre en pause (voir disp_audio_go).
static int disp_audio_init_ex(FILE *src, size_t (*rd)(void *, void *, size_t),
                              double (*pts0_fn)(void *), void *ctx, double offset_ms)
{
	AudState *a = &g_aud;
	memset(a, 0, sizeof(*a));
	a->src = src;
	a->rd = rd; a->pts0_fn = pts0_fn; a->rd_ctx = ctx;
	a->offset_s = offset_ms / 1000.0;

	WAVEFORMATEX wf;
	memset(&wf, 0, sizeof(wf));
	wf.wFormatTag      = WAVE_FORMAT_PCM;
	wf.nChannels       = AUD_CH;
	wf.nSamplesPerSec  = AUD_RATE;
	wf.wBitsPerSample  = 16;
	wf.nBlockAlign     = (WORD)(AUD_CH * 2);
	wf.nAvgBytesPerSec = AUD_RATE * AUD_CH * 2;
	if(waveOutOpen(&a->hwo, WAVE_MAPPER, &wf, 0, 0, CALLBACK_NULL) != MMSYSERR_NOERROR) {
		fprintf(stderr, "Audio: cannot open the audio device\n");
		return -1;
	}
	waveOutPause(a->hwo);
	for(int i = 0; i < AUD_BUFS; i++) {
		a->mem[i] = (char *)malloc(AUD_BUF_BYTES);
		if(!a->mem[i]) return -1;
	}
	a->run = 1;
	a->active = 1;
	pthread_create(&a->th, NULL, aud_thread, a);
	pthread_detach(a->th);
	return 0;
}

// Heure de lecture audio (s), ou -1 si indisponible.
static double disp_audio_time(void)
{
	MMTIME mmt;
	mmt.wType = TIME_SAMPLES;
	if(waveOutGetPosition(g_aud.hwo, &mmt, sizeof(mmt)) != MMSYSERR_NOERROR || mmt.wType != TIME_SAMPLES)
		return -1.0;
	return (double)mmt.u.sample / (double)AUD_RATE;
}

static void disp_audio_close(void)
{
	AudState *a = &g_aud;
	if(!a->active) return;
	a->run = 0;
	Sleep(10);
	waveOutReset(a->hwo);
	for(int i = 0; i < AUD_BUFS; i++)
		if(a->hdr[i].dwFlags & WHDR_PREPARED) waveOutUnprepareHeader(a->hwo, &a->hdr[i], sizeof(a->hdr[i]));
	waveOutClose(a->hwo);
	a->active = 0;
}

static void disp_toggle_fullscreen(void)
{
	DispState *d = &g_disp;
	if(!d->fullscreen) {
		d->win_style = GetWindowLong(d->hwnd, GWL_STYLE);
		GetWindowRect(d->hwnd, &d->win_rect);
		MONITORINFO mi; mi.cbSize = sizeof(mi);
		GetMonitorInfo(MonitorFromWindow(d->hwnd, MONITOR_DEFAULTTONEAREST), &mi);
		SetWindowLong(d->hwnd, GWL_STYLE, WS_POPUP | WS_VISIBLE);
		SetWindowPos(d->hwnd, HWND_TOP, mi.rcMonitor.left, mi.rcMonitor.top,
			mi.rcMonitor.right - mi.rcMonitor.left,
			mi.rcMonitor.bottom - mi.rcMonitor.top, SWP_FRAMECHANGED);
		d->fullscreen = 1;
	} else {
		SetWindowLong(d->hwnd, GWL_STYLE, d->win_style);
		SetWindowPos(d->hwnd, HWND_NOTOPMOST, d->win_rect.left, d->win_rect.top,
			d->win_rect.right - d->win_rect.left,
			d->win_rect.bottom - d->win_rect.top, SWP_FRAMECHANGED);
		d->fullscreen = 0;
	}
}

static LRESULT CALLBACK disp_wndproc(HWND h, UINT m, WPARAM w, LPARAM l)
{
	switch(m) {
		case WM_CLOSE:
		case WM_DESTROY:
			g_disp.quit = 1;
			return 0;
		case WM_KEYDOWN:
			if(w == VK_ESCAPE) g_disp.quit = 1;
			else if(w == VK_F11 || w == 'F') disp_toggle_fullscreen();
			else if(w == VK_SPACE) {
				g_disp.paused = !g_disp.paused;
				if(g_aud.active && g_aud.started) {
					if(g_disp.paused) waveOutPause(g_aud.hwo); else waveOutRestart(g_aud.hwo);
				}
			}
			return 0;
		case WM_ERASEBKGND:
			return 1;
	}
	return DefWindowProcA(h, m, w, l);
}

// Traite les messages de la fenêtre ; renvoie 0 si l'utilisateur a quitté.
static int disp_pump(void)
{
	MSG msg;
	while(PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE)) {
		TranslateMessage(&msg);
		DispatchMessageA(&msg);
	}
	return !g_disp.quit;
}

// fps : cadence d'affichage des frames de sortie. bps : 1 (RGB) ou 2 (RGB48).
static int disp_init(int w, int h, int bps, double fps)
{
	DispState *d = &g_disp;
	memset(d, 0, sizeof(*d));
	d->w = w; d->h = h; d->bps = bps;
	LARGE_INTEGER f; QueryPerformanceFrequency(&f);
	d->qpc_freq = (double)f.QuadPart;
	d->period   = 1.0 / (fps > 0.0 ? fps : 60.0);
	timeBeginPeriod(1);

	WNDCLASSA wc;
	memset(&wc, 0, sizeof(wc));
	wc.lpfnWndProc   = disp_wndproc;
	wc.hInstance     = GetModuleHandleA(NULL);
	wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
	wc.lpszClassName = "TemporalScalingDisplay";
	if(!RegisterClassA(&wc)) { fprintf(stderr, "Display: RegisterClass failed\n"); return -1; }

	// fenêtre initiale : 1280 px de large au plus, même rapport d'aspect
	int ww = w > 1280 ? 1280 : w;
	int wh = (int)((double)ww * h / w + 0.5);
	RECT r = { 0, 0, ww, wh };
	AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE);
	d->hwnd = CreateWindowExA(0, wc.lpszClassName, "temporal_scaling",
		WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT,
		r.right - r.left, r.bottom - r.top, NULL, NULL, wc.hInstance, NULL);
	if(!d->hwnd) { fprintf(stderr, "Display: cannot create window\n"); return -1; }
	d->dc = GetDC(d->hwnd);

	PIXELFORMATDESCRIPTOR pfd;
	memset(&pfd, 0, sizeof(pfd));
	pfd.nSize      = sizeof(pfd);
	pfd.nVersion   = 1;
	pfd.dwFlags    = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
	pfd.iPixelType = PFD_TYPE_RGBA;
	pfd.cColorBits = 24;
	int pf = ChoosePixelFormat(d->dc, &pfd);
	if(!pf || !SetPixelFormat(d->dc, pf, &pfd)) { fprintf(stderr, "Display: no OpenGL pixel format\n"); return -1; }
	d->rc = wglCreateContext(d->dc);
	if(!d->rc || !wglMakeCurrent(d->dc, d->rc)) { fprintf(stderr, "Display: cannot create OpenGL context\n"); return -1; }

	// vsync si disponible (évite le tearing)
	typedef BOOL (WINAPI *swapint_t)(int);
	swapint_t swapint = (swapint_t)wglGetProcAddress("wglSwapIntervalEXT");
	if(swapint) swapint(1);

	GLint maxtex = 0;
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxtex);
	if(maxtex > 0 && (w > maxtex || h > maxtex)) {
		fprintf(stderr, "Display: frame %dx%d exceeds the GPU texture limit (%d)\n", w, h, (int)maxtex);
		return -1;
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &d->tex);
	glBindTexture(GL_TEXTURE_2D, d->tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB,
		bps == 1 ? GL_UNSIGNED_BYTE : GL_UNSIGNED_SHORT, NULL);
	if(glGetError() != GL_NO_ERROR) { fprintf(stderr, "Display: cannot allocate the %dx%d texture\n", w, h); return -1; }

	d->next_t = d->t_start = disp_now();
	return 0;
}

// Raccourci : flux PCM brut (FILE*) lu par fread.
static int disp_audio_init(FILE *src, double offset_ms)
{
	return disp_audio_init_ex(src, NULL, NULL, NULL, offset_ms);
}

// Horloge média (s) = temps de la vidéo actuellement "jouée" : celui de l'audio
// quand il joue, sinon l'horloge système ancrée sur la première image.
// Renvoie -1 tant que l'affichage n'a pas démarré. Appelable depuis un autre thread.
static double disp_media_clock(void)
{
	DispState *d = &g_disp;
	if(!d->started) return -1.0;
	if(g_aud.active && g_aud.started && !g_aud.eof) {
		double at = disp_audio_time();
		if(at >= 0.0) {
			double p0 = g_aud.pts0_fn ? g_aud.pts0_fn(g_aud.rd_ctx) : -1.0;
			if(p0 < 0.0) p0 = d->first_pts;
			d->last_clk = p0 + at;
			return d->last_clk;
		}
	}
	if(d->audio_was) {              // l'audio vient de finir : on continue à l'horloge système
		d->first_pts = d->last_clk >= 0.0 ? d->last_clk : d->first_pts;
		d->t_start = disp_now();
		d->audio_was = 0;
	}
	return d->first_pts + (disp_now() - d->t_start);
}

// Affiche une frame (w*h*3*bps octets). Pendant la pause, la fenêtre reste réactive.
// pts >= 0 : horodatage de la frame en secondes (mode lecteur : calée sur l'horloge
//            média, frames en retard de plus de period_out abandonnées).
// pts <  0 : mode flux : cadence fixe 1/fps (+ audio maître si présent).
// Renvoie 0 si l'utilisateur a quitté.
static int disp_present(const void *pix, double pts, double period_out)
{
	DispState *d = &g_disp;
if(d->t_begin == 0.0) d->t_begin = disp_now();

	if(pts >= 0.0)
	{
		const double due = pts + (g_aud.active ? g_aud.offset_s : 0.0);
		if(!d->started) {
			d->started = 1;
			d->first_pts = pts;
			d->last_clk = -1.0;
			d->t_start = disp_now();
			if(g_aud.active && !g_aud.started) {
				waveOutRestart(g_aud.hwo);
				g_aud.started = 1;
				d->audio_was = 1;
			}
		}
		double clk = 0.0;
		for(;;) {
			if(!disp_pump()) return 0;
			if(d->paused) {
				if(d->pause_mark == 0.0) d->pause_mark = disp_now();
				Sleep(10);
				continue;
			}
			if(d->pause_mark != 0.0) { d->t_start += disp_now() - d->pause_mark; d->pause_mark = 0.0; }
			clk = disp_media_clock();
			double dt = due - clk;
			if(dt <= 0.0) break;
			if(dt > 0.002) Sleep((DWORD)((dt - 0.001) * 1000.0)); else Sleep(0);
		}
{ double lt = clk - due; if(lt < 0.0) lt = 0.0; if(lt > 0.5) lt = 0.5; g_disp_late = 0.9 * g_disp_late + 0.1 * lt; }
		if(clk - due > (period_out > 0.10 ? period_out : 0.10)) { d->dropped++; return 1; }   // trop en retard : on saute
	}
	else

	if(g_aud.active && !g_aud.eof)
	{
		// ---- horloge maître = audio ----
		// L'audio ne démarre qu'à la première frame (il est préchargé, en pause).
		if(!g_aud.started) {
			waveOutRestart(g_aud.hwo);
			g_aud.started = 1;
			d->t_start = disp_now();
		}
		// heure d'affichage de cette frame sur la ligne de temps de l'audio
		const double due = (double)d->vframe * d->period + g_aud.offset_s;
		d->vframe++;
		double clock = 0.0;
		for(;;) {
			if(!disp_pump()) return 0;
			if(d->paused) { Sleep(10); continue; }
			if(g_aud.eof) break;                 // audio fini : on passe à l'horloge système
			clock = disp_audio_time();
			if(clock < 0.0) break;
			double dt = due - clock;
			if(dt <= 0.0) break;
			if(dt > 0.002) Sleep((DWORD)((dt - 0.001) * 1000.0)); else Sleep(0);
		}
		if(!g_aud.eof && clock >= 0.0 && clock - due > d->period) {
			d->dropped++;                        // en retard : on saute l'image pour rester synchrone
			return 1;
		}
		if(g_aud.eof) d->next_t = disp_now();    // re-ancrage pour la suite
	}
	else
	{
		// cadence : attendre la date prévue (sommeil grossier puis attente fine)
		for(;;) {
			if(!disp_pump()) return 0;
			if(d->paused) { Sleep(10); d->next_t = disp_now(); continue; }
			double dt = d->next_t - disp_now();
			if(dt <= 0.0) break;
			if(dt > 0.002) Sleep((DWORD)((dt - 0.001) * 1000.0));
		}
		double now = disp_now();
		// en retard de plus d'une frame : on repart de maintenant (pas de rattrapage)
		d->next_t = (now - d->next_t > d->period) ? now + d->period : d->next_t + d->period;
	}

	RECT rc; GetClientRect(d->hwnd, &rc);
	int cw = rc.right - rc.left, ch = rc.bottom - rc.top;
	if(cw < 1 || ch < 1) return 1;

	glBindTexture(GL_TEXTURE_2D, d->tex);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, d->w, d->h, GL_RGB,
		d->bps == 1 ? GL_UNSIGNED_BYTE : GL_UNSIGNED_SHORT, pix);

	// rapport d'aspect conservé (bandes noires)
	double sa = (double)d->w / d->h, wa = (double)cw / ch;
	int vw = cw, vh = ch;
	if(wa > sa) vw = (int)(ch * sa + 0.5); else vh = (int)(cw / sa + 0.5);
	glViewport(0, 0, cw, ch);
	glClearColor(0.f, 0.f, 0.f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT);
	glViewport((cw - vw) / 2, (ch - vh) / 2, vw, vh);
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
		glTexCoord2f(0.f, 0.f); glVertex2f(-1.f,  1.f);
		glTexCoord2f(1.f, 0.f); glVertex2f( 1.f,  1.f);
		glTexCoord2f(1.f, 1.f); glVertex2f( 1.f, -1.f);
		glTexCoord2f(0.f, 1.f); glVertex2f(-1.f, -1.f);
	glEnd();
	SwapBuffers(d->dc);
	d->frames_shown++;
	return 1;
}

static void disp_close(void)
{
	DispState *d = &g_disp;
	if(d->frames_shown > 0) {
		double t = disp_now() - d->t_begin;
		fprintf(stderr, "Display: %d frames in %.1f s (%.1f fps)", d->frames_shown, t, t > 0 ? d->frames_shown / t : 0.0);
		if(g_aud.active) fprintf(stderr, ", %d dropped to stay in sync", d->dropped);
		fprintf(stderr, "\n");
	}
	disp_audio_close();
	if(d->rc) { wglMakeCurrent(NULL, NULL); wglDeleteContext(d->rc); d->rc = NULL; }
	if(d->hwnd) { ReleaseDC(d->hwnd, d->dc); DestroyWindow(d->hwnd); d->hwnd = NULL; }
	timeEndPeriod(1);
}

#endif // _WIN32
#endif // TEMPORAL_SCALING_DISPLAY_H
