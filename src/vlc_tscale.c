/*****************************************************************************
 * vlc_tscale.c : sortie vidéo VLC 3.0 "tscale" — mise à l'échelle sur GPU (OpenCL)
 *
 * VLC reste le lecteur (interface, codecs, audio, synchronisation). Ce module de
 * sortie vidéo reçoit les images décodées, les met à l'échelle avec le moteur de
 * temporal_scaling (noyaux Kaiser, anisotrope, cadence multipliée, fondu temporel)
 * sur le GPU choisi, puis les affiche lui-même (OpenGL) dans la fenêtre de VLC.
 *
 * Pourquoi une sortie vidéo et pas un filtre vidéo ? Dans VLC 3.0 le cœur de la
 * sortie vidéo ramène toujours la sortie d'un filtre à la taille d'origine (il
 * insère un swscale "pour compenser les changements de format"). Une sortie vidéo
 * n'a pas cette contrainte.
 *
 * Sélection : Outils > Préférences > Vidéo > Sortie > "Temporal scaling (GPU)", ou
 *   vlc --vout=tscale --tscale-height=2160 --tscale-mode=5 fichier
 *****************************************************************************/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <vlc_common.h>
#include <vlc_plugin.h>
#include <vlc_vout_display.h>
#include <vlc_picture_pool.h>
#include <vlc_actions.h>
#include <vlc_mouse.h>

#include <windows.h>
#include <GL/gl.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ts_core.h"

#ifndef N_
# define N_(str) (str)
#endif

#ifndef GL_CLAMP_TO_EDGE
# define GL_CLAMP_TO_EDGE 0x812F
#endif

static int  Open (vlc_object_t *);
static void Close(vlc_object_t *);

#define HEIGHT_TEXT     N_("Output height")
#define HEIGHT_LONGTEXT N_("Vertical resolution the picture is scaled to on the GPU " \
    "(the width follows the display aspect ratio). 0 = height of the screen. " \
    "The scaler only enlarges: it must be larger than the video.")
#define MODE_TEXT       N_("Scaling algorithm")
#define MODE_LONGTEXT   N_("0 nearest, 1 bilinear on duplicated pixels, 2 bilinear with temporal " \
    "alternation, 3 linear weighted 2D, 4 exponential weighted 2D, 5 Kaiser 6x6 EWA, " \
    "6 Kaiser anisotropic edge-directed, 7 Kaiser anisotropic + selective softening.")
#define RATIO_TEXT      N_("Frame rate multiplier")
#define RATIO_LONGTEXT  N_("Number of output pictures shown per decoded picture (2 = 30 fps shown at 60 fps).")
#define TMP_TEXT        N_("Temporal blending mode")
#define TMP_LONGTEXT    N_("0 duplicate pictures, 1 blend 50% on the upper phases, 2 linear, " \
    "3 exponential. Anything but 0 delays the picture by one frame.")
#define CURVE_TEXT      N_("Curve base")
#define CURVE_LONGTEXT  N_("Base of the exponential curves (modes 4-7, temporal mode 3). Must be > 1.")
#define BETA_TEXT       N_("Kaiser beta")
#define BETA_LONGTEXT   N_("Kaiser window beta (modes 5-7): low = sharp with ringing, high = smooth.")
#define ANISO_TEXT      N_("Anisotropy strength")
#define ANISO_LONGTEXT  N_("Edge-directed strength for modes 6 and 7 (0 = isotropic, 1 = full).")
#define DEVICE_TEXT     N_("OpenCL device")
#define DEVICE_LONGTEXT N_("Index of the OpenCL device (see temporal_scaling.exe --ocl-list). -1 = best GPU.")

vlc_module_begin()
    set_shortname("tscale")
    set_description(N_("Temporal scaling (OpenCL GPU)"))
    set_category(CAT_VIDEO)
    set_subcategory(SUBCAT_VIDEO_VOUT)
    set_capability("vout display", 0)
    add_shortcut("tscale")
    add_integer("tscale-height", 2160, HEIGHT_TEXT, HEIGHT_LONGTEXT, false)
    add_integer_with_range("tscale-mode", 5, 0, 7, MODE_TEXT, MODE_LONGTEXT, false)
    add_integer_with_range("tscale-ratio", 1, 1, 8, RATIO_TEXT, RATIO_LONGTEXT, false)
    add_integer_with_range("tscale-tmp-mode", 0, 0, 3, TMP_TEXT, TMP_LONGTEXT, true)
    add_float("tscale-curve", 2.0, CURVE_TEXT, CURVE_LONGTEXT, true)
    add_float("tscale-beta", 2.5, BETA_TEXT, BETA_LONGTEXT, true)
    add_float("tscale-aniso", 1.0, ANISO_TEXT, ANISO_LONGTEXT, true)
    add_integer("tscale-device", -1, DEVICE_TEXT, DEVICE_LONGTEXT, true)
    set_callbacks(Open, Close)
vlc_module_end()

typedef int (WINAPI *swapint_t)(int);

struct vout_display_sys_t {
    TsCtx          *ts;
    TsParams        params;         /* paramètres du moteur (ouvert à la première image) */
    int             ts_failed;
    int             ratio;
    int             swap_uv;        /* YV12 : plans U et V inversés */
    int             x_off, y_off;   /* recadrage éventuel de la source */
    int             out_w, out_h;
    /* anneau d'images prêtes (ratio images RGB24 chacune) : le thread vidéo de VLC écrit,
     * le thread de rendu présente */
#define NSLOT 3
    struct { uint8_t *data; mtime_t date; int n; } slot[NSLOT];
    int             wr, rd, count, stop;
    CRITICAL_SECTION   cs;
    CONDITION_VARIABLE cv;
    HANDLE          render, gl_event;
    int             gl_failed;
    size_t          frame_bytes;
    int             frames_valid;
    picture_pool_t *pool;

    vout_window_t  *window;         /* fenêtre VOUT_WINDOW_TYPE_HWND fournie par VLC */
    HWND            hwnd;
    /* fenêtre propre (quand VLC n'en fournit pas, ex. sans interface Qt) */
    int             own_window;
    HANDLE          thread, ready;
    DWORD           thread_id;
    int             fullscreen;
    RECT            win_rect;
    LONG            win_style;
    vout_display_t *vd;
    HDC             dc;
    HGLRC           rc;
    GLuint          tex;

    mtime_t         last_date;      /* pour estimer la durée d'une image */
    mtime_t         period;
};


/* ---- fenêtre propre : créée et pompée par un thread dédié (comme le module Direct3D11) ---- */

static int VlcKeyFromMessage(WPARAM w)
{
    switch (w) {
        case VK_ESCAPE: return KEY_ESC;
        case VK_LEFT:   return KEY_LEFT;
        case VK_RIGHT:  return KEY_RIGHT;
        case VK_UP:     return KEY_UP;
        case VK_DOWN:   return KEY_DOWN;
        case VK_F11:    return KEY_F11;
        case VK_HOME:   return KEY_HOME;
        case VK_BACK:   return KEY_BACKSPACE;
        case VK_TAB:    return KEY_TAB;
        case VK_RETURN: return KEY_ENTER;
    }
    return 0;
}

static LRESULT CALLBACK OwnWndProc(HWND h, UINT m, WPARAM w, LPARAM l)
{
    if (m == WM_NCCREATE)
        SetWindowLongPtr(h, GWLP_USERDATA, (LONG_PTR)((CREATESTRUCT *)l)->lpCreateParams);
    vout_display_sys_t *sys = (vout_display_sys_t *)GetWindowLongPtr(h, GWLP_USERDATA);
    vout_display_t *vd = sys ? sys->vd : NULL;

    switch (m) {
    case WM_APP:                       /* demande de fermeture venant de Close() */
        DestroyWindow(h);
        return 0;
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    case WM_CLOSE:
        if (vd) vout_display_SendEventClose(vd);
        return 0;
    case WM_SIZE:
        if (vd && w != SIZE_MINIMIZED)
            vout_display_SendEventDisplaySize(vd, LOWORD(l), HIWORD(l));
        break;
    case WM_KEYDOWN:
    case WM_SYSKEYDOWN: {
        int key = VlcKeyFromMessage(w);
        if (key && vd) {
            if (GetKeyState(VK_CONTROL) & 0x8000) key |= KEY_MODIFIER_CTRL;
            vout_display_SendEventKey(vd, key);
        }
        break;
    }
    case WM_CHAR:
        if (vd && w >= 0x20 && w < 0x7f) {
            int key = (int)w;
            if (key >= 'A' && key <= 'Z') key += 'a' - 'A';
            vout_display_SendEventKey(vd, key);
        }
        return 0;
    case WM_MOUSEWHEEL:
        if (vd) vout_display_SendEventKey(vd, (short)HIWORD(w) > 0 ? KEY_MOUSEWHEELUP : KEY_MOUSEWHEELDOWN);
        return 0;
    case WM_LBUTTONDBLCLK:
        if (vd) vout_display_SendEventMouseDoubleClick(vd);
        return 0;
    case WM_LBUTTONDOWN:
        if (vd) vout_display_SendEventMousePressed(vd, MOUSE_BUTTON_LEFT);
        return 0;
    case WM_LBUTTONUP:
        if (vd) vout_display_SendEventMouseReleased(vd, MOUSE_BUTTON_LEFT);
        return 0;
    case WM_MOUSEMOVE:
        if (vd) vout_display_SendEventMouseMoved(vd, (short)LOWORD(l), (short)HIWORD(l));
        return 0;
    case WM_ERASEBKGND:
        return 1;
    }
    return DefWindowProc(h, m, w, l);
}

static DWORD WINAPI OwnWindowThread(LPVOID arg)
{
    vout_display_sys_t *sys = (vout_display_sys_t *)arg;
    WNDCLASS wc;
    memset(&wc, 0, sizeof(wc));
    wc.style         = CS_DBLCLKS | CS_OWNDC;
    wc.lpfnWndProc   = OwnWndProc;
    wc.hInstance     = GetModuleHandle(NULL);
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wc.lpszClassName = TEXT("VLCTemporalScaling");
    RegisterClass(&wc);

    int ww = sys->out_w > 1280 ? 1280 : sys->out_w;
    int wh = (int)((double)ww * sys->out_h / sys->out_w + 0.5);
    RECT r = { 0, 0, ww, wh };
    AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE);
    sys->hwnd = CreateWindowEx(0, wc.lpszClassName, TEXT("VLC (temporal scaling)"),
                               WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT,
                               r.right - r.left, r.bottom - r.top, NULL, NULL, wc.hInstance, sys);
    SetEvent(sys->ready);
    if (!sys->hwnd) return 0;

    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    UnregisterClass(wc.lpszClassName, wc.hInstance);
    return 0;
}

static int OwnWindowOpen(vout_display_sys_t *sys)
{
    sys->ready  = CreateEvent(NULL, TRUE, FALSE, NULL);
    sys->thread = CreateThread(NULL, 0, OwnWindowThread, sys, 0, &sys->thread_id);
    if (!sys->thread || !sys->ready) return -1;
    WaitForSingleObject(sys->ready, 5000);
    return sys->hwnd ? 0 : -1;
}

static void OwnWindowClose(vout_display_sys_t *sys)
{
    if (sys->thread) {
        if (sys->hwnd) PostMessage(sys->hwnd, WM_APP, 0, 0);
        WaitForSingleObject(sys->thread, 3000);
        CloseHandle(sys->thread);
        sys->thread = NULL;
    }
    if (sys->ready) { CloseHandle(sys->ready); sys->ready = NULL; }
}

static void OwnFullscreen(vout_display_sys_t *sys, int fs)
{
    if (!sys->own_window || !sys->hwnd || fs == sys->fullscreen) return;
    if (fs) {
        sys->win_style = GetWindowLong(sys->hwnd, GWL_STYLE);
        GetWindowRect(sys->hwnd, &sys->win_rect);
        MONITORINFO mi; mi.cbSize = sizeof(mi);
        GetMonitorInfo(MonitorFromWindow(sys->hwnd, MONITOR_DEFAULTTONEAREST), &mi);
        SetWindowLong(sys->hwnd, GWL_STYLE, WS_POPUP | WS_VISIBLE);
        SetWindowPos(sys->hwnd, HWND_TOP, mi.rcMonitor.left, mi.rcMonitor.top,
                     mi.rcMonitor.right - mi.rcMonitor.left,
                     mi.rcMonitor.bottom - mi.rcMonitor.top, SWP_FRAMECHANGED);
    } else {
        SetWindowLong(sys->hwnd, GWL_STYLE, sys->win_style);
        SetWindowPos(sys->hwnd, HWND_NOTOPMOST, sys->win_rect.left, sys->win_rect.top,
                     sys->win_rect.right - sys->win_rect.left,
                     sys->win_rect.bottom - sys->win_rect.top, SWP_FRAMECHANGED);
    }
    sys->fullscreen = fs;
}
static picture_pool_t *Pool(vout_display_t *vd, unsigned count)
{
    vout_display_sys_t *sys = vd->sys;
    if (!sys->pool)
        sys->pool = picture_pool_NewFromFormat(&vd->fmt, count);
    return sys->pool;
}

/* ---- OpenGL ---- */


/* Le contexte OpenGL est créé ici mais utilisé par le thread de rendu (RenderThread) :
 * toutes les présentations (avec attente de la synchro verticale) s'y font, pour ne
 * jamais bloquer le thread de sortie vidéo de VLC. */
static int GlOpen(vout_display_t *vd, vout_display_sys_t *sys)
{
    /* sys->hwnd est déjà renseigné (fenêtre VLC embarquée ou fenêtre propre) */
    sys->dc   = GetDC(sys->hwnd);
    if (!sys->dc) return -1;

    PIXELFORMATDESCRIPTOR pfd;
    memset(&pfd, 0, sizeof(pfd));
    pfd.nSize      = sizeof(pfd);
    pfd.nVersion   = 1;
    pfd.dwFlags    = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    int pf = ChoosePixelFormat(sys->dc, &pfd);
    if (!pf || !SetPixelFormat(sys->dc, pf, &pfd)) {
        msg_Err(vd, "no OpenGL pixel format for the video window");
        return -1;
    }
    sys->rc = wglCreateContext(sys->dc);
    if (!sys->rc) {
        msg_Err(vd, "cannot create the OpenGL context");
        return -1;
    }
    return 0;
}

static void GlPresent(vout_display_sys_t *sys, const uint8_t *rgb)
{
    RECT rc;
    GetClientRect(sys->hwnd, &rc);
    int cw = rc.right - rc.left, ch = rc.bottom - rc.top;
    if (cw < 1 || ch < 1) return;

    glBindTexture(GL_TEXTURE_2D, sys->tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, sys->out_w, sys->out_h, GL_RGB, GL_UNSIGNED_BYTE, rgb);

    /* rapport d'aspect conservé (bandes noires) */
    double sa = (double)sys->out_w / sys->out_h, wa = (double)cw / ch;
    int vw = cw, vh = ch;
    if (wa > sa) vw = (int)(ch * sa + 0.5); else vh = (int)(cw / sa + 0.5);
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
    /* débogage : TSCALE_SNAPSHOT=fichier.ppm enregistre la 30e image telle que dessinée dans la fenêtre */
    {
        static int count;
        const char *snap = getenv("TSCALE_SNAPSHOT");
        if (snap && ++count == 30) {
            uint8_t *px = malloc((size_t)cw * ch * 3);
            if (px) {
                glReadBuffer(GL_BACK);
                glReadPixels(0, 0, cw, ch, GL_RGB, GL_UNSIGNED_BYTE, px);
                FILE *f = fopen(snap, "wb");
                if (f) {
                    fprintf(f, "P6\n%d %d\n255\n", cw, ch);
                    for (int y = ch - 1; y >= 0; y--) fwrite(px + (size_t)y * cw * 3, 3, cw, f);
                    fclose(f);
                }
                free(px);
            }
        }
    }
    SwapBuffers(sys->dc);
}

/* ---- thread de rendu ---- */

static int RenderInitGl(vout_display_sys_t *sys)
{
    if (!wglMakeCurrent(sys->dc, sys->rc)) return -1;
    swapint_t swapint = (swapint_t)wglGetProcAddress("wglSwapIntervalEXT");
    if (swapint) swapint(1);

    GLint maxtex = 0;
    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxtex);
    if (maxtex > 0 && (sys->out_w > maxtex || sys->out_h > maxtex)) return -2;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &sys->tex);
    glBindTexture(GL_TEXTURE_2D, sys->tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, sys->out_w, sys->out_h, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    return glGetError() == GL_NO_ERROR ? 0 : -3;
}

static DWORD WINAPI RenderThread(LPVOID arg)
{
    vout_display_sys_t *sys = (vout_display_sys_t *)arg;
    sys->gl_failed = RenderInitGl(sys);
    SetEvent(sys->gl_event);
    if (sys->gl_failed) return 0;

    for (;;) {
        EnterCriticalSection(&sys->cs);
        while (sys->count == 0 && !sys->stop)
            SleepConditionVariableCS(&sys->cv, &sys->cs, INFINITE);
        if (sys->stop) { LeaveCriticalSection(&sys->cs); break; }
        const int i = sys->rd;
        const mtime_t date = sys->slot[i].date;
        const int n = sys->slot[i].n;
        const mtime_t period = sys->period;
        /* rattrapage : image en retard alors qu'une plus récente attend déjà -> on la saute */
        const int skip = (sys->count > 1 && date > 0 && mdate() > date + period);
        LeaveCriticalSection(&sys->cs);

        if (!skip && n > 0) {
            for (int k = 0; k < n; k++) {
                if (k > 0)   /* 6 ms d'avance : l'échange se cale ensuite sur la synchro verticale */
                    mwait(date + (mtime_t)((int64_t)period * k / n) - 6000);
                GlPresent(sys, sys->slot[i].data + (size_t)k * sys->frame_bytes);
            }
        }

        EnterCriticalSection(&sys->cs);
        sys->rd = (sys->rd + 1) % NSLOT;
        sys->count--;
        WakeAllConditionVariable(&sys->cv);
        LeaveCriticalSection(&sys->cs);
    }

    if (sys->tex) glDeleteTextures(1, &sys->tex);
    wglMakeCurrent(NULL, NULL);
    return 0;
}

/* ---- sortie vidéo ---- */

static void Prepare(vout_display_t *vd, picture_t *pic, subpicture_t *subpic)
{
    vout_display_sys_t *sys = vd->sys;
    (void)subpic;

    if (!sys->ts) {
        if (sys->ts_failed) return;
        char err[256] = "";
        sys->ts = ts_open(&sys->params, err, sizeof(err));
        if (!sys->ts) {
            msg_Err(vd, "temporal scaling failed to start: %s", err);
            sys->ts_failed = 1;
            sys->frames_valid = 0;
            return;
        }
        msg_Info(vd, "%dx%d -> %dx%d, mode %d, x%d frames on %s", sys->params.in_w, sys->params.in_h,
                 sys->params.out_w, sys->params.out_h, sys->params.mode, sys->ratio, ts_device_name());
    }

    /* durée d'une image : fréquence de la source si connue, sinon écart des dates */
    if (vd->source.i_frame_rate > 0 && vd->source.i_frame_rate_base > 0)
        sys->period = (mtime_t)CLOCK_FREQ * vd->source.i_frame_rate_base / vd->source.i_frame_rate;
    else if (sys->last_date > 0 && pic->date > sys->last_date) {
        mtime_t d = pic->date - sys->last_date;
        if (d > 5000 && d < 200000) sys->period = d;
    }
    sys->last_date = pic->date;

    const uint8_t *pl[3];
    int pt[3];
    pl[0] = pic->p[0].p_pixels + (size_t)sys->y_off * pic->p[0].i_pitch + sys->x_off;
    pt[0] = pic->p[0].i_pitch;
    int u = sys->swap_uv ? 2 : 1, v = sys->swap_uv ? 1 : 2;
    pl[1] = pic->p[u].p_pixels;  pt[1] = pic->p[u].i_pitch;
    pl[2] = pic->p[v].p_pixels;  pt[2] = pic->p[v].i_pitch;
    /* recadrage des plans de chrominance */
    if (sys->x_off || sys->y_off) {
        const int cx = (pic->p[1].i_visible_pitch && pic->p[0].i_visible_pitch)
                       ? sys->x_off * pic->p[1].i_visible_pitch / pic->p[0].i_visible_pitch : sys->x_off / 2;
        const int cy = (pic->p[1].i_visible_lines && pic->p[0].i_visible_lines)
                       ? sys->y_off * pic->p[1].i_visible_lines / pic->p[0].i_visible_lines : sys->y_off / 2;
        pl[1] += (size_t)cy * pt[1] + cx;
        pl[2] += (size_t)cy * pt[2] + cx;
    }

    /* place libre dans l'anneau (le thread de rendu libère les emplacements) */
    EnterCriticalSection(&sys->cs);
    while (sys->count == NSLOT && !sys->stop)
        SleepConditionVariableCS(&sys->cv, &sys->cs, 100);
    LeaveCriticalSection(&sys->cs);
    sys->frames_valid = (ts_push(sys->ts, pl, pt, sys->slot[sys->wr].data) == 0);
    if (!sys->frames_valid)
        msg_Err(vd, "GPU scaling failed on this picture");
}

static void Display(vout_display_t *vd, picture_t *pic, subpicture_t *subpic)
{
    vout_display_sys_t *sys = vd->sys;

    /* on confie les images au thread de rendu et on rend la main tout de suite */
    EnterCriticalSection(&sys->cs);
    if (!sys->stop && sys->count < NSLOT) {
        sys->slot[sys->wr].date = pic->date;
        sys->slot[sys->wr].n    = sys->frames_valid ? sys->ratio : 0;
        sys->wr = (sys->wr + 1) % NSLOT;
        sys->count++;
        WakeAllConditionVariable(&sys->cv);
    }
    LeaveCriticalSection(&sys->cs);

    picture_Release(pic);
    if (subpic)
        subpicture_Delete(subpic);
}

static int Control(vout_display_t *vd, int query, va_list args)
{
    vout_display_sys_t *sys = vd->sys;
    switch (query) {
    case VOUT_DISPLAY_CHANGE_FULLSCREEN: {            /* fenêtre propre : plein écran sans bordure */
        const vout_display_cfg_t *cfg = va_arg(args, const vout_display_cfg_t *);
        OwnFullscreen(sys, cfg->is_fullscreen);
        return VLC_SUCCESS;
    }
    case VOUT_DISPLAY_CHANGE_WINDOW_STATE:
        return VLC_SUCCESS;
    case VOUT_DISPLAY_CHANGE_DISPLAY_SIZE:
    case VOUT_DISPLAY_CHANGE_DISPLAY_FILLED:
    case VOUT_DISPLAY_CHANGE_ZOOM:
    case VOUT_DISPLAY_CHANGE_SOURCE_ASPECT:
    case VOUT_DISPLAY_CHANGE_SOURCE_CROP:
    case VOUT_DISPLAY_CHANGE_VIEWPOINT:
        return VLC_SUCCESS;     /* la fenêtre est relue à chaque image */
    default:
        return VLC_EGENERIC;
    }
}

static int Open(vlc_object_t *obj)
{
    vout_display_t *vd = (vout_display_t *)obj;

    vout_display_sys_t *sys = calloc(1, sizeof(*sys));
    if (!sys) return VLC_ENOMEM;
    vd->sys = sys;

    sys->vd = vd;

    /* ---- format d'entrée : YUV planaire 8 bits, sinon VLC convertit en I420 ---- */
    video_format_t fmt = vd->fmt;
    TsParams p;
    memset(&p, 0, sizeof(p));
    switch (fmt.i_chroma) {
        case VLC_CODEC_I420: p.chroma = TS_CHROMA_420; break;
        case VLC_CODEC_J420: p.chroma = TS_CHROMA_420; p.yuv_full = 1; break;
        case VLC_CODEC_YV12: p.chroma = TS_CHROMA_420; sys->swap_uv = 1; break;
        case VLC_CODEC_I422: p.chroma = TS_CHROMA_422; break;
        case VLC_CODEC_J422: p.chroma = TS_CHROMA_422; p.yuv_full = 1; break;
        case VLC_CODEC_I444: p.chroma = TS_CHROMA_444; break;
        case VLC_CODEC_J444: p.chroma = TS_CHROMA_444; p.yuv_full = 1; break;
        default:
            fmt.i_chroma = VLC_CODEC_I420;
            p.chroma = TS_CHROMA_420;
            break;
    }
    if (fmt.b_color_range_full) p.yuv_full = 1;
    p.in_w = fmt.i_visible_width;
    p.in_h = fmt.i_visible_height;
    sys->x_off = fmt.i_x_offset;
    sys->y_off = fmt.i_y_offset;
    if (p.chroma != TS_CHROMA_444) { p.in_w &= ~1; sys->x_off &= ~1; }
    if (p.chroma == TS_CHROMA_420) { p.in_h &= ~1; sys->y_off &= ~1; }

    p.yuv_709 = (fmt.space == COLOR_SPACE_BT709 || fmt.space == COLOR_SPACE_BT2020)
                ? 1 : (fmt.space == COLOR_SPACE_UNDEF && p.in_h >= 720);

    /* ---- options ---- */
    p.mode           = var_InheritInteger(vd, "tscale-mode");
    sys->ratio       = p.ratio = var_InheritInteger(vd, "tscale-ratio");
    p.tmp_mode       = var_InheritInteger(vd, "tscale-tmp-mode");
    p.curve_base     = var_InheritFloat(vd, "tscale-curve");
    p.kaiser_beta    = var_InheritFloat(vd, "tscale-beta");
    p.aniso_strength = var_InheritFloat(vd, "tscale-aniso");
    p.ocl_device     = var_InheritInteger(vd, "tscale-device");
    int out_h        = var_InheritInteger(vd, "tscale-height");
    if (out_h <= 0) out_h = GetSystemMetrics(SM_CYSCREEN);

    /* largeur : on garde le rapport d'aspect d'affichage (DAR) de la source */
    double sar = (vd->source.i_sar_num > 0 && vd->source.i_sar_den > 0)
                 ? (double)vd->source.i_sar_num / vd->source.i_sar_den : 1.0;
    double dar = (double)p.in_w / p.in_h * sar;
    p.out_h = (out_h + 1) & ~1;
    p.out_w = ((int)(p.out_h * dar + 0.5) + 1) & ~1;
    sys->out_w = p.out_w;
    sys->out_h = p.out_h;
    sys->period = CLOCK_FREQ / 25;

    /* vérifications bon marché tout de suite ; l'initialisation lourde d'OpenCL
     * (compilation des noyaux, tables) est faite à la première image : VLC ouvre et
     * referme parfois la sortie plusieurs fois au démarrage (essais de décodage matériel). */
    char err[256] = "";
    if (ts_check(&p, err, sizeof(err)) != 0) {
        msg_Err(vd, "temporal scaling unavailable: %s", err);
        goto error;
    }
    sys->params = p;
    sys->frame_bytes = (size_t)p.out_w * p.out_h * 3;
    for (int i = 0; i < NSLOT; i++) {
        sys->slot[i].data = malloc(sys->frame_bytes * (size_t)sys->ratio);
        if (!sys->slot[i].data) goto error;
    }
    InitializeCriticalSection(&sys->cs);
    InitializeConditionVariable(&sys->cv);
    sys->gl_event = CreateEvent(NULL, FALSE, FALSE, NULL);

    /* fenêtre : celle de VLC si elle en fournit une (interface Qt), sinon la nôtre */
    sys->window = vout_display_NewWindow(vd, VOUT_WINDOW_TYPE_HWND);
    if (sys->window && sys->window->type == VOUT_WINDOW_TYPE_HWND) {
        sys->hwnd = (HWND)sys->window->handle.hwnd;
    } else {
        if (sys->window) { vout_display_DeleteWindow(vd, sys->window); sys->window = NULL; }
        sys->own_window = 1;
        if (OwnWindowOpen(sys) != 0) { msg_Err(vd, "cannot create the video window"); goto error; }
    }
    if (GlOpen(vd, sys) != 0) goto error;
    sys->render = CreateThread(NULL, 0, RenderThread, sys, 0, NULL);
    if (!sys->render) goto error;
    WaitForSingleObject(sys->gl_event, 10000);
    if (sys->gl_failed) {
        msg_Err(vd, "OpenGL initialisation failed (%d): texture %dx%d", sys->gl_failed, sys->out_w, sys->out_h);
        goto error;
    }

    msg_Info(vd, "tscale ready: %dx%d -> %dx%d (DAR %.4f), mode %d, x%d frames (GPU starts on the first picture)",
             p.in_w, p.in_h, p.out_w, p.out_h, dar, p.mode, p.ratio);

    vd->fmt     = fmt;
    vd->pool    = Pool;
    vd->prepare = Prepare;
    vd->display = Display;
    vd->control = Control;
    return VLC_SUCCESS;

error:
    Close(obj);
    return VLC_EGENERIC;
}

static void Close(vlc_object_t *obj)
{
    vout_display_t *vd = (vout_display_t *)obj;
    vout_display_sys_t *sys = vd->sys;
    if (!sys) return;

    /* arrêt du thread de rendu, puis libération d'OpenGL et de la fenêtre */
    if (sys->render) {
        EnterCriticalSection(&sys->cs);
        sys->stop = 1;
        WakeAllConditionVariable(&sys->cv);
        LeaveCriticalSection(&sys->cs);
        WaitForSingleObject(sys->render, 5000);
        CloseHandle(sys->render);
    }
    if (sys->gl_event) CloseHandle(sys->gl_event);
    if (sys->rc) { wglMakeCurrent(NULL, NULL); wglDeleteContext(sys->rc); sys->rc = NULL; }
    if (sys->dc) { ReleaseDC(sys->hwnd, sys->dc); sys->dc = NULL; }
    if (sys->window) vout_display_DeleteWindow(vd, sys->window);
    OwnWindowClose(sys);
    if (sys->pool) picture_pool_Release(sys->pool);
    if (sys->ts)   ts_close(sys->ts);
    for (int i = 0; i < NSLOT; i++) free(sys->slot[i].data);
    free(sys);
    vd->sys = NULL;
}
