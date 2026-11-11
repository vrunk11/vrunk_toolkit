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

#define HEIGHT_TEXT     N_("Output height (-1 = automatic)")
#define HEIGHT_LONGTEXT N_("Vertical resolution the picture is scaled to on the GPU (the width follows the " \
    "display aspect ratio). -1 = automatic: follows the size of the picture shown in the window (and VLC " \
    "zoom, crop and aspect ratio); the scaler is rebuilt about half a second after a resize. " \
    "0 = height of the screen. Any other value = fixed. The scaler only enlarges: the output is " \
    "always larger than the video.")
#define MODE_TEXT       N_("Scaling algorithm")
#define MODE_LONGTEXT   N_("0 nearest, 1 bilinear on duplicated pixels, 2 bilinear with temporal " \
    "alternation, 3 linear weighted 2D, 4 exponential weighted 2D, 5 Kaiser 6x6 EWA, " \
    "6 Kaiser anisotropic edge-directed, 7 Kaiser anisotropic + selective softening.")
#define RATIO_TEXT      N_("Frame rate multiplier (0 = automatic)")
#define RATIO_LONGTEXT  N_("Number of output pictures shown per decoded picture (2 = 30 fps shown at 60 fps). " \
    "0 = automatic: the display rate (see 'Target display rate') is reached exactly by alternating " \
    "between two multipliers (pulldown: 24 fps on 60 Hz = 2,3,2,3...).")
#define FPS_TEXT        N_("Target display rate in Hz (0 = screen refresh rate)")
#define FPS_LONGTEXT    N_("Used when the frame rate multiplier is 0. 0 = refresh rate of the screen the window " \
    "is on (follows the window or full screen to another screen); or force e.g. 59.94.")
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
    add_integer("tscale-height", -1, HEIGHT_TEXT, HEIGHT_LONGTEXT, false)
    add_integer_with_range("tscale-mode", 5, 0, 7, MODE_TEXT, MODE_LONGTEXT, false)
    add_integer_with_range("tscale-ratio", 0, 0, 8, RATIO_TEXT, RATIO_LONGTEXT, false)
    add_float("tscale-fps", 0.0, FPS_TEXT, FPS_LONGTEXT, false)
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
    int             ratio;          /* multiplicateur en cours */
    int             opt_ratio;      /* option tscale-ratio : 0 = automatique */
    int             refresh_hz;     /* fréquence de l'écran de la fenêtre (mise en cache) */
    mtime_t         refresh_t;
    int             swap_uv;        /* YV12 : plans U et V inversés */
    int             x_off, y_off;   /* recadrage éventuel de la source */
    int             out_w, out_h;   /* taille de sortie actuelle du scaler */
    double          dar;
    int             opt_height;     /* option tscale-height : -1 auto, 0 écran, >0 fixe */
    int             zoom_num, zoom_den;
    int             tex_w, tex_h, max_tex;
    /* demande de changement de géométrie en attente (temporisée) */
    int             pending, pend_h, pend_x, pend_y, pend_w, pend_hh, pend_ratio;
    double          pend_dar;
    mtime_t         pend_since;
    /* anneau d'images prêtes (ratio images RGB24 chacune) : le thread vidéo de VLC écrit,
     * le thread de rendu présente */
#define NSLOT 3
    /* n images affichées parmi les R produites (pulldown), l'une après l'autre toutes les `step` */
    struct { uint8_t *data; mtime_t date; mtime_t step; int n, R; int w, h; } slot[NSLOT];
    int             n_cur;          /* images à afficher pour l'image en cours (pulldown) */
    double          acc;            /* accumulateur du pulldown */
    double          opt_fps;        /* option tscale-fps : 0 = fréquence de l'écran */
    double          target_fps;     /* cadence d'affichage visée */
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
    int             period_set, period_chg;
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

/* Texture : (ré)allouée par le thread de rendu quand la taille de sortie change. */
static void GlEnsureTexture(vout_display_sys_t *sys, int w, int h)
{
    if (sys->tex_w == w && sys->tex_h == h) return;
    glBindTexture(GL_TEXTURE_2D, sys->tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    sys->tex_w = w;
    sys->tex_h = h;
}

static void GlPresent(vout_display_sys_t *sys, const uint8_t *rgb, int w, int h)
{
    RECT rc;
    GetClientRect(sys->hwnd, &rc);
    int cw = rc.right - rc.left, ch = rc.bottom - rc.top;
    if (cw < 1 || ch < 1) return;

    GlEnsureTexture(sys, w, h);
    glBindTexture(GL_TEXTURE_2D, sys->tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, rgb);

    /* rapport d'aspect conservé (bandes noires), zoom de VLC appliqué */
    double sa = (double)w / h, wa = (double)cw / ch;
    double fw = cw, fh = ch;
    if (wa > sa) fw = ch * sa; else fh = cw / sa;
    if (sys->zoom_den > 0 && sys->zoom_num > 0) {
        fw = fw * sys->zoom_num / sys->zoom_den;
        fh = fh * sys->zoom_num / sys->zoom_den;
    }
    int vw = (int)(fw + 0.5), vh = (int)(fh + 0.5);
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
    sys->max_tex = maxtex > 0 ? maxtex : 8192;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &sys->tex);
    glBindTexture(GL_TEXTURE_2D, sys->tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
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
        const int n = sys->slot[i].n, R = sys->slot[i].R;
        const mtime_t step = sys->slot[i].step;
        const int w = sys->slot[i].w, h = sys->slot[i].h;
        const mtime_t period = sys->period;
        /* rattrapage : image en retard alors qu'une plus récente attend déjà -> on la saute */
        const int skip = (sys->count > 1 && date > 0 && mdate() > date + period);
        LeaveCriticalSection(&sys->cs);

        if (!skip && n > 0) {
            const size_t fb = (size_t)w * h * 3;
            for (int k = 0; k < n; k++) {
                if (k > 0)   /* 6 ms d'avance : l'échange se cale ensuite sur la synchro verticale */
                    mwait(date + step * k - 6000);
                /* n phases régulièrement réparties parmi les R calculées */
                const int idx = (int)((int64_t)k * R / n);
                GlPresent(sys, sys->slot[i].data + (size_t)idx * fb, w, h);
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

/* ---- géométrie : zone source (recadrage), aspect, taille de sortie ---- */

typedef struct {
    int    x, y, in_w, in_h;    /* zone source utilisée (recadrage de VLC appliqué) */
    int    out_w, out_h;
    double dar;
    int    ratio;               /* images calculées par image reçue (R) */
    double mult;                /* images affichées par image reçue, en moyenne (<= R) */
    double target;              /* cadence d'affichage visée (Hz) */
} Geom;

/* Fréquence de l'écran sur lequel se trouve la fenêtre (Hz), relue toutes les 500 ms : suit le
 * déplacement vers un autre écran et le plein écran. */
static int ScreenRefreshHz(vout_display_sys_t *sys)
{
    const mtime_t now = mdate();
    if (sys->refresh_hz > 0 && now - sys->refresh_t < 500000) return sys->refresh_hz;
    sys->refresh_t = now;
    int hz = 0;
    MONITORINFOEXA mi;
    memset(&mi, 0, sizeof(mi));
    mi.cbSize = sizeof(mi);
    if (sys->hwnd && GetMonitorInfoA(MonitorFromWindow(sys->hwnd, MONITOR_DEFAULTTONEAREST), (MONITORINFO *)&mi)) {
        DEVMODEA dm;
        memset(&dm, 0, sizeof(dm));
        dm.dmSize = sizeof(dm);
        if (EnumDisplaySettingsA(mi.szDevice, ENUM_CURRENT_SETTINGS, &dm) && dm.dmDisplayFrequency > 1)
            hz = (int)dm.dmDisplayFrequency;
    }
    if (hz <= 1) hz = 60;
    sys->refresh_hz = hz;
    return hz;
}

/* Multiplicateur moyen voulu m = cadence d'affichage visée / cadence reçue (1..8). Le moteur calcule
 * R = ceil(m) phases par image (le jeu de cartes « au-dessus ») ; quand m n'est pas entier on n'en
 * affiche que n = floor(m) ou floor(m)+1 selon un accumulateur (pulldown, ex. 24 -> 60 : 2,3,2,3),
 * en choisissant des phases régulièrement réparties (le jeu « en dessous »). Si m est à 0,03 d'un
 * entier (29,97 -> 60 Hz) il est ramené à cet entier : pas de pulldown. */
static void AutoMult(vout_display_sys_t *sys, int *R, double *m, double *target)
{
    const double fps = sys->period > 0 ? (double)CLOCK_FREQ / sys->period : 25.0;
    const double tgt = sys->opt_fps > 0 ? sys->opt_fps : (double)ScreenRefreshHz(sys);
    double mm = tgt / fps;
    if (mm < 1.0) mm = 1.0;
    if (mm > 8.0) mm = 8.0;
    const double r = floor(mm + 0.5);
    if (fabs(mm - r) < 0.03) { mm = r; *R = (int)r; }
    else                     { *R = (int)ceil(mm); }
    *m = mm;
    *target = tgt;
}

/* Calcule la géométrie voulue maintenant : recadrage et aspect de VLC (vd->source), taille de la
 * fenêtre, zoom. Hauteur de sortie : option tscale-height (-1 = suit l'image affichée dans la fenêtre,
 * 0 = hauteur de l'écran, >0 = fixe). Le scaler ne sait qu'agrandir : la sortie reste plus grande
 * que la source dans les deux dimensions. */
static void ComputeGeom(vout_display_t *vd, vout_display_sys_t *sys, Geom *g)
{
    const video_format_t *s = &vd->source;
    int x = s->i_x_offset, y = s->i_y_offset;
    int w = s->i_visible_width, h = s->i_visible_height;
    if (sys->params.chroma != TS_CHROMA_444) { x &= ~1; w &= ~1; }
    if (sys->params.chroma == TS_CHROMA_420) { y &= ~1; h &= ~1; }
    if (w < 2) w = 2;
    if (h < 2) h = 2;

    const double sar = (s->i_sar_num > 0 && s->i_sar_den > 0) ? (double)s->i_sar_num / s->i_sar_den : 1.0;
    const double dar = (double)w / h * sar;

    int target;
    if (sys->opt_height > 0) {
        target = sys->opt_height;
    } else if (sys->opt_height == 0) {
        target = GetSystemMetrics(SM_CYSCREEN);
    } else {
        RECT rc;
        GetClientRect(sys->hwnd, &rc);
        const double ww = rc.right - rc.left, wh = rc.bottom - rc.top;
        if (ww < 16 || wh < 16) {                       /* fenêtre réduite : on garde la taille actuelle */
            target = sys->out_h > 0 ? sys->out_h : 1080;
        } else {
            double disp_h = (ww / wh > dar) ? wh : ww / dar;     /* hauteur de l'image ajustée à la fenêtre */
            if (sys->zoom_den > 0 && sys->zoom_num > 0) disp_h = disp_h * sys->zoom_num / sys->zoom_den;
            target = (int)(disp_h + 0.5);
        }
    }

    int min_h = h + 2;
    const int min_h_for_w = (int)((w + 2) / dar) + 2;
    if (min_h_for_w > min_h) min_h = min_h_for_w;
    int out_h = target > min_h ? target : min_h;
    const int cap = sys->max_tex > 0 ? sys->max_tex : 16384;
    if (out_h > cap) out_h = cap;
    out_h = (out_h + 1) & ~1;
    int out_w = ((int)(out_h * dar + 0.5) + 1) & ~1;
    if (out_w > cap) { out_w = cap & ~1; out_h = ((int)(out_w / dar + 0.5) + 1) & ~1; }

    g->x = x; g->y = y; g->in_w = w; g->in_h = h;
    g->out_w = out_w; g->out_h = out_h;
    g->dar = dar;
    if (sys->opt_ratio > 0) {                                   /* ratio fixe : pas de pulldown */
        g->ratio = sys->opt_ratio; g->mult = sys->opt_ratio; g->target = 0;
    } else {
        AutoMult(sys, &g->ratio, &g->mult, &g->target);
    }
}

static int GeomChangedEnough(const vout_display_sys_t *sys, const Geom *g)
{
    if (g->x != sys->x_off || g->y != sys->y_off || g->in_w != sys->params.in_w || g->in_h != sys->params.in_h)
        return 1;                                        /* recadrage */
    if (fabs(g->dar - sys->dar) > 0.002)
        return 1;                                        /* rapport d'aspect */
    if (g->ratio != sys->ratio)
        return 1;                                        /* cadence / écran */
    const int dh = abs(g->out_h - sys->out_h);
    return dh * 100 > sys->out_h * 5;                    /* taille : au moins 5 % d'écart */
}

/* (Re)construit le scaler pour la géométrie g : vide la file de présentation, réalloue les
 * tampons et rouvre le moteur (quelques centaines de ms). Renvoie 0 si OK. */
static int Rebuild(vout_display_t *vd, vout_display_sys_t *sys, const Geom *g)
{
    /* 1) plus aucune image en cours de présentation */
    EnterCriticalSection(&sys->cs);
    for (int t = 0; t < 40 && sys->count > 0 && !sys->stop; t++)
        SleepConditionVariableCS(&sys->cv, &sys->cs, 50);
    LeaveCriticalSection(&sys->cs);

    /* 2) ancien moteur et tampons */
    if (sys->ts) { ts_close(sys->ts); sys->ts = NULL; }
    for (int i = 0; i < NSLOT; i++) { free(sys->slot[i].data); sys->slot[i].data = NULL; }

    /* 3) nouveaux paramètres */
    sys->params.in_w = g->in_w;
    sys->params.in_h = g->in_h;
    sys->params.out_w = g->out_w;
    sys->params.out_h = g->out_h;
    sys->x_off = g->x;
    sys->y_off = g->y;
    sys->out_w = g->out_w;
    sys->out_h = g->out_h;
    sys->dar   = g->dar;
    sys->ratio = sys->params.ratio = g->ratio;
    sys->frame_bytes = (size_t)g->out_w * g->out_h * 3;
    for (int i = 0; i < NSLOT; i++) {
        sys->slot[i].data = malloc(sys->frame_bytes * (size_t)sys->ratio);
        if (!sys->slot[i].data) return -1;
    }

    char err[256] = "";
    sys->ts = ts_open(&sys->params, err, sizeof(err));
    if (!sys->ts) {
        msg_Err(vd, "temporal scaling failed to start: %s", err);
        return -1;
    }
    msg_Info(vd, "%dx%d (crop %d,%d) -> %dx%d (DAR %.4f), mode %d, x%d frames computed (%.2f fps in, target %.2f Hz, x%.3f shown on average%s) on %s",
             g->in_w, g->in_h, g->x, g->y, g->out_w, g->out_h, g->dar,
             sys->params.mode, sys->ratio, sys->period > 0 ? (double)CLOCK_FREQ / sys->period : 0.0,
             g->target, g->mult, sys->opt_ratio > 0 ? ", fixed ratio" : ", pulldown", ts_device_name());
    return 0;
}

/* ---- sortie vidéo ---- */

static void Prepare(vout_display_t *vd, picture_t *pic, subpicture_t *subpic)
{
    vout_display_sys_t *sys = vd->sys;
    (void)subpic;
    sys->frames_valid = 0;
    if (sys->ts_failed) return;

    /* durée d'une image : fréquence de la source au départ, puis écart réel des dates (lissé) ;
     * nécessaire avec un désentrelacement VLC à double cadence (yadif2x, bob...) : la cadence
     * annoncée par la source (25 i/s) n'est alors plus celle des images reçues (50 i/s) */
    if (!sys->period_set && vd->source.i_frame_rate > 0 && vd->source.i_frame_rate_base > 0)
        sys->period = (mtime_t)CLOCK_FREQ * vd->source.i_frame_rate_base / vd->source.i_frame_rate;
    if (sys->last_date > 0 && pic->date > sys->last_date) {
        const mtime_t d = pic->date - sys->last_date;
        if (d > 5000 && d < 200000) {
            if (!sys->period_set) {
                sys->period = d; sys->period_set = 1;
            } else if (d * 4 >= sys->period * 3 && d * 3 <= sys->period * 4) {
                sys->period = (sys->period * 7 + d) / 8;        /* même cadence : lissage */
                sys->period_chg = 0;
            } else if (++sys->period_chg >= 3) {                /* cadence changée : 3 écarts de suite */
                sys->period = d; sys->period_chg = 0;
            }
        }
    }
    sys->last_date = pic->date;

    /* géométrie voulue (recadrage / aspect / taille de la fenêtre / multiplicateur de cadence selon
     * l'écran) ; le moteur n'est reconstruit qu'une fois la demande stable depuis 400 ms (pas
     * pendant un redimensionnement à la souris) */
    Geom g;
    ComputeGeom(vd, sys, &g);
    if (!sys->ts) {
        if (Rebuild(vd, sys, &g) != 0) { sys->ts_failed = 1; return; }
    } else if (GeomChangedEnough(sys, &g)) {
        const mtime_t now = mdate();
        if (!sys->pending || g.out_h != sys->pend_h || g.x != sys->pend_x || g.y != sys->pend_y ||
            g.in_w != sys->pend_w || g.in_h != sys->pend_hh || fabs(g.dar - sys->pend_dar) > 0.002 ||
            g.ratio != sys->pend_ratio) {
            sys->pending = 1;
            sys->pend_ratio = g.ratio;
            sys->pend_since = now;
            sys->pend_h = g.out_h; sys->pend_x = g.x; sys->pend_y = g.y;
            sys->pend_w = g.in_w;  sys->pend_hh = g.in_h; sys->pend_dar = g.dar;
        } else if (now - sys->pend_since >= 400000) {
            sys->pending = 0;
            if (Rebuild(vd, sys, &g) != 0) { sys->ts_failed = 1; return; }
        }
    } else {
        sys->pending = 0;
    }

    /* pulldown : nombre d'images à afficher pour cette image (moyenne g.mult, au plus R) */
    {
        double m = g.mult > sys->ratio ? sys->ratio : g.mult;
        int base = (int)floor(m + 1e-9);
        if (base < 1) base = 1;
        int n = base;
        const double fr = m - base;
        if (fr > 1e-4) {
            sys->acc += fr;
            if (sys->acc >= 1.0 - 1e-9) { n = base + 1; sys->acc -= 1.0; }
        } else {
            sys->acc = 0;
        }
        sys->n_cur = n > sys->ratio ? sys->ratio : n;
        sys->target_fps = g.target;
    }

    const uint8_t *pl[3];
    int pt[3];
    pl[0] = pic->p[0].p_pixels + (size_t)sys->y_off * pic->p[0].i_pitch + sys->x_off;
    pt[0] = pic->p[0].i_pitch;
    int u = sys->swap_uv ? 2 : 1, v = sys->swap_uv ? 1 : 2;
    pl[1] = pic->p[u].p_pixels;  pt[1] = pic->p[u].i_pitch;
    pl[2] = pic->p[v].p_pixels;  pt[2] = pic->p[v].i_pitch;
    /* recadrage des plans de chrominance (4:2:0 : moitié, 4:2:2 : moitié en x, 4:4:4 : rien) */
    if (sys->x_off || sys->y_off) {
        const int cx = (sys->params.chroma == TS_CHROMA_444) ? sys->x_off : sys->x_off / 2;
        const int cy = (sys->params.chroma == TS_CHROMA_420) ? sys->y_off / 2 : sys->y_off;
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
    if (!sys->stop && sys->count < NSLOT && sys->ts) {
        sys->slot[sys->wr].date = pic->date;
        sys->slot[sys->wr].n    = sys->frames_valid ? sys->n_cur : 0;
        sys->slot[sys->wr].R    = sys->ratio;
        /* écart entre deux images affichées : une période d'écran en mode automatique (la synchro
         * verticale fait le reste), sinon la durée de l'image répartie sur les n images */
        sys->slot[sys->wr].step = sys->target_fps > 0 ? (mtime_t)(CLOCK_FREQ / sys->target_fps)
                                : sys->period / (sys->n_cur > 0 ? sys->n_cur : 1);
        sys->slot[sys->wr].w    = sys->out_w;
        sys->slot[sys->wr].h    = sys->out_h;
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
    case VOUT_DISPLAY_CHANGE_ZOOM: {
        const vout_display_cfg_t *cfg = va_arg(args, const vout_display_cfg_t *);
        sys->zoom_num = cfg->zoom.num;
        sys->zoom_den = cfg->zoom.den;
        return VLC_SUCCESS;
    }
    case VOUT_DISPLAY_CHANGE_WINDOW_STATE:
    case VOUT_DISPLAY_CHANGE_DISPLAY_SIZE:      /* la fenêtre et la géométrie sont relues à chaque image */
    case VOUT_DISPLAY_CHANGE_DISPLAY_FILLED:
    case VOUT_DISPLAY_CHANGE_SOURCE_ASPECT:     /* vd->source est déjà mis à jour par VLC */
    case VOUT_DISPLAY_CHANGE_SOURCE_CROP:
    case VOUT_DISPLAY_CHANGE_VIEWPOINT:
        return VLC_SUCCESS;
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
    TsParams *p = &sys->params;
    switch (fmt.i_chroma) {
        case VLC_CODEC_I420: p->chroma = TS_CHROMA_420; break;
        case VLC_CODEC_J420: p->chroma = TS_CHROMA_420; p->yuv_full = 1; break;
        case VLC_CODEC_YV12: p->chroma = TS_CHROMA_420; sys->swap_uv = 1; break;
        case VLC_CODEC_I422: p->chroma = TS_CHROMA_422; break;
        case VLC_CODEC_J422: p->chroma = TS_CHROMA_422; p->yuv_full = 1; break;
        case VLC_CODEC_I444: p->chroma = TS_CHROMA_444; break;
        case VLC_CODEC_J444: p->chroma = TS_CHROMA_444; p->yuv_full = 1; break;
        case VLC_CODEC_D3D9_OPAQUE: case VLC_CODEC_D3D9_OPAQUE_10B:
        case VLC_CODEC_D3D11_OPAQUE: case VLC_CODEC_D3D11_OPAQUE_10B:
            /* surfaces du décodage matériel : le désentrelacement de VLC ne peut pas s'y appliquer
             * (VLC retire alors tout le filtre -> image en peigne). On refuse : le décodeur
             * retombe sur le décodage logiciel, où le désentrelacement et les sous-titres marchent. */
            free(sys);
            vd->sys = NULL;
            return VLC_EGENERIC;
        default:
            fmt.i_chroma = VLC_CODEC_I420;
            p->chroma = TS_CHROMA_420;
            break;
    }
    if (fmt.b_color_range_full) p->yuv_full = 1;
    p->yuv_709 = (fmt.space == COLOR_SPACE_BT709 || fmt.space == COLOR_SPACE_BT2020)
                 ? 1 : (fmt.space == COLOR_SPACE_UNDEF && (int)fmt.i_visible_height >= 720);

    /* ---- options ---- */
    p->mode           = var_InheritInteger(vd, "tscale-mode");
    sys->opt_fps      = var_InheritFloat(vd, "tscale-fps");
    sys->opt_ratio    = var_InheritInteger(vd, "tscale-ratio");     /* 0 = automatique (écran) */
    sys->ratio        = p->ratio = sys->opt_ratio > 0 ? sys->opt_ratio : 1;
    p->tmp_mode       = var_InheritInteger(vd, "tscale-tmp-mode");
    p->curve_base     = var_InheritFloat(vd, "tscale-curve");
    p->kaiser_beta    = var_InheritFloat(vd, "tscale-beta");
    p->aniso_strength = var_InheritFloat(vd, "tscale-aniso");
    p->ocl_device     = var_InheritInteger(vd, "tscale-device");
    sys->opt_height   = var_InheritInteger(vd, "tscale-height");
    sys->zoom_num     = vd->cfg->zoom.num;
    sys->zoom_den     = vd->cfg->zoom.den;
    sys->period       = CLOCK_FREQ / 25;
    sys->max_tex      = 16384;

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
        sys->out_w = fmt.i_visible_width * 2;       /* taille initiale de la fenêtre propre */
        sys->out_h = fmt.i_visible_height * 2;
        if (OwnWindowOpen(sys) != 0) { msg_Err(vd, "cannot create the video window"); goto error; }
        sys->out_w = sys->out_h = 0;
    }

    /* vérifications bon marché tout de suite ; l'initialisation lourde d'OpenCL (noyaux, tables) est
     * faite à la première image : VLC ouvre et referme parfois la sortie plusieurs fois au démarrage
     * (essais de décodage matériel) */
    {
        Geom g;
        ComputeGeom(vd, sys, &g);
        TsParams chk = *p;
        chk.in_w = g.in_w; chk.in_h = g.in_h; chk.out_w = g.out_w; chk.out_h = g.out_h;
        chk.ratio = g.ratio;
        char err[256] = "";
        if (ts_check(&chk, err, sizeof(err)) != 0) {
            msg_Err(vd, "temporal scaling unavailable: %s", err);
            goto error;
        }
        msg_Info(vd, "tscale ready: %dx%d -> %dx%d (DAR %.4f), mode %d, x%d frames, height %s",
                 g.in_w, g.in_h, g.out_w, g.out_h, g.dar, p->mode, g.ratio,
                 sys->opt_height < 0 ? "automatic" : sys->opt_height == 0 ? "= screen" : "fixed");
    }

    if (GlOpen(vd, sys) != 0) goto error;
    sys->render = CreateThread(NULL, 0, RenderThread, sys, 0, NULL);
    if (!sys->render) goto error;
    WaitForSingleObject(sys->gl_event, 10000);
    if (sys->gl_failed) {
        msg_Err(vd, "OpenGL initialisation failed (%d)", sys->gl_failed);
        goto error;
    }

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
