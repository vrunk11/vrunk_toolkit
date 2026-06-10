#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include <unistd.h>

#ifndef _WIN32
	#include <unistd.h>
	#define sleep_ms(ms)	usleep(ms*1000)
	#else
	#include <windows.h>
	#include <io.h>
	#include <fcntl.h>
	#define sleep_ms(ms)	Sleep(ms)
#endif

#ifdef _WIN32
    #define ALIGNED_MALLOC(size, align) _aligned_malloc(size, align)
    #define ALIGNED_FREE(ptr)           _aligned_free(ptr)
#else
    #define ALIGNED_MALLOC(size, align) aligned_alloc(align, size)
    #define ALIGNED_FREE(ptr)           free(ptr)
#endif

#define _FILE_OFFSET_BITS 64

#ifdef _WIN64
#define FSEEK fseeko64
#else
#define FSEEK fseeko
#endif

#define NB_THREADS 16

// Pondération en virgule fixe : poids stockés sur 0..WEIGHT_MAX
#define WEIGHT_SHIFT 8
#define WEIGHT_MAX   (1 << WEIGHT_SHIFT)
#define CLAMP255(v)  ((v) < 0 ? 0 : (v) > 255 ? 255 : (v))  // ← ajouter ici

typedef enum {
    PIX_RGB    = 0,
    PIX_444P   = 1,
    PIX_422P   = 2,
    PIX_420P   = 3,
    PIX_411P   = 4,
    PIX_410P   = 5,
    PIX_RGB48  = 6,
    PIX_444P16 = 7,
    PIX_422P16 = 8,
    PIX_420P16 = 9,
} PixFmt;

typedef struct BresenhamMap BresenhamMap;
struct BresenhamMap {
    int *w_src;     // index source
    int *w_flag;    // -1, 0, +1 : direction d'interpolation
    int *w_weight;  // 0..WEIGHT_MAX : poids du voisin (mode 3)
    int *w_heavy;   // 0/1 : pixel dans la moitié haute du run (mode 3)
	int *w_pos;   // position continue src*WEIGHT_MAX + frac, pour 2D isotrope
    int *h_src;
    int *h_flag;
    int *h_weight;
    int *h_heavy;
	int *h_pos;
    int w_length;
    int h_length;
};

typedef struct {
    int16_t gx;   // gradient horizontal
    int16_t gy;   // gradient vertical
    int16_t mag;  // confiance 0..WEIGHT_MAX (magnitude normalisée)
    int16_t conf;
} GradPixel;

typedef struct ThreadArgs ThreadArgs;
struct ThreadArgs
{
	unsigned char *in_buf;
	unsigned char *next_buf;   // frame d'entrée suivante pour blend temporel
	unsigned char *out_buf;
	unsigned char *tmp_buf;
	unsigned char *blend_buf;  // résultat du lerp temporel (taille in_w*in_h*3)
	BresenhamMap *map;
	BresenhamMap *map_uv;
	const int *tmp_weights;    // poids temporel par phase (taille = ratio)
	GradPixel *grad;
	int16_t   *dist_map;
	int        blur_artifact;
	int aniso_strength_fp;  // 0..WEIGHT_MAX
	int mode_1d;
	PixFmt in_pix_fmt;
	PixFmt out_pix_fmt;
	int phase_v;
	int bps_in;
	int bps_out;
	int in_w, in_h, out_w, out_h, ratio, mode;
};

typedef struct {
    pthread_t      tid;
    pthread_mutex_t mutex;
    pthread_cond_t  cond_work;
    pthread_cond_t  cond_done;
    ThreadArgs      args;
    int             ready;  // 1=travail, -1=exit, 0=repos
    int             done;
} Worker;

static const int sobel5_x[5][5] = {
    { -1, -2, 0,  2,  1 },
    { -4, -8, 0,  8,  4 },
    { -6,-12, 0, 12,  6 },
    { -4, -8, 0,  8,  4 },
    { -1, -2, 0,  2,  1 }
};

// ============================================================
// RECONSTRUCTION DE TRANSITION
//
// Sur une bande de pixels source flaggés artefact (≥2 adjacents dans
// un axe), on reconstruit une rampe smoothstep entre les deux pixels
// fiables (non flaggés) qui bordent la bande. Les pixels libres de
// l'upscale déroulent la transition en continu. Cardinalité préservée :
// seuls les pixels flaggés (déjà non-cardinaux) sont réécrits.
// ============================================================

// LUT smoothstep : index = t en virgule fixe [0..WEIGHT_MAX],
// valeur = 3t² − 2t³ en virgule fixe [0..WEIGHT_MAX].
static int smoothstep_lut[WEIGHT_MAX + 1];

void compute_smoothstep_lut(void)
{
    for(int i = 0; i <= WEIGHT_MAX; i++) {
        double t = (double)i / (double)WEIGHT_MAX;
        double s = t * t * (3.0 - 2.0 * t);
        smoothstep_lut[i] = (int)(s * (double)WEIGHT_MAX + 0.5);
    }
}

static inline int gsample(const unsigned char *src, int x, int y,
                          int w, int h, int bps)
{
    x = x < 0 ? 0 : x >= w ? w-1 : x;
    y = y < 0 ? 0 : y >= h ? h-1 : y;
    return (bps == 1) ? (int)src[y*w + x]
                      : (int)((const uint16_t*)src)[y*w + x] >> 8;
}

//    LUT GLOBALE — 4 poids kaiser pour chaque position fractionnaire
//    Taille : (WEIGHT_MAX+1) * 4 * sizeof(int) = 4KB, négligeable.
//    Précalculée une fois au démarrage dans main().
 
static int kaiser4_lut[WEIGHT_MAX + 1][4];
 
// sinc normalisé : sinc(0)=1, sinc(x)=sin(πx)/(πx)
static double sinc_norm(double x)
{
    if(x == 0.0) return 1.0;
    double px = M_PI * x;
    return sin(px) / px;
}
 
// Fenêtre de Kaiser : I0(beta * sqrt(1 - (x/a)²)) / I0(beta)
// Approximation I0 via série de Taylor (suffisante pour beta<=10)
static double bessel_i0(double x)
{
    double sum = 1.0, term = 1.0;
    double x2 = x * x * 0.25;
    for(int k = 1; k <= 20; k++) {
        term *= x2 / ((double)(k * k));
        sum  += term;
        if(term < 1e-12 * sum) break;
    }
    return sum;
}
 
// kaiser_sinc(x, lobes, beta) :
//   x     : position en pixels
//   lobes : nombre de lobes (2 = lanczos-2, ici on utilise 2)
//   beta  : paramètre kaiser (6..12 typique, 6=doux, 12=agressif)
static double kaiser_sinc(double x, double lobes, double beta)
{
    double ax = fabs(x);
    if(ax >= lobes) return 0.0;
    double w = bessel_i0(beta * sqrt(1.0 - (x / lobes) * (x / lobes)))
             / bessel_i0(beta);
    return sinc_norm(x) * w;
}
 
// Précalcule la LUT : pour chaque t dans [0..WEIGHT_MAX],
// calcule les 4 poids normalisés (en virgule fixe 0..WEIGHT_MAX)
// correspondant à l'interpolation à la position fractionnaire
// t / WEIGHT_MAX entre src et src+1.
//
// Convention :
//   k[0] = poids pour src-1   (x = t+1 par rapport à src-1)
//   k[1] = poids pour src     (x = t)
//   k[2] = poids pour src+1   (x = t-1)
//   k[3] = poids pour src+2   (x = t-2)
//
// À t=0 : k[]={0,1,0,0} → pixel source pur (dégénère en nearest)
// À t=1 : k[]={0,0,1,0} → pixel src+1 pur
//
// beta=6.0 est un bon défaut pour du contenu vidéo compressé :
// assez doux pour ne pas amplifier les artefacts DCT.
void compute_kaiser4_lut(double beta)
{
    const double lobes = 2.0;
    for(int ti = 0; ti <= WEIGHT_MAX; ti++)
    {
        double t = (double)ti / (double)WEIGHT_MAX;
 
        // distances de chaque tap source par rapport à la position cible
        double xs[4] = { t + 1.0, t, t - 1.0, t - 2.0 };
 
        double w[4], sum = 0.0;
        for(int i = 0; i < 4; i++) {
            w[i] = kaiser_sinc(xs[i], lobes, beta);
            sum += w[i];
        }
 
        // normalisation + virgule fixe
        int isum = 0;
        for(int i = 0; i < 4; i++) {
            kaiser4_lut[ti][i] = (int)(w[i] / sum * (double)WEIGHT_MAX + 0.5);
            isum += kaiser4_lut[ti][i];
        }
 
        // correction d'arrondi : ajuster k[1] pour que la somme == WEIGHT_MAX
        kaiser4_lut[ti][1] += WEIGHT_MAX - isum;
    }
}

// LUT radiale pour kaiser 2D isotrope.
// r dans [0 .. sqrt(2)*2] ≈ 2.83, indexé par r * WEIGHT_MAX
#define KAISER_RAD_SIZE (5 * WEIGHT_MAX + 1)   // couvre sqrt(2)*3 ≈ 4.24
static int kaiser_rad_lut[KAISER_RAD_SIZE];

#define CURVE_RAD_SIZE  (5 * WEIGHT_MAX + 1)
static int curve_rad_lut[CURVE_RAD_SIZE];

#define KERNEL_MAX_TAPS 24
#define NQUANT_SHIFT 1   // log2(NQUANT / WEIGHT_MAX)
                          // 0 → NQUANT=256, 1 → NQUANT=512, 2 → NQUANT=1024

#define NQUANT (WEIGHT_MAX << NQUANT_SHIFT)

#define FQ(frac) ((frac) << NQUANT_SHIFT)

#define NQUANT_ANGLE 64       // orientations quantifiées sur [0, pi)
#define ANISO_K      2.5      // raideur de l'atténuation perpendiculaire
// (plus grand = bord plus net mais plus d'aliasing)

#define ANISO_RATIO   2.5f
#define ANISO_LOW     20
#define ANISO_HIGH    120

// max_mag : magnitude au-delà de laquelle conf sature à WEIGHT_MAX.
// Réglé empiriquement ; plus bas = anisotropie qui se déclenche plus tôt.
#define GRAD_MAX_MAG 2048

#define SHARP_HI       200   // borne haute rampe / seuil contour franc (mode 7)
#define NARROW_DROP      2    // narrow si  mag_voisin * NARROW_DROP < mag_centre
#define CONT_KEEP        2    // cont   si  mag_voisin * CONT_KEEP  >= mag_centre
 
typedef struct {
    int16_t w [KERNEL_MAX_TAPS];   // poids normalisés (somme = WEIGHT_MAX)
    int16_t ox[KERNEL_MAX_TAPS];   // offset pixel linéaire : dy * in_w + dx
    int8_t  dx[KERNEL_MAX_TAPS];   // delta x (-2..+3) pour slow path clamp
    int8_t  dy[KERNEL_MAX_TAPS];   // delta y (-2..+3) pour slow path clamp
    int8_t  n;                      // nombre de taps non-nuls
    int8_t  pad[7];
} __attribute__((aligned(64))) Kernel2D;
 
static Kernel2D kaiser_k2d   [NQUANT * NQUANT];
static Kernel2D kaiser_k2d_uv[NQUANT * NQUANT];
static Kernel2D curve_k2d    [NQUANT * NQUANT];
static Kernel2D curve_k2d_uv [NQUANT * NQUANT];

#define MID_RAD_SIZE  CURVE_RAD_SIZE          // même taille d'index que les autres
#define MID_LOBE_GAIN 0.5                       // atténuation des lobes négatifs
 
static int mid_rad_lut[MID_RAD_SIZE];
 
static Kernel2D mid_k2d   [NQUANT * NQUANT];
static Kernel2D mid_k2d_uv[NQUANT * NQUANT];
 
void compute_mid_rad_lut(double beta)
{
    const double radius = 2.0;
    const double inv_i0 = 1.0 / bessel_i0(beta);
 
    for(int ri = 0; ri < MID_RAD_SIZE; ri++)
    {
        double r = (double)ri / (double)WEIGHT_MAX;
 
        if(r >= radius) { mid_rad_lut[ri] = 0; continue; }
 
        // sinc radial cardinal
        double s;
        if(r == 0.0) s = 1.0;
        else { double px = M_PI * r; s = sin(px) / px; }
 
        // fenêtre kaiser sur le support radius 2
        double k = bessel_i0(beta * sqrt(1.0 - (r/radius)*(r/radius))) * inv_i0;
 
        double val = s * k;
 
        // atténue uniquement les lobes négatifs
        if(val < 0.0) val *= MID_LOBE_GAIN;
 
        int vi = (val >= 0.0)
               ?  (int)( val * (double)WEIGHT_MAX + 0.5)
               : -(int)(-val * (double)WEIGHT_MAX + 0.5);
        mid_rad_lut[ri] = vi;
    }
}

// lecture clampée de mag dans la map de gradient
static inline int mag_at(const GradPixel *grad, int x, int y, int in_w, int in_h)
{
    x = x < 0 ? 0 : x >= in_w ? in_w - 1 : x;
    y = y < 0 ? 0 : y >= in_h ? in_h - 1 : y;
    return (int)grad[y * in_w + x].mag;
}
 
// détecteur séparable d'artefact de compression / bruit axial.
// Lit mag au centre + 4 voisins cardinaux. Pas de trigonométrie.
static inline int is_compression_artifact(const GradPixel *grad,
                                          int sx, int sy, int mag_c,
                                          int in_w, int in_h)
{
    int mx_l = mag_at(grad, sx - 1, sy,     in_w, in_h);
    int mx_r = mag_at(grad, sx + 1, sy,     in_w, in_h);
    int my_u = mag_at(grad, sx,     sy - 1, in_w, in_h);
    int my_d = mag_at(grad, sx,     sy + 1, in_w, in_h);
 
    // étroit sur un axe : le centre domine et retombe fort des deux côtés
    int narrow_x = (mx_l * NARROW_DROP < mag_c) && (mx_r * NARROW_DROP < mag_c);
    int narrow_y = (my_u * NARROW_DROP < mag_c) && (my_d * NARROW_DROP < mag_c);
 
    // continu sur un axe : le gradient reste soutenu des deux côtés
    int cont_x = (mx_l * CONT_KEEP >= mag_c) && (mx_r * CONT_KEEP >= mag_c);
    int cont_y = (my_u * CONT_KEEP >= mag_c) && (my_d * CONT_KEEP >= mag_c);
 
    // 1) croix = étroit dans les deux axes (intersection de grille DCT)
    int is_cross = narrow_x && narrow_y;
 
    // 2) arête axiale = étroit sur un axe, continu sur l'autre
    int is_axial_edge = (narrow_x && cont_y) || (narrow_y && cont_x);
 
    return is_cross || is_axial_edge;
}


static void build_k2d(Kernel2D *out, const int *rad_lut, int rad_size, int stride)
{
    for(int fy_q = 0; fy_q < NQUANT; fy_q++)
    for(int fx_q = 0; fx_q < NQUANT; fx_q++)
    {
        int fx = fx_q * WEIGHT_MAX / NQUANT;
        int fy = fy_q * WEIGHT_MAX / NQUANT;

        int tw[36], tdx[36], tdy[36];
        int wsum = 0, n = 0;

        for(int dj = -2; dj <= 3; dj++)
        for(int di = -2; di <= 3; di++)
        {
            float dxf = (float)fx / WEIGHT_MAX - di;
            float dyf = (float)fy / WEIGHT_MAX - dj;
            int ri = (int)(sqrtf(dxf*dxf + dyf*dyf) * WEIGHT_MAX + 0.5f);
            int w  = (ri < rad_size) ? rad_lut[ri] : 0;
            if(w != 0) {
                tw [n] = w;
                tdx[n] = di;
                tdy[n] = dj;
                wsum += w;
                n++;
            }
        }

        int idx = fy_q * NQUANT + fx_q;
        Kernel2D *k = &out[idx];
        k->n = (int8_t)(n > KERNEL_MAX_TAPS ? KERNEL_MAX_TAPS : n);
        memset(k->w,  0, sizeof(k->w));
        memset(k->ox, 0, sizeof(k->ox));
        memset(k->dx, 0, sizeof(k->dx));
        memset(k->dy, 0, sizeof(k->dy));

        for(int i = 0; i < k->n; i++) {
            long long num = (long long)tw[i] * WEIGHT_MAX;
            k->w [i] = (int16_t)((num + (num >= 0 ? wsum/2 : -wsum/2)) / wsum);
            k->ox[i] = (int16_t)(tdy[i] * stride + tdx[i]);
            k->dx[i] = (int8_t)tdx[i];
            k->dy[i] = (int8_t)tdy[i];
        }

        int wcheck = 0;
        for(int i = 0; i < k->n; i++) wcheck += k->w[i];
        if(k->n > 0) k->w[0] += (int16_t)(WEIGHT_MAX - wcheck);
    }
}

static void compute_gradient_plane(const unsigned char *src,
                                   GradPixel *grad,
                                   int w, int h, int bps)
{
    // bordure : avec clamp (4 bandes de 2 pixels)
    // bulk : accès direct, zéro branche
    for(int y = 0; y < h; y++)
    {
        const int border = (y < 2 || y >= h-2);
        for(int x = 0; x < w; x++)
        {
            int gx = 0, gy = 0;
            if(__builtin_expect(border || x < 2 || x >= w-2, 0))
            {
                // slow path : clamp
                for(int dy = -2; dy <= 2; dy++)
                for(int dx = -2; dx <= 2; dx++) {
                    int v = gsample(src, x+dx, y+dy, w, h, bps);
                    gx += sobel5_x[dy+2][dx+2] * v;
                    gy += sobel5_x[dx+2][dy+2] * v;
                }
            }
            else
            {
                // fast path : accès direct, vectorisable
                if(bps == 1) {
                    const unsigned char *row = src + (y-2)*w + (x-2);
                    for(int dy = 0; dy < 5; dy++, row += w) {
                        gx += sobel5_x[dy][0]*row[0] + sobel5_x[dy][1]*row[1]
                            + sobel5_x[dy][3]*row[3] + sobel5_x[dy][4]*row[4];
                        gy += sobel5_x[0][dy]*row[0] + sobel5_x[1][dy]*row[1]
                            + sobel5_x[3][dy]*row[3] + sobel5_x[4][dy]*row[4];
                    }
                } else {
                    const uint16_t *row = (const uint16_t*)src + (y-2)*w + (x-2);
                    for(int dy = 0; dy < 5; dy++, row += w) {
                        gx += sobel5_x[dy][0]*((int)row[0]>>8) + sobel5_x[dy][1]*((int)row[1]>>8)
                            + sobel5_x[dy][3]*((int)row[3]>>8) + sobel5_x[dy][4]*((int)row[4]>>8);
                        gy += sobel5_x[0][dy]*((int)row[0]>>8) + sobel5_x[1][dy]*((int)row[1]>>8)
                            + sobel5_x[3][dy]*((int)row[3]>>8) + sobel5_x[4][dy]*((int)row[4]>>8);
                    }
                }
            }
            gx >>= 6; gy >>= 6;
            grad[y*w + x].gx  = (int16_t)gx;
            grad[y*w + x].gy  = (int16_t)gy;
            grad[y*w + x].mag = (int16_t)(abs(gx) + abs(gy));
        }
    }
}

// luma = (77*R + 150*G + 29*B) >> 8  (approx Rec.601)
static void compute_gradient_rgb(const unsigned char *src,
                                 GradPixel *grad,
                                 int w, int h, int bps)
{
    const int bps3 = 3 * bps;

    for(int y = 0; y < h; y++)
    {
        const int border = (y < 2 || y >= h-2);
        for(int x = 0; x < w; x++)
        {
            int gx = 0, gy = 0;
            if(__builtin_expect(border || x < 2 || x >= w-2, 0))
            {
                // slow path : clamp
                for(int dy = -2; dy <= 2; dy++)
                for(int dx = -2; dx <= 2; dx++) {
                    int cx = x+dx < 0 ? 0 : x+dx >= w ? w-1 : x+dx;
                    int cy = y+dy < 0 ? 0 : y+dy >= h ? h-1 : y+dy;
                    const unsigned char *p = src + (cy*w + cx) * bps3;
                    int v;
                    if(bps == 1)
                        v = (77*p[0] + 150*p[1] + 29*p[2]) >> 8;
                    else {
                        const uint16_t *q = (const uint16_t*)p;
                        v = (77*(q[0]>>8) + 150*(q[1]>>8) + 29*(q[2]>>8)) >> 8;
                    }
                    gx += sobel5_x[dy+2][dx+2] * v;
                    gy += sobel5_x[dx+2][dy+2] * v;
                }
            }
            else
            {
                // fast path : accès direct, colonne centrale skippée (=0)
                if(bps == 1) {
                    const unsigned char *row = src + ((y-2)*w + (x-2)) * bps3;
                    const int stride = w * bps3;
                    for(int dy = 0; dy < 5; dy++, row += stride) {
                        // luma des 5 pixels de la ligne (skip dx=2, coef=0)
                        int v0 = (77*row[0]        + 150*row[1]        + 29*row[2])        >> 8;
                        int v1 = (77*row[3]        + 150*row[4]        + 29*row[5])        >> 8;
                        // v2 = row[6..8] : colonne centrale, coef=0, skip
                        int v3 = (77*row[9]        + 150*row[10]       + 29*row[11])       >> 8;
                        int v4 = (77*row[12]       + 150*row[13]       + 29*row[14])       >> 8;
                        gx += sobel5_x[dy][0]*v0 + sobel5_x[dy][1]*v1
                            + sobel5_x[dy][3]*v3 + sobel5_x[dy][4]*v4;
                        gy += sobel5_x[0][dy]*v0 + sobel5_x[1][dy]*v1
                            + sobel5_x[3][dy]*v3 + sobel5_x[4][dy]*v4;
                    }
                } else {
                    const uint16_t *row = (const uint16_t*)src + ((y-2)*w + (x-2)) * 3;
                    const int stride = w * 3;
                    for(int dy = 0; dy < 5; dy++, row += stride) {
                        int v0 = (77*(row[0]>>8) + 150*(row[1]>>8) + 29*(row[2]>>8)) >> 8;
                        int v1 = (77*(row[3]>>8) + 150*(row[4]>>8) + 29*(row[5]>>8)) >> 8;
                        int v3 = (77*(row[9]>>8) + 150*(row[10]>>8)+ 29*(row[11]>>8))>> 8;
                        int v4 = (77*(row[12]>>8)+ 150*(row[13]>>8)+ 29*(row[14]>>8))>> 8;
                        gx += sobel5_x[dy][0]*v0 + sobel5_x[dy][1]*v1
                            + sobel5_x[dy][3]*v3 + sobel5_x[dy][4]*v4;
                        gy += sobel5_x[0][dy]*v0 + sobel5_x[1][dy]*v1
                            + sobel5_x[3][dy]*v3 + sobel5_x[4][dy]*v4;
                    }
                }
            }
            gx >>= 6; gy >>= 6;
            grad[y*w + x].gx  = (int16_t)gx;
            grad[y*w + x].gy  = (int16_t)gy;
            grad[y*w + x].mag = (int16_t)(abs(gx) + abs(gy));
        }
    }
}

static void compute_conf_map(GradPixel *grad, int in_w, int in_h,
                             int mode, int aniso_strength_fp)
{
    for(int y = 0; y < in_h; y++)
    {
        for(int x = 0; x < in_w; x++)
        {
            const int i   = y * in_w + x;
            const int mag = grad[i].mag;
 
            int conf;
            if(mode == 7) {
                if(mag <= ANISO_LOW
                   || is_compression_artifact(grad, x, y, mag, in_w, in_h)) {
                    conf = 0;
                }
                else if(mag >= ANISO_HIGH) {
                    conf = WEIGHT_MAX;
                }
                else {
                    conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
                }
            }
            else {  // mode 6
                if(mag <= ANISO_LOW)       conf = 0;
                else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
                else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
            }
            conf = (conf * aniso_strength_fp) >> WEIGHT_SHIFT;
 
            grad[i].conf = (int16_t)conf;
        }
    }
}

#define COUTURE_N    2                              // largeur bande, px source
#define COUTURE_FAR  ((COUTURE_N + 1) * WEIGHT_MAX) // sentinelle "plein doux"
 
static void compute_couture_dist(const GradPixel * restrict grad,
                                 int16_t * restrict dist_map,
                                 int in_w, int in_h)
{
    const int R   = COUTURE_N;
    const int R2  = R * R;
 
    for(int y = 0; y < in_h; y++)
    {
        for(int x = 0; x < in_w; x++)
        {
            const int i = y * in_w + x;
 
            // pixel net : distance nulle, plein net
            if(grad[i].conf > 0) {
                dist_map[i] = 0;
                continue;
            }
 
            // pixel doux : cherche le net le plus proche dans le disque r<=R
            int best2 = R2 + 1;   // carré de la meilleure distance trouvée
            for(int dj = -R; dj <= R; dj++)
            {
                const int yy = y + dj;
                if(yy < 0 || yy >= in_h) continue;
                const int dj2 = dj * dj;
                for(int di = -R; di <= R; di++)
                {
                    const int d2 = di * di + dj2;
                    if(d2 == 0 || d2 > R2) continue;   // hors disque ou centre
                    if(d2 >= best2)        continue;   // déjà pire
                    const int xx = x + di;
                    if(xx < 0 || xx >= in_w) continue;
                    if(grad[yy * in_w + xx].conf > 0)
                        best2 = d2;
                }
            }
 
            if(best2 > R2) {
                dist_map[i] = (int16_t)COUTURE_FAR;    // aucun net proche
            } else {
                float d = sqrtf((float)best2);
                int df = (int)(d * (float)WEIGHT_MAX + 0.5f);
                dist_map[i] = (int16_t)df;
            }
        }
    }
}



static inline int aniso_pixel_plane(const unsigned char *src,
                                    int sx, int sy, int fx, int fy,
                                    int gx, int gy, int conf,
                                    int in_w, int in_h, int bps)
{
    float glen = sqrtf((float)gx*gx + (float)gy*gy);
    if(glen < 0.5f) glen = 0.5f;
    float nx = gx / glen, ny = gy / glen;
    float ratio = 1.0f + (ANISO_RATIO - 1.0f) * (conf / (float)WEIGHT_MAX);
    const float ffx = fx / (float)WEIGHT_MAX;
    const float ffy = fy / (float)WEIGHT_MAX;
 
    int acc = 0, wsum = 0;
    int ws[36], vs[36];
    int n = 0;
 
    for(int dj = -2; dj <= 3; dj++)
    for(int di = -2; di <= 3; di++)
    {
        float ddx = di - ffx;
        float ddy = dj - ffy;
        float d_across = ddx*nx + ddy*ny;
        float d_along  = -ddx*ny + ddy*nx;
        float r = sqrtf(d_across*d_across + (d_along / ratio) * (d_along / ratio));
 
        int ri = (int)(r * WEIGHT_MAX + 0.5f);
        int ww = (ri < KAISER_RAD_SIZE) ? kaiser_rad_lut[ri] : 0;
        if(ww == 0) continue;
 
        int lx = sx + di, ly = sy + dj;
        lx = lx < 0 ? 0 : lx >= in_w ? in_w-1 : lx;
        ly = ly < 0 ? 0 : ly >= in_h ? in_h-1 : ly;
 
        int v = (bps == 1)
            ? ((int)src[ly*in_w + lx] << 8)
            : (int)((const uint16_t*)src)[ly*in_w + lx];
 
        ws[n] = ww;
        vs[n] = v;
        wsum += ww;
        n++;
    }
 
    if(wsum <= 0) return -1;

    long long acc_ll = 0;
    for(int i = 0; i < n; i++)
        acc_ll += (long long)vs[i] * ws[i];
    return (int)(acc_ll / wsum);
}
 
static inline void aniso_pixel_rgb(const unsigned char *src,
                                   int sx, int sy, int fx, int fy,
                                   int gx, int gy, int conf,
                                   int in_w, int in_h, int bps,
                                   int *out0, int *out1, int *out2)
{
    float glen = sqrtf((float)gx*gx + (float)gy*gy);
    if(glen < 0.5f) glen = 0.5f;
    float nx = gx / glen, ny = gy / glen;
    float ratio = 1.0f + (ANISO_RATIO - 1.0f) * (conf / (float)WEIGHT_MAX);
 
    const float ffx = fx / (float)WEIGHT_MAX;
    const float ffy = fy / (float)WEIGHT_MAX;
    const int bps3 = 3 * bps;
 
    int acc0 = 0, acc1 = 0, acc2 = 0, wsum = 0;
    for(int dj = -2; dj <= 3; dj++)
    for(int di = -2; di <= 3; di++)
    {
        float ddx = di - ffx, ddy = dj - ffy;
        float d_across = ddx*nx + ddy*ny;
        float d_along  = -ddx*ny + ddy*nx;
        float r = sqrtf(d_across*d_across + (d_along / ratio) * (d_along / ratio));
 
        int ri = (int)(r * WEIGHT_MAX + 0.5f);
        int w;
        w = (ri < KAISER_RAD_SIZE) ? kaiser_rad_lut[ri] : 0;
        if(w == 0) continue;
 
        int lx = sx + di, ly = sy + dj;
        lx = lx < 0 ? 0 : lx >= in_w ? in_w-1 : lx;
        ly = ly < 0 ? 0 : ly >= in_h ? in_h-1 : ly;
        const unsigned char *p = src + (ly*in_w + lx) * bps3;
 
        if(bps == 1) {
            acc0 += ((int)p[0] << 8)*w; acc1 += ((int)p[1] << 8)*w; acc2 += ((int)p[2] << 8)*w;
        } else {
            const uint16_t *q = (const uint16_t*)p;
            acc0 += (int)q[0]*w; acc1 += (int)q[1]*w; acc2 += (int)q[2]*w;
        }
        wsum += w;
    }
    if(wsum > 0) {
        *out0 = acc0 / wsum; *out1 = acc1 / wsum; *out2 = acc2 / wsum;
    } else {
        *out0 = -1;   // signale le fallback à l'appelant
    }
}
 
void compute_kaiser_rad_lut(double beta, double base)
{
    const double radius = 3.0;
    const double inv_i0 = 1.0 / bessel_i0(beta);
    const double denom  = pow(base, radius) - 1.0;

    for(int ri = 0; ri < KAISER_RAD_SIZE; ri++)
    {
        double r = (double)ri / (double)WEIGHT_MAX;

        if(r >= radius) {
            kaiser_rad_lut[ri] = 0;
            continue;
        }

        // base cardinale : sinc radial (zéros aux entiers)
        double s;
        if(r == 0.0) s = 1.0;
        else { double px = M_PI * r; s = sin(px) / px; }

        // fenêtre kaiser
        double k = bessel_i0(beta * sqrt(1.0 - (r/radius)*(r/radius))) * inv_i0;

        // courbe exponentielle
        double curve = (pow(base, radius - r) - 1.0) / denom;

        double val = s * k * curve;

        int vi = (val >= 0.0)
               ?  (int)( val * (double)WEIGHT_MAX + 0.5)
               : -(int)(-val * (double)WEIGHT_MAX + 0.5);
        kaiser_rad_lut[ri] = vi;
    }
}


void compute_curve_rad_lut(int mode, double base)
{
    // radius = 1.0 : le noyau ne touche que les pixels dans le cercle
    // unité. À r=1 le poids est 0 → cardinal (pixel source préservé).
    // Pas de lobes négatifs : décroissance monotone de 1 (r=0) à 0 (r=1).
    // C'est le "bilinéaire radial" — pendant 2D du bilinéaire pondéré 1D.
    const double radius = 1.0;
    const double denom_exp = pow(base, radius) - 1.0;
 
    for(int ri = 0; ri < CURVE_RAD_SIZE; ri++)
    {
        double r = (double)ri / (double)WEIGHT_MAX;
 
        double v;
        if(r >= radius) {
            v = 0.0;
        } else if(mode == 3) {
            // tent radiale : décroissance linéaire
            v = 1.0 - r / radius;
        } else {
            // version exponentielle de la tent
            v = (pow(base, radius - r) - 1.0) / denom_exp;
        }
 
        int vi = (int)(v * (double)WEIGHT_MAX + 0.5);
        curve_rad_lut[ri] = vi < 0 ? 0 : vi;
    }
}

static inline int luma_at(const unsigned char *s,int x,int y,int w,int h,int bps,int rgb){
    x=x<0?0:x>=w?w-1:x; y=y<0?0:y>=h?h-1:y;
    if(rgb){ const unsigned char *p=s+(y*w+x)*3*bps;
        if(bps==1) return (77*p[0]+150*p[1]+29*p[2])>>8;
        const uint16_t *q=(const uint16_t*)p; return (77*(q[0]>>8)+150*(q[1]>>8)+29*(q[2]>>8))>>8; }
    return (bps==1)? (int)s[y*w+x] : (int)((const uint16_t*)s)[y*w+x]>>8;
}

// Détecte la phase de grille verticale (frontières DCT) sur une frame.
// Retourne la phase [0..7], ou -1 si aucune grille fiable détectée.
static int compute_grid_phase(const unsigned char *src,
                              int w, int h, int bps, int rgb)
{
    #define GRID_SHARP_MIN 1.30

    long *profil_v = calloc(w, sizeof(long));
    if(!profil_v) return -1;

    #define LU(X,Y) luma_at(src,(X),(Y),w,h,bps,rgb)
    // rupture de courbure horizontale par pixel, accumulée par colonne
    for(int y = 0; y < h; y++) {
        for(int x = 2; x < w - 2; x++) {
            int gL = LU(x,y)   - LU(x-1,y);
            int gR = LU(x+1,y) - LU(x,y);
            int rupture = gR - gL;
            if(rupture < 0) rupture = -rupture;
            profil_v[x] += rupture;
        }
    }
    #undef LU

    // détection de phase : 8 décalages, on garde le max
    long best = -1, sumall = 0; int phase = 0;
    for(int phi = 0; phi < 8; phi++) {
        long s = 0;
        for(int x = phi; x < w; x += 8) s += profil_v[x];
        sumall += s;
        if(s > best) { best = s; phase = phi; }
    }
    double mean = (double)sumall / 8.0;
    double sharp = (mean > 0) ? (double)best / mean : 0.0;

    free(profil_v);

    return (sharp >= GRID_SHARP_MIN) ? phase : -1;
    #undef GRID_SHARP_MIN
}

// Post-passe de déblocking : floute horizontalement les colonnes de grille
// dans l'image de sortie. Purement horizontal (même ligne uniquement),
// donc le vertical interpolé par le scaling reste intact.
// - estompage vers les bords de zone (pas de bordure nette)
// - radius réduit sur fort contraste (pas de moyenne aberrante au milieu)
static void post_deblock_rgb(unsigned char *out_buf,
                             const BresenhamMap *map,
                             int out_w, int out_h, int ratio,
                             int bps_out, int phase_v)
{
    if(phase_v < 0) return;

    const int bps3 = 3 * bps_out;
    #define SRC_SPAN    2     // pixels source de chaque côté de la frontière
    #define MAXOUT      8     // radius max du flou (faible contraste)
    #define MAXOUT_MIN  2     // radius min du flou (fort contraste)
    #define CONTRAST_LO 20    // sous ce contraste : radius max
    #define CONTRAST_HI 80    // au-dessus : radius min
    static const int exp_w[MAXOUT+1] = {256,128,64,32,16,8,4,2,1};

    int *tmp = malloc(out_w * 3 * sizeof(int));
    if(!tmp) return;

    for(int f = 0; f < ratio; f++)
    {
        const int map_w_base = f * out_w;
        unsigned char *frame = out_buf + (size_t)f * out_w * out_h * 3 * bps_out;

        for(int y = 0; y < out_h; y++)
        {
            unsigned char *row = frame + (size_t)y * out_w * bps3;

            for(int x = 0; x < out_w; x++)
            {
                int sx = map->w_src[map_w_base + x];
                if(((sx - phase_v) % 8 + 8) % 8 != 0) continue;  // pas sur la grille

                // étend de SRC_SPAN pixels source de chaque côté via changements de w_src
                int lo = x, cnt = 0, prev = sx;
                while(lo > 0 && cnt < SRC_SPAN) {
                    lo--;
                    int s = map->w_src[map_w_base + lo];
                    if(s != prev) { cnt++; prev = s; }
                }
                int hi = x; cnt = 0; prev = sx;
                while(hi < out_w-1 && cnt < SRC_SPAN) {
                    hi++;
                    int s = map->w_src[map_w_base + hi];
                    if(s != prev) { cnt++; prev = s; }
                }

                // contraste à travers la zone (luma bord gauche vs bord droit)
                int lum_lo, lum_hi;
                {
                    const unsigned char *pl = row + lo*bps3;
                    const unsigned char *ph = row + hi*bps3;
                    if(bps_out==1){
                        lum_lo=(77*pl[0]+150*pl[1]+29*pl[2])>>8;
                        lum_hi=(77*ph[0]+150*ph[1]+29*ph[2])>>8;
                    } else {
                        const uint16_t*ql=(const uint16_t*)pl;
                        const uint16_t*qh=(const uint16_t*)ph;
                        lum_lo=(77*(ql[0]>>8)+150*(ql[1]>>8)+29*(ql[2]>>8))>>8;
                        lum_hi=(77*(qh[0]>>8)+150*(qh[1]>>8)+29*(qh[2]>>8))>>8;
                    }
                }
                int contrast = lum_hi - lum_lo; if(contrast<0) contrast=-contrast;

                // radius effectif : grand si contraste faible, petit si fort
                int radius;
                if(contrast <= CONTRAST_LO) radius = MAXOUT;
                else if(contrast >= CONTRAST_HI) radius = MAXOUT_MIN;
                else radius = MAXOUT_MIN +
                     (MAXOUT - MAXOUT_MIN) * (CONTRAST_HI - contrast) / (CONTRAST_HI - CONTRAST_LO);

                // floute chaque pixel de [lo..hi], avec estompage vers les bords
                for(int px = lo; px <= hi; px++)
                {
                    int a0=0,a1=0,a2=0,wsum=0;
                    for(int d = -radius; d <= radius; d++) {
                        int lx = px + d;
                        if(lx < lo) lx = lo; if(lx > hi) lx = hi;
                        int wt = exp_w[d<0?-d:d];
                        const unsigned char *p = row + lx*bps3;
                        if(bps_out == 1) { a0+=p[0]*wt; a1+=p[1]*wt; a2+=p[2]*wt; }
                        else { const uint16_t*q=(const uint16_t*)p; a0+=q[0]*wt; a1+=q[1]*wt; a2+=q[2]*wt; }
                        wsum += wt;
                    }
                    int b0=a0/wsum, b1=a1/wsum, b2=a2/wsum;

                    // estompage : 1 au centre (colonne de grille x), 0 aux bords
                    int dc = px - x; if(dc < 0) dc = -dc;
                    int span = (px <= x) ? (x - lo) : (hi - x);
                    if(span < 1) span = 1;
                    int aw = WEIGHT_MAX - (dc * WEIGHT_MAX / span);
                    if(aw < 0) aw = 0; if(aw > WEIGHT_MAX) aw = WEIGHT_MAX;
                    aw = smoothstep_lut[aw];
                    int iaw = WEIGHT_MAX - aw;

                    const unsigned char *po = row + px*bps3;
                    int o0,o1,o2;
                    if(bps_out==1){o0=po[0];o1=po[1];o2=po[2];}
                    else{const uint16_t*q=(const uint16_t*)po;o0=q[0];o1=q[1];o2=q[2];}

                    tmp[px*3+0]=(b0*aw+o0*iaw)>>WEIGHT_SHIFT;
                    tmp[px*3+1]=(b1*aw+o1*iaw)>>WEIGHT_SHIFT;
                    tmp[px*3+2]=(b2*aw+o2*iaw)>>WEIGHT_SHIFT;
                }

                // recopie séparée (ne pas flouter à partir de pixels déjà modifiés)
                for(int px = lo; px <= hi; px++) {
                    unsigned char *p = row + px*bps3;
                    if(bps_out == 1) { p[0]=tmp[px*3]; p[1]=tmp[px*3+1]; p[2]=tmp[px*3+2]; }
                    else { uint16_t*q=(uint16_t*)p; q[0]=tmp[px*3]; q[1]=tmp[px*3+1]; q[2]=tmp[px*3+2]; }
                }

                x = hi;
            }
        }
    }
    free(tmp);
    #undef SRC_SPAN
    #undef MAXOUT
    #undef MAXOUT_MIN
    #undef CONTRAST_LO
    #undef CONTRAST_HI
}

// Post-passe de déblocking pour un plan (1 canal). Floute horizontalement
// les colonnes de grille. Purement horizontal → vertical intact.
static void post_deblock_plane(unsigned char *out_plane,
                               const BresenhamMap *map,
                               int out_w, int out_h, int ratio,
                               int bps_out, int phase_v)
{
    if(phase_v < 0) return;

    #define SRC_SPAN    2
    #define MAXOUT      8
    #define MAXOUT_MIN  2
    #define CONTRAST_LO 20
    #define CONTRAST_HI 80
    static const int exp_w[MAXOUT+1] = {256,128,64,32,16,8,4,2,1};

    int *tmp = malloc(out_w * sizeof(int));
    if(!tmp) return;

    for(int f = 0; f < ratio; f++)
    {
        const int map_w_base = f * out_w;
        unsigned char *frame = out_plane + (size_t)f * out_w * out_h * bps_out;

        for(int y = 0; y < out_h; y++)
        {
            unsigned char *row = frame + (size_t)y * out_w * bps_out;

            #define RDP(px) ((bps_out==1) ? (int)row[px] : (int)((uint16_t*)row)[px])

            for(int x = 0; x < out_w; x++)
            {
                int sx = map->w_src[map_w_base + x];
                if(((sx - phase_v) % 8 + 8) % 8 != 0) continue;

                int lo = x, cnt = 0, prev = sx;
                while(lo > 0 && cnt < SRC_SPAN) {
                    lo--;
                    int s = map->w_src[map_w_base + lo];
                    if(s != prev) { cnt++; prev = s; }
                }
                int hi = x; cnt = 0; prev = sx;
                while(hi < out_w-1 && cnt < SRC_SPAN) {
                    hi++;
                    int s = map->w_src[map_w_base + hi];
                    if(s != prev) { cnt++; prev = s; }
                }

                int contrast = RDP(hi) - RDP(lo); if(contrast<0) contrast=-contrast;
                // note : contraste sur le plan, en échelle 16-bit si bps==2
                int contrast8 = (bps_out==1) ? contrast : (contrast >> 8);

                int radius;
                if(contrast8 <= CONTRAST_LO) radius = MAXOUT;
                else if(contrast8 >= CONTRAST_HI) radius = MAXOUT_MIN;
                else radius = MAXOUT_MIN +
                     (MAXOUT - MAXOUT_MIN) * (CONTRAST_HI - contrast8) / (CONTRAST_HI - CONTRAST_LO);

                for(int px = lo; px <= hi; px++)
                {
                    int acc=0, wsum=0;
                    for(int d = -radius; d <= radius; d++) {
                        int lx = px + d;
                        if(lx < lo) lx = lo; if(lx > hi) lx = hi;
                        int wt = exp_w[d<0?-d:d];
                        acc += RDP(lx)*wt;
                        wsum += wt;
                    }
                    int b = acc/wsum;

                    int dc = px - x; if(dc < 0) dc = -dc;
                    int span = (px <= x) ? (x - lo) : (hi - x);
                    if(span < 1) span = 1;
                    int aw = WEIGHT_MAX - (dc * WEIGHT_MAX / span);
                    if(aw < 0) aw = 0; if(aw > WEIGHT_MAX) aw = WEIGHT_MAX;
                    aw = smoothstep_lut[aw];
                    int iaw = WEIGHT_MAX - aw;

                    int o = RDP(px);
                    tmp[px] = (b*aw + o*iaw) >> WEIGHT_SHIFT;
                }

                for(int px = lo; px <= hi; px++) {
                    if(bps_out==1) row[px] = (unsigned char)tmp[px];
                    else ((uint16_t*)row)[px] = (uint16_t)tmp[px];
                }

                x = hi;
            }
            #undef RDP
        }
    }
    free(tmp);
    #undef SRC_SPAN
    #undef MAXOUT
    #undef MAXOUT_MIN
    #undef CONTRAST_LO
    #undef CONTRAST_HI
}

void usage(void)
{
    fprintf(stderr,
        "temporal_scaling - temporal pixel art upscaler\n\n"
        "Usage:\n"
        "\t-i filename        input raw video file (use '-' for stdin)\n\n"
        "Options:\n"
        "\t--in-width  INT    input width  (default: 720)\n"
        "\t--in-height INT    input height (default: 480)\n"
        "\t--out-width INT    output width  (default: 1920)\n"
        "\t--out-height INT   output height (default: 1080)\n"
        "\t--fps-ratio INT    output/input fps multiplier (default: 2)\n"
        "\t--mode INT         scaling mode (default: 0)\n"
		"\t                     0 = nearest neighbor\n"
		"\t                     1 = bilinear on duplicated pixels\n"
		"\t                     2 = bilinear with temporal alternation\n"
		"\t                     3 = linear weighted 2D with temporal alternation\n"
		"\t                     4 = exponential weighted 2D with temporal alternation\n"
		"\t                     5 = kaiser 6x6 EWA isotropic\n"
		"\t                     6 = kaiser 6x6 EWA anisotropic edge-directed\n"
		"\t                     7 = kaiser 6x6 EWA aniso + selective softening\n"
        "\t--tmp-mode INT     temporal blending mode (default: 0)\n"
        "\t                     0 = duplicate frames (no blending)\n"
        "\t                     1 = blend 50%% on the upper minority of phases\n"
        "\t                     2 = linear weighted blending\n"
        "\t                     3 = exponential weighted blending (uses --curve-base)\n"
        "\t--curve-base FLOAT exponential curve base for mode 4 / tmp-mode 3 (default: 2.0, must be > 1.0)\n\n"
		"\t--kaiser-beta FLOAT  beta parameter of the kaiser kernel for modes 5/6 (default: 2.5)\n"
		"\t                     low = sharp with ringing, high = smooth without ringing\n"
		"\t                     note: --1d pipeline uses a separable 4-tap kaiser instead\n"
		"\t--aniso-strength FLOAT  anisotropy strength for mode 6 (default: 1.0)\n"
		"\t                     0.0 = isotropic (identical to mode 5)\n"
		"\t                     1.0 = full anisotropy on strong edges\n"
		"\t--1d               use separable 1D pipeline (default: 2D)\n"
		"\t--in-pix-fmt  STR  RGB, 444p, 422p, 420p, 411p, 410p,\n"
		"\t                   RGB48, 444p16, 422p16, 420p16\n"
		"\t--out-pix-fmt STR  output pixel format : RGB, 444p, 422p, 420p, 411p, 410p,\n"
		"\t                   RGB48, 444p16, 422p16, 420p16\n"
		"                     RGB<->YUV conversion not supported\n"
    );
    exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

void adjust_phase_direction(BresenhamMap *map, int in_w, int in_h, int out_w, int out_h, int ratio) {
    const int do_w = (out_w * 2 >= in_w * 5);
    const int do_h = (out_h * 2 >= in_h * 5);

    for (int r = 0; r < ratio; r++) {
        const int is_odd_frame = (r % 2 == 1);  // Vrai si frame impaire

        if (do_w) {
            int base = r * out_w;
            int px = 0;
            while (px < out_w) {
                int src = map->w_src[base + px];
                int run_start = px;
                // Détecte la fin du run de duplication (même src)
                while (px + 1 < out_w && map->w_src[base + px + 1] == src) {
                    px++;
                }
                int run_length = px - run_start + 1;

                // On ne traite que les runs de duplication (longueur >= 2)
                if (run_length >= 2) {
                    int half = run_length / 2;
                    if (is_odd_frame) {
                        // Frame impaire : inverse les half premiers pixels du run
                        for (int i = 0; i < half; i++) {
                            map->w_flag[base + run_start + i] = -1;
                        }
                        for (int i = half; i < run_length; i++) {
                            map->w_flag[base + run_start + i] = 0;
                        }
                    } else {
                        // Frame paire : met les half derniers pixels du run à 1
                        for (int i = 0; i < run_length - half; i++) {
                            map->w_flag[base + run_start + i] = 0;
                        }
                        for (int i = run_length - half; i < run_length; i++) {
                            map->w_flag[base + run_start + i] = 1;
                        }
                    }
                }
                px++;  // Passe au pixel suivant
            }
        }

        if (do_h) {
            int base = r * out_h;
            int px = 0;
            while (px < out_h) {
                int src = map->h_src[base + px];
                int run_start = px;
                // Détecte la fin du run de duplication (même src)
                while (px + 1 < out_h && map->h_src[base + px + 1] == src) {
                    px++;
                }
                int run_length = px - run_start + 1;

                // On ne traite que les runs de duplication (longueur >= 2)
                if (run_length >= 2) {
                    int half = run_length / 2;
                    if (is_odd_frame) {
                        // Frame impaire : inverse les half premiers pixels du run
                        for (int i = 0; i < half; i++) {
                            map->h_flag[base + run_start + i] = -1;
                        }
                        for (int i = half; i < run_length; i++) {
                            map->h_flag[base + run_start + i] = 0;
                        }
                    } else {
                        // Frame paire : met les half derniers pixels du run à 1
                        for (int i = 0; i < run_length - half; i++) {
                            map->h_flag[base + run_start + i] = 0;
                        }
                        for (int i = run_length - half; i < run_length; i++) {
                            map->h_flag[base + run_start + i] = 1;
                        }
                    }
                }
                px++;  // Passe au pixel suivant
            }
        }
    }
}

void adjust_phase_direction_weighted(BresenhamMap *map, int in_w, int in_h, int out_w, int out_h, int ratio)
{
    const int do_w = (out_w * 2 >= in_w * 5);
    const int do_h = (out_h * 2 >= in_h * 5);

    for(int r = 1; r < ratio; r += 2)
    {
        if(do_w)
        {
            int base = r * out_w;
            int px = 0;
            while(px < out_w)
            {
                int src = map->w_src[base + px];
                int run_end = px + 1;
                while(run_end < out_w && map->w_src[base + run_end] == src)
                    run_end++;

                const int N = run_end - px;
                if(N > 1)
                {
                    int tmp_w[N], tmp_h[N];
                    for(int i = 0; i < N; i++) {
                        tmp_w[i] = map->w_weight[base + px + i];
                        tmp_h[i] = map->w_heavy [base + px + i];
                    }
                    for(int i = 0; i < N; i++) {
                        int mi = N - 1 - i;
                        map->w_weight[base + px + i] = tmp_w[mi];
                        map->w_heavy [base + px + i] = tmp_h[mi];
                        map->w_flag  [base + px + i] = (tmp_w[mi] > 0) ? -1 : 0;
                    }
                }
				else  // N == 1 : pixel isolé dans une phase impaire
				{
					if(map->w_weight[base + px] > 0)
						map->w_flag[base + px] = -1;
				}
                px = run_end;
            }
        }

        if(do_h)
        {
            int base = r * out_h;
            int px = 0;
            while(px < out_h)
            {
                int src = map->h_src[base + px];
                int run_end = px + 1;
                while(run_end < out_h && map->h_src[base + run_end] == src)
                    run_end++;

                const int N = run_end - px;
                if(N > 1)
                {
                    int tmp_w[N], tmp_h[N];
                    for(int i = 0; i < N; i++) {
                        tmp_w[i] = map->h_weight[base + px + i];
                        tmp_h[i] = map->h_heavy [base + px + i];
                    }
                    for(int i = 0; i < N; i++) {
                        int mi = N - 1 - i;
                        map->h_weight[base + px + i] = tmp_w[mi];
                        map->h_heavy [base + px + i] = tmp_h[mi];
                        map->h_flag  [base + px + i] = (tmp_w[mi] > 0) ? -1 : 0;
                    }
                }
				else  // N == 1
				{
					if(map->h_weight[base + px] > 0)
						map->h_flag[base + px] = -1;
				}
                px = run_end;
            }
        }
    }
}

void compute_scale_map(BresenhamMap *map, int in_w, int in_h, int out_w, int out_h, int ratio)
{
    for(int r = 0; r < ratio; r++)
    {
        int err = (int)(((long long)(2*r + 1) * out_w) / (2 * ratio));
        if(err >= out_w) err -= out_w;
        int src = 0;
        int dst = r * out_w;
        for(int px = 0; px < out_w; px++)
        {
            map->w_src [dst] = (src < in_w) ? src : in_w - 1;
            map->w_flag[dst] = (2 * err >= out_w) ? 1 : 0;
			int clamped_src = (src < in_w) ? src : in_w - 1;
			map->w_pos[dst] = clamped_src * WEIGHT_MAX + (err * WEIGHT_MAX) / out_w;
            dst++;
            err += in_w;
            if(err >= out_w) { err -= out_w; src++; }
        }
    }

    for(int r = 0; r < ratio; r++)
    {
        int err = (int)(((long long)(2*r + 1) * out_h) / (2 * ratio));
        if(err >= out_h) err -= out_h;
        int src = 0;
        int dst = r * out_h;
        for(int px = 0; px < out_h; px++)
        {
            map->h_src [dst] = (src < in_h) ? src : in_h - 1;
            map->h_flag[dst] = (2 * err >= out_h) ? 1 : 0;
			int clamped_src = (src < in_h) ? src : in_h - 1;
			map->h_pos[dst] = clamped_src * WEIGHT_MAX + (err * WEIGHT_MAX) / out_h;
            dst++;
            err += in_h;
            if(err >= out_h) { err -= out_h; src++; }
        }
    }
}

typedef enum {
    CURVE_EXPONENTIAL = 0,  // mode 3
    CURVE_LINEAR      = 1   // mode 4
} CurveType;

typedef enum {
    TMP_DUPLICATE   = 0,  // pas de blend, comportement par défaut
    TMP_HALF        = 1,  // blend 50% sur la minorité supérieure des phases
    TMP_LINEAR      = 2,  // poids linéaire i/R
    TMP_EXPONENTIAL = 3   // poids exponentiel (base^i - 1) / base^(R-1)
} TmpMode;

// Pré-calcule les poids temporels (en virgule fixe 0..WEIGHT_MAX) pour
// chaque phase r ∈ [0, R-1]. w=0 => pas de blend (frame source pure).
//
// Logique miroir du spatial :
//   TMP_HALF        : seules les floor(R/2) dernières phases sont blendées
//                     à 50%, le reste reste pur. Garantit qu'au moins la
//                     moitié des frames produites sont identiques à la
//                     frame source originale (plus de poids à l'original
//                     qu'à l'interpolation).
//   TMP_LINEAR      : w_r = r / R
//   TMP_EXPONENTIAL : w_r = (base^r - 1) / base^(R-1), 1ère moitié quasi
//                     identique à la source, cassure sur la fin.
//
// Exemples :
//   R=3 mode 1 : 0%, 0%, 50%
//   R=3 mode 2 : 0%, 33%, 67%
//   R=3 mode 3 base=2 : 0%, 25%, 75%
//   R=4 mode 1 : 0%, 0%, 50%, 50%
//   R=4 mode 2 : 0%, 25%, 50%, 75%
//   R=4 mode 3 base=2 : 0%, 12.5%, 37.5%, 87.5%
void compute_temporal_weights(int *tmp_weights, int ratio, TmpMode mode, double base)
{
    if(ratio <= 1)
    {
        if(ratio == 1) tmp_weights[0] = 0;
        return;
    }

    if(mode == TMP_DUPLICATE)
    {
        for(int r = 0; r < ratio; r++) tmp_weights[r] = 0;
        return;
    }

    if(mode == TMP_HALF)
    {
        const int threshold = (ratio + 1) / 2;  // ceil(R/2)
        for(int r = 0; r < ratio; r++)
            tmp_weights[r] = (r >= threshold) ? (WEIGHT_MAX / 2) : 0;
        return;
    }

    // TMP_LINEAR ou TMP_EXPONENTIAL
    const double denom = (mode == TMP_EXPONENTIAL) ? pow(base, (double)(ratio - 1)) : 1.0;
    for(int r = 0; r < ratio; r++)
    {
        double w;
        if(mode == TMP_LINEAR)
            w = (double)r / (double)ratio;
        else  // TMP_EXPONENTIAL
            w = (pow(base, (double)r) - 1.0) / denom;

        int wi = (int)(w * WEIGHT_MAX + 0.5);
        if(wi < 0)          wi = 0;
        if(wi > WEIGHT_MAX) wi = WEIGHT_MAX;
        tmp_weights[r] = wi;
    }
}

// Mode 3/4 : pondération par run de pixels dupliqués.
//
// Pour un run de N pixels partageant la même source, le i-ème pixel reçoit :
//   CURVE_EXPONENTIAL : w_i = (base^i - 1) / base^(N-1)
//   CURVE_LINEAR      : w_i = i / N
//
// Les deux formules garantissent w_0 = 0 (1er pixel = source pure) et ne
// touchent jamais 100% pour le dernier pixel du run, ce qui préserve
// toujours un peu d'info source à la frontière entre runs consécutifs.
//
// Le 1er pixel du run (i=0) a flag=0 (pas d'interpolation). Les autres
// interpolent vers le voisin de droite par défaut (flag=+1).
//
// On marque w_heavy/h_heavy = 1 pour les floor(N/2) pixels du run qui
// portent les poids les plus forts (les derniers indices, par construction
// monotone croissante de la courbe). C'est le critère utilisé par
// adjust_phase_direction_weighted pour décider d'inverser la direction.
//
// Heavy par run (i >= ceil(N/2), donc floor(N/2) pixels) :
//   N=2 : i=1                (1 pixel)
//   N=3 : i=2                (1 pixel)
//   N=4 : i=2,3              (2 pixels)
//   N=5 : i=3,4              (2 pixels)
//
// Cette fonction écrase les valeurs de w_flag/h_flag posées par
// compute_scale_map (logique différente du mode 2).
//
// Exemples de poids (en %) :
//   linéaire (i/N)         exponentiel base=2
//   N=2 : 0, 50            0, 50
//   N=3 : 0, 33, 67        0, 25, 75
//   N=4 : 0, 25, 50, 75    0, 12.5, 37.5, 87.5
//   N=5 : 0, 20, 40, 60, 80    0, 6.25, 18.75, 43.75, 93.75
void compute_curve_weights(BresenhamMap *map, int out_w, int out_h, int ratio,
                            CurveType curve, double base)
{
    // axe horizontal
    for(int r = 0; r < ratio; r++)
    {
        const int row = r * out_w;
        int run_start = 0;
        while(run_start < out_w)
        {
            const int src = map->w_src[row + run_start];
            int run_end = run_start + 1;
            while(run_end < out_w && map->w_src[row + run_end] == src)
                run_end++;

            const int N = run_end - run_start;

            if(N == 1)
            {
                // Cas spécial : run d'un seul pixel (typique des ratios
                // proches de 1, ex 1080->1440). Pas de "position dans le run"
                // à exploiter, donc on conserve le flag Bresenham déjà posé
                // par compute_scale_map et on applique 50% comme en mode 2.
                // Sinon ces pixels seraient laissés purs source et le rendu
                // serait perceptuellement plus pixelisé que le mode 2.
                const int f = map->w_flag[row + run_start];
                map->w_weight[row + run_start] = (f != 0) ? (WEIGHT_MAX / 2) : 0;
                map->w_heavy [row + run_start] = (f != 0) ? 1 : 0;
            }
            else
            {
                const int heavy_threshold = (N + 1) / 2;  // ceil(N/2)
                const double denom = (curve == CURVE_EXPONENTIAL)
                                      ? pow(base, (double)(N - 1)) : 1.0;

                for(int i = 0; i < N; i++)
                {
                    double w;
                    if(curve == CURVE_LINEAR)
                        w = (double)i / (double)N;
                    else  // CURVE_EXPONENTIAL
                        w = (pow(base, (double)i) - 1.0) / denom;

                    int wi = (int)(w * WEIGHT_MAX + 0.5);
                    if(wi < 0)          wi = 0;
                    if(wi > WEIGHT_MAX) wi = WEIGHT_MAX;
                    map->w_weight[row + run_start + i] = wi;
                    map->w_flag  [row + run_start + i] = (i == 0) ? 0 : 1;
                    map->w_heavy [row + run_start + i] = (i >= heavy_threshold) ? 1 : 0;
                }
            }
            run_start = run_end;
        }
    }

    // axe vertical
    for(int r = 0; r < ratio; r++)
    {
        const int row = r * out_h;
        int run_start = 0;
        while(run_start < out_h)
        {
            const int src = map->h_src[row + run_start];
            int run_end = run_start + 1;
            while(run_end < out_h && map->h_src[row + run_end] == src)
                run_end++;

            const int N = run_end - run_start;

            if(N == 1)
            {
                const int f = map->h_flag[row + run_start];
                map->h_weight[row + run_start] = (f != 0) ? (WEIGHT_MAX / 2) : 0;
                map->h_heavy [row + run_start] = (f != 0) ? 1 : 0;
            }
            else
            {
                const int heavy_threshold = (N + 1) / 2;
                const double denom = (curve == CURVE_EXPONENTIAL)
                                      ? pow(base, (double)(N - 1)) : 1.0;

                for(int i = 0; i < N; i++)
                {
                    double w;
                    if(curve == CURVE_LINEAR)
                        w = (double)i / (double)N;
                    else  // CURVE_EXPONENTIAL
                        w = (pow(base, (double)i) - 1.0) / denom;

                    int wi = (int)(w * WEIGHT_MAX + 0.5);
                    if(wi < 0)          wi = 0;
                    if(wi > WEIGHT_MAX) wi = WEIGHT_MAX;
                    map->h_weight[row + run_start + i] = wi;
                    map->h_flag  [row + run_start + i] = (i == 0) ? 0 : 1;
                    map->h_heavy [row + run_start + i] = (i >= heavy_threshold) ? 1 : 0;
                }
            }
            run_start = run_end;
        }
    }
}

// ============================================================
// process_files_plane — version 1D, 1 canal planar
// Remplace process_files pour les formats YUV planar.
// Appelé 3 fois (Y, U, V) avec des pointeurs différents.
// ============================================================

void process_files_plane(const unsigned char * restrict in_plane,
                         const unsigned char * restrict next_plane,
                         unsigned char * restrict blend_plane,
                         unsigned char * restrict out_plane,
                         unsigned char * restrict tmp_plane,
                         const BresenhamMap * restrict map,
                         const int *tmp_weights,
                         const int in_w, const int in_h,
                         const int out_w, const int out_h,
                         const int ratio, const int mode,
                         const int bps_in, const int bps_out)
{
    const int in_size  = in_w * in_h;  // en pixels, pas en bytes
    const int bilinear = (mode >= 1 && mode != 5);
    const int weighted = (mode == 3 || mode == 4);
 
    // Lecture : normalise vers l'espace de calcul interne (toujours 0..65535)
    // pour que les interpolations soient cohérentes quelle que soit la profondeur source
    #define RD(ptr, x) \
        ((bps_in == 1) ? ((int)(ptr)[x] << 8) \
                       : (int)((const uint16_t*)(ptr))[x])
	
	#define RD_TMP(ptr, x) \
    ((bps_out == 1) ? ((int)(ptr)[x] << 8) \
                    : (int)((const uint16_t*)(ptr))[x])
 
    // Écriture : dénormalise depuis l'espace 16-bit vers la profondeur cible
    #define WR(ptr, x, v) do { \
        int _v = (v); \
        _v = _v < 0 ? 0 : _v > 65535 ? 65535 : _v; \
        if(bps_out == 1) (ptr)[x]               = (unsigned char)(_v >> 8); \
        else ((uint16_t*)(ptr))[x]              = (uint16_t)_v; \
    } while(0)
 
    unsigned int offset = 0;
    for(int r = 0; r < ratio; r++)
    {
        const int tw = tmp_weights[r];
        const unsigned char *src_plane;
        if(tw > 0 && next_plane != in_plane)
        {
            const int itw = WEIGHT_MAX - tw;
            if(bps_in == 1) {
                for(int i = 0; i < in_size; i++)
                    blend_plane[i] = (in_plane[i]*itw + next_plane[i]*tw) >> WEIGHT_SHIFT;
            } else {
                const uint16_t *in16   = (const uint16_t *)in_plane;
                const uint16_t *next16 = (const uint16_t *)next_plane;
                      uint16_t *bl16   = (uint16_t *)blend_plane;
                for(int i = 0; i < in_size; i++)
                    bl16[i] = (in16[i]*itw + next16[i]*tw) >> WEIGHT_SHIFT;
            }
            src_plane = blend_plane;
        }
        else
        {
            src_plane = in_plane;
        }
 
        const int map_w_base = r * out_w;
        const int map_h_base = r * out_h;
 
        // passe horizontale : in_plane → tmp_plane
        for(int h = 0; h < in_h; h++)
        {
            const unsigned char *src_row = src_plane + h * in_w  * bps_in;
			unsigned char       *dst_row = tmp_plane + h * out_w * bps_out;
 
            for(int px = 0; px < out_w; px++)
            {
                const int src_px      = map->w_src [map_w_base + px];
                const int is_bilinear = map->w_flag[map_w_base + px];
 
                if(mode == 5 && is_bilinear != 0)
                {
                    const int w = map->w_weight[map_w_base + px];
                    if(w == 0) {
                        WR(dst_row, px, RD(src_row, src_px));
                    } else {
                        const int kh0 = kaiser4_lut[w][0];
                        const int kh1 = kaiser4_lut[w][1];
                        const int kh2 = kaiser4_lut[w][2];
                        const int kh3 = kaiser4_lut[w][3];
                        int p1, p2, p3, p4;
                        if(is_bilinear == 1) {
                            p1 = __builtin_expect(src_px >= 1,       1) ? src_px-1 : src_px;
                            p2 = src_px;
                            p3 = __builtin_expect(src_px+1 < in_w,   1) ? src_px+1 : src_px;
                            p4 = __builtin_expect(src_px+2 < in_w,   1) ? src_px+2 : p3;
                        } else {
                            p1 = __builtin_expect(src_px+1 < in_w,   1) ? src_px+1 : src_px;
                            p2 = src_px;
                            p3 = __builtin_expect(src_px >= 1,        1) ? src_px-1 : src_px;
                            p4 = __builtin_expect(src_px >= 2,        1) ? src_px-2 : p3;
                        }
                        int v = (RD(src_row,p1)*kh0 + RD(src_row,p2)*kh1
                               + RD(src_row,p3)*kh2 + RD(src_row,p4)*kh3) >> WEIGHT_SHIFT;
                        WR(dst_row, px, v);
                    }
                }
                else if(weighted && is_bilinear == 1 && src_px+1 < in_w)
                {
                    const int w  = map->w_weight[map_w_base + px];
                    const int iw = WEIGHT_MAX - w;
                    WR(dst_row, px, (RD(src_row,src_px)*iw + RD(src_row,src_px+1)*w) >> WEIGHT_SHIFT);
                }
                else if(weighted && is_bilinear == -1 && src_px >= 1)
                {
                    const int w  = map->w_weight[map_w_base + px];
                    const int iw = WEIGHT_MAX - w;
                    WR(dst_row, px, (RD(src_row,src_px)*iw + RD(src_row,src_px-1)*w) >> WEIGHT_SHIFT);
                }
                else if(bilinear && is_bilinear == 1 && src_px+1 < in_w)
                {
                    WR(dst_row, px, (RD(src_row,src_px) + RD(src_row,src_px+1)) >> 1);
                }
                else if(bilinear && is_bilinear == -1 && src_px >= 1)
                {
                    WR(dst_row, px, (RD(src_row,src_px-1) + RD(src_row,src_px)) >> 1);
                }
                else
                {
                    WR(dst_row, px, RD(src_row, src_px));
                }
            }
        }
 
        // passe verticale : tmp_plane → out_plane
        for(int px = 0; px < out_h; px++)
        {
            const int src_line    = map->h_src [map_h_base + px];
            const int is_bilinear = map->h_flag[map_h_base + px];
 
            const unsigned char *tmp_row0 = tmp_plane + src_line * out_w * bps_out;
			unsigned char *dst_row = out_plane + offset + px * out_w * bps_out;
 
            if(mode == 5 && is_bilinear != 0)
            {
                const int w = map->h_weight[map_h_base + px];
                if(w == 0) {
                    memcpy(dst_row, tmp_row0, out_w * bps_out);
                } else {
                    const int kv0 = kaiser4_lut[w][0];
                    const int kv1 = kaiser4_lut[w][1];
                    const int kv2 = kaiser4_lut[w][2];
                    const int kv3 = kaiser4_lut[w][3];
                    int l1, l2, l3, l4;
                    if(is_bilinear == 1) {
                        l1 = __builtin_expect(src_line >= 1,       1) ? src_line-1 : src_line;
                        l2 = src_line;
                        l3 = __builtin_expect(src_line+1 < in_h,   1) ? src_line+1 : src_line;
                        l4 = __builtin_expect(src_line+2 < in_h,   1) ? src_line+2 : l3;
                    } else {
                        l1 = __builtin_expect(src_line+1 < in_h,   1) ? src_line+1 : src_line;
                        l2 = src_line;
                        l3 = __builtin_expect(src_line >= 1,        1) ? src_line-1 : src_line;
                        l4 = __builtin_expect(src_line >= 2,        1) ? src_line-2 : l3;
                    }
                    const unsigned char *r1 = tmp_plane + l1 * out_w * bps_out;
                    const unsigned char *r2 = tmp_plane + l2 * out_w * bps_out;
                    const unsigned char *r3 = tmp_plane + l3 * out_w * bps_out;
                    const unsigned char *r4 = tmp_plane + l4 * out_w * bps_out;
                    for(int i = 0; i < out_w; i++) {
                        int v = (RD_TMP(r1,i)*kv0 + RD_TMP(r2,i)*kv1
                               + RD_TMP(r3,i)*kv2 + RD_TMP(r4,i)*kv3) >> WEIGHT_SHIFT;
                        WR(dst_row, i, v);
                    }
                }
            }
            else if(weighted && is_bilinear == 1 && src_line+1 < in_h)
            {
                const int w = map->h_weight[map_h_base + px];
                if(w == 0) {
                    memcpy(dst_row, tmp_row0, out_w * bps_out);
                } else {
                    const int iw = WEIGHT_MAX - w;
                    const unsigned char *tmp_row1 = tmp_row0 + out_w * bps_out;
                    for(int i = 0; i < out_w; i++)
                        WR(dst_row, i, (RD_TMP(tmp_row0,i)*iw + RD_TMP(tmp_row1,i)*w) >> WEIGHT_SHIFT);
                }
            }
            else if(weighted && is_bilinear == -1 && src_line >= 1)
            {
                const int w = map->h_weight[map_h_base + px];
                if(w == 0) {
                    memcpy(dst_row, tmp_row0, out_w * bps_out);
                } else {
                    const int iw = WEIGHT_MAX - w;
                    const unsigned char *tmp_row_1 = tmp_row0 - out_w * bps_out;
                    for(int i = 0; i < out_w; i++)
                        WR(dst_row, i, (RD_TMP(tmp_row0,i)*iw + RD_TMP(tmp_row_1,i)*w) >> WEIGHT_SHIFT);
                }
            }
            else if(bilinear && is_bilinear == 1 && src_line+1 < in_h)
            {
                const unsigned char *tmp_row1 = tmp_row0 + out_w * bps_out;
                for(int i = 0; i < out_w; i++)
                    WR(dst_row, i, (RD_TMP(tmp_row0,i) + RD_TMP(tmp_row1,i)) >> 1);
            }
            else if(bilinear && is_bilinear == -1 && src_line >= 1)
            {
                const unsigned char *tmp_row_1 = tmp_row0 - out_w * bps_out;
                for(int i = 0; i < out_w; i++)
                    WR(dst_row, i, (RD_TMP(tmp_row_1,i) + RD_TMP(tmp_row0,i)) >> 1);
            }
            else
            {
                memcpy(dst_row, tmp_row0, out_w * bps_out);
            }
        }
        offset += out_w * out_h * bps_out;
    }
    #undef RD
	#undef RD_TMP
    #undef WR
}

// ============================================================
// process_files_2d_plane — version 2D, 1 canal planar
// ============================================================

__attribute__((optimize("O3,unroll-loops")))
void process_files_2d_plane(const unsigned char * restrict in_plane,
                             const unsigned char * restrict next_plane,
                             unsigned char * restrict blend_plane,
                             unsigned char * restrict out_plane,
                             const BresenhamMap * restrict map,
                             const Kernel2D * restrict k2d,
                             const Kernel2D * restrict k2d_soft,
                             const Kernel2D * restrict k2d_mid,
                             const int *tmp_weights,
                             const int in_w, const int in_h,
                             const int out_w, const int out_h,
                             const int ratio, const int mode,
                             const int bps_in, const int bps_out,
                             GradPixel * restrict grad,
                             int16_t * restrict dist_map,
                             int aniso_strength_fp)
{
    const int in_size       = in_w * in_h;
    const int use_2d_kernel = (mode >= 3);
    const int is_next_plane = (next_plane != in_plane);
 
    #define WR(x, v) do { \
        int _v = (v); \
        _v = _v < 0 ? 0 : _v > 65535 ? 65535 : _v; \
        if(bps_out == 1) dst_row_8[x] = (unsigned char)(_v >> 8); \
        else             dst16    [x] = (uint16_t)_v; \
    } while(0)
    #define RD(row, x) \
        ((bps_in == 1) ? ((int)(row##_8)[x] << 8) \
                       : (int)((row##_16)[x]))
 
	if(mode == 6 || mode == 7)
	{
		compute_gradient_plane(in_plane, grad, in_w, in_h, bps_in);
		compute_conf_map(grad, in_w, in_h, mode, aniso_strength_fp);
		compute_couture_dist(grad, dist_map, in_w, in_h);
	}
 
    unsigned int offset = 0;
    for(int r = 0; r < ratio; r++)
    {
        const int tw = tmp_weights[r];
        const unsigned char *src_plane;
        if(tw > 0 && is_next_plane)
        {
            const int itw = WEIGHT_MAX - tw;
            if(bps_in == 1) {
                for(int i = 0; i < in_size; i++)
                    blend_plane[i] = (in_plane[i]*itw + next_plane[i]*tw) >> WEIGHT_SHIFT;
            } else {
                const uint16_t *in16   = (const uint16_t *)in_plane;
                const uint16_t *next16 = (const uint16_t *)next_plane;
                      uint16_t *bl16   = (uint16_t *)blend_plane;
                for(int i = 0; i < in_size; i++)
                    bl16[i] = (in16[i]*itw + next16[i]*tw) >> WEIGHT_SHIFT;
            }
            src_plane = blend_plane;
        }
        else { src_plane = in_plane; }
 
        const int map_w_base = r * out_w;
        const int map_h_base = r * out_h;
 
        for(int y = 0; y < out_h; y++)
        {
            // --- résolution verticale (flag + weight) ---
            int sy, fy;
			if(use_2d_kernel) {
				const int hpos = map->h_pos[map_h_base + y];
				sy = hpos >> WEIGHT_SHIFT;
				fy = hpos & (WEIGHT_MAX - 1);
			} else {
				sy = map->h_src[map_h_base + y];
				fy = 0;
			}
 
            uint16_t      *dst16     = (uint16_t *)(out_plane + offset + y * out_w * bps_out);
            unsigned char *dst_row_8 =              out_plane + offset + y * out_w * bps_out;
            const unsigned char *row2_8  = src_plane + sy * in_w * bps_in;
            const uint16_t      *row2_16 = (const uint16_t *)row2_8;
 
            if(use_2d_kernel)
            {
                const Kernel2D *k_base = &k2d[FQ(fy) * NQUANT];
				const Kernel2D *k_base_soft = (mode == 7) ? &k2d_soft[FQ(fy) * NQUANT] : k_base;
				const Kernel2D *k_base_mid  = (mode == 7) ? &k2d_mid [FQ(fy)*NQUANT] : k_base;
 
                for(int x = 0; x < out_w; x++)
                {
                    // --- résolution horizontale (flag + weight) ---
                    const int wpos = map->w_pos[map_w_base + x];
					const int sx   = wpos >> WEIGHT_SHIFT;
					const int fx = wpos & (WEIGHT_MAX - 1);
 
                    if(mode == 6 || mode == 7)
                    {
                        // gradient interpolé bilinéairement
                        const int sx1 = sx + 1 < in_w ? sx + 1 : sx;
                        const int sy1 = sy + 1 < in_h ? sy + 1 : sy;
                        const GradPixel *g00 = &grad[sy  * in_w + sx ];
                        const GradPixel *g10 = &grad[sy  * in_w + sx1];
                        const GradPixel *g01 = &grad[sy1 * in_w + sx ];
                        const GradPixel *g11 = &grad[sy1 * in_w + sx1];
                        const int gx_t = (g00->gx*(WEIGHT_MAX-fx) + g10->gx*fx) >> WEIGHT_SHIFT;
                        const int gy_t = (g00->gy*(WEIGHT_MAX-fx) + g10->gy*fx) >> WEIGHT_SHIFT;
                        const int mg_t = (g00->mag*(WEIGHT_MAX-fx)+ g10->mag*fx) >> WEIGHT_SHIFT;
                        const int gx_b = (g01->gx*(WEIGHT_MAX-fx) + g11->gx*fx) >> WEIGHT_SHIFT;
                        const int gy_b = (g01->gy*(WEIGHT_MAX-fx) + g11->gy*fx) >> WEIGHT_SHIFT;
                        const int mg_b = (g01->mag*(WEIGHT_MAX-fx)+ g11->mag*fx) >> WEIGHT_SHIFT;
                        const int igx = (gx_t*(WEIGHT_MAX-fy) + gx_b*fy) >> WEIGHT_SHIFT;
                        const int igy = (gy_t*(WEIGHT_MAX-fy) + gy_b*fy) >> WEIGHT_SHIFT;
                        const int mag = (mg_t*(WEIGHT_MAX-fy) + mg_b*fy) >> WEIGHT_SHIFT;
 
                        int conf;
						if(mode == 7) {
							if(mag <= ANISO_LOW
							   || is_compression_artifact(grad, sx, sy, mag, in_w, in_h)) {
								conf = 0;                                   // doux
							}
							else if(mag >= ANISO_HIGH) {
								conf = WEIGHT_MAX;
							}
							else {
								conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
							}
						}
						else {  // mode 6
							if(mag <= ANISO_LOW)       conf = 0;
							else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
							else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
						}
						conf = (conf * aniso_strength_fp) >> WEIGHT_SHIFT;
 
                        int v;
                        if(conf == 0) {
							// distance à la couture, interpolée bilinéairement
							const int d00 = dist_map[sy  * in_w + sx ];
							const int d10 = dist_map[sy  * in_w + sx1];
							const int d01 = dist_map[sy1 * in_w + sx ];
							const int d11 = dist_map[sy1 * in_w + sx1];
							const int d_t = (d00*(WEIGHT_MAX-fx) + d10*fx) >> WEIGHT_SHIFT;
							const int d_b = (d01*(WEIGHT_MAX-fx) + d11*fx) >> WEIGHT_SHIFT;
							const int d_fp = (d_t*(WEIGHT_MAX-fy) + d_b*fy) >> WEIGHT_SHIFT;

							const int dmax = (COUTURE_N + 1) * WEIGHT_MAX;
							// alpha = part de "mid" : 1 à la frontière (d=0), 0 au fond du soft
							int alpha = (d_fp >= dmax) ? 0 : (WEIGHT_MAX - (d_fp * WEIGHT_MAX / dmax));

							// valeur soft (curve)
							int v_soft;
							{
								const Kernel2D *k = &k_base_soft[FQ(fx)];
								const int ntaps = (int)k->n;
								int acc = 0;
								if(sx>=2 && sx+3<in_w && sy>=2 && sy+3<in_h) {
									if(bps_in == 1) {
										const unsigned char *c = src_plane + sy*in_w + sx;
										for(int t=0;t<ntaps;t++) acc += ((int)c[k->ox[t]]<<8)*(int)k->w[t];
									} else {
										const uint16_t *c = (const uint16_t*)src_plane + sy*in_w + sx;
										for(int t=0;t<ntaps;t++) acc += (int)c[k->ox[t]]*(int)k->w[t];
									}
								} else {
									for(int t=0;t<ntaps;t++){
										int lx=sx+k->dx[t], ly=sy+k->dy[t];
										lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
										int pv=(bps_in==1)?((int)src_plane[ly*in_w+lx]<<8)
														  :(int)((const uint16_t*)src_plane)[ly*in_w+lx];
										acc += pv*(int)k->w[t];
									}
								}
								v_soft = acc >> WEIGHT_SHIFT;
							}

							if(alpha == 0) {
								v = v_soft;
							} else {
								// valeur mid (sinc r2 lobes réduits)
								int v_mid;
								{
									const Kernel2D *k = &k_base_mid[FQ(fx)];
									const int ntaps = (int)k->n;
									int acc = 0;
									if(sx>=2 && sx+3<in_w && sy>=2 && sy+3<in_h) {
										if(bps_in == 1) {
											const unsigned char *c = src_plane + sy*in_w + sx;
											for(int t=0;t<ntaps;t++) acc += ((int)c[k->ox[t]]<<8)*(int)k->w[t];
										} else {
											const uint16_t *c = (const uint16_t*)src_plane + sy*in_w + sx;
											for(int t=0;t<ntaps;t++) acc += (int)c[k->ox[t]]*(int)k->w[t];
										}
									} else {
										for(int t=0;t<ntaps;t++){
											int lx=sx+k->dx[t], ly=sy+k->dy[t];
											lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
											int pv=(bps_in==1)?((int)src_plane[ly*in_w+lx]<<8)
															  :(int)((const uint16_t*)src_plane)[ly*in_w+lx];
											acc += pv*(int)k->w[t];
										}
									}
									v_mid = acc >> WEIGHT_SHIFT;
								}
								v = (v_mid * alpha + v_soft * (WEIGHT_MAX - alpha)) >> WEIGHT_SHIFT;
							}
						} else {
                                v = aniso_pixel_plane(src_plane, sx, sy, fx, fy, igx, igy, conf,
													  in_w, in_h, bps_in);
                            if(v < 0) {
                                const Kernel2D *k = &k_base[FQ(fx)];
                                const int ntaps = (int)k->n;
                                int acc = 0;
                                for(int t=0;t<ntaps;t++){
                                    int lx=sx+k->dx[t], ly=sy+k->dy[t];
                                    lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
                                    int pv=(bps_in==1)?((int)src_plane[ly*in_w+lx]<<8)
                                                      :(int)((const uint16_t*)src_plane)[ly*in_w+lx];
                                    acc += pv*(int)k->w[t];
                                }
                                v = acc >> WEIGHT_SHIFT;
                            }
                        }
                        WR(x, v);
                    }
                    else  // modes 3/4/5
                    {
                        const Kernel2D *k = &k_base[FQ(fx)];
                        const int ntaps = (int)k->n;
                        int acc = 0;
                        if(sx>=2 && sx+3<in_w && sy>=2 && sy+3<in_h) {
                            if(bps_in == 1) {
                                const unsigned char *c = src_plane + sy*in_w + sx;
                                for(int t=0;t<ntaps;t++) acc += ((int)c[k->ox[t]]<<8)*(int)k->w[t];
                            } else {
                                const uint16_t *c = (const uint16_t*)src_plane + sy*in_w + sx;
                                for(int t=0;t<ntaps;t++) acc += (int)c[k->ox[t]]*(int)k->w[t];
                            }
                        } else {
                            for(int t=0;t<ntaps;t++){
                                int lx=sx+k->dx[t], ly=sy+k->dy[t];
                                lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
                                int pv=(bps_in==1)?((int)src_plane[ly*in_w+lx]<<8)
                                                  :(int)((const uint16_t*)src_plane)[ly*in_w+lx];
                                acc += pv*(int)k->w[t];
                            }
                        }
                        WR(x, acc >> WEIGHT_SHIFT);
                    }
                }
            }
            else
            {
                // modes 0/1/2 — inchangé (utilise déjà flag)
                const int hl = (mode >= 1) ? map->h_flag[map_h_base + y] : 0;
                int sy_v2 = sy;
                if(hl ==  1) sy_v2 = (sy + 1 < in_h) ? sy + 1 : sy;
                if(hl == -1) sy_v2 = (sy - 1 >= 0)   ? sy - 1 : sy;
                const unsigned char *row3_8  = src_plane + sy_v2 * in_w * bps_in;
                const uint16_t      *row3_16 = (const uint16_t *)row3_8;
 
                if(hl == 0) {
                    for(int x = 0; x < out_w; x++) {
                        const int sx = map->w_src [map_w_base + x];
                        const int wl = (mode >= 1) ? map->w_flag[map_w_base + x] : 0;
                        if(wl == 0) WR(x, RD(row2, sx));
                        else {
                            const int sx_b = (wl==1)?(sx+1<in_w?sx+1:sx):(sx-1>=0?sx-1:sx);
                            WR(x, (RD(row2,sx)+RD(row2,sx_b))>>1);
                        }
                    }
                } else {
                    for(int x = 0; x < out_w; x++) {
                        const int sx = map->w_src [map_w_base + x];
                        const int wl = (mode >= 1) ? map->w_flag[map_w_base + x] : 0;
                        const int sx_b = (wl!=0)?((wl==1)?(sx+1<in_w?sx+1:sx):(sx-1>=0?sx-1:sx)):sx;
                        if(wl == 0) WR(x, (RD(row2,sx)+RD(row3,sx))>>1);
                        else WR(x, (RD(row2,sx)+RD(row2,sx_b)+RD(row3,sx)+RD(row3,sx_b))>>2);
                    }
                }
            }
        }
        offset += out_w * out_h * bps_out;
    }
    #undef RD
    #undef WR
}

// ============================================================
// process_files_2d — version optimisée
// ============================================================

__attribute__((optimize("O3,unroll-loops")))
void process_files_2d(const unsigned char * restrict in_buf,
                      const unsigned char * restrict next_buf,
                      unsigned char * restrict blend_buf,
                      unsigned char * restrict out_buf,
                      unsigned char * restrict tmp_buf,
                      const BresenhamMap * restrict map,
                      const int *tmp_weights,
                      const int in_w, const int in_h,
                      const int out_w, const int out_h,
                      const int ratio, const int mode,
                      const int bps_in, const int bps_out,
                      GradPixel * restrict grad,
					  int16_t * restrict dist_map,
                      int aniso_strength_fp)
{
    const int in_w3         = in_w * 3;
    const int out_w3        = out_w * 3;
    const int in_size       = in_w * in_h * 3;
    const int bps3_in       = 3 * bps_in;
    const int use_2d_kernel = (mode >= 3);
    const int is_next_buf   = (next_buf != in_buf);
 
    #define RD3(ptr, i) \
        ((bps_in == 1) ? ((int)(ptr)[i] << 8) \
                       : (int)((const uint16_t*)(ptr))[i])
    #define WR3(ptr, i, v) do { \
        int _v = (v); \
        _v = _v < 0 ? 0 : _v > 65535 ? 65535 : _v; \
        if(bps_out == 1) (ptr)[i]      = (unsigned char)(_v >> 8); \
        else ((uint16_t*)(ptr))[i]     = (uint16_t)_v; \
    } while(0)
 
	if(mode == 6 || mode == 7)
	{
		compute_gradient_rgb(in_buf, grad, in_w, in_h, bps_in);
		compute_conf_map(grad, in_w, in_h, mode, aniso_strength_fp);
		compute_couture_dist(grad, dist_map, in_w, in_h);
	}
 
    unsigned int offset = 0;
    for(int r = 0; r < ratio; r++)
    {
        const int tw = tmp_weights[r];
        const unsigned char *src_buf_ptr;
        if(tw > 0 && is_next_buf)
        {
            const int itw = WEIGHT_MAX - tw;
            if(bps_in == 1) {
                for(int i = 0; i < in_size; i++)
                    blend_buf[i] = (in_buf[i]*itw + next_buf[i]*tw) >> WEIGHT_SHIFT;
            } else {
                const uint16_t *in16   = (const uint16_t *)in_buf;
                const uint16_t *next16 = (const uint16_t *)next_buf;
                      uint16_t *bl16   = (uint16_t *)blend_buf;
                for(int i = 0; i < in_size; i++)
                    bl16[i] = (in16[i]*itw + next16[i]*tw) >> WEIGHT_SHIFT;
            }
            src_buf_ptr = blend_buf;
        }
        else { src_buf_ptr = in_buf; }
 
        const int map_w_base = r * out_w;
        const int map_h_base = r * out_h;
        const Kernel2D *k2d = (mode >= 5) ? kaiser_k2d : curve_k2d;
 
        for(int y = 0; y < out_h; y++)
        {
            int sy, fy;
			if(use_2d_kernel) {
				const int hpos = map->h_pos[map_h_base + y];
				sy = hpos >> WEIGHT_SHIFT;
				fy = hpos & (WEIGHT_MAX - 1);
			} else {
				sy = map->h_src[map_h_base + y];
				fy = 0;
			}
 
            unsigned char *dst_row = out_buf + offset + y * out_w3 * bps_out;
            const unsigned char *row2 = src_buf_ptr + sy * in_w3 * bps_in;
 
            if(use_2d_kernel)
            {
                const Kernel2D *k_base = &k2d[FQ(fy) * NQUANT];
				const Kernel2D *k_base_soft = (mode == 7) ? &curve_k2d[FQ(fy) * NQUANT] : k_base;
				const Kernel2D *k_base_mid  = (mode == 7) ? &mid_k2d  [FQ(fy) * NQUANT] : k_base;
                int x3 = 0;
                for(int x = 0; x < out_w; x++, x3 += 3)
                {
                    const int wpos = map->w_pos[map_w_base + x];
					const int sx   = wpos >> WEIGHT_SHIFT;
					const int fx = wpos & (WEIGHT_MAX - 1);
 
                    if(mode == 6 || mode == 7)
                    {
                        const int sx1 = sx + 1 < in_w ? sx + 1 : sx;
                        const int sy1 = sy + 1 < in_h ? sy + 1 : sy;
                        const GradPixel *g00 = &grad[sy  * in_w + sx ];
                        const GradPixel *g10 = &grad[sy  * in_w + sx1];
                        const GradPixel *g01 = &grad[sy1 * in_w + sx ];
                        const GradPixel *g11 = &grad[sy1 * in_w + sx1];
                        const int gx_t = (g00->gx*(WEIGHT_MAX-fx) + g10->gx*fx) >> WEIGHT_SHIFT;
                        const int gy_t = (g00->gy*(WEIGHT_MAX-fx) + g10->gy*fx) >> WEIGHT_SHIFT;
                        const int mg_t = (g00->mag*(WEIGHT_MAX-fx)+ g10->mag*fx) >> WEIGHT_SHIFT;
                        const int gx_b = (g01->gx*(WEIGHT_MAX-fx) + g11->gx*fx) >> WEIGHT_SHIFT;
                        const int gy_b = (g01->gy*(WEIGHT_MAX-fx) + g11->gy*fx) >> WEIGHT_SHIFT;
                        const int mg_b = (g01->mag*(WEIGHT_MAX-fx)+ g11->mag*fx) >> WEIGHT_SHIFT;
                        const int igx = (gx_t*(WEIGHT_MAX-fy) + gx_b*fy) >> WEIGHT_SHIFT;
                        const int igy = (gy_t*(WEIGHT_MAX-fy) + gy_b*fy) >> WEIGHT_SHIFT;
                        const int mag = (mg_t*(WEIGHT_MAX-fy) + mg_b*fy) >> WEIGHT_SHIFT;
 
                        int conf;
						if(mode == 7) {
							if(mag <= ANISO_LOW
							   || is_compression_artifact(grad, sx, sy, mag, in_w, in_h)) {
								conf = 0;                                   // doux
							}
							else if(mag >= ANISO_HIGH) {
								conf = WEIGHT_MAX;
							}
							else {
								conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
							}
						}
						else {  // mode 6
							if(mag <= ANISO_LOW)       conf = 0;
							else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
							else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
						}
						conf = (conf * aniso_strength_fp) >> WEIGHT_SHIFT;
 
                        int v0, v1, v2;
                        if(conf == 0) {
							const int d00 = dist_map[sy  * in_w + sx ];
							const int d10 = dist_map[sy  * in_w + sx1];
							const int d01 = dist_map[sy1 * in_w + sx ];
							const int d11 = dist_map[sy1 * in_w + sx1];
							const int d_t = (d00*(WEIGHT_MAX-fx) + d10*fx) >> WEIGHT_SHIFT;
							const int d_b = (d01*(WEIGHT_MAX-fx) + d11*fx) >> WEIGHT_SHIFT;
							const int d_fp = (d_t*(WEIGHT_MAX-fy) + d_b*fy) >> WEIGHT_SHIFT;

							const int dmax = (COUTURE_N + 1) * WEIGHT_MAX;
							int alpha = (d_fp >= dmax) ? 0 : (WEIGHT_MAX - (d_fp * WEIGHT_MAX / dmax));

							// --- soft (curve), 3 canaux ---
							int s0,s1,s2;
							{
								const Kernel2D *k = &k_base_soft[FQ(fx)];
								const int ntaps = (int)k->n;
								int a0=0,a1=0,a2=0;
								if(sx>=2 && sx+3<in_w && sy>=2 && sy+3<in_h) {
									const unsigned char *c = src_buf_ptr + sy*in_w3 + sx*bps3_in;
									for(int t=0;t<ntaps;t++){
										const unsigned char *p=c+k->ox[t]*bps3_in; int wt=(int)k->w[t];
										a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
									}
								} else {
									for(int t=0;t<ntaps;t++){
										int lx=sx+k->dx[t],ly=sy+k->dy[t];
										lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
										const unsigned char *p=src_buf_ptr+(ly*in_w+lx)*bps3_in; int wt=(int)k->w[t];
										a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
									}
								}
								s0=a0>>WEIGHT_SHIFT; s1=a1>>WEIGHT_SHIFT; s2=a2>>WEIGHT_SHIFT;
							}

							if(alpha == 0) {
								v0=s0; v1=s1; v2=s2;
							} else {
								// --- mid (sinc r2 lobes réduits), 3 canaux ---
								int m0,m1,m2;
								{
									const Kernel2D *k = &k_base_mid[FQ(fx)];
									const int ntaps = (int)k->n;
									int a0=0,a1=0,a2=0;
									if(sx>=2 && sx+3<in_w && sy>=2 && sy+3<in_h) {
										const unsigned char *c = src_buf_ptr + sy*in_w3 + sx*bps3_in;
										for(int t=0;t<ntaps;t++){
											const unsigned char *p=c+k->ox[t]*bps3_in; int wt=(int)k->w[t];
											a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
										}
									} else {
										for(int t=0;t<ntaps;t++){
											int lx=sx+k->dx[t],ly=sy+k->dy[t];
											lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
											const unsigned char *p=src_buf_ptr+(ly*in_w+lx)*bps3_in; int wt=(int)k->w[t];
											a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
										}
									}
									m0=a0>>WEIGHT_SHIFT; m1=a1>>WEIGHT_SHIFT; m2=a2>>WEIGHT_SHIFT;
								}
								const int ia = WEIGHT_MAX - alpha;
								v0 = (m0*alpha + s0*ia) >> WEIGHT_SHIFT;
								v1 = (m1*alpha + s1*ia) >> WEIGHT_SHIFT;
								v2 = (m2*alpha + s2*ia) >> WEIGHT_SHIFT;
							}
						} else {
                                aniso_pixel_rgb(src_buf_ptr, sx, sy, fx, fy, igx, igy, conf,
												in_w, in_h, bps_in, &v0, &v1, &v2);
                            if(v0 < 0) {
                                const Kernel2D *k = &k_base[FQ(fx)];
                                const int ntaps = (int)k->n;
                                int a0=0,a1=0,a2=0;
                                for(int t=0;t<ntaps;t++){
                                    int lx=sx+k->dx[t],ly=sy+k->dy[t];
                                    lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
                                    const unsigned char *p=src_buf_ptr+(ly*in_w+lx)*bps3_in; int wt=(int)k->w[t];
                                    a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
                                }
                                v0=a0>>WEIGHT_SHIFT; v1=a1>>WEIGHT_SHIFT; v2=a2>>WEIGHT_SHIFT;
                            }
                        }
                        WR3(dst_row, x3,   v0);
                        WR3(dst_row, x3+1, v1);
                        WR3(dst_row, x3+2, v2);
                    }
                    else  // modes 3/4/5
                    {
                        const Kernel2D *k = &k_base[FQ(fx)];
                        const int ntaps = (int)k->n;
                        int a0=0,a1=0,a2=0;
                        if(sx>=2 && sx+3<in_w && sy>=2 && sy+3<in_h) {
                            const unsigned char *c = src_buf_ptr + sy*in_w3 + sx*bps3_in;
                            for(int t=0;t<ntaps;t++){
                                const unsigned char *p=c+k->ox[t]*bps3_in; int wt=(int)k->w[t];
                                a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
                            }
                        } else {
                            for(int t=0;t<ntaps;t++){
                                int lx=sx+k->dx[t],ly=sy+k->dy[t];
                                lx=lx<0?0:lx>=in_w?in_w-1:lx; ly=ly<0?0:ly>=in_h?in_h-1:ly;
                                const unsigned char *p=src_buf_ptr+(ly*in_w+lx)*bps3_in; int wt=(int)k->w[t];
                                a0+=RD3(p,0)*wt; a1+=RD3(p,1)*wt; a2+=RD3(p,2)*wt;
                            }
                        }
                        WR3(dst_row,x3,  a0>>WEIGHT_SHIFT);
                        WR3(dst_row,x3+1,a1>>WEIGHT_SHIFT);
                        WR3(dst_row,x3+2,a2>>WEIGHT_SHIFT);
                    }
                }
            }
            else
            {
                // modes 0/1/2 — inchangé
                const int hl   = (mode >= 1) ? map->h_flag[map_h_base + y] : 0;
                int sy_v2 = sy;
                if(hl ==  1) sy_v2 = (sy + 1 < in_h) ? sy + 1 : sy;
                if(hl == -1) sy_v2 = (sy - 1 >= 0)   ? sy - 1 : sy;
                const unsigned char *row3 = src_buf_ptr + sy_v2 * in_w3 * bps_in;
 
                if(hl == 0) {
                    int x3 = 0;
                    for(int x = 0; x < out_w; x++, x3 += 3) {
                        const int sx  = map->w_src [map_w_base + x];
                        const int wl  = (mode >= 1) ? map->w_flag[map_w_base + x] : 0;
                        const int sx3 = sx * 3;
                        if(wl == 0) {
                            WR3(dst_row,x3,RD3(row2,sx3)); WR3(dst_row,x3+1,RD3(row2,sx3+1)); WR3(dst_row,x3+2,RD3(row2,sx3+2));
                        } else {
                            const int sx3b=(wl==1)?(sx3+3<in_w3?sx3+3:sx3):(sx3-3>=0?sx3-3:sx3);
                            for(int c=0;c<3;c++) WR3(dst_row,x3+c,(RD3(row2,sx3+c)+RD3(row2,sx3b+c))>>1);
                        }
                    }
                } else {
                    int x3 = 0;
                    for(int x = 0; x < out_w; x++, x3 += 3) {
                        const int sx  = map->w_src [map_w_base + x];
                        const int wl  = (mode >= 1) ? map->w_flag[map_w_base + x] : 0;
                        const int sx3 = sx * 3;
                        if(wl == 0) {
                            for(int c=0;c<3;c++) WR3(dst_row,x3+c,(RD3(row2,sx3+c)+RD3(row3,sx3+c))>>1);
                        } else {
                            const int sx3b=(wl==1)?(sx3+3<in_w3?sx3+3:sx3):(sx3-3>=0?sx3-3:sx3);
                            for(int c=0;c<3;c++) WR3(dst_row,x3+c,(RD3(row2,sx3+c)+RD3(row2,sx3b+c)+RD3(row3,sx3+c)+RD3(row3,sx3b+c))>>2);
                        }
                    }
                }
            }
        }
        offset += out_w * out_h * 3 * bps_out;
    }
    #undef RD3
    #undef WR3
}

void process_files(const unsigned char * restrict in_buf,
                   const unsigned char * restrict next_buf,
                   unsigned char * restrict blend_buf,
                   unsigned char * restrict out_buf,
                   unsigned char * restrict tmp_buf,
                   const BresenhamMap * restrict map,
                   const int *tmp_weights,
                   const int in_w, const int in_h,
                   const int out_w, const int out_h,
                   const int ratio, const int mode,
                   const int bps_in, const int bps_out)
{
    const int in_w3    = in_w * 3;   // en pixels, stride en canaux
    const int out_w3   = out_w * 3;
    const int in_size  = in_w * in_h * 3;
    const int bilinear = (mode >= 1 && mode != 5);
    const int weighted = (mode == 3 || mode == 4);
	 
		#define RD3(ptr, i) \
			((bps_in == 1) ? ((int)(ptr)[i] << 8) \
						   : (int)((const uint16_t*)(ptr))[i])
	 
		#define WR3(ptr, i, v) do { \
			int _v = (v); \
			_v = _v < 0 ? 0 : _v > 65535 ? 65535 : _v; \
			if(bps_out == 1) (ptr)[i]               = (unsigned char)(_v >> 8); \
			else ((uint16_t*)(ptr))[i]              = (uint16_t)_v; \
		} while(0)
 
    unsigned int offset = 0;
    for(int r = 0; r < ratio; r++)
    {
        const int tw = tmp_weights[r];
        const unsigned char *src_buf_ptr;
        if(tw > 0 && next_buf != in_buf)
        {
            const int itw = WEIGHT_MAX - tw;
            if(bps_in == 1) {
                for(int i = 0; i < in_size; i++)
                    blend_buf[i] = (in_buf[i]*itw + next_buf[i]*tw) >> WEIGHT_SHIFT;
            } else {
                const uint16_t *in16   = (const uint16_t *)in_buf;
                const uint16_t *next16 = (const uint16_t *)next_buf;
                      uint16_t *bl16   = (uint16_t *)blend_buf;
                for(int i = 0; i < in_size; i++)
                    bl16[i] = (in16[i]*itw + next16[i]*tw) >> WEIGHT_SHIFT;
            }
            src_buf_ptr = blend_buf;
        }
        else
        {
            src_buf_ptr = in_buf;
        }
 
        const int map_w_base = r * out_w;
 
        for(int h = 0; h < in_h; h++)
        {
            const unsigned char *src_row = src_buf_ptr + h * in_w3 * bps_in;
            unsigned char       *dst_row = tmp_buf     + h * out_w3 * bps_out;
 
            int px3 = 0;
            for(int px = 0; px < out_w; px++, px3 += 3)
            {
                const int src_px      = map->w_src [map_w_base + px] * 3;
                const int is_bilinear = map->w_flag[map_w_base + px];
 
                if(mode == 5 && is_bilinear != 0)
                {
                    const int w = map->w_weight[map_w_base + px];
                    if(w == 0) {
                        WR3(dst_row, px3,   RD3(src_row, src_px));
                        WR3(dst_row, px3+1, RD3(src_row, src_px+1));
                        WR3(dst_row, px3+2, RD3(src_row, src_px+2));
                    } else {
                        const int *k = kaiser4_lut[w];
                        int p1, p2, p3, p4;
                        if(is_bilinear == 1) {
                            p1 = (src_px >= 3)         ? src_px-3 : src_px;
                            p2 = src_px;
                            p3 = (src_px+3 < in_w3)    ? src_px+3 : src_px;
                            p4 = (src_px+6 < in_w3)    ? src_px+6 : p3;
                        } else {
                            p1 = (src_px+3 < in_w3)    ? src_px+3 : src_px;
                            p2 = src_px;
                            p3 = (src_px >= 3)          ? src_px-3 : src_px;
                            p4 = (src_px >= 6)          ? src_px-6 : p3;
                        }
                        WR3(dst_row, px3,   (RD3(src_row,p1)*k[0]+RD3(src_row,p2)*k[1]+RD3(src_row,p3)*k[2]+RD3(src_row,p4)*k[3]) >> WEIGHT_SHIFT);
                        WR3(dst_row, px3+1, (RD3(src_row,p1+1)*k[0]+RD3(src_row,p2+1)*k[1]+RD3(src_row,p3+1)*k[2]+RD3(src_row,p4+1)*k[3]) >> WEIGHT_SHIFT);
                        WR3(dst_row, px3+2, (RD3(src_row,p1+2)*k[0]+RD3(src_row,p2+2)*k[1]+RD3(src_row,p3+2)*k[2]+RD3(src_row,p4+2)*k[3]) >> WEIGHT_SHIFT);
                    }
                }
                else if(weighted && is_bilinear == 1 && src_px+3 < in_w3)
                {
                    const int w = map->w_weight[map_w_base + px];
                    const int iw = WEIGHT_MAX - w;
                    WR3(dst_row, px3,   (RD3(src_row,src_px)  *iw + RD3(src_row,src_px+3)*w) >> WEIGHT_SHIFT);
                    WR3(dst_row, px3+1, (RD3(src_row,src_px+1)*iw + RD3(src_row,src_px+4)*w) >> WEIGHT_SHIFT);
                    WR3(dst_row, px3+2, (RD3(src_row,src_px+2)*iw + RD3(src_row,src_px+5)*w) >> WEIGHT_SHIFT);
                }
                else if(weighted && is_bilinear == -1 && src_px >= 3)
                {
                    const int w = map->w_weight[map_w_base + px];
                    const int iw = WEIGHT_MAX - w;
                    WR3(dst_row, px3,   (RD3(src_row,src_px)  *iw + RD3(src_row,src_px-3)*w) >> WEIGHT_SHIFT);
                    WR3(dst_row, px3+1, (RD3(src_row,src_px+1)*iw + RD3(src_row,src_px-2)*w) >> WEIGHT_SHIFT);
                    WR3(dst_row, px3+2, (RD3(src_row,src_px+2)*iw + RD3(src_row,src_px-1)*w) >> WEIGHT_SHIFT);
                }
                else if(bilinear && !weighted && is_bilinear == 1 && src_px+3 < in_w3)
                {
                    WR3(dst_row, px3,   (RD3(src_row,src_px)   + RD3(src_row,src_px+3)) >> 1);
                    WR3(dst_row, px3+1, (RD3(src_row,src_px+1) + RD3(src_row,src_px+4)) >> 1);
                    WR3(dst_row, px3+2, (RD3(src_row,src_px+2) + RD3(src_row,src_px+5)) >> 1);
                }
                else if(bilinear && !weighted && is_bilinear == -1 && src_px >= 3)
                {
                    WR3(dst_row, px3,   (RD3(src_row,src_px-3) + RD3(src_row,src_px))   >> 1);
                    WR3(dst_row, px3+1, (RD3(src_row,src_px-2) + RD3(src_row,src_px+1)) >> 1);
                    WR3(dst_row, px3+2, (RD3(src_row,src_px-1) + RD3(src_row,src_px+2)) >> 1);
                }
                else
                {
                    WR3(dst_row, px3,   RD3(src_row, src_px));
                    WR3(dst_row, px3+1, RD3(src_row, src_px+1));
                    WR3(dst_row, px3+2, RD3(src_row, src_px+2));
                }
            }
        }
 
        const int map_h_base = r * out_h;
        for(int px = 0; px < out_h; px++)
        {
            const int src_line    = map->h_src [map_h_base + px];
            const int is_bilinear = map->h_flag[map_h_base + px];
 
            unsigned char       *dst_row  = out_buf + offset + px * out_w3 * bps_out;
            const unsigned char *tmp_row0 = tmp_buf + src_line * out_w3 * bps_out;
            const int row_bytes = out_w3 * bps_out;
 
            if(mode == 5 && is_bilinear != 0)
            {
                const int w = map->h_weight[map_h_base + px];
                if(w == 0) {
                    memcpy(dst_row, tmp_row0, row_bytes);
                } else {
                    const int *k = kaiser4_lut[w];
                    int l1, l2, l3, l4;
                    if(is_bilinear == 1) {
                        l1 = (src_line >= 1)      ? src_line-1 : src_line;
                        l2 = src_line;
                        l3 = (src_line+1 < in_h)  ? src_line+1 : src_line;
                        l4 = (src_line+2 < in_h)  ? src_line+2 : l3;
                    } else {
                        l1 = (src_line+1 < in_h)  ? src_line+1 : src_line;
                        l2 = src_line;
                        l3 = (src_line >= 1)       ? src_line-1 : src_line;
                        l4 = (src_line >= 2)       ? src_line-2 : l3;
                    }
                    const unsigned char *r1 = tmp_buf + l1 * row_bytes;
                    const unsigned char *r2 = tmp_buf + l2 * row_bytes;
                    const unsigned char *r3 = tmp_buf + l3 * row_bytes;
                    const unsigned char *r4 = tmp_buf + l4 * row_bytes;
                    for(int i = 0; i < out_w3; i++) {
                        int v = (RD3(r1,i)*k[0]+RD3(r2,i)*k[1]+RD3(r3,i)*k[2]+RD3(r4,i)*k[3]) >> WEIGHT_SHIFT;
                        WR3(dst_row, i, v);
                    }
                }
            }
            else if(weighted && is_bilinear == 1 && src_line+1 < in_h)
            {
                const int w = map->h_weight[map_h_base + px];
                if(w == 0) {
                    memcpy(dst_row, tmp_row0, row_bytes);
                } else {
                    const int iw = WEIGHT_MAX - w;
                    const unsigned char *tmp_row1 = tmp_row0 + row_bytes;
                    for(int i = 0; i < out_w3; i++)
                        WR3(dst_row, i, (RD3(tmp_row0,i)*iw + RD3(tmp_row1,i)*w) >> WEIGHT_SHIFT);
                }
            }
            else if(weighted && is_bilinear == -1 && src_line >= 1)
            {
                const int w = map->h_weight[map_h_base + px];
                if(w == 0) {
                    memcpy(dst_row, tmp_row0, row_bytes);
                } else {
                    const int iw = WEIGHT_MAX - w;
                    const unsigned char *tmp_row_1 = tmp_row0 - row_bytes;
                    for(int i = 0; i < out_w3; i++)
                        WR3(dst_row, i, (RD3(tmp_row0,i)*iw + RD3(tmp_row_1,i)*w) >> WEIGHT_SHIFT);
                }
            }
            else if(bilinear && !weighted && is_bilinear == 1 && src_line+1 < in_h)
            {
                const unsigned char *tmp_row1 = tmp_row0 + row_bytes;
                for(int i = 0; i < out_w3; i++)
                    WR3(dst_row, i, (RD3(tmp_row0,i) + RD3(tmp_row1,i)) >> 1);
            }
            else if(bilinear && !weighted && is_bilinear == -1 && src_line >= 1)
            {
                const unsigned char *tmp_row_1 = tmp_row0 - row_bytes;
                for(int i = 0; i < out_w3; i++)
                    WR3(dst_row, i, (RD3(tmp_row_1,i) + RD3(tmp_row0,i)) >> 1);
            }
            else
            {
                memcpy(dst_row, tmp_row0, row_bytes);
            }
        }
        offset += out_w * out_h * 3 * bps_out;
    }
    #undef RD3
    #undef WR3
}

// Layout mémoire YUV444p dans les buffers existants :
//   in_buf  : [Y plane : in_w*in_h] [U plane : in_w*in_h] [V plane : in_w*in_h]
//   out_buf : [Y plane : out_w*out_h*ratio] [U plane] [V plane]
//   tmp_buf : [Y plane : out_w*in_h] [U plane] [V plane]
//   blend_buf: [Y plane : in_w*in_h] [U plane] [V plane]
 
void *worker_loop(void *arg) {
    Worker *w = (Worker *)arg;
    while (1) {
        pthread_mutex_lock(&w->mutex);
        while (w->ready == 0)
            pthread_cond_wait(&w->cond_work, &w->mutex);
        int cmd = w->ready;
        w->ready = 0;
        pthread_mutex_unlock(&w->mutex);

        if (cmd == -1) break;

        ThreadArgs *a = &w->args;

        const int uv_in_w  = (a->in_pix_fmt  == PIX_444P || a->in_pix_fmt  == PIX_444P16) ? a->in_w    :
							 (a->in_pix_fmt  == PIX_411P  || a->in_pix_fmt  == PIX_410P)  ? a->in_w/4  : a->in_w/2;
		const int uv_in_h  = (a->in_pix_fmt  == PIX_420P  || a->in_pix_fmt  == PIX_420P16) ? a->in_h/2 :
							 (a->in_pix_fmt  == PIX_410P)                                  ? a->in_h/4 : a->in_h;
		const int uv_out_w = (a->out_pix_fmt == PIX_444P || a->out_pix_fmt == PIX_444P16) ? a->out_w   :
							 (a->out_pix_fmt == PIX_411P  || a->out_pix_fmt == PIX_410P)  ? a->out_w/4 : a->out_w/2;
		const int uv_out_h = (a->out_pix_fmt == PIX_420P  || a->out_pix_fmt == PIX_420P16) ? a->out_h/2 :
							 (a->out_pix_fmt == PIX_410P)                                  ? a->out_h/4 : a->out_h;

		const int in_y_size   = a->in_w  * a->in_h  * a->bps_in;
		const int out_y_size  = a->out_w * a->out_h * a->ratio * a->bps_out;
		const int tmp_y_size  = a->out_w * a->in_h  * a->bps_out;
		const int in_uv_size  = uv_in_w  * uv_in_h  * a->bps_in;
		const int out_uv_size = uv_out_w * uv_out_h * a->ratio * a->bps_out;
		const int tmp_uv_size = uv_out_w * uv_in_h  * a->bps_out;

        switch(a->mode_1d) {

            case 1:  // pipeline 1D séparable
                switch(a->in_pix_fmt) {
                    case PIX_RGB:
					case PIX_RGB48:
                        process_files(
                            a->in_buf, a->next_buf, a->blend_buf,
                            a->out_buf, a->tmp_buf,
                            a->map, a->tmp_weights,
                            a->in_w, a->in_h, a->out_w, a->out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
                        break;

                    case PIX_444P:
					case PIX_444P16:
					case PIX_422P:
					case PIX_422P16:
					case PIX_420P:
					case PIX_420P16:
					case PIX_411P:
					case PIX_410P:
                        // Y — toujours pleine résolution
                        process_files_plane(
                            a->in_buf,
                            a->next_buf,
                            a->blend_buf,
                            a->out_buf,
                            a->tmp_buf,
                            a->map, a->tmp_weights,
                            a->in_w, a->in_h, a->out_w, a->out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
                        // U
                        process_files_plane(
                            a->in_buf    + in_y_size,
                            a->next_buf  + in_y_size,
                            a->blend_buf + in_y_size,
                            a->out_buf   + out_y_size,
                            a->tmp_buf   + tmp_y_size,
                            a->map_uv, a->tmp_weights,
                            uv_in_w, uv_in_h, uv_out_w, uv_out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
                        // V
                        process_files_plane(
                            a->in_buf    + in_y_size + in_uv_size,
                            a->next_buf  + in_y_size + in_uv_size,
                            a->blend_buf + in_y_size + in_uv_size,
                            a->out_buf   + out_y_size + out_uv_size,
                            a->tmp_buf   + tmp_y_size + tmp_uv_size,
                            a->map_uv, a->tmp_weights,
                            uv_in_w, uv_in_h, uv_out_w, uv_out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
                        break;

                    default:
                        break;
                }
                break;

            default:  // pipeline 2D direct
                switch(a->in_pix_fmt) {
                    case PIX_RGB:
					case PIX_RGB48:
						process_files_2d(
							a->in_buf, a->next_buf, a->blend_buf,
							a->out_buf, a->tmp_buf,
							a->map, a->tmp_weights,
							a->in_w, a->in_h, a->out_w, a->out_h,
							a->ratio, a->mode, a->bps_in, a->bps_out,
							a->grad, a->dist_map,
							a->aniso_strength_fp);
						if(a->blur_artifact) {
							const int phase = compute_grid_phase(a->in_buf, a->in_w, a->in_h, a->bps_in, 1);
							post_deblock_rgb(a->out_buf, a->map, a->out_w, a->out_h, a->ratio, a->bps_out, phase);
						}
						break;

                    case PIX_444P:
					case PIX_444P16:
					case PIX_422P:
					case PIX_422P16:
					case PIX_420P:
					case PIX_420P16:
					case PIX_411P:
					case PIX_410P:
                        // Y
						process_files_2d_plane(
							a->in_buf, a->next_buf, a->blend_buf, a->out_buf,
							a->map,
							(a->mode >= 5) ? kaiser_k2d : curve_k2d,
							curve_k2d,
							mid_k2d,
							a->tmp_weights,
							a->in_w, a->in_h, a->out_w, a->out_h,
							a->ratio, a->mode, a->bps_in, a->bps_out,
							a->grad, a->dist_map,
							a->aniso_strength_fp);

						// U
						process_files_2d_plane(
							a->in_buf    + in_y_size,
							a->next_buf  + in_y_size,
							a->blend_buf + in_y_size,
							a->out_buf   + out_y_size,
							a->map_uv,
							(a->mode >= 5) ? kaiser_k2d_uv : curve_k2d_uv,
							curve_k2d_uv,
							mid_k2d_uv,
							a->tmp_weights,
							uv_in_w, uv_in_h, uv_out_w, uv_out_h,
							a->ratio, a->mode, a->bps_in, a->bps_out,
							a->grad, a->dist_map,
							a->aniso_strength_fp);

						// V — même table que U
						process_files_2d_plane(
							a->in_buf    + in_y_size + in_uv_size,
							a->next_buf  + in_y_size + in_uv_size,
							a->blend_buf + in_y_size + in_uv_size,
							a->out_buf   + out_y_size + out_uv_size,
							a->map_uv,
							(a->mode >= 5) ? kaiser_k2d_uv : curve_k2d_uv,
							curve_k2d_uv,
							mid_k2d_uv,
							a->tmp_weights,
							uv_in_w, uv_in_h, uv_out_w, uv_out_h,
							a->ratio, a->mode, a->bps_in, a->bps_out,
							a->grad, a->dist_map,
							a->aniso_strength_fp);
						
						//luma only
						if(a->blur_artifact) {
						    int phase = compute_grid_phase(a->in_buf, a->in_w, a->in_h, a->bps_in, 0);
						    post_deblock_plane(a->out_buf, a->map, a->out_w, a->out_h, a->ratio, a->bps_out, phase);
						}

                        break;

                    default:
                        break;
                }
                break;
        }

        pthread_mutex_lock(&w->mutex);
        w->done = 1;
        pthread_cond_signal(&w->cond_done);
        pthread_mutex_unlock(&w->mutex);
    }
    return NULL;
}


int main(int argc, char **argv)
{
#ifdef _WIN32 || _WIN64
	_setmode(_fileno(stdout), O_BINARY);
	_setmode(_fileno(stdin), O_BINARY);	
#endif

	/*setvbuf(stdin,  NULL, _IOFBF, 64 * 1024 * 1024);  // 1MB
	setvbuf(stdout, NULL, _IOFBF, 128 * 1024 * 1024);  // 4MB*/

	int opt;

	//file adress
	char *input_name = NULL;
	
	//input file
	FILE *input;
	
	//default parameters
	int in_w = 720;
	int in_h = 480;
	int out_w = 1920;
	int out_h = 1080;
	int frames_ratio = 2;
	int mode = 0;
	int tmp_mode = 0;         // mode de blending temporel (0 = duplication)
	double curve_base = 2.0;  // base de la courbe exponentielle (modes 4 / tmp 3)
	double kaiser_beta = 2.5;  // défaut agressif, netteté maximale
	int mode_1d = 0;  // défaut : 2D direct
	int *tmp_weights = NULL;  // poids temporel par phase (taille = frames_ratio)
	double aniso_strength = 1.0;  // défaut : anisotropie pleine
	int blur_artifact = 0;
	
	PixFmt in_pix_fmt  = PIX_RGB;
	PixFmt out_pix_fmt = PIX_RGB;
	
	//buffer : in_buf a 1 slot supplémentaire pour le look-ahead temporel.
	//in_buf[NB_THREADS] est la frame "next" du dernier worker, qui devient
	//in_buf[0] (current du worker 0) du batch suivant.
	unsigned char *in_buf[NB_THREADS + 1];
	unsigned char *out_buf[NB_THREADS];
	unsigned char *tmp_buf[NB_THREADS];
	unsigned char *blend_buf[NB_THREADS];  // buffer du lerp temporel par worker
	GradPixel *grad_buf[NB_THREADS];
	int16_t   *dist_buf[NB_THREADS];
		
	for(int t=0; t<NB_THREADS; t++)
	{
		in_buf[t]    = NULL;
		out_buf[t]   = NULL;
		tmp_buf[t]   = NULL;
		blend_buf[t] = NULL;
		grad_buf[t]  = NULL;
		dist_buf[t]  = NULL;
	}
	in_buf[NB_THREADS] = NULL;
	
	//buffer length
	unsigned int in_buf_length = 0;
	unsigned int out_buf_length = 0;
	
	BresenhamMap map;
	map.w_src    = NULL;
	map.w_flag   = NULL;
	map.w_weight = NULL;
	map.w_heavy  = NULL;
	map.w_pos    = NULL;
	map.h_src    = NULL;
	map.h_flag   = NULL;
	map.h_weight = NULL;
	map.h_heavy  = NULL;
	map.h_pos    = NULL;
	map.w_length = 0;
	map.h_length = 0;
	
	BresenhamMap map_uv;
	map_uv.w_src    = NULL;
	map_uv.w_flag   = NULL;
	map_uv.w_weight = NULL;
	map_uv.w_heavy  = NULL;
	map_uv.w_pos    = NULL;
	map_uv.h_src    = NULL;
	map_uv.h_flag   = NULL;
	map_uv.h_weight = NULL;
	map_uv.h_heavy  = NULL;
	map_uv.h_pos    = NULL;
	map_uv.w_length = 0;
	map_uv.h_length = 0;
	
	int option_index = 0;
	static struct option long_options[] = {
		{"in-width",   1, 0, 1},
		{"in-height",  1, 0, 2},
		{"out-width",  1, 0, 3},
		{"out-height", 1, 0, 4},
		{"fps-ratio",  1, 0, 5},
		{"mode",       1, 0, 6},
		{"curve-base", 1, 0, 7},
		{"tmp-mode",   1, 0, 8},
		{"kaiser-beta", 1, 0, 9},
		{"1d", 0, 0, 10},
		{"in-pix-fmt",  1, 0, 11},
		{"out-pix-fmt", 1, 0, 12},
		{"aniso-strength", 1, 0, 13},
		{"blur-artifact", 0, 0, 14},
		{0, 0, 0, 0}//reminder : letter value are from 65 to 122
	};

	while ((opt = getopt_long_only(argc, argv, "i:",long_options, &option_index)) != -1) {
		switch (opt) {
		case 'i':
			input_name = optarg;
			break;
		case 1:
			in_w = atoi(optarg);
			break;
		case 2:
			in_h = atoi(optarg);
			break;
		case 3:
			out_w = atoi(optarg);
			break;
		case 4:
			out_h = atoi(optarg);
			break;
		case 5:
			frames_ratio = atoi(optarg);
			break;
		case 6:
			mode = atoi(optarg);
			break;
		case 7:
			curve_base = atof(optarg);
			break;
		case 8:
			tmp_mode = atoi(optarg);
			break;
		case 9:
			kaiser_beta = atof(optarg);
			break;
		case 10:
			mode_1d = 1;
			break;
		case 11:
			if     (strcmp(optarg, "rgb24") == 0 || strcmp(optarg, "RGB24") == 0 || strcmp(optarg, "rgb") == 0 || strcmp(optarg, "RGB") == 0) in_pix_fmt = PIX_RGB;
			else if(strcmp(optarg, "444p") == 0 || strcmp(optarg, "444P") == 0 || strcmp(optarg, "444") == 0) in_pix_fmt = PIX_444P;
			else if(strcmp(optarg, "422p") == 0 || strcmp(optarg, "422P") == 0 || strcmp(optarg, "422") == 0) in_pix_fmt = PIX_422P;
			else if(strcmp(optarg, "420p") == 0 || strcmp(optarg, "420P") == 0 || strcmp(optarg, "420") == 0) in_pix_fmt = PIX_420P;
			else if(strcmp(optarg, "411p") == 0 || strcmp(optarg, "411P") == 0 || strcmp(optarg, "411") == 0) in_pix_fmt = PIX_411P;
			else if(strcmp(optarg, "410p") == 0 || strcmp(optarg, "410P") == 0 || strcmp(optarg, "410") == 0) in_pix_fmt = PIX_410P;
			else if(strcmp(optarg, "rgb48") == 0|| strcmp(optarg, "RGB48") == 0) in_pix_fmt = PIX_RGB48;
			else if(strcmp(optarg, "444p16") == 0 || strcmp(optarg, "444P16") == 0) in_pix_fmt = PIX_444P16;
			else if(strcmp(optarg, "422p16") == 0 || strcmp(optarg, "422P16") == 0) in_pix_fmt = PIX_422P16;
			else if(strcmp(optarg, "420p16") == 0 || strcmp(optarg, "420P16") == 0) in_pix_fmt = PIX_420P16;
			else { fprintf(stderr, "Error: unknown in-pix-fmt '%s'\n", optarg); return -1; }
			break;
		case 12:
			if     (strcmp(optarg, "rgb24") == 0 || strcmp(optarg, "RGB24") == 0 || strcmp(optarg, "rgb") == 0 || strcmp(optarg, "RGB") == 0) out_pix_fmt = PIX_RGB;
			else if(strcmp(optarg, "444p") == 0 || strcmp(optarg, "444P") == 0 || strcmp(optarg, "444") == 0) out_pix_fmt = PIX_444P;
			else if(strcmp(optarg, "422p") == 0 || strcmp(optarg, "422P") == 0 || strcmp(optarg, "422") == 0) out_pix_fmt = PIX_422P;
			else if(strcmp(optarg, "420p") == 0 || strcmp(optarg, "420P") == 0 || strcmp(optarg, "420") == 0) out_pix_fmt = PIX_420P;
			else if(strcmp(optarg, "411p") == 0 || strcmp(optarg, "411P") == 0 || strcmp(optarg, "411") == 0) out_pix_fmt = PIX_411P;
			else if(strcmp(optarg, "410p") == 0 || strcmp(optarg, "410P") == 0 || strcmp(optarg, "410") == 0) out_pix_fmt = PIX_410P;
			else if(strcmp(optarg, "rgb48") == 0 || strcmp(optarg, "RGB48") == 0) out_pix_fmt = PIX_RGB48;
			else if(strcmp(optarg, "444p16") == 0 || strcmp(optarg, "444P16") == 0) out_pix_fmt = PIX_444P16;
			else if(strcmp(optarg, "422p16") == 0 || strcmp(optarg, "422P16") == 0) out_pix_fmt = PIX_422P16;
			else if(strcmp(optarg, "420p16") == 0 || strcmp(optarg, "420P16") == 0) out_pix_fmt = PIX_420P16;
			else { fprintf(stderr, "Error: unknown out-pix-fmt '%s'\n", optarg); return -1; }
			break;
		case 13:
			aniso_strength = atof(optarg);
			break;
		case 14:
			blur_artifact = 1;
			break;
		default:
			usage();
			break;
		}
	}
	
	if(input_name == NULL)
	{
		fprintf(stderr, "no input file provided\n");
		return -1;
	}
	
	if(mode < 0 || mode > 7)
	{
		fprintf(stderr, "Error : mode '%d' is not supported\n", mode);
		return -1;
	}
	
	if(mode >= 6 && mode_1d)
	{
		fprintf(stderr, "Error: mode 6 requires 2D pipeline, --1d not supported\n");
		return -1;
	}
	
	if(blur_artifact && (mode < 3 || mode_1d)) {
		fprintf(stderr, "Error: --blur-artifact requires mode >= 3 and the 2D pipeline (not --1d)\n");
		return -1;
	}
	
	if(tmp_mode < 0 || tmp_mode > 3)
	{
		fprintf(stderr, "Error : tmp-mode '%d' is not supported\n", tmp_mode);
		return -1;
	}
	
	if((mode == 4 || tmp_mode == 3) && curve_base <= 1.0)
	{
		fprintf(stderr, "Error : --curve-base must be > 1.0 (got %f)\n", curve_base);
		return -1;
	}

	if(mode == 3 && curve_base != 2.0)
	{
		fprintf(stderr, "Warning : --curve-base is ignored in mode 3 (linear curve)\n");
	}
	
	//reading file
	if (strcmp(input_name, "-") == 0)// Read samples from stdin
	{
		input = stdin;
	}
	else
	{
		input = fopen(input_name, "rb");
		if (!input) {
			fprintf(stderr, "Failed to open %s\n", input_name);
			return -ENOENT;
		}
	}
	
	if((in_pix_fmt == PIX_RGB || in_pix_fmt == PIX_RGB48) != (out_pix_fmt == PIX_RGB || out_pix_fmt == PIX_RGB48))
	{
		fprintf(stderr, "Error: RGB<->YUV conversion not supported\n");
		return -1;
	}
	
	if(in_pix_fmt == PIX_422P || in_pix_fmt == PIX_420P)
	{
		if(in_w % 2 != 0)
		{
			fprintf(stderr, "Error: in-width %d must be even for 422/420\n", in_w);
			return -1;
		}
	}
	if(in_pix_fmt == PIX_420P)
	{
		if(in_h % 2 != 0)
		{
			fprintf(stderr, "Error: in-height %d must be even for 420\n", in_h);
			return -1;
		}
	}
	if(out_pix_fmt == PIX_422P || out_pix_fmt == PIX_420P)
	{
		if(out_w % 2 != 0)
		{
			fprintf(stderr, "Error: out-width %d must be even for 422/420\n", out_w);
			return -1;
		}
	}
	if(out_pix_fmt == PIX_420P)
	{
		if(out_h % 2 != 0)
		{
			fprintf(stderr, "Error: out-height %d must be even for 420\n", out_h);
			return -1;
		}
	}
	
	if(in_pix_fmt == PIX_411P || in_pix_fmt == PIX_410P)
	{
		if(in_w % 4 != 0)
		{
			fprintf(stderr, "Error: in-width must be multiple of 4 for 411/410\n"); return -1;
		}
	}
	if(in_pix_fmt == PIX_410P)
	{
		if(in_h % 4 != 0)
		{
			fprintf(stderr, "Error: in-height must be multiple of 4 for 410\n"); return -1;
		}
	}
	
	if(out_pix_fmt == PIX_411P || out_pix_fmt == PIX_410P)
	{
		if(out_w % 4 != 0)
		{
			fprintf(stderr, "Error: in-width must be multiple of 4 for 411/410\n"); return -1;
		}
	}
	if(out_pix_fmt == PIX_410P)
	{
		if(out_h % 4 != 0)
		{
			fprintf(stderr, "Error: in-height must be multiple of 4 for 410\n"); return -1;
		}
	}
	
	const int is_rgb_in  = (in_pix_fmt  == PIX_RGB  || in_pix_fmt  == PIX_RGB48);
	const int is_rgb_out = (out_pix_fmt == PIX_RGB  || out_pix_fmt == PIX_RGB48);
	const int is_444_in  = (in_pix_fmt  == PIX_444P || in_pix_fmt  == PIX_444P16);
	const int is_444_out = (out_pix_fmt == PIX_444P || out_pix_fmt == PIX_444P16);
	 
	const int need_map_uv = !is_rgb_in && !is_rgb_out && !(is_444_in && is_444_out);
	
	const int bps_in  = (in_pix_fmt  >= PIX_RGB48) ? 2 : 1;
	const int bps_out = (out_pix_fmt >= PIX_RGB48) ? 2 : 1;
	
	if(out_h > in_h && out_w > in_w && in_h > 0 && in_w > 0)
	{
		switch(in_pix_fmt) {
			case PIX_RGB:	 in_buf_length = in_w * in_h * 3;              break;
			case PIX_RGB48:	 in_buf_length = in_w * in_h * 3 * 2;          break;
			case PIX_444P:
			case PIX_444P16: in_buf_length = in_w * in_h * 3 * bps_in;     break;
			case PIX_422P:
			case PIX_422P16: in_buf_length = in_w * in_h * 2 * bps_in;     break;
			case PIX_420P:
			case PIX_420P16: in_buf_length = in_w * in_h * 3 / 2 * bps_in; break;
			case PIX_411P:	 in_buf_length = in_w * in_h * 3 / 2;          break;
			case PIX_410P:   in_buf_length = in_w * in_h * 9 / 8;          break;
		}
		switch(out_pix_fmt) {
			case PIX_RGB:	 out_buf_length = out_w * out_h * 3     * frames_ratio;           break;
			case PIX_RGB48:  out_buf_length = out_w * out_h * 3 * 2 * frames_ratio;           break;
			case PIX_444P:
			case PIX_444P16: out_buf_length = out_w * out_h * 3 * bps_out * frames_ratio;     break;
			case PIX_422P:
			case PIX_422P16: out_buf_length = out_w * out_h * 2 * bps_out * frames_ratio;     break;
			case PIX_420P:
			case PIX_420P16: out_buf_length = out_w * out_h * 3 / 2 * bps_out * frames_ratio; break;
			case PIX_411P:	 out_buf_length = out_w * out_h * 3 / 2 * frames_ratio;           break;
			case PIX_410P:   out_buf_length = out_w * out_h * 9 / 8 * frames_ratio;           break;
		}
		
		map.w_length = out_w * frames_ratio;
		map.h_length = out_h * frames_ratio;
		
		map.w_src    = ALIGNED_MALLOC(map.w_length * sizeof(int), 64);
		map.w_flag   = ALIGNED_MALLOC(map.w_length * sizeof(int), 64);
		map.w_weight = ALIGNED_MALLOC(map.w_length * sizeof(int), 64);
		map.w_heavy  = ALIGNED_MALLOC(map.w_length * sizeof(int), 64);
		map.w_pos    = ALIGNED_MALLOC(map.w_length * sizeof(int), 64);
		map.h_src    = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_flag   = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_weight = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_heavy  = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_pos    = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);

		if(need_map_uv)
		{
			    const int uv_in_w  = (in_pix_fmt  == PIX_444P  || in_pix_fmt  == PIX_444P16) ? in_w    :
									 (in_pix_fmt  == PIX_411P  || in_pix_fmt  == PIX_410P)   ? in_w/4  : in_w/2;
				const int uv_in_h  = (in_pix_fmt  == PIX_420P  || in_pix_fmt  == PIX_420P16) ? in_h/2  :
									 (in_pix_fmt  == PIX_410P)                               ? in_h/4  : in_h;
				const int uv_out_w = (out_pix_fmt == PIX_444P  || out_pix_fmt == PIX_444P16) ? out_w   :
									 (out_pix_fmt == PIX_411P  || out_pix_fmt == PIX_410P)   ? out_w/4 : out_w/2;
				const int uv_out_h = (out_pix_fmt == PIX_420P  || out_pix_fmt == PIX_420P16) ? out_h/2 :
									 (out_pix_fmt == PIX_410P)                               ? out_h/4 : out_h;

			map_uv.w_length = uv_out_w * frames_ratio;
			map_uv.h_length = uv_out_h * frames_ratio;

			map_uv.w_src    = ALIGNED_MALLOC(map_uv.w_length * sizeof(int), 64);
			map_uv.w_flag   = ALIGNED_MALLOC(map_uv.w_length * sizeof(int), 64);
			map_uv.w_weight = ALIGNED_MALLOC(map_uv.w_length * sizeof(int), 64);
			map_uv.w_heavy  = ALIGNED_MALLOC(map_uv.w_length * sizeof(int), 64);
			map_uv.w_pos    = ALIGNED_MALLOC(map_uv.w_length * sizeof(int), 64);
			map_uv.h_src    = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_flag   = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_weight = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_heavy  = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_pos    = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);

			if(!map_uv.w_src || !map_uv.w_flag || !map_uv.w_weight || !map_uv.w_heavy ||
			   !map_uv.h_src || !map_uv.h_flag || !map_uv.h_weight || !map_uv.h_heavy)
			{
				ALIGNED_FREE(map_uv.w_src);   ALIGNED_FREE(map_uv.w_flag);
				ALIGNED_FREE(map_uv.w_weight);ALIGNED_FREE(map_uv.w_heavy);
				ALIGNED_FREE(map_uv.w_pos);
				ALIGNED_FREE(map_uv.h_src);   ALIGNED_FREE(map_uv.h_flag);
				ALIGNED_FREE(map_uv.h_weight);ALIGNED_FREE(map_uv.h_heavy);
				ALIGNED_FREE(map_uv.h_pos);
				fprintf(stderr, "ALIGNED_MALLOC error (map_uv)\n");
				return -1;
			}
		}
		
		int alloc_err = 0;
		for(int t=0; t<NB_THREADS; t++)
		{
			in_buf[t]    = ALIGNED_MALLOC(in_buf_length, 64);
			out_buf[t]   = ALIGNED_MALLOC(out_buf_length, 64);
			grad_buf[t]  = ALIGNED_MALLOC(in_w * in_h * sizeof(GradPixel), 64);
			dist_buf[t]  = ALIGNED_MALLOC(in_w * in_h * sizeof(int16_t), 64);
			
			const int uv_out_w_alloc = (out_pix_fmt == PIX_444P || out_pix_fmt == PIX_444P16) ? out_w  :
									   (out_pix_fmt == PIX_411P || out_pix_fmt == PIX_410P)   ? out_w/4 : out_w/2;
			const int uv_in_h_alloc  = (in_pix_fmt  == PIX_420P || in_pix_fmt  == PIX_420P16) ? in_h/2 :
									   (in_pix_fmt  == PIX_410P)                              ? in_h/4 : in_h;
		 
			switch(in_pix_fmt) {
				case PIX_RGB:
				case PIX_RGB48:
					tmp_buf[t] = ALIGNED_MALLOC(out_w * in_h * 3 * bps_out, 64);
					break;
				default:  // tous les formats planar
					tmp_buf[t] = ALIGNED_MALLOC((out_w * in_h + uv_out_w_alloc * uv_in_h_alloc * 2)	* bps_out, 64);
					break;
			}

			blend_buf[t] = ALIGNED_MALLOC(in_buf_length, 64);
			if(in_buf[t] == NULL || out_buf[t] == NULL || tmp_buf[t] == NULL || blend_buf[t] == NULL || grad_buf[t] == NULL || dist_buf[t] == NULL)
				alloc_err = 1;
		}
		// Slot supplémentaire pour le look-ahead temporel
		in_buf[NB_THREADS] = ALIGNED_MALLOC(in_buf_length, 64);
		if(in_buf[NB_THREADS] == NULL) alloc_err = 1;

		// Table des poids temporels (1 entry par phase)
		tmp_weights = ALIGNED_MALLOC(frames_ratio * sizeof(int), 64);
		if(tmp_weights == NULL) alloc_err = 1;
		
		if(!map.w_src || !map.w_flag || !map.w_weight || !map.w_heavy ||
		   !map.h_src || !map.h_flag || !map.h_weight || !map.h_heavy || alloc_err)//check allocation error
		{
			ALIGNED_FREE(map.w_src);
			ALIGNED_FREE(map.w_flag);
			ALIGNED_FREE(map.w_weight);
			ALIGNED_FREE(map.w_heavy);
			ALIGNED_FREE(map.w_pos);
			ALIGNED_FREE(map.h_src);
			ALIGNED_FREE(map.h_flag);
			ALIGNED_FREE(map.h_weight);
			ALIGNED_FREE(map.h_heavy);
			ALIGNED_FREE(tmp_weights);
			ALIGNED_FREE(map.h_pos);			
			
			for(int t=0; t<NB_THREADS; t++)
			{
				ALIGNED_FREE(in_buf[t]);
				ALIGNED_FREE(out_buf[t]);
				ALIGNED_FREE(tmp_buf[t]);
				ALIGNED_FREE(blend_buf[t]);
				ALIGNED_FREE(grad_buf[t]);
				ALIGNED_FREE(dist_buf[t]);
			}
			ALIGNED_FREE(in_buf[NB_THREADS]);
			fprintf(stderr, "ALIGNED_MALLOC error (buf)\n");
			return -1;
		}
	}
	else//error
	{
		fprintf(stderr, "invalid size error\n");
		return -1;
	}
	
	//compute the maping used for scaling the input image into an output one
	compute_scale_map(&map, in_w, in_h, out_w, out_h, frames_ratio);
	
	compute_smoothstep_lut();
	
	for(int r = 0; r < frames_ratio; r++) {
		fprintf(stderr, "phase %d w_src : ", r);
		for(int px = 0; px < 16; px++)
			fprintf(stderr, "%d ", map.w_src[r * out_w + px]);
		fprintf(stderr, "\n");

		fprintf(stderr, "phase %d w_flag: ", r);
		for(int px = 0; px < 16; px++)
			fprintf(stderr, "%d ", map.w_flag[r * out_w + px]);
		fprintf(stderr, "\n");
	}
	
	const int uv_in_w  = (in_pix_fmt  == PIX_444P  || in_pix_fmt  == PIX_444P16) ? in_w    :
						 (in_pix_fmt  == PIX_411P  || in_pix_fmt  == PIX_410P)   ? in_w/4  : in_w/2;
	const int uv_in_h  = (in_pix_fmt  == PIX_420P  || in_pix_fmt  == PIX_420P16) ? in_h/2  :
						 (in_pix_fmt  == PIX_410P)                               ? in_h/4  : in_h;
	const int uv_out_w = (out_pix_fmt == PIX_444P  || out_pix_fmt == PIX_444P16) ? out_w   :
						 (out_pix_fmt == PIX_411P  || out_pix_fmt == PIX_410P)   ? out_w/4 : out_w/2;
	const int uv_out_h = (out_pix_fmt == PIX_420P  || out_pix_fmt == PIX_420P16) ? out_h/2 :
						 (out_pix_fmt == PIX_410P)                               ? out_h/4 : out_h;
	
	if(mode == 3)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_LINEAR, curve_base);
		adjust_phase_direction_weighted(&map, in_w, in_h, out_w, out_h, frames_ratio);
		compute_curve_rad_lut(3, curve_base);
		build_k2d(curve_k2d, curve_rad_lut, CURVE_RAD_SIZE, in_w);
		if(need_map_uv)
			build_k2d(curve_k2d_uv, curve_rad_lut, CURVE_RAD_SIZE, uv_in_w);
		}
	else if(mode == 4)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_EXPONENTIAL, curve_base);
		adjust_phase_direction_weighted(&map, in_w, in_h, out_w, out_h, frames_ratio);
		compute_curve_rad_lut(4, curve_base);
		build_k2d(curve_k2d, curve_rad_lut, CURVE_RAD_SIZE, in_w);
		if(need_map_uv)
			build_k2d(curve_k2d_uv, curve_rad_lut, CURVE_RAD_SIZE, uv_in_w);
		}
	else if(mode == 2)
	{
		adjust_phase_direction(&map, in_w, in_h, out_w, out_h, frames_ratio);
		for(int r = 0; r < frames_ratio; r++) {
			int cw = 0, ch = 0;
			for(int px = 0; px < out_w; px++)
				if(map.w_flag[r * out_w + px] != 0) cw++;
			for(int px = 0; px < out_h; px++)
				if(map.h_flag[r * out_h + px] != 0) ch++;
			fprintf(stderr, "phase %d: w=%d/%d (%.1f%%)  h=%d/%d (%.1f%%)\n",
					r, cw, out_w, 100.0*cw/out_w, ch, out_h, 100.0*ch/out_h);
		}
	}
	else if(mode == 5 || mode == 6 || mode == 7)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio,
							  CURVE_EXPONENTIAL, curve_base);
		adjust_phase_direction_weighted(&map, in_w, in_h, out_w, out_h, frames_ratio);
		compute_kaiser4_lut(kaiser_beta);
		compute_kaiser_rad_lut(kaiser_beta, curve_base);
		build_k2d(kaiser_k2d, kaiser_rad_lut, KAISER_RAD_SIZE, in_w);
		if(need_map_uv)
			build_k2d(kaiser_k2d_uv, kaiser_rad_lut, KAISER_RAD_SIZE, uv_in_w);
	 
		if(mode == 7) {
			// table douce : courbe exponentielle radius=1 (mode 4 style),
			// cardinale, sans lobes → utilisée quand conf==0
			compute_curve_rad_lut(4, curve_base);
			compute_mid_rad_lut(kaiser_beta);
			build_k2d(curve_k2d, curve_rad_lut, CURVE_RAD_SIZE, in_w);
			build_k2d(mid_k2d, mid_rad_lut, MID_RAD_SIZE, in_w);
			if(need_map_uv)
			{
				build_k2d(curve_k2d_uv, curve_rad_lut, CURVE_RAD_SIZE, uv_in_w);
				build_k2d(mid_k2d_uv, mid_rad_lut, MID_RAD_SIZE, uv_in_w);
			}
		}
	}
	
	if(need_map_uv)
	{
		compute_scale_map(&map_uv, uv_in_w, uv_in_h, uv_out_w, uv_out_h, frames_ratio);

		if(mode == 3) {
			compute_curve_weights(&map_uv, uv_out_w, uv_out_h, frames_ratio, CURVE_LINEAR, curve_base);
			adjust_phase_direction_weighted(&map_uv, uv_in_w, uv_in_h, uv_out_w, uv_out_h, frames_ratio);
		} else if(mode == 4 || mode == 5 || mode == 6 || mode == 7) {
			compute_curve_weights(&map_uv, uv_out_w, uv_out_h, frames_ratio, CURVE_EXPONENTIAL, curve_base);
			adjust_phase_direction_weighted(&map_uv, uv_in_w, uv_in_h, uv_out_w, uv_out_h, frames_ratio);
		} else if(mode == 2) {
			adjust_phase_direction_weighted(&map_uv, uv_in_w, uv_in_h, uv_out_w, uv_out_h, frames_ratio);
		}
	}

	// Pré-calcul des poids temporels (une fois pour toute la vidéo)
	compute_temporal_weights(tmp_weights, frames_ratio, (TmpMode)tmp_mode, curve_base);
	
	ThreadArgs args[NB_THREADS];
	for(int t=0; t<NB_THREADS; t++)
	{
		args[t].map         = &map;
		args[t].tmp_weights = tmp_weights;
		args[t].in_w        = in_w;
		args[t].in_h        = in_h;
		args[t].out_w       = out_w;
		args[t].out_h       = out_h;
		args[t].ratio       = frames_ratio;
		args[t].mode        = mode;
		args[t].mode_1d     = mode_1d;
		args[t].in_pix_fmt  = in_pix_fmt;
		args[t].out_pix_fmt = out_pix_fmt;
		args[t].bps_in  = (in_pix_fmt  >= PIX_RGB48) ? 2 : 1;
		args[t].bps_out = (out_pix_fmt >= PIX_RGB48) ? 2 : 1;
		args[t].map_uv = need_map_uv ? &map_uv : &map;
		args[t].grad        = grad_buf[t];
		args[t].dist_map    = dist_buf[t];
		args[t].blur_artifact = blur_artifact;
		args[t].aniso_strength_fp = (int)(aniso_strength * WEIGHT_MAX + 0.5);
	}
	
	// --- init thread pool ---
	Worker workers[NB_THREADS];
	for (int t = 0; t < NB_THREADS; t++) {
		workers[t].ready = 0;
		workers[t].done  = 0;
		pthread_mutex_init(&workers[t].mutex, NULL);
		pthread_cond_init(&workers[t].cond_work, NULL);
		pthread_cond_init(&workers[t].cond_done, NULL);
		pthread_create(&workers[t].tid, NULL, worker_loop, &workers[t]);
	}

	// Pré-lecture de la frame 0 dans in_buf[0]. Si vide, on saute la boucle.
	int has_current = (fread(in_buf[0], in_buf_length, 1, input) == 1);

	while (has_current)
	{
		// On a déjà 1 frame valide dans in_buf[0]. Lire jusqu'à NB_THREADS
		// frames supplémentaires dans in_buf[1..NB_THREADS]. Le dernier slot
		// (NB_THREADS) sert de "next" pour le worker NB_THREADS-1 et sera
		// recyclé comme "current" du worker 0 du batch suivant.
		int frames_loaded = 1;  // in_buf[0] déjà chargé
		for (int t = 1; t <= NB_THREADS; t++) {
			if (fread(in_buf[t], in_buf_length, 1, input) == 1)
				frames_loaded++;
			else
				break;
		}

		// Nombre de workers à dispatcher : on a frames_loaded frames valides,
		// on en traite jusqu'à NB_THREADS. Si frames_loaded == NB_THREADS+1
		// le dernier worker a un vrai "next" ; sinon le dernier worker
		// utilise son current comme next (blend trivial = current).
		int frames_to_process = (frames_loaded > NB_THREADS) ? NB_THREADS : frames_loaded;
		int has_full_lookahead = (frames_loaded == NB_THREADS + 1);

		for (int t = 0; t < frames_to_process; t++) {
			workers[t].args           = args[t];
			workers[t].args.in_buf    = in_buf[t];
			workers[t].args.out_buf   = out_buf[t];
			workers[t].args.tmp_buf   = tmp_buf[t];
			workers[t].args.blend_buf = blend_buf[t];
			// "next" = frame suivante si dispo, sinon current (= no-op blend)
			if (t < frames_to_process - 1)
				workers[t].args.next_buf = in_buf[t + 1];
			else if (has_full_lookahead)
				workers[t].args.next_buf = in_buf[NB_THREADS];
			else
				workers[t].args.next_buf = in_buf[t];  // dernière frame du fichier
		}

		// signal aux workers
		for (int t = 0; t < frames_to_process; t++) {
			pthread_mutex_lock(&workers[t].mutex);
			workers[t].done  = 0;
			workers[t].ready = 1;
			pthread_cond_signal(&workers[t].cond_work);
			pthread_mutex_unlock(&workers[t].mutex);
		}

		// attente
		for (int t = 0; t < frames_to_process; t++) {
			pthread_mutex_lock(&workers[t].mutex);
			while (!workers[t].done)
				pthread_cond_wait(&workers[t].cond_done, &workers[t].mutex);
			pthread_mutex_unlock(&workers[t].mutex);
		}

		// écriture
		if (!isatty(STDOUT_FILENO)) {
			for (int t = 0; t < frames_to_process; t++)
			{
				if(out_pix_fmt == PIX_RGB || out_pix_fmt == PIX_RGB48)
				{
					fwrite(out_buf[t], out_buf_length, 1, stdout);
				}
				else
				{
					const int frame_y_size = out_w * out_h * bps_out;
					const int uv_out_w_f = (out_pix_fmt == PIX_444P  || out_pix_fmt == PIX_444P16) ? out_w   :
										   (out_pix_fmt == PIX_411P  || out_pix_fmt == PIX_410P)   ? out_w/4 : out_w/2;
					const int uv_out_h_f = (out_pix_fmt == PIX_420P  || out_pix_fmt == PIX_420P16) ? out_h/2 :
										   (out_pix_fmt == PIX_410P)                               ? out_h/4 : out_h;
					const int frame_uv_size = uv_out_w_f * uv_out_h_f * bps_out;
					const int out_y_plane   = frame_y_size  * frames_ratio;
					const int out_uv_plane  = frame_uv_size * frames_ratio;
				 
					for(int f = 0; f < frames_ratio; f++)
					{
						fwrite(out_buf[t] + f * frame_y_size,                               frame_y_size,  1, stdout);
						fwrite(out_buf[t] + out_y_plane + f * frame_uv_size,                frame_uv_size, 1, stdout);
						fwrite(out_buf[t] + out_y_plane + out_uv_plane + f * frame_uv_size, frame_uv_size, 1, stdout);
					}
				}
			}
			fflush(stdout);
		}

		// Préparer le batch suivant : si on avait un look-ahead complet, la
		// frame in_buf[NB_THREADS] devient le current du worker 0 du prochain
		// batch (swap des pointeurs pour éviter la copie).
		if (has_full_lookahead) {
			unsigned char *swap = in_buf[0];
			in_buf[0] = in_buf[NB_THREADS];
			in_buf[NB_THREADS] = swap;
			has_current = 1;
		} else {
			has_current = 0;  // EOF déjà atteint, on sort
		}
	}

	// --- arrêt thread pool ---
	for (int t = 0; t < NB_THREADS; t++) {
		pthread_mutex_lock(&workers[t].mutex);
		workers[t].ready = -1;
		pthread_cond_signal(&workers[t].cond_work);
		pthread_mutex_unlock(&workers[t].mutex);
		pthread_join(workers[t].tid, NULL);
		pthread_mutex_destroy(&workers[t].mutex);
		pthread_cond_destroy(&workers[t].cond_work);
		pthread_cond_destroy(&workers[t].cond_done);
	}

////ending of the program
	
	ALIGNED_FREE(map.w_src);
	ALIGNED_FREE(map.w_flag);
	ALIGNED_FREE(map.w_weight);
	ALIGNED_FREE(map.w_heavy);
	ALIGNED_FREE(map.w_pos);
	ALIGNED_FREE(map.h_src);
	ALIGNED_FREE(map.h_flag);
	ALIGNED_FREE(map.h_weight);
	ALIGNED_FREE(map.h_heavy);
	ALIGNED_FREE(tmp_weights);
	ALIGNED_FREE(map.h_pos);
	
	if(need_map_uv)
	{
		ALIGNED_FREE(map_uv.w_src);
		ALIGNED_FREE(map_uv.w_flag);
		ALIGNED_FREE(map_uv.w_weight);
		ALIGNED_FREE(map_uv.w_heavy);
		ALIGNED_FREE(map_uv.w_pos);
		ALIGNED_FREE(map_uv.h_src);
		ALIGNED_FREE(map_uv.h_flag);
		ALIGNED_FREE(map_uv.h_weight);
		ALIGNED_FREE(map_uv.h_heavy);
		ALIGNED_FREE(map_uv.h_pos);
	}

	
	for(int t=0; t<NB_THREADS; t++)
	{
		ALIGNED_FREE(in_buf[t]);
		ALIGNED_FREE(out_buf[t]);
		ALIGNED_FREE(tmp_buf[t]);
		ALIGNED_FREE(blend_buf[t]);
		ALIGNED_FREE(grad_buf[t]);
		ALIGNED_FREE(dist_buf[t]);
	}
	ALIGNED_FREE(in_buf[NB_THREADS]);
	
	//Close file
	if (input && (input != stdin))
	{
		fclose(input);
	}

	return 0;
}