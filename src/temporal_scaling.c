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
    int *h_src;
    int *h_flag;
    int *h_weight;
    int *h_heavy;
    int w_length;
    int h_length;
};

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
	int mode_1d;
	PixFmt in_pix_fmt;
	PixFmt out_pix_fmt;
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
        "\t                     3 = linear weighted bilinear with temporal alternation\n"
        "\t                     4 = exponential weighted bilinear with temporal alternation\n"
		"\t                     5 = kaiser 4-tap avec poids exponentiels\n"
        "\t--tmp-mode INT     temporal blending mode (default: 0)\n"
        "\t                     0 = duplicate frames (no blending)\n"
        "\t                     1 = blend 50%% on the upper minority of phases\n"
        "\t                     2 = linear weighted blending\n"
        "\t                     3 = exponential weighted blending (uses --curve-base)\n"
        "\t--curve-base FLOAT exponential curve base for mode 4 / tmp-mode 3 (default: 2.0, must be > 1.0)\n\n"
		"\t--kaiser-beta FLOAT  beta du noyau kaiser pour mode 5 (default: 2.5)\n"
		"\t                     low = sharp/ringing, high = smooth/no ringing\n"
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
            map->h_flag[dst] = (2 * err >= out_h) ? 1 : 0;  // h_flag pas w_flag
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

void process_files_2d_plane(const unsigned char * restrict in_plane,
                             const unsigned char * restrict next_plane,
                             unsigned char * restrict blend_plane,
                             unsigned char * restrict out_plane,
                             const BresenhamMap * restrict map,
                             const int *tmp_weights,
                             const int in_w, const int in_h,
                             const int out_w, const int out_h,
                             const int ratio, const int mode,
                             const int bps_in, const int bps_out)
{
    const int in_size  = in_w * in_h;
    const int weighted = (mode == 3 || mode == 4 || mode == 5);
    const int mode5    = (mode == 5);
		 
	// Lecture : normalise vers l'espace de calcul interne (toujours 0..65535)
	// pour que les interpolations soient cohérentes quelle que soit la profondeur source
	#define RD(row, x) \
		((bps_in == 1) ? ((int)(row##_8)[x] << 8) \
					   : (int)((row##_16)[x]))
		 
	// Écriture : dénormalise depuis l'espace 16-bit vers la profondeur cible
	#define WR(x, v) do { \
		int _v = (v); \
		_v = _v < 0 ? 0 : _v > 65535 ? 65535 : _v; \
		if(bps_out == 1) dst_row_8[x]  = (unsigned char)(_v >> 8); \
		else             dst16    [x]  = (uint16_t)_v; \
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
 
        for(int y = 0; y < out_h; y++)
        {
            const int sy = map->h_src   [map_h_base + y];
            const int hl = (mode >= 1) ? map->h_flag  [map_h_base + y] : 0;
            const int wy = (mode >= 1) ? map->h_weight[map_h_base + y] : 0;
 
            // pointeurs typés pour les lignes verticales
            uint16_t       *dst16 = (uint16_t *)(out_plane + offset + y * out_w * bps_out);
 
            int l1, l2, l3, l4;
            if(hl == 1) {
                l1 = __builtin_expect(sy >= 1,       1) ? sy - 1 : sy;
                l2 = sy;
                l3 = __builtin_expect(sy + 1 < in_h, 1) ? sy + 1 : sy;
                l4 = __builtin_expect(sy + 2 < in_h, 1) ? sy + 2 : l3;
            } else if(hl == -1) {
                l1 = __builtin_expect(sy + 1 < in_h, 1) ? sy + 1 : sy;
                l2 = sy;
                l3 = __builtin_expect(sy >= 1,       1) ? sy - 1 : sy;
                l4 = __builtin_expect(sy >= 2,       1) ? sy - 2 : l3;
            } else {
                l1 = l2 = l3 = l4 = sy;
            }
 
            // pointeurs lignes — un seul calcul, utilisés pour 8 ET 16-bit
            const unsigned char *row1_8  = src_plane + l1 * in_w * bps_in;
            const unsigned char *row2_8  = src_plane + l2 * in_w * bps_in;
            const unsigned char *row3_8  = src_plane + l3 * in_w * bps_in;
            const unsigned char *row4_8  = src_plane + l4 * in_w * bps_in;
            const uint16_t      *row1_16 = (const uint16_t *)row1_8;
            const uint16_t      *row2_16 = (const uint16_t *)row2_8;
            const uint16_t      *row3_16 = (const uint16_t *)row3_8;
            const uint16_t      *row4_16 = (const uint16_t *)row4_8;
 
            unsigned char *dst_row_8 = out_plane + offset + y * out_w * bps_out;
 
            int kv0 = 0, kv1 = 0, kv2 = 0, kv3 = 0;
            int iwy = 0;
            if(mode5) {
                kv0 = kaiser4_lut[wy][0];
                kv1 = kaiser4_lut[wy][1];
                kv2 = kaiser4_lut[wy][2];
                kv3 = kaiser4_lut[wy][3];
            } else if(weighted) {
                iwy = WEIGHT_MAX - wy;
            }
 
            if(hl == 0)
            {
                for(int x = 0; x < out_w; x++)
                {
                    const int sx = map->w_src   [map_w_base + x];
                    const int wl = (mode >= 1) ? map->w_flag  [map_w_base + x] : 0;
                    const int wx = (mode >= 1) ? map->w_weight[map_w_base + x] : 0;
 
                    if(wl == 0)
                    {
                        WR(x, RD(row2, sx));
                    }
                    else if(mode5)
                    {
                        const int kh0 = kaiser4_lut[wx][0];
                        const int kh1 = kaiser4_lut[wx][1];
                        const int kh2 = kaiser4_lut[wx][2];
                        const int kh3 = kaiser4_lut[wx][3];
                        int p1, p2, p3, p4;
                        if(wl == 1) {
                            p1 = __builtin_expect(sx >= 1,       1) ? sx - 1 : sx;
                            p2 = sx;
                            p3 = __builtin_expect(sx + 1 < in_w, 1) ? sx + 1 : sx;
                            p4 = __builtin_expect(sx + 2 < in_w, 1) ? sx + 2 : p3;
                        } else {
                            p1 = __builtin_expect(sx + 1 < in_w, 1) ? sx + 1 : sx;
                            p2 = sx;
                            p3 = __builtin_expect(sx >= 1,       1) ? sx - 1 : sx;
                            p4 = __builtin_expect(sx >= 2,       1) ? sx - 2 : p3;
                        }
                        int v = (RD(row2,p1)*kh0 + RD(row2,p2)*kh1
                               + RD(row2,p3)*kh2 + RD(row2,p4)*kh3) >> WEIGHT_SHIFT;
                        WR(x, v);
                    }
                    else if(weighted)
                    {
                        const int iw = WEIGHT_MAX - wx;
                        int sx1 = (wl == 1) ? (__builtin_expect(sx+1 < in_w,1) ? sx+1 : sx)
                                            : (__builtin_expect(sx >= 1,    1) ? sx-1 : sx);
                        WR(x, (RD(row2,sx)*iw + RD(row2,sx1)*wx) >> WEIGHT_SHIFT);
                    }
                    else
                    {
                        int sx1 = (wl == 1) ? (__builtin_expect(sx+1 < in_w,1) ? sx+1 : sx)
                                            : (__builtin_expect(sx >= 1,    1) ? sx-1 : sx);
                        WR(x, (RD(row2,sx) + RD(row2,sx1)) >> 1);
                    }
                }
            }
            else
            {
                for(int x = 0; x < out_w; x++)
                {
                    const int sx = map->w_src   [map_w_base + x];
                    const int wl = (mode >= 1) ? map->w_flag  [map_w_base + x] : 0;
                    const int wx = (mode >= 1) ? map->w_weight[map_w_base + x] : 0;
 
                    int p1, p2, p3, p4;
                    if(wl == 1) {
                        p1 = __builtin_expect(sx >= 1,       1) ? sx - 1 : sx;
                        p2 = sx;
                        p3 = __builtin_expect(sx + 1 < in_w, 1) ? sx + 1 : sx;
                        p4 = __builtin_expect(sx + 2 < in_w, 1) ? sx + 2 : p3;
                    } else if(wl == -1) {
                        p1 = __builtin_expect(sx + 1 < in_w, 1) ? sx + 1 : sx;
                        p2 = sx;
                        p3 = __builtin_expect(sx >= 1,       1) ? sx - 1 : sx;
                        p4 = __builtin_expect(sx >= 2,       1) ? sx - 2 : p3;
                    } else {
                        p1 = p2 = p3 = p4 = sx;
                    }
 
                    if(wl == 0)
                    {
                        if(mode5)
                        {
                            int v = (RD(row1,sx)*kv0 + RD(row2,sx)*kv1
                                   + RD(row3,sx)*kv2 + RD(row4,sx)*kv3) >> WEIGHT_SHIFT;
                            WR(x, v);
                        }
                        else if(weighted)
                        {
                            WR(x, (RD(row2,sx)*iwy + RD(row3,sx)*wy) >> WEIGHT_SHIFT);
                        }
                        else
                        {
                            WR(x, (RD(row2,sx) + RD(row3,sx)) >> 1);
                        }
                    }
                    else
                    {
                        if(mode5)
                        {
                            const int kh0 = kaiser4_lut[wx][0];
                            const int kh1 = kaiser4_lut[wx][1];
                            const int kh2 = kaiser4_lut[wx][2];
                            const int kh3 = kaiser4_lut[wx][3];
                            int h1 = (RD(row1,p1)*kh0 + RD(row1,p2)*kh1
                                    + RD(row1,p3)*kh2 + RD(row1,p4)*kh3) >> WEIGHT_SHIFT;
                            int h2 = (RD(row2,p1)*kh0 + RD(row2,p2)*kh1
                                    + RD(row2,p3)*kh2 + RD(row2,p4)*kh3) >> WEIGHT_SHIFT;
                            int h3 = (RD(row3,p1)*kh0 + RD(row3,p2)*kh1
                                    + RD(row3,p3)*kh2 + RD(row3,p4)*kh3) >> WEIGHT_SHIFT;
                            int h4 = (RD(row4,p1)*kh0 + RD(row4,p2)*kh1
                                    + RD(row4,p3)*kh2 + RD(row4,p4)*kh3) >> WEIGHT_SHIFT;
                            // clamp intermédiaire
                            h1 = h1 < 0 ? 0 : h1 > 65535 ? 65535 : h1;
							h2 = h2 < 0 ? 0 : h2 > 65535 ? 65535 : h2;
							h3 = h3 < 0 ? 0 : h3 > 65535 ? 65535 : h3;
							h4 = h4 < 0 ? 0 : h4 > 65535 ? 65535 : h4;
                            int v = (h1*kv0 + h2*kv1 + h3*kv2 + h4*kv3) >> WEIGHT_SHIFT;
                            WR(x, v);
                        }
                        else if(weighted)
                        {
                            const int iw = WEIGHT_MAX - wx;
                            int h2 = (RD(row2,sx)*iw + RD(row2,p3)*wx) >> WEIGHT_SHIFT;
                            int h3 = (RD(row3,sx)*iw + RD(row3,p3)*wx) >> WEIGHT_SHIFT;
                            WR(x, (h2*iwy + h3*wy) >> WEIGHT_SHIFT);
                        }
                        else
                        {
                            int h2 = (RD(row2,sx) + RD(row2,p3)) >> 1;
                            int h3 = (RD(row3,sx) + RD(row3,p3)) >> 1;
                            WR(x, (h2 + h3) >> 1);
                        }
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
                      const int bps_in, const int bps_out)
{
    const int in_w3    = in_w * 3;
    const int out_w3   = out_w * 3;
    const int in_size  = in_w * in_h * 3;
    const int weighted = (mode == 3 || mode == 4 || mode == 5);
    const int mode5    = (mode == 5);
 
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
        else { src_buf_ptr = in_buf; }
 
        const int map_w_base = r * out_w;
        const int map_h_base = r * out_h;
 
        for(int y = 0; y < out_h; y++)
        {
            const int sy = map->h_src   [map_h_base + y];
            const int hl = (mode >= 1) ? map->h_flag  [map_h_base + y] : 0;
            const int wy = (mode >= 1) ? map->h_weight[map_h_base + y] : 0;
 
            unsigned char *dst_row = out_buf + offset + y * out_w3 * bps_out;
 
            int l1, l2, l3, l4;
            if(hl == 1) {
                l1 = __builtin_expect(sy >= 1,       1) ? sy-1 : sy;
                l2 = sy;
                l3 = __builtin_expect(sy+1 < in_h,   1) ? sy+1 : sy;
                l4 = __builtin_expect(sy+2 < in_h,   1) ? sy+2 : l3;
            } else if(hl == -1) {
                l1 = __builtin_expect(sy+1 < in_h,   1) ? sy+1 : sy;
                l2 = sy;
                l3 = __builtin_expect(sy >= 1,        1) ? sy-1 : sy;
                l4 = __builtin_expect(sy >= 2,        1) ? sy-2 : l3;
            } else { l1 = l2 = l3 = l4 = sy; }
 
            const unsigned char *row1 = src_buf_ptr + l1 * in_w3 * bps_in;
            const unsigned char *row2 = src_buf_ptr + l2 * in_w3 * bps_in;
            const unsigned char *row3 = src_buf_ptr + l3 * in_w3 * bps_in;
            const unsigned char *row4 = src_buf_ptr + l4 * in_w3 * bps_in;
 
            int kv0 = 0, kv1 = 0, kv2 = 0, kv3 = 0, iwy = 0;
            if(mode5) {
                kv0 = kaiser4_lut[wy][0]; kv1 = kaiser4_lut[wy][1];
                kv2 = kaiser4_lut[wy][2]; kv3 = kaiser4_lut[wy][3];
            } else if(weighted) { iwy = WEIGHT_MAX - wy; }
 
            if(hl == 0)
            {
                int x3 = 0;
                for(int x = 0; x < out_w; x++, x3 += 3)
                {
                    const int sx  = map->w_src   [map_w_base + x];
                    const int wl  = (mode >= 1) ? map->w_flag  [map_w_base + x] : 0;
                    const int wx  = (mode >= 1) ? map->w_weight[map_w_base + x] : 0;
                    const int sx3 = sx * 3;
 
                    if(wl == 0) {
                        WR3(dst_row, x3,   RD3(row2, sx3));
                        WR3(dst_row, x3+1, RD3(row2, sx3+1));
                        WR3(dst_row, x3+2, RD3(row2, sx3+2));
                    } else if(mode5) {
                        const int kh0=kaiser4_lut[wx][0], kh1=kaiser4_lut[wx][1];
                        const int kh2=kaiser4_lut[wx][2], kh3=kaiser4_lut[wx][3];
                        int p1,p2,p3,p4;
                        if(wl==1){p1=(sx3>=3)?sx3-3:sx3;p2=sx3;p3=(sx3+3<in_w3)?sx3+3:sx3;p4=(sx3+6<in_w3)?sx3+6:p3;}
                        else     {p1=(sx3+3<in_w3)?sx3+3:sx3;p2=sx3;p3=(sx3>=3)?sx3-3:sx3;p4=(sx3>=6)?sx3-6:p3;}
                        for(int c=0;c<3;c++)
                            WR3(dst_row, x3+c, (RD3(row2,p1+c)*kh0+RD3(row2,p2+c)*kh1+RD3(row2,p3+c)*kh2+RD3(row2,p4+c)*kh3)>>WEIGHT_SHIFT);
                    } else if(weighted) {
                        const int iw=WEIGHT_MAX-wx;
                        int sx3b=(wl==1)?((sx3+3<in_w3)?sx3+3:sx3):((sx3>=3)?sx3-3:sx3);
                        for(int c=0;c<3;c++)
                            WR3(dst_row, x3+c, (RD3(row2,sx3+c)*iw+RD3(row2,sx3b+c)*wx)>>WEIGHT_SHIFT);
                    } else {
                        int sx3b=(wl==1)?((sx3+3<in_w3)?sx3+3:sx3):((sx3>=3)?sx3-3:sx3);
                        for(int c=0;c<3;c++)
                            WR3(dst_row, x3+c, (RD3(row2,sx3+c)+RD3(row2,sx3b+c))>>1);
                    }
                }
            }
            else
            {
                int x3 = 0;
                for(int x = 0; x < out_w; x++, x3 += 3)
                {
                    const int sx  = map->w_src   [map_w_base + x];
                    const int wl  = map->w_flag  [map_w_base + x];
                    const int wx  = map->w_weight[map_w_base + x];
                    const int sx3 = sx * 3;
                    int p1,p2,p3,p4;
                    if(wl==1){p1=(sx3>=3)?sx3-3:sx3;p2=sx3;p3=(sx3+3<in_w3)?sx3+3:sx3;p4=(sx3+6<in_w3)?sx3+6:p3;}
                    else if(wl==-1){p1=(sx3+3<in_w3)?sx3+3:sx3;p2=sx3;p3=(sx3>=3)?sx3-3:sx3;p4=(sx3>=6)?sx3-6:p3;}
                    else{p1=p2=p3=p4=sx3;}
 
                    if(wl == 0) {
                        if(mode5) {
                            for(int c=0;c<3;c++){
                                int v=(RD3(row1,sx3+c)*kv0+RD3(row2,sx3+c)*kv1+RD3(row3,sx3+c)*kv2+RD3(row4,sx3+c)*kv3)>>WEIGHT_SHIFT;
                                WR3(dst_row,x3+c,v);
                            }
                        } else if(weighted) {
                            for(int c=0;c<3;c++)
                                WR3(dst_row,x3+c,(RD3(row2,sx3+c)*iwy+RD3(row3,sx3+c)*wy)>>WEIGHT_SHIFT);
                        } else {
                            for(int c=0;c<3;c++)
                                WR3(dst_row,x3+c,(RD3(row2,sx3+c)+RD3(row3,sx3+c))>>1);
                        }
                    } else {
                        if(mode5) {
                            const int kh0=kaiser4_lut[wx][0],kh1=kaiser4_lut[wx][1];
                            const int kh2=kaiser4_lut[wx][2],kh3=kaiser4_lut[wx][3];
                            for(int c=0;c<3;c++){
                                int h1=(RD3(row1,p1+c)*kh0+RD3(row1,p2+c)*kh1+RD3(row1,p3+c)*kh2+RD3(row1,p4+c)*kh3)>>WEIGHT_SHIFT;
                                int h2=(RD3(row2,p1+c)*kh0+RD3(row2,p2+c)*kh1+RD3(row2,p3+c)*kh2+RD3(row2,p4+c)*kh3)>>WEIGHT_SHIFT;
                                int h3=(RD3(row3,p1+c)*kh0+RD3(row3,p2+c)*kh1+RD3(row3,p3+c)*kh2+RD3(row3,p4+c)*kh3)>>WEIGHT_SHIFT;
                                int h4=(RD3(row4,p1+c)*kh0+RD3(row4,p2+c)*kh1+RD3(row4,p3+c)*kh2+RD3(row4,p4+c)*kh3)>>WEIGHT_SHIFT;
                                h1=h1<0?0:h1>65535?65535:h1; h2=h2<0?0:h2>65535?65535:h2;
                                h3=h3<0?0:h3>65535?65535:h3; h4=h4<0?0:h4>65535?65535:h4;
                                WR3(dst_row,x3+c,(h1*kv0+h2*kv1+h3*kv2+h4*kv3)>>WEIGHT_SHIFT);
                            }
                        } else if(weighted) {
                            const int iw=WEIGHT_MAX-wx;
                            for(int c=0;c<3;c++){
                                int h2=(RD3(row2,sx3+c)*iw+RD3(row2,p3+c)*wx)>>WEIGHT_SHIFT;
                                int h3=(RD3(row3,sx3+c)*iw+RD3(row3,p3+c)*wx)>>WEIGHT_SHIFT;
                                WR3(dst_row,x3+c,(h2*iwy+h3*wy)>>WEIGHT_SHIFT);
                            }
                        } else {
                            for(int c=0;c<3;c++){
                                int h2=(RD3(row2,sx3+c)+RD3(row2,p3+c))>>1;
                                int h3=(RD3(row3,sx3+c)+RD3(row3,p3+c))>>1;
                                WR3(dst_row,x3+c,(h2+h3)>>1);
                            }
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
                        process_files_2d(
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
                        process_files_2d_plane(
                            a->in_buf,
                            a->next_buf,
                            a->blend_buf,
                            a->out_buf,
                            a->map, a->tmp_weights,
                            a->in_w, a->in_h, a->out_w, a->out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
                        // U
                        process_files_2d_plane(
                            a->in_buf    + in_y_size,
                            a->next_buf  + in_y_size,
                            a->blend_buf + in_y_size,
                            a->out_buf   + out_y_size,
                            a->map_uv, a->tmp_weights,
                            uv_in_w, uv_in_h, uv_out_w, uv_out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
                        // V
                        process_files_2d_plane(
                            a->in_buf    + in_y_size + in_uv_size,
                            a->next_buf  + in_y_size + in_uv_size,
                            a->blend_buf + in_y_size + in_uv_size,
                            a->out_buf   + out_y_size + out_uv_size,
                            a->map_uv, a->tmp_weights,
                            uv_in_w, uv_in_h, uv_out_w, uv_out_h,
                            a->ratio, a->mode, a->bps_in, a->bps_out);
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
	
	PixFmt in_pix_fmt  = PIX_RGB;
	PixFmt out_pix_fmt = PIX_RGB;
	
	//buffer : in_buf a 1 slot supplémentaire pour le look-ahead temporel.
	//in_buf[NB_THREADS] est la frame "next" du dernier worker, qui devient
	//in_buf[0] (current du worker 0) du batch suivant.
	unsigned char *in_buf[NB_THREADS + 1];
	unsigned char *out_buf[NB_THREADS];
	unsigned char *tmp_buf[NB_THREADS];
	unsigned char *blend_buf[NB_THREADS];  // buffer du lerp temporel par worker
	for(int t=0; t<NB_THREADS; t++)
	{
		in_buf[t]    = NULL;
		out_buf[t]   = NULL;
		tmp_buf[t]   = NULL;
		blend_buf[t] = NULL;
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
	map.h_src    = NULL;
	map.h_flag   = NULL;
	map.h_weight = NULL;
	map.h_heavy  = NULL;
	map.w_length = 0;
	map.h_length = 0;
	
	BresenhamMap map_uv;
	map_uv.w_src    = NULL;
	map_uv.w_flag   = NULL;
	map_uv.w_weight = NULL;
	map_uv.w_heavy  = NULL;
	map_uv.h_src    = NULL;
	map_uv.h_flag   = NULL;
	map_uv.h_weight = NULL;
	map_uv.h_heavy  = NULL;
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
	
	if(mode < 0 || mode > 5)
	{
		fprintf(stderr, "Error : mode '%d' is not supported\n", mode);
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
	
	if((in_pix_fmt == PIX_RGB) != (out_pix_fmt == PIX_RGB))
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
		map.h_src    = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_flag   = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_weight = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);
		map.h_heavy  = ALIGNED_MALLOC(map.h_length * sizeof(int), 64);

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
			map_uv.h_src    = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_flag   = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_weight = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);
			map_uv.h_heavy  = ALIGNED_MALLOC(map_uv.h_length * sizeof(int), 64);

			if(!map_uv.w_src || !map_uv.w_flag || !map_uv.w_weight || !map_uv.w_heavy ||
			   !map_uv.h_src || !map_uv.h_flag || !map_uv.h_weight || !map_uv.h_heavy)
			{
				ALIGNED_FREE(map_uv.w_src);   ALIGNED_FREE(map_uv.w_flag);
				ALIGNED_FREE(map_uv.w_weight);ALIGNED_FREE(map_uv.w_heavy);
				ALIGNED_FREE(map_uv.h_src);   ALIGNED_FREE(map_uv.h_flag);
				ALIGNED_FREE(map_uv.h_weight);ALIGNED_FREE(map_uv.h_heavy);
				fprintf(stderr, "ALIGNED_MALLOC error (map_uv)\n");
				return -1;
			}
		}
		
		int alloc_err = 0;
		for(int t=0; t<NB_THREADS; t++)
		{
			in_buf[t]    = ALIGNED_MALLOC(in_buf_length, 64);
			out_buf[t]   = ALIGNED_MALLOC(out_buf_length, 64);
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
			if(in_buf[t] == NULL || out_buf[t] == NULL || tmp_buf[t] == NULL || blend_buf[t] == NULL)
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
			ALIGNED_FREE(map.h_src);
			ALIGNED_FREE(map.h_flag);
			ALIGNED_FREE(map.h_weight);
			ALIGNED_FREE(map.h_heavy);
			ALIGNED_FREE(tmp_weights);
			
			for(int t=0; t<NB_THREADS; t++)
			{
				ALIGNED_FREE(in_buf[t]);
				ALIGNED_FREE(out_buf[t]);
				ALIGNED_FREE(tmp_buf[t]);
				ALIGNED_FREE(blend_buf[t]);
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
	
	if(mode == 3)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_LINEAR, curve_base);
		adjust_phase_direction_weighted(&map, in_w, in_h, out_w, out_h, frames_ratio);
	}
	else if(mode == 4)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_EXPONENTIAL, curve_base);
		adjust_phase_direction_weighted(&map, in_w, in_h, out_w, out_h, frames_ratio);
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
	// --- 4b. init map + LUT (après les blocs mode 2/3/4 existants) ---
	else if(mode == 5)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_EXPONENTIAL, curve_base);
		adjust_phase_direction_weighted(&map, in_w, in_h, out_w, out_h, frames_ratio);
		compute_kaiser4_lut(kaiser_beta);
	}
	
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

		compute_scale_map(&map_uv, uv_in_w, uv_in_h, uv_out_w, uv_out_h, frames_ratio);

		if(mode == 3) {
			compute_curve_weights(&map_uv, uv_out_w, uv_out_h, frames_ratio, CURVE_LINEAR, curve_base);
			adjust_phase_direction_weighted(&map_uv, uv_in_w, uv_in_h, uv_out_w, uv_out_h, frames_ratio);
		} else if(mode == 4 || mode == 5) {
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
				if(out_pix_fmt == PIX_RGB)
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
	ALIGNED_FREE(map.h_src);
	ALIGNED_FREE(map.h_flag);
	ALIGNED_FREE(map.h_weight);
	ALIGNED_FREE(map.h_heavy);
	ALIGNED_FREE(tmp_weights);
	
	if(need_map_uv)
	{
		ALIGNED_FREE(map_uv.w_src);
		ALIGNED_FREE(map_uv.w_flag);
		ALIGNED_FREE(map_uv.w_weight);
		ALIGNED_FREE(map_uv.w_heavy);
		ALIGNED_FREE(map_uv.h_src);
		ALIGNED_FREE(map_uv.h_flag);
		ALIGNED_FREE(map_uv.h_weight);
		ALIGNED_FREE(map_uv.h_heavy);
	}

	
	for(int t=0; t<NB_THREADS; t++)
	{
		ALIGNED_FREE(in_buf[t]);
		ALIGNED_FREE(out_buf[t]);
		ALIGNED_FREE(tmp_buf[t]);
		ALIGNED_FREE(blend_buf[t]);
	}
	ALIGNED_FREE(in_buf[NB_THREADS]);
	
	//Close file
	if (input && (input != stdin))
	{
		fclose(input);
	}

	return 0;
}