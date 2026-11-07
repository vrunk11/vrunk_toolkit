#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>

#ifndef _WIN32
    #include <unistd.h>
#else
    #include <windows.h>
    #include <io.h>
    #include <fcntl.h>
#endif

#ifdef _WIN32
    #define ALIGNED_MALLOC(size, align) _aligned_malloc(size, align)
    #define ALIGNED_FREE(ptr)           _aligned_free(ptr)
#else
    #define ALIGNED_MALLOC(size, align) aligned_alloc(align, size)
    #define ALIGNED_FREE(ptr)           free(ptr)
#endif

#define _FILE_OFFSET_BITS 64

#define WEIGHT_SHIFT 8
#define WEIGHT_MAX   (1 << WEIGHT_SHIFT)

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

typedef struct {
    int *w_src;
    int *w_flag;
    int *w_weight;
    int *w_heavy;
    int *w_pos;
    int *h_src;
    int *h_flag;
    int *h_weight;
    int *h_heavy;
    int *h_pos;
    int w_length;
    int h_length;
} BresenhamMap;

// ============================================================
// Reprise exacte des fonctions de l'upscaler
// ============================================================

void compute_scale_map(BresenhamMap *map,
                       int in_w, int in_h,
                       int out_w, int out_h,
                       int ratio)
{
    for(int r = 0; r < ratio; r++)
    {
        int err = (int)(((long long)(2*r + 1) * out_w) / (2 * ratio));
        if(err >= out_w) err -= out_w;
        int src = 0, dst = r * out_w;
        for(int px = 0; px < out_w; px++)
        {
            map->w_src [dst] = (src < in_w) ? src : in_w - 1;
            map->w_flag[dst] = (2 * err >= out_w) ? 1 : 0;
            map->w_pos [dst] = ((src < in_w) ? src : in_w-1) * WEIGHT_MAX
                               + (err * WEIGHT_MAX) / out_w;
            dst++;
            err += in_w;
            if(err >= out_w) { err -= out_w; src++; }
        }
    }
    for(int r = 0; r < ratio; r++)
    {
        int err = (int)(((long long)(2*r + 1) * out_h) / (2 * ratio));
        if(err >= out_h) err -= out_h;
        int src = 0, dst = r * out_h;
        for(int px = 0; px < out_h; px++)
        {
            map->h_src [dst] = (src < in_h) ? src : in_h - 1;
            map->h_flag[dst] = (2 * err >= out_h) ? 1 : 0;
            map->h_pos [dst] = ((src < in_h) ? src : in_h-1) * WEIGHT_MAX
                               + (err * WEIGHT_MAX) / out_h;
            dst++;
            err += in_h;
            if(err >= out_h) { err -= out_h; src++; }
        }
    }
}

// Reprise exacte de compute_curve_weights (les deux modes)
// Nécessaire pour que flag=0 + weight=0 soit le vrai critère de pixel pur.
typedef enum { CURVE_EXPONENTIAL = 0, CURVE_LINEAR = 1 } CurveType;

void compute_curve_weights(BresenhamMap *map,
                           int out_w, int out_h, int ratio,
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
                const int f = map->w_flag[row + run_start];
                map->w_weight[row + run_start] = (f != 0) ? (WEIGHT_MAX / 2) : 0;
                map->w_heavy [row + run_start] = (f != 0) ? 1 : 0;
            }
            else
            {
                const int heavy_threshold = (N + 1) / 2;
                const double denom = (curve == CURVE_EXPONENTIAL)
                                   ? pow(base, (double)(N-1)) : 1.0;
                for(int i = 0; i < N; i++)
                {
                    double w = (curve == CURVE_LINEAR)
                             ? (double)i / (double)N
                             : (pow(base, (double)i) - 1.0) / denom;
                    int wi = (int)(w * WEIGHT_MAX + 0.5);
                    wi = wi < 0 ? 0 : wi > WEIGHT_MAX ? WEIGHT_MAX : wi;
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
                                   ? pow(base, (double)(N-1)) : 1.0;
                for(int i = 0; i < N; i++)
                {
                    double w = (curve == CURVE_LINEAR)
                             ? (double)i / (double)N
                             : (pow(base, (double)i) - 1.0) / denom;
                    int wi = (int)(w * WEIGHT_MAX + 0.5);
                    wi = wi < 0 ? 0 : wi > WEIGHT_MAX ? WEIGHT_MAX : wi;
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
// Table inverse : src → out_px
// Critère pixel pur : flag==0 ET weight==0 (après curve_weights)
// ============================================================
static int *build_inverse_map_1d(const int *src_map,
                                 const int *flag_map,
                                 const int *weight_map,
                                 int out_len,
                                 int in_len,
                                 int *missing_out)
{
    int *inv = malloc(in_len * sizeof(int));
    if(!inv) return NULL;
    for(int i = 0; i < in_len; i++) inv[i] = -1;

    for(int ox = 0; ox < out_len; ox++)
    {
        int sx = src_map[ox];
        // pixel pur = flag==0 ET weight==0
        if(flag_map[ox] == 0 && weight_map[ox] == 0 && inv[sx] == -1)
            inv[sx] = ox;
    }

    int missing = 0;
    for(int sx = 0; sx < in_len; sx++)
        if(inv[sx] == -1) missing++;
    *missing_out = missing;
    return inv;
}

// ============================================================
// Reconstruction d'un plan planar
// ============================================================
static int reconstruct_plane(
    const unsigned char *upscale_plane,
    unsigned char       *out_plane,
    int in_w, int in_h,
    int out_w,
    const int *inv_w, const int *inv_h,
    int bps)
{
    int errors = 0;
    for(int sy = 0; sy < in_h; sy++)
    {
        int oy = inv_h[sy];
        for(int sx = 0; sx < in_w; sx++)
        {
            int ox = inv_w[sx];
            if(oy < 0 || ox < 0) {
                if(bps == 1) out_plane[sy*in_w + sx] = 0;
                else ((uint16_t*)out_plane)[sy*in_w + sx] = 0;
                errors++;
            } else {
                if(bps == 1) out_plane[sy*in_w + sx] = upscale_plane[oy*out_w + ox];
                else ((uint16_t*)out_plane)[sy*in_w + sx] =
                         ((const uint16_t*)upscale_plane)[oy*out_w + ox];
            }
        }
    }
    return errors;
}

static int reconstruct_rgb(
    const unsigned char *upscale_buf,
    unsigned char       *out_buf,
    int in_w, int in_h,
    int out_w,
    const int *inv_w, const int *inv_h,
    int bps)
{
    int errors = 0;
    const int bps3 = 3 * bps;
    for(int sy = 0; sy < in_h; sy++)
    {
        int oy = inv_h[sy];
        for(int sx = 0; sx < in_w; sx++)
        {
            int ox = inv_w[sx];
            unsigned char       *dst = out_buf     + (sy*in_w  + sx)*bps3;
            const unsigned char *src = upscale_buf + (oy*out_w + ox)*bps3;
            if(oy < 0 || ox < 0) { memset(dst, 0, bps3); errors++; }
            else                   memcpy(dst, src, bps3);
        }
    }
    return errors;
}

// ============================================================
// Aide
// ============================================================
void usage(void)
{
    fprintf(stderr,
        "reverse_scale - reconstruit la source depuis un upscale\n\n"
        "Usage:\n"
        "\t-i filename        fichier upscalé (ou '-' pour stdin)\n\n"
        "Options (identiques à l'upscaler) :\n"
        "\t--in-width  INT    largeur  SOURCE  (défaut: 720)\n"
        "\t--in-height INT    hauteur  SOURCE  (défaut: 480)\n"
        "\t--out-width INT    largeur  UPSCALE (défaut: 1920)\n"
        "\t--out-height INT   hauteur  UPSCALE (défaut: 1080)\n"
        "\t--fps-ratio INT    multiplicateur fps (défaut: 2)\n"
        "\t--mode INT         mode upscaler utilisé (défaut: 7)\n"
        "\t                   important : doit correspondre exactement\n"
        "\t--curve-base FLOAT base exponentielle (défaut: 2.0)\n"
        "\t--pix-fmt   STR    RGB, RGB48, 444p, 422p, 420p, 411p, 410p,\n"
        "\t                   444p16, 422p16, 420p16\n"
    );
    exit(1);
}

// ============================================================
// main
// ============================================================
int main(int argc, char **argv)
{
#if defined(_WIN32) || defined(_WIN64)
    _setmode(_fileno(stdout), O_BINARY);
    _setmode(_fileno(stdin),  O_BINARY);
#endif

    char  *input_name   = NULL;
    FILE  *input        = NULL;

    int    in_w         = 720;
    int    in_h         = 480;
    int    out_w        = 1920;
    int    out_h        = 1080;
    int    frames_ratio = 2;
    int    mode         = 7;
    double curve_base   = 2.0;
    PixFmt pix_fmt      = PIX_RGB;

    static struct option long_options[] = {
        {"in-width",   1, 0, 1},
        {"in-height",  1, 0, 2},
        {"out-width",  1, 0, 3},
        {"out-height", 1, 0, 4},
        {"fps-ratio",  1, 0, 5},
        {"pix-fmt",    1, 0, 6},
        {"mode",       1, 0, 7},
        {"curve-base", 1, 0, 8},
        {0, 0, 0, 0}
    };
    int opt, option_index = 0;
    while((opt = getopt_long_only(argc, argv, "i:", long_options, &option_index)) != -1)
    {
        switch(opt) {
        case 'i': input_name  = optarg;        break;
        case 1:   in_w        = atoi(optarg);  break;
        case 2:   in_h        = atoi(optarg);  break;
        case 3:   out_w       = atoi(optarg);  break;
        case 4:   out_h       = atoi(optarg);  break;
        case 5:   frames_ratio = atoi(optarg); break;
        case 6:
            if     (!strcasecmp(optarg,"rgb")   || !strcasecmp(optarg,"rgb24")) pix_fmt = PIX_RGB;
            else if(!strcasecmp(optarg,"rgb48"))                                 pix_fmt = PIX_RGB48;
            else if(!strcasecmp(optarg,"444p")  || !strcasecmp(optarg,"444"))   pix_fmt = PIX_444P;
            else if(!strcasecmp(optarg,"422p")  || !strcasecmp(optarg,"422"))   pix_fmt = PIX_422P;
            else if(!strcasecmp(optarg,"420p")  || !strcasecmp(optarg,"420"))   pix_fmt = PIX_420P;
            else if(!strcasecmp(optarg,"411p")  || !strcasecmp(optarg,"411"))   pix_fmt = PIX_411P;
            else if(!strcasecmp(optarg,"410p")  || !strcasecmp(optarg,"410"))   pix_fmt = PIX_410P;
            else if(!strcasecmp(optarg,"444p16"))                                pix_fmt = PIX_444P16;
            else if(!strcasecmp(optarg,"422p16"))                                pix_fmt = PIX_422P16;
            else if(!strcasecmp(optarg,"420p16"))                                pix_fmt = PIX_420P16;
            else { fprintf(stderr, "format inconnu : %s\n", optarg); return -1; }
            break;
        case 7:   mode        = atoi(optarg);  break;
        case 8:   curve_base  = atof(optarg);  break;
        default:  usage(); break;
        }
    }

    if(!input_name) { fprintf(stderr, "pas de fichier d'entrée\n"); return -1; }
    if(out_w <= in_w || out_h <= in_h) {
        fprintf(stderr, "erreur : out doit être > in\n"); return -1;
    }

    const int bps = (pix_fmt >= PIX_RGB48) ? 2 : 1;

    const int uv_in_w  = (pix_fmt==PIX_444P||pix_fmt==PIX_444P16) ? in_w    :
                         (pix_fmt==PIX_411P ||pix_fmt==PIX_410P)  ? in_w/4  : in_w/2;
    const int uv_in_h  = (pix_fmt==PIX_420P||pix_fmt==PIX_420P16) ? in_h/2  :
                         (pix_fmt==PIX_410P)                      ? in_h/4  : in_h;
    const int uv_out_w = (pix_fmt==PIX_444P||pix_fmt==PIX_444P16) ? out_w   :
                         (pix_fmt==PIX_411P ||pix_fmt==PIX_410P)  ? out_w/4 : out_w/2;
    const int uv_out_h = (pix_fmt==PIX_420P||pix_fmt==PIX_420P16) ? out_h/2 :
                         (pix_fmt==PIX_410P)                      ? out_h/4 : out_h;

    // ---- tailles buffers ----
    unsigned int upscale_frame_size = 0, source_frame_size = 0;
    switch(pix_fmt) {
        case PIX_RGB:    upscale_frame_size = out_w*out_h*3;           break;
        case PIX_RGB48:  upscale_frame_size = out_w*out_h*3*2;         break;
        case PIX_444P:
        case PIX_444P16: upscale_frame_size = out_w*out_h*3*bps;       break;
        case PIX_422P:
        case PIX_422P16: upscale_frame_size = out_w*out_h*2*bps;       break;
        case PIX_420P:
        case PIX_420P16: upscale_frame_size = out_w*out_h*3/2*bps;     break;
        case PIX_411P:   upscale_frame_size = out_w*out_h*3/2;         break;
        case PIX_410P:   upscale_frame_size = out_w*out_h*9/8;         break;
    }
    switch(pix_fmt) {
        case PIX_RGB:    source_frame_size = in_w*in_h*3;              break;
        case PIX_RGB48:  source_frame_size = in_w*in_h*3*2;            break;
        case PIX_444P:
        case PIX_444P16: source_frame_size = in_w*in_h*3*bps;          break;
        case PIX_422P:
        case PIX_422P16: source_frame_size = in_w*in_h*2*bps;          break;
        case PIX_420P:
        case PIX_420P16: source_frame_size = in_w*in_h*3/2*bps;        break;
        case PIX_411P:   source_frame_size = in_w*in_h*3/2;            break;
        case PIX_410P:   source_frame_size = in_w*in_h*9/8;            break;
    }

    const unsigned int upscale_group_size = upscale_frame_size * frames_ratio;

    unsigned char *upscale_buf = ALIGNED_MALLOC(upscale_group_size, 64);
    unsigned char *source_buf  = ALIGNED_MALLOC(source_frame_size,  64);
    if(!upscale_buf || !source_buf) { fprintf(stderr,"erreur malloc\n"); return -1; }

    // ---- construction de la map Y ----
    BresenhamMap map;
    memset(&map, 0, sizeof(map));
    map.w_length = out_w * frames_ratio;
    map.h_length = out_h * frames_ratio;
    map.w_src    = ALIGNED_MALLOC(map.w_length*sizeof(int), 64);
    map.w_flag   = ALIGNED_MALLOC(map.w_length*sizeof(int), 64);
    map.w_weight = ALIGNED_MALLOC(map.w_length*sizeof(int), 64);
    map.w_heavy  = ALIGNED_MALLOC(map.w_length*sizeof(int), 64);
    map.w_pos    = ALIGNED_MALLOC(map.w_length*sizeof(int), 64);
    map.h_src    = ALIGNED_MALLOC(map.h_length*sizeof(int), 64);
    map.h_flag   = ALIGNED_MALLOC(map.h_length*sizeof(int), 64);
    map.h_weight = ALIGNED_MALLOC(map.h_length*sizeof(int), 64);
    map.h_heavy  = ALIGNED_MALLOC(map.h_length*sizeof(int), 64);
    map.h_pos    = ALIGNED_MALLOC(map.h_length*sizeof(int), 64);
    if(!map.w_src||!map.w_flag||!map.w_weight||!map.w_heavy||!map.w_pos||
       !map.h_src||!map.h_flag||!map.h_weight||!map.h_heavy||!map.h_pos)
    { fprintf(stderr,"erreur malloc (map)\n"); return -1; }

    compute_scale_map(&map, in_w, in_h, out_w, out_h, frames_ratio);

    // Appliquer compute_curve_weights exactement comme l'upscaler le fait
    // pour les modes 3/4/5/6/7 (tous sauf 0/1/2 qui n'appellent pas curve_weights)
    if(mode >= 3) {
        CurveType ct = (mode == 3) ? CURVE_LINEAR : CURVE_EXPONENTIAL;
        compute_curve_weights(&map, out_w, out_h, frames_ratio, ct, curve_base);
    }
    // modes 0/1/2 : pas de curve_weights → flag seul suffit (weight reste 0)

    // ---- map UV ----
    const int is_rgb = (pix_fmt==PIX_RGB || pix_fmt==PIX_RGB48);
    const int is_444 = (pix_fmt==PIX_444P|| pix_fmt==PIX_444P16);

    BresenhamMap map_uv;
    memset(&map_uv, 0, sizeof(map_uv));
    int *inv_w_uv = NULL, *inv_h_uv = NULL;

    if(!is_rgb)
    {
        const int uw = is_444 ? out_w : uv_out_w;
        const int uh = is_444 ? out_h : uv_out_h;
        const int sw = is_444 ? in_w  : uv_in_w;
        const int sh = is_444 ? in_h  : uv_in_h;

        if(!is_444) {
            map_uv.w_length = uw * frames_ratio;
            map_uv.h_length = uh * frames_ratio;
            map_uv.w_src    = ALIGNED_MALLOC(map_uv.w_length*sizeof(int), 64);
            map_uv.w_flag   = ALIGNED_MALLOC(map_uv.w_length*sizeof(int), 64);
            map_uv.w_weight = ALIGNED_MALLOC(map_uv.w_length*sizeof(int), 64);
            map_uv.w_heavy  = ALIGNED_MALLOC(map_uv.w_length*sizeof(int), 64);
            map_uv.w_pos    = ALIGNED_MALLOC(map_uv.w_length*sizeof(int), 64);
            map_uv.h_src    = ALIGNED_MALLOC(map_uv.h_length*sizeof(int), 64);
            map_uv.h_flag   = ALIGNED_MALLOC(map_uv.h_length*sizeof(int), 64);
            map_uv.h_weight = ALIGNED_MALLOC(map_uv.h_length*sizeof(int), 64);
            map_uv.h_heavy  = ALIGNED_MALLOC(map_uv.h_length*sizeof(int), 64);
            map_uv.h_pos    = ALIGNED_MALLOC(map_uv.h_length*sizeof(int), 64);
            if(!map_uv.w_src||!map_uv.w_flag||!map_uv.w_weight||!map_uv.w_heavy||
               !map_uv.w_pos||!map_uv.h_src||!map_uv.h_flag||!map_uv.h_weight||
               !map_uv.h_heavy||!map_uv.h_pos)
            { fprintf(stderr,"erreur malloc (map_uv)\n"); return -1; }

            compute_scale_map(&map_uv, sw, sh, uw, uh, frames_ratio);
            if(mode >= 3) {
                CurveType ct = (mode==3) ? CURVE_LINEAR : CURVE_EXPONENTIAL;
                compute_curve_weights(&map_uv, uw, uh, frames_ratio, ct, curve_base);
            }
        }

        // Tables inverses UV
        const BresenhamMap *m = is_444 ? &map : &map_uv;
        const int fw = is_444 ? out_w : uw;
        const int fh = is_444 ? out_h : uh; (void)fh;
        const int fsw = is_444 ? in_w : sw;
        const int fsh = is_444 ? in_h : sh; (void)fsh;

        int mw=0, mh=0;
        inv_w_uv = build_inverse_map_1d(m->w_src, m->w_flag, m->w_weight, fw, fsw, &mw);
        inv_h_uv = build_inverse_map_1d(m->h_src, m->h_flag, m->h_weight,
                                        is_444?out_h:uh, fsh, &mh);
        if(!inv_w_uv||!inv_h_uv) { fprintf(stderr,"erreur malloc (inv_uv)\n"); return -1; }
        if(mw>0||mh>0)
            fprintf(stderr,"ERREUR FATALE (chroma) : %d col. %d lig. sans pixel intact\n",mw,mh);
    }

    // ---- tables inverses Y / RGB ----
    int missing_w=0, missing_h=0;
    int *inv_w = build_inverse_map_1d(map.w_src, map.w_flag, map.w_weight,
                                      out_w, in_w, &missing_w);
    int *inv_h = build_inverse_map_1d(map.h_src, map.h_flag, map.h_weight,
                                      out_h, in_h, &missing_h);
    if(!inv_w||!inv_h) { fprintf(stderr,"erreur malloc (inv)\n"); return -1; }

    if(missing_w>0||missing_h>0)
        fprintf(stderr,
            "ERREUR FATALE : %d col. et %d lig. source sans pixel intact.\n"
            "Ces pixels seront noirs.\n", missing_w, missing_h);

    fprintf(stderr,
        "Map inverse (Y/RGB) : %d/%d colonnes intactes, %d/%d lignes intactes\n",
        in_w-missing_w, in_w, in_h-missing_h, in_h);

    // ---- ouverture fichier ----
    if(strcmp(input_name,"-")==0) input=stdin;
    else {
        input=fopen(input_name,"rb");
        if(!input) { fprintf(stderr,"impossible d'ouvrir %s\n",input_name); return -ENOENT; }
    }

    // ---- boucle principale ----
    int total_frames=0, total_errors=0, fatal_reported=0;

    while(fread(upscale_buf, upscale_group_size, 1, input)==1)
    {
        // Phase 0 = premiers upscale_frame_size octets du groupe
        const unsigned char *phase0 = upscale_buf;
        int frame_errors = 0;

        if(is_rgb)
        {
            frame_errors += reconstruct_rgb(phase0, source_buf,
                                            in_w, in_h, out_w,
                                            inv_w, inv_h, bps);
        }
        else
        {
            const int out_y_size = out_w * out_h * bps;
            const int in_y_size  = in_w  * in_h  * bps;
            const int uv_ow = is_444 ? out_w : uv_out_w;
            const int uv_iw = is_444 ? in_w  : uv_in_w;
            const int uv_ih = is_444 ? in_h  : uv_in_h;
            const int out_uv_size = uv_ow*(is_444?out_h:uv_out_h)*bps;
            const int in_uv_size  = uv_iw*uv_ih*bps;
            const int *iw_uv = inv_w_uv ? inv_w_uv : inv_w;
            const int *ih_uv = inv_h_uv ? inv_h_uv : inv_h;

            frame_errors += reconstruct_plane(phase0, source_buf,
                                              in_w,in_h,out_w, inv_w,inv_h, bps);
            frame_errors += reconstruct_plane(phase0+out_y_size,
                                              source_buf+in_y_size,
                                              uv_iw,uv_ih,uv_ow, iw_uv,ih_uv, bps);
            frame_errors += reconstruct_plane(phase0+out_y_size+out_uv_size,
                                              source_buf+in_y_size+in_uv_size,
                                              uv_iw,uv_ih,uv_ow, iw_uv,ih_uv, bps);
        }

        if(frame_errors>0 && !fatal_reported) {
            fprintf(stderr,"ERREUR FATALE frame %d : %d pixel(s) noirs\n",
                    total_frames, frame_errors);
            fatal_reported=1;
        }
        total_errors += frame_errors;

        if(!isatty(STDOUT_FILENO))
            fwrite(source_buf, source_frame_size, 1, stdout);

        total_frames++;
    }

    fprintf(stderr, "%d frame(s) reconstruite(s), %d erreur(s) totale(s).\n",
            total_frames, total_errors);
    if(total_errors>0)
        fprintf(stderr,"=> NON lossless : %d pixels source sans exemplaire intact.\n",
                total_errors);
    else
        fprintf(stderr,"=> Lossless confirmé.\n");

    // ---- nettoyage ----
    ALIGNED_FREE(upscale_buf); ALIGNED_FREE(source_buf);
    free(inv_w); free(inv_h); free(inv_w_uv); free(inv_h_uv);
    ALIGNED_FREE(map.w_src);   ALIGNED_FREE(map.w_flag);
    ALIGNED_FREE(map.w_weight);ALIGNED_FREE(map.w_heavy); ALIGNED_FREE(map.w_pos);
    ALIGNED_FREE(map.h_src);   ALIGNED_FREE(map.h_flag);
    ALIGNED_FREE(map.h_weight);ALIGNED_FREE(map.h_heavy); ALIGNED_FREE(map.h_pos);
    if(!is_rgb && !is_444) {
        ALIGNED_FREE(map_uv.w_src);   ALIGNED_FREE(map_uv.w_flag);
        ALIGNED_FREE(map_uv.w_weight);ALIGNED_FREE(map_uv.w_heavy);ALIGNED_FREE(map_uv.w_pos);
        ALIGNED_FREE(map_uv.h_src);   ALIGNED_FREE(map_uv.h_flag);
        ALIGNED_FREE(map_uv.h_weight);ALIGNED_FREE(map_uv.h_heavy);ALIGNED_FREE(map_uv.h_pos);
    }
    if(input && input!=stdin) fclose(input);
    return (total_errors>0) ? 1 : 0;
}