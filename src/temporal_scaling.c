#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
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
	unsigned char *out_buf;
	unsigned char *tmp_buf;
	BresenhamMap *map;
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
        "\t--curve-base FLOAT exponential curve base for mode 4 (default: 2.0, must be > 1.0)\n\n"
    );
    exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

void adjust_bilinear_phases(BresenhamMap *map, int out_w, int out_h, int in_w, int in_h, int ratio)
{
    // ajustement axe horizontal
    for(int r=1; r<ratio; r++)
    {
        // etape 1 : compte le chevauchement global
        int overlap = 0;
        for(int px=0; px<out_w; px++)
            if(map->w_flag[r*out_w+px] == 1 && map->w_flag[(r-1)*out_w+px] == 1)
                overlap++;
        
        // etape 2 : si trop de chevauchement on inverse tous les 1 en -1
        if(overlap > out_w/2)
            for(int px=0; px<out_w; px++)
                if(map->w_flag[r*out_w+px] == 1)
                    map->w_flag[r*out_w+px] = -1;
    }

    // ajustement axe vertical
    for(int r=1; r<ratio; r++)
    {
        // etape 1 : compte le chevauchement global
        int overlap = 0;
        for(int px=0; px<out_h; px++)
            if(map->h_flag[r*out_h+px] == 1 && map->h_flag[(r-1)*out_h+px] == 1)
                overlap++;
        
        // etape 2 : si trop de chevauchement on inverse tous les 1 en -1
        if(overlap > out_h/2)
            for(int px=0; px<out_h; px++)
                if(map->h_flag[r*out_h+px] == 1)
                    map->h_flag[r*out_h+px] = -1;
    }
}

void compute_scale_map(BresenhamMap *map, int in_w, int in_h, int out_w, int out_h, int ratio)
{
    for(int r = 0; r < ratio; r++)
    {
        int err = (r * in_w) / ratio;  // phase en unités source, pas dest
        int src = 0;
        int dst = r * out_w;
        for(int px = 0; px < out_w; px++)
        {
			map->w_src[dst]  = src;
			map->w_flag[dst] = (err > out_w / 2) ? 1 : 0;
			dst++;
            err += in_w;
            if(err >= out_w) { err -= out_w; src++; }
        }
    }

    for(int r = 0; r < ratio; r++)
    {
        int err = (r * in_h) / ratio;
        int src = 0;
        int dst = r * out_h;
        for(int px = 0; px < out_h; px++)
        {
            map->h_src[dst]  = src;
			map->h_flag[dst] = (err > out_h / 2) ? 1 : 0;
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
// adjust_bilinear_phases_weighted pour décider d'inverser la direction.
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

// Variante de adjust_bilinear_phases pour le mode 3 :
// le critère de chevauchement utilise w_heavy/h_heavy (pixels dans la
// moitié haute du run, ie. ceux fortement interpolés) au lieu de w_flag.
// On veut basculer la direction quand les "vraies" interpolations s'alignent
// entre phases consécutives, pas quand les pixels juste légèrement
// interpolés s'alignent (sinon on flip presque toujours).
void adjust_bilinear_phases_weighted(BresenhamMap *map, int out_w, int out_h, int ratio)
{
    // axe horizontal
    for(int r = 1; r < ratio; r++)
    {
        int overlap = 0;
        for(int px = 0; px < out_w; px++)
            if(map->w_heavy[r*out_w + px] && map->w_heavy[(r-1)*out_w + px])
                overlap++;

        if(overlap > out_w / 2)
            for(int px = 0; px < out_w; px++)
                if(map->w_flag[r*out_w + px] == 1)
                    map->w_flag[r*out_w + px] = -1;
    }

    // axe vertical
    for(int r = 1; r < ratio; r++)
    {
        int overlap = 0;
        for(int px = 0; px < out_h; px++)
            if(map->h_heavy[r*out_h + px] && map->h_heavy[(r-1)*out_h + px])
                overlap++;

        if(overlap > out_h / 2)
            for(int px = 0; px < out_h; px++)
                if(map->h_flag[r*out_h + px] == 1)
                    map->h_flag[r*out_h + px] = -1;
    }
}

void process_files(const unsigned char * restrict in_buf, unsigned char * restrict out_buf, unsigned char * restrict tmp_buf,const BresenhamMap *  restrict map,const int in_w,const int in_h,const int out_w,const int out_h,const int ratio,const int mode)
{
	const int in_w3    = in_w * 3;
	const int out_w3   = out_w * 3;
	const int bilinear = (mode >= 1);
	const int weighted = (mode == 3 || mode == 4);
    unsigned int offset = 0;
    for(int r=0; r<ratio; r++)
    {
        //scale the width line by line
        const int map_w_base = r * out_w;  // hissé hors des boucles

		for(int h=0; h<in_h; h++)
		{
			const unsigned char * restrict src_row = in_buf  + h * in_w3;   // calculé une fois par ligne
			unsigned char       *dst_row = tmp_buf + h * out_w3;

			int px3 = 0;
			for(int px=0; px<out_w; px++, px3+=3)  // px3 remplace px/3
			{
				const int src_px      = map->w_src [map_w_base + px] * 3;
				const int is_bilinear = map->w_flag[map_w_base + px];

				if(weighted && is_bilinear == 1 && src_px + 3 < in_w3)
				{
					const int w  = map->w_weight[map_w_base + px];
					const int iw = WEIGHT_MAX - w;
					dst_row[px3]   = (src_row[src_px]   * iw + src_row[src_px+3] * w) >> WEIGHT_SHIFT;
					dst_row[px3+1] = (src_row[src_px+1] * iw + src_row[src_px+4] * w) >> WEIGHT_SHIFT;
					dst_row[px3+2] = (src_row[src_px+2] * iw + src_row[src_px+5] * w) >> WEIGHT_SHIFT;
				}
				else if(weighted && is_bilinear == -1 && src_px >= 3)
				{
					const int w  = map->w_weight[map_w_base + px];
					const int iw = WEIGHT_MAX - w;
					dst_row[px3]   = (src_row[src_px]   * iw + src_row[src_px-3] * w) >> WEIGHT_SHIFT;
					dst_row[px3+1] = (src_row[src_px+1] * iw + src_row[src_px-2] * w) >> WEIGHT_SHIFT;
					dst_row[px3+2] = (src_row[src_px+2] * iw + src_row[src_px-1] * w) >> WEIGHT_SHIFT;
				}
				else if(bilinear && !weighted && is_bilinear == 1 && src_px + 3 < in_w3)
				{
					dst_row[px3]   = (src_row[src_px]   + src_row[src_px+3]) >> 1;
					dst_row[px3+1] = (src_row[src_px+1] + src_row[src_px+4]) >> 1;
					dst_row[px3+2] = (src_row[src_px+2] + src_row[src_px+5]) >> 1;
				}
				else if(bilinear && !weighted && is_bilinear == -1 && src_px >= 3)
				{
					dst_row[px3]   = (src_row[src_px-3] + src_row[src_px])   >> 1;
					dst_row[px3+1] = (src_row[src_px-2] + src_row[src_px+1]) >> 1;
					dst_row[px3+2] = (src_row[src_px-1] + src_row[src_px+2]) >> 1;
				}
				else
				{
					dst_row[px3]   = src_row[src_px];
					dst_row[px3+1] = src_row[src_px+1];
					dst_row[px3+2] = src_row[src_px+2];
				}
			}
		}

        //scale the height line by line
		const int map_h_base = r * out_h;

		for(int px=0; px<out_h; px++)
		{
			const int src_line    = map->h_src [map_h_base + px];
			const int is_bilinear = map->h_flag[map_h_base + px];

			unsigned char       *dst_row  = out_buf + offset + px * out_w3;
			const unsigned char * restrict tmp_row0 = tmp_buf + src_line * out_w3;

			if(weighted && is_bilinear == 1 && src_line + 1 < in_h)
			{
				const int w = map->h_weight[map_h_base + px];
				if(w == 0)
				{
					memcpy(dst_row, tmp_row0, out_w3);
				}
				else
				{
					const int iw = WEIGHT_MAX - w;
					const unsigned char * restrict tmp_row1 = tmp_row0 + out_w3;
					for(int i=0; i<out_w3; i++)
						dst_row[i] = (tmp_row0[i] * iw + tmp_row1[i] * w) >> WEIGHT_SHIFT;
				}
			}
			else if(weighted && is_bilinear == -1 && src_line >= 1)
			{
				const int w = map->h_weight[map_h_base + px];
				if(w == 0)
				{
					memcpy(dst_row, tmp_row0, out_w3);
				}
				else
				{
					const int iw = WEIGHT_MAX - w;
					const unsigned char * restrict tmp_row_1 = tmp_row0 - out_w3;
					for(int i=0; i<out_w3; i++)
						dst_row[i] = (tmp_row0[i] * iw + tmp_row_1[i] * w) >> WEIGHT_SHIFT;
				}
			}
			else if(bilinear && !weighted && is_bilinear == 1 && src_line + 1 < in_h)
			{
				const unsigned char * restrict tmp_row1 = tmp_row0 + out_w3;
				for(int i=0; i<out_w3; i++)
					dst_row[i] = (tmp_row0[i] + tmp_row1[i]) >> 1;  // ← SIMD-vectorisable
			}
			else if(bilinear && !weighted && is_bilinear == -1 && src_line >= 1)
			{
				const unsigned char *tmp_row_1 = tmp_row0 - out_w3;
				for(int i=0; i<out_w3; i++)
					dst_row[i] = (tmp_row_1[i] + tmp_row0[i]) >> 1;
			}
			else
			{
				memcpy(dst_row, tmp_row0, out_w3);  // ← bien plus rapide que pixel par pixel
			}
		}
        offset += out_w*out_h*3;
    }
}

void *worker_loop(void *arg) {
    Worker *w = (Worker *)arg;
    while (1) {
        pthread_mutex_lock(&w->mutex);
        while (w->ready == 0)                    // attend un signal
            pthread_cond_wait(&w->cond_work, &w->mutex);
        int cmd = w->ready;
        w->ready = 0;                            // ← reset ICI sous le mutex
        pthread_mutex_unlock(&w->mutex);

        if (cmd == -1) break;

        ThreadArgs *a = &w->args;
        process_files(a->in_buf, a->out_buf, a->tmp_buf, a->map,
                      a->in_w, a->in_h, a->out_w, a->out_h, a->ratio, a->mode);

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
	double curve_base = 2.0;  // base de la courbe exponentielle (mode 3)
	
	//buffer
	unsigned char *in_buf[NB_THREADS];
	unsigned char *out_buf[NB_THREADS];
	unsigned char *tmp_buf[NB_THREADS];
	for(int t=0; t<NB_THREADS; t++)
	{
		in_buf[t] = NULL;
		out_buf[t] = NULL;
		tmp_buf[t] = NULL;
	}
	
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
	
	int option_index = 0;
	static struct option long_options[] = {
		{"in-width",   1, 0, 1},
		{"in-height",  1, 0, 2},
		{"out-width",  1, 0, 3},
		{"out-height", 1, 0, 4},
		{"fps-ratio",  1, 0, 5},
		{"mode",       1, 0, 6},
		{"curve-base", 1, 0, 7},
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
	
	if(mode < 0 || mode > 4)
	{
		fprintf(stderr, "Error : mode '%d' is not supported\n", mode);
		return -1;
	}
	
	if(mode == 4 && curve_base <= 1.0)
	{
		fprintf(stderr, "Error : --curve-base must be > 1.0 (got %f)\n", curve_base);
		return -1;
	}

	if(mode == 3 && curve_base != 2.0)
		fprintf(stderr, "Warning : --curve-base is ignored in mode 3 (linear curve)\n");
	
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
	
	if(out_h > in_h && out_w > in_w && in_h > 0 && in_w > 0)
	{
		in_buf_length  = in_w * in_h * 3;
		out_buf_length = out_w * out_h * 3 * frames_ratio;
		
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
		
		int alloc_err = 0;
		for(int t=0; t<NB_THREADS; t++)
		{
			in_buf[t]  = ALIGNED_MALLOC(in_buf_length, 64);
			out_buf[t] = ALIGNED_MALLOC(out_buf_length, 64);
			tmp_buf[t] = ALIGNED_MALLOC(out_w * in_h * 3, 64);
			if(in_buf[t] == NULL || out_buf[t] == NULL || tmp_buf[t] == NULL)
				alloc_err = 1;
		}
		
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
			
			for(int t=0; t<NB_THREADS; t++)
			{
				ALIGNED_FREE(in_buf[t]);
				ALIGNED_FREE(out_buf[t]);
				ALIGNED_FREE(tmp_buf[t]);
			}
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
	
	if(mode == 3)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_LINEAR, curve_base);
		adjust_bilinear_phases_weighted(&map, out_w, out_h, frames_ratio);
	}
	else if(mode == 4)
	{
		compute_curve_weights(&map, out_w, out_h, frames_ratio, CURVE_EXPONENTIAL, curve_base);
		adjust_bilinear_phases_weighted(&map, out_w, out_h, frames_ratio);
	}
	else if(mode == 2)
	{
		adjust_bilinear_phases(&map, out_w, out_h, in_w, in_h, frames_ratio);
	}
	
	ThreadArgs args[NB_THREADS];
	for(int t=0; t<NB_THREADS; t++)
	{
		args[t].map   = &map;
		args[t].in_w  = in_w;
		args[t].in_h  = in_h;
		args[t].out_w = out_w;
		args[t].out_h = out_h;
		args[t].ratio = frames_ratio;
		args[t].mode  = mode;
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

	while (!feof(input))
	{
		int frames_read = 0;
		for (int t = 0; t < NB_THREADS; t++) {
			if (fread(in_buf[t], in_buf_length, 1, input) == 1) {
				workers[t].args         = args[t];
				workers[t].args.in_buf  = in_buf[t];
				workers[t].args.out_buf = out_buf[t];
				workers[t].args.tmp_buf = tmp_buf[t];
				frames_read++;
			} else break;
		}

		// signal aux workers
		for (int t = 0; t < frames_read; t++) {
			pthread_mutex_lock(&workers[t].mutex);
			workers[t].done  = 0;
			workers[t].ready = 1;
			pthread_cond_signal(&workers[t].cond_work);
			pthread_mutex_unlock(&workers[t].mutex);
		}

		// attente
		for (int t = 0; t < frames_read; t++) {
			pthread_mutex_lock(&workers[t].mutex);
			while (!workers[t].done)
				pthread_cond_wait(&workers[t].cond_done, &workers[t].mutex);
			pthread_mutex_unlock(&workers[t].mutex);
		}

		// écriture
		if (!isatty(STDOUT_FILENO)) {
			for (int t = 0; t < frames_read; t++)
				fwrite(out_buf[t], out_buf_length, 1, stdout);
			fflush(stdout);
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
	for(int t=0; t<NB_THREADS; t++)
	{
		ALIGNED_FREE(in_buf[t]);
		ALIGNED_FREE(out_buf[t]);
		ALIGNED_FREE(tmp_buf[t]);
	}
	
	//Close file
	if (input && (input != stdin))
	{
		fclose(input);
	}

	return 0;
}