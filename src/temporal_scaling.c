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

#define _FILE_OFFSET_BITS 64

#ifdef _WIN64
#define FSEEK fseeko64
#else
#define FSEEK fseeko
#endif

#define NB_THREADS 4

typedef struct BresenhamMap BresenhamMap;
struct BresenhamMap {
    int *w_src;   // index source
    int *w_flag;  // flag bilinear
    int *h_src;
    int *h_flag;
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
        "\t                     2 = bilinear with temporal alternation\n\n"
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

void process_files(unsigned char *in_buf, unsigned char *out_buf, unsigned char *tmp_buf, BresenhamMap *map, int in_w, int in_h, int out_w, int out_h, int ratio, int mode)
{
	const int in_w3    = in_w * 3;
	const int out_w3   = out_w * 3;
	const int bilinear = (mode >= 1);
    unsigned int offset = 0;
    for(int r=0; r<ratio; r++)
    {
        //scale the width line by line
        const int map_w_base = r * out_w;  // hissé hors des boucles

		for(int h=0; h<in_h; h++)
		{
			const unsigned char *src_row = in_buf  + h * in_w3;   // calculé une fois par ligne
			unsigned char       *dst_row = tmp_buf + h * out_w3;

			int px3 = 0;
			for(int px=0; px<out_w; px++, px3+=3)  // px3 remplace px/3
			{
				int src_px     = map->w_src[map_w_base + px] * 3;
				int is_bilinear = map->w_flag[map_w_base + px];

				if(bilinear && is_bilinear == 1 && src_px + 3 < in_w3)
				{
					dst_row[px3]   = (src_row[src_px]   + src_row[src_px+3]) >> 1;
					dst_row[px3+1] = (src_row[src_px+1] + src_row[src_px+4]) >> 1;
					dst_row[px3+2] = (src_row[src_px+2] + src_row[src_px+5]) >> 1;
				}
				else if(bilinear && is_bilinear == -1 && src_px >= 3)
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
			int src_line   = map->h_src[map_h_base + px];
			int is_bilinear = map->h_flag[map_h_base + px];

			unsigned char       *dst_row  = out_buf + offset + px * out_w3;
			const unsigned char *tmp_row0 = tmp_buf + src_line * out_w3;

			if(bilinear && is_bilinear == 1 && src_line + 1 < in_h)
			{
				const unsigned char *tmp_row1 = tmp_row0 + out_w3;
				for(int i=0; i<out_w3; i++)
					dst_row[i] = (tmp_row0[i] + tmp_row1[i]) >> 1;  // ← SIMD-vectorisable
			}
			else if(bilinear && is_bilinear == -1 && src_line >= 1)
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
	map.w_src = NULL;
	map.w_flag = NULL;
	map.h_src = NULL;
	map.h_flag = NULL;
	map.w_length = 0;
	map.h_length = 0;
	
	int option_index = 0;
	static struct option long_options[] = {
		{"in-width", 1, 0, 1},
		{"in-height", 1, 0, 2},
		{"out-width", 1, 0, 3},
		{"out-height", 1, 0, 4},
		{"fps-ratio", 1, 0, 5},
		{"mode", 1, 0, 6},
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
		in_buf_length = in_w*in_h*3;
		out_buf_length = out_w*out_h*3*frames_ratio;
		
		map.w_length = out_w * frames_ratio;
		map.h_length = out_h * frames_ratio;
		
		map.w_src  = malloc(map.w_length * sizeof(int));
		map.w_flag = malloc(map.w_length * sizeof(int));
		map.h_src  = malloc(map.h_length * sizeof(int));
		map.h_flag = malloc(map.h_length * sizeof(int));
		
		int alloc_err = 0;
		for(int t=0; t<NB_THREADS; t++)
		{
			in_buf[t] = malloc(in_buf_length);
			out_buf[t] = malloc(out_buf_length);
			tmp_buf[t] = malloc(out_w * in_h * 3);
			if(in_buf[t] == NULL || out_buf[t] == NULL || tmp_buf[t] == NULL)
				alloc_err = 1;
		}
		
		if(!map.w_src || !map.w_flag || !map.h_src || !map.h_flag || alloc_err)//check allocation error
		{
			free(map.w_src);
			free(map.w_flag);
			free(map.h_src);
			free(map.h_flag);
			
			for(int t=0; t<NB_THREADS; t++)
			{
				free(in_buf[t]);
				free(out_buf[t]);
				free(tmp_buf[t]);
			}
			fprintf(stderr, "malloc error (buf)\n");
			return -1;
		}
	}
	else//error
	{
		fprintf(stderr, "invalid size error\n");
		return -1;
	}
	
	if(mode < 0 || mode > 2)
	{
		fprintf(stderr, "Error : mode '%d' is not supported\n", mode);
		return -1;
	}
	
	//compute the maping used for scaling the input image into an output one
	compute_scale_map(&map, in_w, in_h, out_w, out_h, frames_ratio);
	
	if(mode == 2)
	{
		adjust_bilinear_phases(&map, out_w, out_h, in_w, in_h, frames_ratio);
	}
	
	ThreadArgs args[NB_THREADS];
	for(int t=0; t<NB_THREADS; t++)
	{
		args[t].map = &map;
		args[t].in_w = in_w;
		args[t].in_h = in_h;
		args[t].out_w = out_w;
		args[t].out_h = out_h;
		args[t].ratio = frames_ratio;
		args[t].mode = mode;
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
	
	free(map.w_src);
	free(map.w_flag);
	free(map.h_src);
	free(map.h_flag);
	for(int t=0; t<NB_THREADS; t++)
	{
		free(in_buf[t]);
		free(out_buf[t]);
		free(tmp_buf[t]);
	}
	
	//Close file
	if (input && (input != stdin))
	{
		fclose(input);
	}

	return 0;
}