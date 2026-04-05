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
struct BresenhamMap
{
	int *w;
	int *h;
	unsigned int w_length;
	unsigned int h_length;
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

void usage(void)
{
	fprintf(stderr,
		"Vrunk11 toolkit, a simple program for process video\n\n"
		"Usage:\n"
		"\t[-i filename (use '-' to read from stdin)\n"
	);
	exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
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
            map->w[dst++] = src;
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
            map->h[dst++] = src;
            err += in_h;
            if(err >= out_h) { err -= out_h; src++; }
        }
    }
}

void process_files(unsigned char *in_buf, unsigned char *out_buf, unsigned char *tmp_buf, BresenhamMap *map, int in_w, int in_h, int out_w, int out_h, int ratio, int mode)
{
    unsigned int offset = 0;

    for(int r=0; r<ratio; r++)
    {
        //scale the width line by line
        for(int h=0; h<in_h; h++)
        {
            for(int px=0; px<(out_w*3); px+=3)
            {
                int src_px = map->w[r*out_w + px/3] * 3;
                tmp_buf[(h*out_w*3)+px]   = in_buf[(h*in_w*3)+src_px];
                tmp_buf[(h*out_w*3)+px+1] = in_buf[(h*in_w*3)+src_px+1];
                tmp_buf[(h*out_w*3)+px+2] = in_buf[(h*in_w*3)+src_px+2];
            }
        }

        //scale the height column by column
        for(int w=0; w<out_w; w++)
        {
            for(int px=0; px<(out_h*3); px+=3)
            {
                int src_line = map->h[r*out_h + px/3];
                out_buf[offset+(px/3*out_w+w)*3]   = tmp_buf[(src_line*out_w+w)*3];
                out_buf[(offset+(px/3*out_w+w)*3)+1] = tmp_buf[((src_line*out_w+w)*3)+1];
                out_buf[(offset+(px/3*out_w+w)*3)+2] = tmp_buf[((src_line*out_w+w)*3)+2];
            }
        }
        offset += out_w*out_h*3;
    }
}

void *thread_func(void *arg)
{
	ThreadArgs *a = (ThreadArgs *)arg;
	process_files(a->in_buf, a->out_buf, a->tmp_buf, a->map, a->in_w, a->in_h, a->out_w, a->out_h, a->ratio, a->mode);
	return NULL;
}

int main(int argc, char **argv)
{
#ifdef _WIN32 || _WIN64
	_setmode(_fileno(stdout), O_BINARY);
	_setmode(_fileno(stdin), O_BINARY);	
#endif

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
	map.w = NULL;
	map.h = NULL;
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
		
		map.w = malloc(map.w_length * sizeof(int));
		map.h = malloc(map.h_length * sizeof(int));
		
		int alloc_err = 0;
		for(int t=0; t<NB_THREADS; t++)
		{
			in_buf[t] = malloc(in_buf_length);
			out_buf[t] = malloc(out_buf_length);
			tmp_buf[t] = malloc(out_w * in_h * 3);
			if(in_buf[t] == NULL || out_buf[t] == NULL || tmp_buf[t] == NULL)
				alloc_err = 1;
		}
		
		if(map.w == NULL || map.h == NULL || alloc_err)//check allocation error
		{
			free(map.w);
			free(map.h);
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
	
	if(mode < 0 || mode > 1)
	{
		fprintf(stderr, "Error : mode '%d' is not supported\n", mode);
		return -1;
	}
	
	//compute the maping used for scaling the input image into an output one
	compute_scale_map(&map, in_w, in_h, out_w, out_h, frames_ratio);
	
	pthread_t thread_1;
	pthread_t thread_2;
	pthread_t thread_3;
	pthread_t thread_4;
	
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
	
	while(!feof(input))
	{
		//read NB_THREADS frames
		int frames_read = 0;
		for(int t=0; t<NB_THREADS; t++)
		{
			if(fread(in_buf[t], in_buf_length, 1, input) == 1)
			{
				args[t].in_buf = in_buf[t];
				args[t].out_buf = out_buf[t];
				args[t].tmp_buf = tmp_buf[t];
				frames_read++;
			}
			else break;
		}
		
		//launch threads
		pthread_t threads[NB_THREADS];
		for(int t=0; t<frames_read; t++)
			pthread_create(&threads[t], NULL, thread_func, &args[t]);
		
		//wait for all threads
		for(int t=0; t<frames_read; t++)
			pthread_join(threads[t], NULL);
		
		//write to stdout (pipe)
		if(isatty(STDOUT_FILENO) == 0)
		{
			for(int t=0; t<frames_read; t++)
			{
				fwrite(out_buf[t], out_buf_length, 1, stdout);
			}
			fflush(stdout);
		}
	}

////ending of the program
	
	free(map.w);
	free(map.h);
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