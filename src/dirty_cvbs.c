#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

void usage(void)
{
	fprintf(stderr,
		"Vrunk11 toolkit, dirty_cvbs a simple adding phase error into cvbs video\n\n"
		"Usage:\n"
		"\t[-i input file (use '-' to read from stdin)\n"
		"\t[-o output file (use '-' to read from stdin)\n"
		"\t[-b bit depth (8 or 16) default = 16\n"
		"\t[-s video standard (pal / ntsc) default = ntsc\n"
	);
	exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

void process_files(FILE *F1,unsigned short *out_buf,unsigned int buf_size,int bitdepth)
{
	unsigned int i = 0;
	int random = 0;
	
	unsigned short *tmp_r_buf = malloc(buf_size);//cvbs input
	
	//8bit cast
	unsigned char *tmp_r_buf8 = (void*)tmp_r_buf;//8bit luma input
	unsigned char *out_buf8 = (void*)out_buf;//8bit writing buffer
	
	fread(tmp_r_buf,buf_size,1,F1);
	
	if(bitdepth == 16)//16bit
	{
		while(i < (buf_size/2))
		{
			random = (rand() % 20) - 10;
			if(random > 0)
			{
				if(i <= random)
				{
					out_buf[i] = 0;
				}
				else
				{
					out_buf[i] = tmp_r_buf[i - random];
				}
			}
			else
			{
				if(i >= ((buf_size/2) + random))
				{
					out_buf[i] = 0;
				}
				else
				{
					out_buf[i] = tmp_r_buf[i - random];
				}
			}
			i++;
		}
	}
	else//8bit
	{
		while(i < (buf_size))
		{
			random = (rand() % 20) - 10;
			if(random > 0)
			{
				if(i <= random)
				{
					out_buf8[i] = 0;
				}
				else
				{
					out_buf8[i] = tmp_r_buf8[i - random];
				}
			}
			else
			{
				if(i >= ((buf_size/2) + random))
				{
					out_buf8[i] = 0;
				}
				else
				{
					out_buf8[i] = tmp_r_buf8[i - random];
				}
			}
			i++;
		}
	}
	
	free(tmp_r_buf);
}

int main(int argc, char **argv)
{
#ifdef _WIN32 || _WIN64
	_setmode(_fileno(stdout), O_BINARY);
	_setmode(_fileno(stdin), O_BINARY);	
#endif

	int opt;

	//file adress
	char *input_name_1 = NULL;
	char *output_name_f = NULL;
	
	//input file
	FILE *input_1;
	FILE *output_f;
	
	char standard = 'N';//default ntsc
	int bitdepth = 16;//default 16bit
	
	unsigned short *buf = NULL;
	unsigned int buf_length = 0;
	
	srand(1);

	while ((opt = getopt(argc, argv, "s:b:i:o:")) != -1) {
		switch (opt) {
		case 's':
			if((strcmp(optarg, "ntsc" ) == 0) || (strcmp(optarg, "NTSC" ) == 0) || (strcmp(optarg, "Ntsc" ) == 0) || (strcmp(optarg, "N" ) == 0) || (strcmp(optarg, "n" ) == 0))
			{
				standard = 'N';
			}
			else if((strcmp(optarg, "pal" ) == 0) || (strcmp(optarg, "PAL" ) == 0) || (strcmp(optarg, "Pal" ) == 0) || (strcmp(optarg, "P" ) == 0) || (strcmp(optarg, "p" ) == 0))
			{
				standard = 'P';
			}
			else
			{
				standard = 'E';//error
			}
			break;
		case 'b':
			bitdepth = atoi(optarg);
			break;
		case 'i':
			input_name_1 = optarg;
			break;
		case 'o':
			output_name_f = optarg;
			break;
		default:
			usage();
			break;
		}
	}
	
	//check bitdepth
	if(bitdepth != 16 && bitdepth != 8)
	{
		fprintf(stderr, "unknown value '%d' for option -b / supported value = (16,8)\n", bitdepth);
		usage();
	}
	
	//reading file 1
	if (strcmp(input_name_1, "-") == 0)// Read samples from stdin
	{
		input_1 = stdin;
	}
	else
	{
		input_1 = fopen(input_name_1, "rb");
		if (!input_1) {
			fprintf(stderr, "(1) : Failed to open %s\n", input_1);
			return -ENOENT;
		}
	}
	
	//opening output file
	if (strcmp(output_name_f, "-") == 0)// Read samples from stdin
	{
		output_f = stdout;
	}
	else
	{
		output_f = fopen(output_name_f, "wb");
		if (!output_f) {
			fprintf(stderr, "(2) : Failed to open %s\n", output_f);
			return -ENOENT;
		}
	}
	
	if(standard == 'N' || standard == 'P')
	{
		if(standard == 'N')
		{
			buf_length = bitdepth == 16 ? 910*2 : 478660;//(910 * 526)*2 for 16bit / (910 * 526) for 8bit
		}
		else if(standard == 'P')
		{
			buf_length = bitdepth == 16 ? 1135*2 : 710510;//(1135 * 626)*2 for 16bit / (1135 * 626) for 8bit
		}
		
		buf_length = (buf_length *10);
		
		buf = malloc(buf_length);
		if(buf == NULL)//check allocation error
		{
			free(buf);
			fprintf(stderr, "malloc error (buf)\n");
			return -1;
		}
	}
	else//error
	{
		fprintf(stderr, "unknow standard error / supported value = (pal,PAL,Pal,p,P) (ntsc,NTSC,Ntsc,n,N)\n");
		return -1;
	}
	
	while(!feof(input_1))
	{
		process_files(input_1,buf,buf_length,bitdepth);
	
		//write output
		fwrite(buf, buf_length,1,output_f);
		fflush(output_f);
	}

////ending of the program
	
	free(buf);
	
	//Close file 1
	if (input_1 && (input_1 != stdin))
	{
		fclose(input_1);
	}
	
	//Close out file
	if (output_f && (output_f != stdout))
	{
		fclose(output_f);
	}

	return 0;
}
