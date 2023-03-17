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

void usage(void)
{
	fprintf(stderr,
		"Vrunk11 toolkit, yc2cvbs a simple program for merging y/c into cvbs\n\n"
		"Usage:\n"
		"\t[-l luma file (use '-' to read from stdin)\n"
		"\t[-c chroma file (use '-' to read from stdin)\n"
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

void process_files(FILE *F1,FILE *F2,unsigned short *out_buf,unsigned int buf_size,int bitdepth)
{
	unsigned int i = 0;
	long calc = 0;	
	
	unsigned short *tmp_r_buf = malloc(buf_size);//luma input
	unsigned short *tmp_r2_buf = malloc(buf_size);//chroma input
	unsigned short *tmp_w_buf = malloc(buf_size);//writing buffer
	
	/*//16bit cast
	short *tmp_r_buf16 = (void*)tmp_r_buf;//16bit luma input
	short *tmp_r2_buf16 = (void*)tmp_r_buf;//16bit chroma input*/
	
	//8bit cast
	unsigned char *tmp_r_buf8 = (void*)tmp_r_buf;//8bit luma input
	unsigned char *tmp_r2_buf8 = (void*)tmp_r_buf;//8bit chroma input
	unsigned char *tmp_w_buf8 = (void*)tmp_w_buf;//8bit writing buffer
	
	unsigned short value16 = 0;
	unsigned short value16_2 = 0;
	unsigned char value8 = 0;
	unsigned char value8_2 = 0;
	
	short *value16_signed = (void *)&value16;
	short *value16_2_signed = (void *)&value16_2;
	char *value8_signed = (void *)&value8;
	char *value8_2_signed = (void *)&value8_2;
	
	fread(tmp_r_buf,buf_size,1,F1);
	fread(tmp_r2_buf,buf_size,1,F2);
	
	if(bitdepth == 16)//16bit
	{
		while(i < (buf_size/2))
		{
			value16 = tmp_r_buf[i];
			value16_2 = tmp_r2_buf[i];
			/*if((tmp_r_buf16[i] + tmp_r2_buf16[i] + 32768) > 65535)
			{
				tmp_w_buf[i] = 65535;
			}
			else if((tmp_r_buf16[i] + tmp_r2_buf16[i] + 32768) < 1024)
			{
				tmp_w_buf[i] = 1024;
			}*/
			tmp_w_buf[i] = *value16_signed + *value16_2_signed + 32768;//add luma + chroma 
			i++;
		}
	}
	else//8bit
	{
		while(i < (buf_size))
		{
			value8 = tmp_r_buf8[i];
			value8_2 = tmp_r2_buf8[i];
			/*if((tmp_r_buf8[i] + tmp_r2_buf8[i] + 32768) > 255)
			{
				tmp_w_buf8[i] = 255;
			}
			else if((tmp_r_buf8[i] + tmp_r2_buf8[i] + 32768) < 3)
			{
				tmp_w_buf8[i] = 3;
			}*/
			tmp_w_buf8[i] = *value8_signed + *value8_2_signed + 128;//add luma + chroma
			i++;
		}
	}
	
	memcpy(out_buf, tmp_w_buf, buf_size);
	
	free(tmp_r_buf);
	free(tmp_r2_buf);
	free(tmp_w_buf);

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
	char *input_name_2 = NULL;
	
	//input file
	FILE *input_1;
	FILE *input_2;
	
	char standard = 'N';//default ntsc
	int bitdepth = 16;//default 16bit
	
	unsigned short *buf = NULL;
	unsigned int buf_length = 0;

	while ((opt = getopt(argc, argv, "s:b:l:c:")) != -1) {
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
		case 'l':
			input_name_1 = optarg;
			break;
		case 'c':
			input_name_2 = optarg;
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
	
	//reading file 2
	if (strcmp(input_name_2, "-") == 0)// Read samples from stdin
	{
		input_2 = stdin;
	}
	else
	{
		input_2 = fopen(input_name_2, "rb");
		if (!input_2) {
			fprintf(stderr, "(2) : Failed to open %s\n", input_2);
			return -ENOENT;
		}
	}
	
	if(standard == 'N' || standard == 'P')
	{
		if(standard == 'N')
		{
			buf_length = bitdepth == 16 ? 478660*2 : 478660;//(910 * 526)*2 for 16bit / (910 * 526) for 8bit
		}
		else if(standard == 'P')
		{
			buf_length = bitdepth == 16 ? 710510*2 : 710510;//(1135 * 626)*2 for 16bit / (1135 * 626) for 8bit
		}
		
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
	
	while(!feof(input_1) || !feof(input_2))
	{
		process_files(input_1,input_2,buf,buf_length,bitdepth);
	
		//write to stdout (pipe)
		if(isatty(STDOUT_FILENO) == 0)
		{
			fwrite(buf, buf_length,1,stdout);
			fflush(stdout);
		}
	}

////ending of the program
	
	free(buf);
	
	//Close file 1
	if (input_1 && (input_1 != stdin))
	{
		fclose(input_1);
	}
	
	//Close file 2
	if (input_2 && (input_2 != stdin))
	{
		fclose(input_2);
	}

	return 0;
}
