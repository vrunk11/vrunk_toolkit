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
		"Vrunk11 toolkit, a simple program for process video\n\n"
		"Usage:\n"
		"\t[-i1 filename 1 (use '-' to read from stdin)\n"
		"\t[-i2 filename 2 (use '-' to read from stdin)\n"
		"\t[-i3 filename 3 (use '-' to read from stdin)\n"
	);
	exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

void process_files(FILE *F1,FILE *F2,FILE *F3,unsigned short *out_buf,unsigned int buf_size,char standard)
{
	unsigned int y = 0;
	unsigned int i = 0;
	int black = 0;
	int white = 0;
	
	if(standard == "N")//NTSC
	{
		black = 18048;
		white = 51200;
	}
	else//PAL
	{
		black = 16384;
		white = 54016;
	}

	double VideoRange = (white - black);
	double Yscale = (56064.0 / VideoRange);
	double Uscale = (28672 / (0.886 * 0.49211104112248356308804691718185)) / VideoRange;
	double Vscale = (28672 / (0.701 * 0.87728321993817866838972487283129)) / VideoRange;
	
	long calc = 0;	
	
	unsigned short *tmp_r_buf = malloc(buf_size/3);//reading
	//short *tmp_r_buf_u = (void*)tmp_r_buf;//reading
	unsigned short *tmp_w_buf = malloc(buf_size);//writing
	
	fread(tmp_r_buf,buf_size/3,1,F1);
	
	while(i < (buf_size/2)/3)//Y
	{
		if(((tmp_r_buf[i] - 18048) * Yscale) + 4096 >= 65216)//clipping white
		{
			tmp_w_buf[y] = 65216;
		}
		else if(((tmp_r_buf[i] - 18048) * Yscale) + 4096 <= 256)//clipping black
		{
			tmp_w_buf[y] = 256;
		}
		else
		{
			tmp_w_buf[y] = ((tmp_r_buf[i] - 18048) * Yscale) + 4096;
		}
		i++;
		y++;
	}
	i = 0;
	fread(tmp_r_buf,buf_size/3,1,F2);
	 
	while(i < (buf_size/2)/3)//U
	{
		calc = tmp_r_buf[i] - 32767;
		//tmp_r_buf[i] = tmp_r_buf[i] - 32767;
		/*if(((calc * Uscale) + 4096) >= 65216)
		{
			tmp_w_buf[y] = 65216;
		}
		else if(((calc * Uscale) + 4096) <= 32768)
		{
			tmp_w_buf[y] = 32768;
		}
		else
		{*/
			tmp_w_buf[y] = ((calc * Uscale))+32767;
		//}
		i++;
		y++;
	}
	i = 0;
	fread(tmp_r_buf,buf_size/3,1,F3);
	
	while(i < (buf_size/2)/3)//V
	{
		calc = tmp_r_buf[i] - 32767;
		//tmp_r_buf[i] = tmp_r_buf[i] - 32767;
		/*if(((calc * Vscale) + 4096) >= 65216)
		{
			tmp_w_buf[y] = 65216;
		}
		else if(((calc * Vscale) + 4096) <= 32768)
		{
			tmp_w_buf[y] = 32768;
		}
		else
		{*/
			tmp_w_buf[y] = ((calc * Vscale))+32767;
		//}
		i++;
		y++;
	}
	
	memcpy(out_buf, tmp_w_buf, buf_size);
	
	free(tmp_r_buf);
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
	char *input_name_3 = NULL;
	
	//input file
	FILE *input_1;
	FILE *input_2;
	FILE *input_3;
	
	char standard = 'N';//default ntsc
	
	unsigned short *buf = NULL;
	unsigned int buf_length = 0;
	
	int option_index = 0;
	static struct option long_options[] = {
		{"i1", 1, 0, 1},
		{"i2", 1, 0, 2},
		{"i3", 1, 0, 3},
		{0, 0, 0, 0}//reminder : letter value are from 65 to 122
	};

	while ((opt = getopt_long_only(argc, argv, "s:",long_options, &option_index)) != -1) {
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
		case 1:
			input_name_1 = optarg;
			break;
		case 2:
			input_name_2 = optarg;
			break;
		case 3:
			input_name_3 = optarg;
			break;
		default:
			usage();
			break;
		}
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
	
	//reading file 3
	if (strcmp(input_name_3, "-") == 0)// Read samples from stdin 
	{
		input_3 = stdin;
	}
	else
	{
		input_3 = fopen(input_name_3, "rb");
		if (!input_3) {
			fprintf(stderr, "(3) : Failed to open %s\n", input_3);
			return -ENOENT;
		}
	}
	
	if(standard == 'N' || standard == 'P')
	{
		if(standard == 'N')
		{
			buf_length = 2871960;//((910 * 526) * 3)*2
		}
		else if(standard == 'P')
		{
			buf_length = 4263060;//((1135 * 626) * 3)*2
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
		fprintf(stderr, "unknow standard error\n");
		return -1;
	}
	
	while(!feof(input_1) || !feof(input_2) || !feof(input_3))
	{
		process_files(input_1,input_2,input_3,buf,buf_length,standard);
	
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
	
	//Close file 3
	if (input_3 && (input_3 != stdin))
	{
		fclose(input_3);
	}

	return 0;
}
