#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include <unistd.h>
#include <stdint.h>

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

int freez = 0;
int pal_offset = 1;

void usage(void)
{
	fprintf(stderr,
		"Vrunk11 toolkit, cvbs decode DdD is a simple program for align and correcting DdD cvbs sample\n\n"
		"Usage:\n"
		"\t[-i input file (use '-' to read from stdin)\n"
		"\t[-o output file (use '-' to read from stdin)\n"
		"\t[-f freeze the stabilisation\n"
		"\t[-h hsync precision [1 - 525]\n"
		"\t[-s video standard (pal / ntsc) default = ntsc\n"
	);
	exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

void process_files(FILE *F1,unsigned short *out_buf,unsigned int buf_size,char std,int hsync_nb_line)
{
	unsigned int i = 0;
	unsigned int y = 0;
	long calc = 0;	
	
	unsigned short *tmp_r_buf = malloc(buf_size);//luma input
	unsigned short *tmp_w_buf = malloc(buf_size);//writing buffer
	
	unsigned short value16 = 0;
	
	short *value16_signed = (void *)&value16;
	
	int hsync_pos = 0;
	int vsync_pos = 0;
	int line_length = 0;
	int hsync_size = 0;
	int vsync_size = 0;
	int vsync_gap = 0;
	
	uint32_t tmp_sum = 0;           //variable used to sum
	uint32_t tmp_min = 2524971008;//lowest score
	
	if(std == 'N')
	{
		line_length = 910;
	}
	else
	{
		line_length = 1135;
		hsync_size = 80;
		vsync_size = 39;
		vsync_gap = 567;
	}
	
	fread(tmp_r_buf,buf_size,1,F1);

//alignment		
	if(freez != 1)
	{	
		//hsync detection
		while((y + hsync_size) < line_length)
		{
			tmp_sum = 0;
			for(int cnt=0;cnt < hsync_size;cnt++)
			{
				//start position (y) + cnt
				for(int line=0;line < hsync_nb_line;line++)
				{
					tmp_sum += tmp_r_buf[y+cnt+(line_length*line)];
				}
			}
			if(tmp_sum < tmp_min)
			{
				tmp_min = tmp_sum;
				hsync_pos = y;
			}
			y++;
		}
		
		y=0;
		tmp_min = 2524971008;//big value by default
		
		//vsync detection
		while(i < (((buf_size/2)) - line_length*7))//scaning only 1 field
		{
			tmp_sum = 0;
			//summ of all sample from the vsync rectangle
			for(int cnt=0;cnt < vsync_size;cnt++)
			{
				tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*0)];
				tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*1)];
				tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*2)];
				//tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*3)];
				//tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*4)];
				//tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*5)];
				tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*6)];
				tmp_sum += tmp_r_buf[i+hsync_pos+vsync_gap+cnt+(line_length*7)];
			}
			if(tmp_sum < tmp_min)
			{
				tmp_min = tmp_sum;
				vsync_pos = i;
			}
			i += line_length;
		}
		if(freez == 2)
		{
			freez = 1;
		}
	}
//output preparation (lire pour compenser)
	i = 0;//output position
	y = ((vsync_pos + hsync_pos) - 2);
	
	if(std == 'P' && pal_offset == 1)
	{
		y += (line_length*3);	
		pal_offset = 0;
	}
	
	while(y < (buf_size/2))
	{
		tmp_w_buf[i] = tmp_r_buf[y];
		i++;
		y++;
	}
	y = 0;
	if(i < (buf_size/2))
	{
		fread(tmp_r_buf,((buf_size)-(i*2)),1,F1);
		while(i < (buf_size/2))
		{
			tmp_w_buf[i] = tmp_r_buf[y];
			i++;
			y++;
		}
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
	char *output_name_f = NULL;
	
	//input file
	FILE *input_1;
	FILE *output_f;
	
	char standard = 'N';//default ntsc
	int bitdepth = 16;
	int hsync_nb_line = 525;
	
	unsigned short *buf = NULL;
	unsigned int buf_length = 0;

	while ((opt = getopt(argc, argv, "s:h:i:o:f")) != -1) {
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
		case 'h':
			hsync_nb_line = atoi(optarg);
			if(hsync_nb_line > 525)
			{
				hsync_nb_line = 525;
			}
			else if(hsync_nb_line < 1)
			{
				hsync_nb_line = 1;
			}
			break;
		case 'i':
			input_name_1 = optarg;
			break;
		case 'o':
			output_name_f = optarg;
			break;
		case 'f':
			freez = 2;
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
			buf_length = bitdepth == 16 ? 477750*2 : 477750;//(910 * 526)*2 for 16bit / (910 * 526) for 8bit
		}
		else if(standard == 'P')
		{
			buf_length = bitdepth == 16 ? 709375*2 : 709375;//(1135 * 626)*2 for 16bit / (1135 * 626) for 8bit
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
	
	while(!feof(input_1))
	{
		process_files(input_1,buf,buf_length,standard,hsync_nb_line);
	
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
