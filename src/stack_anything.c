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
		"Vrunk11 toolkit, dirty_cvbs a simple adding phase error into cvbs video\n\n"
		"Usage:\n"
		"\t[--Fx input file number x where (0 < X < 33)\n"
		"\t[-o output file (use '-' to read from stdin)\n"
		"\t[-b bit depth (8 or 16) default = 16\n"
		"\t[--mode stacking mode 1: mean 2: median 3: smart\n"
		"\t[-t treshold for the smart mode (0 <-> 128) (median - x) < value < (median + x) \n"
	);
	exit(1);
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

int cmpfunc(const void* a, const void* b)
{
    return (*(int*)a - *(int*)b);
}

// Function for calculating median
double median16(unsigned short a[], int n)
{
    // First we sort the array
    qsort(a, n, sizeof(short), cmpfunc);
 
    // check for even case
    if (n % 2 != 0)
        return a[n / 2];
 
    return (a[(n - 1) / 2] + a[n / 2]) / 2.0;
}

double median8(unsigned char a[], int n)
{
    // First we sort the array
    qsort(a, n, sizeof(char), cmpfunc);
 
    // check for even case
    if (n % 2 != 0)
        return a[n / 2];
 
    return (a[(n - 1) / 2] + a[n / 2]) / 2.0;
}

int process_files(FILE *F1[32], int nb_files, unsigned short *out_buf, unsigned int buf_size, int bitdepth, int mode, int treshold)
{
	unsigned int y = 0;
	unsigned long calc = 0;	
	double median = 0;
	int nb_select = 0;
	
	unsigned short *tmp_w_buf = malloc(buf_size);//writing buffer
	
	//8bit cast
	unsigned char *tmp_w_buf8 = (void*)tmp_w_buf;//8bit writing buffer
	
	unsigned short *value16 = malloc(nb_files*2);
	unsigned char *value8 = malloc(nb_files);
	
	if(bitdepth == 16)//16bit
	{
		while(y < (buf_size/2))
		{
			calc = 0;
			if(mode == 0)
			{
				for(int i = 0;i < nb_files;i++)
				{
					fread(&value16[i],2,1,F1[i]);
					if(feof(F1[i]))
					{
						fprintf(stderr, "Reached end of file %d \n", i);
						return 1;
					}
					else
					{
						calc += value16[i];
					}
				}
				tmp_w_buf[y] = (calc/nb_files);
			}
			else if(mode == 1)
			{
				for(int i = 0;i < nb_files;i++)
				{
					fread(&value16[i],2,1,F1[i]);
					tmp_w_buf[y] = median16(value16,nb_files);
				}
			}
			else if(mode == 2)
			{
				nb_select = 0;
				for(int i = 0;i < nb_files;i++)
				{
					fread(&value16[i],2,1,F1[i]);
					if(feof(F1[i]))
					{
						fprintf(stderr, "Reached end of file %d \n", i);
						return 1;
					}
				}
				median = median16(value16,nb_files);
				for(int i = 0;i < nb_files;i++)
				{
					if(feof(F1[i]))
					{
						fprintf(stderr, "Reached end of file %d \n", i);
						return 1;
					}
					else
					{
						if((value16[i] < (median + (treshold * 256))) && (value16[i] > (median - treshold * 256)))//255 * 4 = 1020 (16bit)
						{
							nb_select++;
							calc += value16[i];
						}
					}
				}
				if(nb_select > 0)
				{
					tmp_w_buf[y] = (calc/nb_select);
				}
				else
				{
					tmp_w_buf[y] = median;
				}
				//fprintf(stderr, "selected %d sample \n", nb_select);
			}
			y++;
		}
	}
	else//8bit
	{
		while(y < (buf_size))
		{
			calc = 0;
			if(mode == 0)
			{
				for(int i = 0;i < nb_files;i++)
				{
					fread(&value8[i],1,1,F1[i]);
					if(feof(F1[i]))
					{
						fprintf(stderr, "Reached end of file %d \n", i);
						return 1;
					}
					else
					{
						calc += value8[i];
					}
				}
				tmp_w_buf[y] = (calc/nb_files);
			}
			else if(mode == 1)
			{
				for(int i = 0;i < nb_files;i++)
				{
					fread(&value8[i],1,1,F1[i]);
					tmp_w_buf[y] = median8(value8,nb_files);
				}
			}
			else if(mode == 2)
			{
				nb_select = 0;
				for(int i = 0;i < nb_files;i++)
				{
					fread(&value8[i],1,1,F1[i]);
					if(feof(F1[i]))
					{
						fprintf(stderr, "Reached end of file %d \n", i);
						return 1;
					}
				}
				median = median8(value8,nb_files);
				for(int i = 0;i < nb_files;i++)
				{
					if(feof(F1[i]))
					{
						fprintf(stderr, "Reached end of file %d \n", i);
						return 1;
					}
					else
					{
						if((value8[i] < (median + treshold)) && (value8[i] > (median - treshold)))//255 * 4 = 1020 (16bit)
						{
							nb_select++;
							calc += value8[i];
						}
					}
				}
				if(nb_select > 0)
				{
					tmp_w_buf[y] = (calc/nb_select);
				}
				else
				{
					tmp_w_buf[y] = median;
				}
				//fprintf(stderr, "selected %d sample \n", nb_select);
			}
			y++;
		}
		
	}
	
	memcpy(out_buf, tmp_w_buf, buf_size);
	free(value16);
	free(value8);
	free(tmp_w_buf);
	return 0;

}

int main(int argc, char **argv)
{
#ifdef _WIN32 || _WIN64
	_setmode(_fileno(stdout), O_BINARY);
	_setmode(_fileno(stdin), O_BINARY);	
#endif

	int opt;

	//file adress
	char *input_name[32] = { NULL };
	char *output_name_f = NULL;
	
	int mode = 0;
	int treshold = 4;
	
	//input file
	FILE *input[32] = { NULL };
	FILE *output_f;
	
	int bitdepth = 16;//default 16bit
	
	int nb_files = 0;
	
	unsigned short *buf = NULL;
	unsigned int buf_length = 0;
	
	int option_index = 0;
	static struct option long_options[] = {
		{"F1", 1, 0, 1},
		{"F2", 1, 0, 2},
		{"F3", 1, 0, 3},
		{"F4", 1, 0, 4},
		{"F5", 1, 0, 5},
		{"F6", 1, 0, 6},
		{"F7", 1, 0, 7},
		{"F8", 1, 0, 8},
		{"F9", 1, 0, 9},
		{"F10", 1, 0, 10},
		{"F11", 1, 0, 11},
		{"F12", 1, 0, 12},
		{"F13", 1, 0, 13},
		{"F14", 1, 0, 14},
		{"F15", 1, 0, 15},
		{"F16", 1, 0, 16},
		{"F17", 1, 0, 17},
		{"F18", 1, 0, 18},
		{"F19", 1, 0, 19},
		{"F20", 1, 0, 20},
		{"F21", 1, 0, 21},
		{"F22", 1, 0, 22},
		{"F23", 1, 0, 23},
		{"F24", 1, 0, 24},
		{"F25", 1, 0, 25},
		{"F26", 1, 0, 26},
		{"F27", 1, 0, 27},
		{"F28", 1, 0, 28},
		{"F29", 1, 0, 29},
		{"F30", 1, 0, 30},
		{"F31", 1, 0, 31},
		{"F32", 1, 0, 32},
		{"mode", 1, 0, 33},
		{0, 0, 0, 0}//reminder : letter value are from 65 to 122
	};

	while ((opt = getopt_long_only(argc, argv, "b:o:t:",long_options, &option_index)) != -1) {
		switch (opt) {
		case 'b':
			bitdepth = atoi(optarg);
			break;
		case 't':
			treshold = atoi(optarg);
			break;
		case 1:
			input_name[0] = optarg;
			nb_files++;
			break;
		case 2:
			input_name[1] = optarg;
			nb_files++;
			break;
		case 3:
			input_name[2] = optarg;
			nb_files++;
			break;
		case 4:
			input_name[3] = optarg;
			nb_files++;
			break;
		case 5:
			input_name[4] = optarg;
			nb_files++;
			break;
		case 6:
			input_name[5] = optarg;
			nb_files++;
			break;
		case 7:
			input_name[6] = optarg;
			nb_files++;
			break;
		case 8:
			input_name[7] = optarg;
			nb_files++;
			break;
		case 9:
			input_name[8] = optarg;
			nb_files++;
			break;
		case 10:
			input_name[9] = optarg;
			nb_files++;
			break;
		case 11:
			input_name[10] = optarg;
			nb_files++;
			break;
		case 12:
			input_name[11] = optarg;
			nb_files++;
			break;
		case 13:
			input_name[12] = optarg;
			nb_files++;
			break;
		case 14:
			input_name[13] = optarg;
			nb_files++;
			break;
		case 15:
			input_name[14] = optarg;
			nb_files++;
			break;
		case 16:
			input_name[15] = optarg;
			nb_files++;
			break;
		case 17:
			input_name[16] = optarg;
			nb_files++;
			break;
		case 18:
			input_name[17] = optarg;
			nb_files++;
			break;
		case 19:
			input_name[18] = optarg;
			nb_files++;
			break;
		case 20:
			input_name[19] = optarg;
			nb_files++;
			break;
		case 21:
			input_name[20] = optarg;
			nb_files++;
			break;
		case 22:
			input_name[21] = optarg;
			nb_files++;
			break;
		case 23:
			input_name[22] = optarg;
			nb_files++;
			break;
		case 24:
			input_name[23] = optarg;
			nb_files++;
			break;
		case 25:
			input_name[24] = optarg;
			nb_files++;
			break;
		case 26:
			input_name[25] = optarg;
			nb_files++;
			break;
		case 27:
			input_name[26] = optarg;
			nb_files++;
			break;
		case 28:
			input_name[27] = optarg;
			nb_files++;
			break;
		case 29:
			input_name[28] = optarg;
			nb_files++;
			break;
		case 30:
			input_name[29] = optarg;
			nb_files++;
			break;
		case 31:
			input_name[30] = optarg;
			nb_files++;
			break;
		case 32:
			input_name[31] = optarg;
			nb_files++;
			break;
		case 33:
			mode = (int)atoi(optarg);
			break;
		case 'o':
			output_name_f = optarg;
			break;
		default:
			usage();
			break;
		}
	}
	
	//check mode
	if(mode < 0 || mode > 2)
	{
		mode = 0;
		fprintf(stderr, "unknown mode '%d' using 0 instead\n", mode);
	}
	
	//check treshold
	if(treshold < 0 || treshold > 128)
	{
		mode = 4;
		fprintf(stderr, "Out of range value '%d' using 4 instead\n", mode);
	}
	
	//check bitdepth
	if(bitdepth != 16 && bitdepth != 8)
	{
		fprintf(stderr, "unknown value '%d' for option -b / supported value = (16,8)\n", bitdepth);
		usage();
	}
	
	int y = 0;
	for(int i=0;i < 32;i++)
	{
		if(input_name[i] != NULL)
		{
			input[y] = fopen(input_name[i], "rb");
			if (!input[i]) {
				fprintf(stderr, "(%d) : Failed to open %s\n", i, input[i]);
				return -ENOENT;
			}
			else
			{
				y++;
			}
		}
	}
	y=0;
	
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
	
	buf_length = bitdepth == 16 ? 20000 : 10000;//*2 for 16bit
	buf = malloc(buf_length);
	if(buf == NULL)//check allocation error
	{
		free(buf);
		fprintf(stderr, "malloc error (buf)\n");
		return -1;
	}
	
	while(!process_files(input,nb_files,buf,buf_length,bitdepth,mode,treshold))
	{
		//write output
		fwrite(buf, buf_length, 1, output_f);
		fflush(output_f);
	}

////ending of the program
	
	free(buf);
	
	//Close files
	for(int i = 0;i < nb_files;i++)
	{
		fclose(input[i]);
	}
	
	//Close out file
	if (output_f && (output_f != stdout))
	{
		fclose(output_f);
	}

	return 0;
}
