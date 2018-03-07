#include <iostream>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "FileUtils.h"
//#include "ConfUtils.h"
#include "DataUtils.h"
#include "filter.h"
#include "genome_size.h"
extern "C"
{
#include "pnorm.h"
}

using namespace std;

void get_local_maxima(_mark_sig* mark_sig,int num,int half_size,int* list, int& list_num)
	// mark_sig - significant mark bins of a single marker 
	// num - number of significant bins
	// half_size - default 2500
	// list - index array for mark_sig
	// list_num - ?
	// all local maximum are stored in list, and the length of list in list_num
{

	int* identical = new int[2*half_size/100+1];
	if(!identical) return;

	list_num = 0;
	for(int i = 0; i < num;i++)
	{
		int pos = mark_sig[i].pos;
		double pos_z = mark_sig[i].z;

		int is_local_max = 1;
		int count = 0;
		memset(identical, 0, sizeof(int)*(2*half_size/100 + 1));

		//search left
		int index2 = i - 1;
		while(index2 > 0)
		{
			int j = mark_sig[index2].pos;
			double j_z = mark_sig[index2].z;
			if(j < pos - half_size) break;
			if(pos_z == j_z)
			{
				identical[count] = j / 100;
				count++;
			}
			if(j_z > pos_z)
			{
				is_local_max = 0;
				break;
			}
			index2--;
		}

		// search right
		if(is_local_max)
		{
			index2 = i + 1;
			while(index2 < num)
			{
				int j = mark_sig[index2].pos;
				double j_z = mark_sig[index2].z;
				if(j > pos + half_size) break;
				if(pos_z == j_z)
				{
					identical[count] = j / 100;
					count++;
				}
				if(j_z > pos_z)
				{
					is_local_max = 0;
					break;
				}
				index2++;
			}
		}

		// in case that z-score is identical, only consider local maximum when the bin is in the middle of the 5k window
		if(is_local_max && count > 0)
		{
			identical[count] = pos / 100;
			count++;
			qsort(identical,count,sizeof(int),cmp_int);
			int idx = count / 2;
			if(identical[idx] != pos/100) is_local_max = 0;
		}
		if(is_local_max)
		{
			list[list_num] = i;
			list_num++;
		}
	}
	delete[] identical;
}

// filter significant region on all mark level, and print them out to "filter.idx"
void print_non_overlapping(_sig* sig,int sig_num,int chr,int half_size,GS* gs)
{
	//string chrom = id2chr(chr);
	string chrom = gs->id_chr_map[chr-1];

	FileUtils file1(4096);
	int total = 0;
	//add an newline
	if(FileLength("filter.idx") > 0)
	{
		char buf[64];
		strcpy(buf,"");
		file1.WriteLine("filter.idx",buf,"a+");
	}

	int* identical = NULL;
	for(int i = 0; i < sig_num;i++)
	{
		//printf("index=%d, mark num=%d, pos=%d, z-value=%lf\n", i, sig[i].marks, sig[i].pos, sig[i].sum);
		if(sig[i].pos == 0) continue;

		int is_local_max = 1;
		int pos = sig[i].pos;
		if(identical==NULL) 
			identical = new int[2*half_size+1];
		int count = 0;

		int min = (pos - half_size) / 100;
		int max = (pos + half_size) / 100;
		if(min < 0) min = 0;
		if(max >= sig_num) max = sig_num;
		// overall marker, check enrich region at 10k
		for(int j = min;j <= max;j++)
		{
			if(pos/100 == j || sig[j].pos == 0) continue;

			if(sig[j].sum == sig[i].sum)
			{
				identical[count] = j;
				count++;
			}
			if(sig[j].sum > sig[i].sum)
			{
				is_local_max = 0;
				break;
			}
		}
		if(is_local_max && count > 0)
		{
			identical[count] = pos / 100;
			count++;
			qsort(identical,count,sizeof(int),cmp_int);
			int idx = count / 2;
			if(identical[idx] != pos/100) is_local_max = 0;
		}
		if(is_local_max)
		{
			char buf[256];
			sprintf(buf,"%s\t%d",chrom.c_str(),pos);
			file1.WriteLine("filter.idx",buf,"a+");
			total++;
		}
	}
	if(identical) delete[] identical;

	file1.LastWrite();
	file1.Release();
	printf("%s finished, total=%d\n",chrom.c_str(),total);
}

// type 2 core function, find significant region, filter overlap region on single mark level and all mark level
//void get_sig_locs(int stat_half_size, int maxima_half_size, int overlap_half_size, double pval_cutoff, _mark **mark, GS* gs)
void get_sig_locs(int stat_half_size, int maxima_half_size, int overlap_half_size, double pval_cutoff, GS* gs)
	// default value: 1000 2500 5000 1e-5
{
	_mark* mark;
	float* mark_buf[MAX_MARK]; // marker signal
	_mark_sig* mark_sig[MAX_MARK]; // significant marks, position and zscore
	int mark_sig_num[MAX_MARK]; // number of significant region for each marks
	int mem_n[MAX_MARK]; // size of mark_sig in unit of 20000

	int mark_num = load_marks("mark_.para",&mark); // load mark mean and std information to mark[]
	if(mark_num < 0 ) return;
	remove("filter.idx");

	int old_chr = -1;
	for(int i=0; i<gs->g_cnt; i++)
	{
		int min = 0; // start bin number for chr i
		int max = gs->g_size[i] / 100; // end bin number for chr i
		//int max = hg18[i] / 100; // end bin number for chr i

		memset(mark_sig_num,0,mark_num*sizeof(int)); // initialization for mark_sig_num

		// initialization for men_n
		// memory allocation for mark_sig
		for(int j = 0; j < mark_num; j++)
		{
			mem_n[j] = 1;
			mark_sig[j] = (_mark_sig*)malloc(sizeof(_mark_sig) * 20000);
		}
		alloc_mem(mark_buf, mark_num,i+1,old_chr,mark,gs); // load all marks of chromosome i+1 to mark_buf
		old_chr = i + 1;

		// find all significant marks for chr i
		for(int j = min; j < max; j++)
		{
			int min1 = j - stat_half_size / 100;
			int max1 = j + stat_half_size / 100;
			if(min1 < 0) min1 = 0;
			if(max1 >= max) max1 = max;
			for(int k = 0; k < mark_num; k++)
			{
				double val = 0;
				// sum overall 2000bp regions (21 bins)
				for(int m = min1; m < max1; m++) 
					val += (double)mark_buf[k][m];

				double z_value= (val - mark[k].mean) / mark[k].std; // z-value to the marker on all encode regions
				double cum = 0,ccum = 0;
				pnorm_both(z_value,&cum,&ccum,0,0); // get p-value for z-value
				double pval = 1.0 - cum;
				if(pval <= pval_cutoff)
				{
					int idx = mark_sig_num[k];
					if(idx >= mem_n[k] * 20000)
					{
						mem_n[k] += 1;
						mark_sig[k]=(_mark_sig*)realloc(mark_sig[k],mem_n[k]*20000*sizeof(_mark_sig));
						if(mark_sig[k] == NULL) exit(-1);
					}
					mark_sig[k][idx].pos = j*100;
					mark_sig[k][idx].z = z_value;
					mark_sig_num[k] += 1;
				}
			}
		}
		//second, find local maxima and remove all significant hits between maxima that are too close
		for(int j = 0; j < mark_num; j++)
		{
			int num = mark_sig_num[j];
			if(num <= 0) continue;
			int* list = new int[num+1];

			int list_num;
			// get local maximum in the 5k window
			get_local_maxima(mark_sig[j], num, maxima_half_size, list, list_num);

			for(int k = 0; k < list_num - 1; k++)
			{
				int index1 = list[k];
				int index2 = list[k+1];

				int loc1 = mark_sig[j][index1].pos;
				int loc2 = mark_sig[j][index2].pos;
				if(loc2 - loc1 <= overlap_half_size) // if distance of the two peaks is less than 5k
				{
					//remove left
					int temp_index = index1 - 1;
					while(temp_index > 0)
					{
						int temp_loc = mark_sig[j][temp_index].pos;
						if(temp_loc >= loc1 - maxima_half_size) mark_sig[j][temp_index].z = 0;
						else break;
						temp_index--;
					}
					//remove middle
					temp_index = index1;
					while(1)
					{
						int temp_loc = mark_sig[j][temp_index].pos;
						if(temp_loc >= loc1 && temp_loc <= loc2) mark_sig[j][temp_index].z = 0;
						else break;
						temp_index++;
					}
					//remove right
					temp_index = index2 + 1;
					while(temp_index < list_num)
					{
						int temp_loc = mark_sig[j][temp_index].pos;
						if(temp_loc <= loc2 + maxima_half_size) mark_sig[j][temp_index].z = 0;
						else break;
						temp_index++;
					}
				}
			}
			delete[] list;
		}
		//third, combined everything remain
		_sig* sig = new  _sig[max + 1];
		memset(sig,0,(max + 1)*sizeof(_sig));
		for(int j = 0; j < mark_num; j++)
		{
			int num = mark_sig_num[j];
			if(num <= 0) continue;

			for(int k = 0; k < num; k++)
			{
				int pos = mark_sig[j][k].pos;
				double z_val = mark_sig[j][k].z;
				if(z_val > 0)
				{
					int idx = pos / 100;
					// statistics for all marks on chr i
					sig[idx].sum += z_val;
					sig[idx].marks++;
					sig[idx].pos = pos;
				}
			}
		}

		// check overlap on all mark level, write to filter.idx
		print_non_overlapping(sig,max,i+1,overlap_half_size,gs);
		for(int j = 0; j < mark_num; j++) free(mark_sig[j]);
		delete[] sig;
	}
	for(int i = 0; i < mark_num; i++) free(mark_buf[i]);
	delete[] mark;
}

// type 1 core function, mean and std of all marks on encode region
//void get_normal_params(int stat_half_window_size, string mark_file, _mark **mark, _encode **encode, GS* gs)
void get_normal_params(int stat_half_window_size, string mark_file, string encode_file, int flag_global_background, GS* gs)
	// default stat_half_window_size=1000
{
	_mark* mark;
	_encode* encode;
	float* mark_buf[MAX_MARK]; // guess, marker infor for each marker
	int mark_num;
	int encode_num;

	mark_num = get_marks(mark_file,&mark); // initial marker information, get_marks(), DataUtils.cpp
	if( flag_global_background )
	{
		encode_num = gs->g_cnt;
		encode = new _encode[encode_num];
		for( int i=0; i<encode_num; i++ )
		{
			encode[i].chrom = i+1;
			encode[i].min = 1;
			encode[i].max = gs->g_size[i];
			encode[i].name = gs->id_chr_map[i];
		}
	}
	else
	{
		encode_num = get_encode("encode_hg18.bed",&encode,gs); // regions in encode pilot project, get encode information, DataUtil.cpp
	}

	if(mark_num < 0 || encode_num < 0) return;

	int old_chr = -1;
	double count = 0;
	// calculate mean of mark value for all encode regions
	for(int i = 0; i < encode_num; i++) // loop over all encode regions
	{
		int chr_idx = encode[i].chrom; // chr id, 1-24 for 1-22+XY
		int min = encode[i].min / 100; // bin number of region start, bin size 100bp
		int max = encode[i].max / 100; // bin number of region end
		//int buf_size = hg18[chr_idx-1]/100; // hg18[], chromosome size, DataUtil.cpp
		int buf_size = gs->g_size[chr_idx-1]/100; // hg18[], chromosome size, DataUtil.cpp

		// allocate and load data to mark_buf
		// if working on the chromosome, do nothing and continue to next step
		alloc_mem(mark_buf, mark_num,chr_idx,old_chr,mark,gs); // DataUtil.cpp
		old_chr = chr_idx;

		for(int j = min; j < max; j++)
		{
			// assign extended region: 1k--encode--1k
			int min1 = j - stat_half_window_size / 100;
			int max1 = j + stat_half_window_size / 100;
			if(min1 < 0) min1 = 0;
			if(max1 >= buf_size) max1 = buf_size;

			for(int k = 0; k < mark_num; k++)
			{
				double val = 0;
				for(int l = min1; l <= max1; l++) val += (double)mark_buf[k][l];
				mark[k].mean += val; // mean value of all encode regions
			}
			count += 1.0; // count for bin number
		}
	}
	for(int i = 0; i < mark_num; i++)
	{
		mark[i].mean /= count; // mean value of all encode regions
		cout<<mark[i].name.c_str()<<": "<<"mean="<<mark[i].mean<<endl;
	}

	// calculate std of mark value for all encode regions
	for(int i = 0; i < encode_num; i++) 
	{
		int chr_idx = encode[i].chrom;
		int min = encode[i].min / 100;
		int max = encode[i].max / 100;
		//int buf_size = hg18[chr_idx-1]/100;
		int buf_size = gs->g_size[chr_idx-1]/100;

		alloc_mem(mark_buf,mark_num,chr_idx,old_chr,mark,gs);
		old_chr = chr_idx;

		for(int j = min; j < max; j++)
		{
			int min1 = j - stat_half_window_size / 100;
			int max1 = j + stat_half_window_size / 100;
			if(min1 < 0) min1 = 0;
			if(max1 >= buf_size) max1 = buf_size;
			for(int k = 0; k < mark_num; k++)
			{
				double val = 0;
				for(int l = min1; l <= max1; l++) val += (double)mark_buf[k][l];
				mark[k].std += (val - mark[k].mean)*(val - mark[k].mean);
			}
		}
	}

	FileUtils file1;
	char rt_buf[1024];
	for(int i = 0; i < mark_num; i++)
	{
		mark[i].std /= (count - 1.0);
		mark[i].std = sqrt(mark[i].std); 
		sprintf(rt_buf,"%d\t%s\t%.15f\t%.15f",i+1,mark[i].name.c_str(),mark[i].mean,mark[i].std);
		file1.WriteLine("mark_.para",rt_buf);
		cout<<mark[i].name.c_str()<<": "<<"std="<<mark[i].std<<endl;
	}
	file1.LastWrite();
	file1.Release();
	for(int i = 0; i < mark_num; i++) free(mark_buf[i]);
	delete[] mark;
	delete[] encode;
}

// type 4 core function, +/- 20k around enriched region as background region, write to background_.para
//void get_background_dist_global(int bg_buf_size, string mark_file, string filter_file,_mark **mark, _filter **filter, GS* gs)
void get_background_dist_global(int bg_buf_size, string mark_file, string filter_file, GS* gs)
	// bg_buf_size: default 20000
	// mark_file: mark_info, "mark.idx", filename of marks
	// filter_file: filter, "filter.idx", enriched region information
{
	_mark* mark;
	_filter* filter;
	float* mark_buf[MAX_MARK];

	int mark_num = get_marks(mark_file,&mark);
	int filter_num = load_filter(filter_file,&filter,gs); // load_filter, DataUtil.cpp
	if(mark_num < 0 || filter_num < 0) return;

	qsort(filter,filter_num,sizeof(_filter),cmp_filter);
	double count = filter_num*(bg_buf_size/50 + 1) * 1.0;

	int old_chr = -1;
	printf("\ncomputing background mean......\n");
	for(int i = 0; i < filter_num; i++)
	{
		int chr_idx = filter[i].en;
		int buf_size = gs->g_size[chr_idx-1]/100;
		//int buf_size = hg18[chr_idx-1]/100;

		alloc_mem(mark_buf, mark_num,chr_idx,old_chr,mark,gs);
		old_chr = chr_idx;

		int center = filter[i].pos;
		int min = (center - bg_buf_size) / 100;
		int max = (center + bg_buf_size) / 100;
		if(min < 0) min = 0;
		if(max >= buf_size) max = buf_size;
		for(int k = 0; k < mark_num; k++)
		{
			double val = 0;
			for(int j = min; j <= max; j++) val += (double)mark_buf[k][j];
			mark[k].mean += val;
		}

	}
	for(int i = 0; i < mark_num; i++)
	{
		mark[i].mean /= count; // mean of signal on all filterd marks +/- 20k region
		cout<<mark[i].name.c_str()<<": "<<"mean="<<mark[i].mean<<endl;
	}

	printf("computing background std......\n");
	for(int i = 0; i < filter_num; i++)
	{
		int chr_idx = filter[i].en;
		int buf_size = gs->g_size[chr_idx-1]/100;
		//int buf_size = hg18[chr_idx-1]/100;

		alloc_mem(mark_buf, mark_num,chr_idx,old_chr,mark,gs);
		old_chr = chr_idx;

		int center = filter[i].pos;
		int min = (center - bg_buf_size) / 100;
		int max = (center + bg_buf_size) / 100;
		if(min < 0) min = 0;
		if(max >= buf_size) max = buf_size;
		for(int k = 0; k < mark_num; k++)
		{
			double val = 0;
			for(int j = min; j <= max; j++) 
			{
				val = (double)mark_buf[k][j];
				mark[k].std += (val - mark[k].mean)*(val - mark[k].mean);
			}
		}
	}

	FileUtils file1;
	char rt_buf[1024];
	for(int i = 0; i < mark_num; i++)
	{
		mark[i].std /= (count - 1.0); // standard deviation of signal on all filterd marks +/- 20k region
		mark[i].std = sqrt(mark[i].std); 
		sprintf(rt_buf,"%d\t%s\t%.15f\t%.15f",i+1,mark[i].name.c_str(),mark[i].mean,mark[i].std);
		file1.WriteLine("background_.para",rt_buf);
		cout<<mark[i].name.c_str()<<": "<<"std="<<mark[i].std<<endl;
	}
	file1.LastWrite();
	file1.Release();
	for(int i = 0; i < mark_num; i++) free(mark_buf[i]);
	delete[] mark;
	delete[] filter;
}
