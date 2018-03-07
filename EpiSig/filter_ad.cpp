#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<string> 
#include<vector>
#include <sys/stat.h>

#include "fast_pmf_cycle.h"
#include "utils.h"
#include "FileUtils.h"
#include "DataUtils.h"
#include "EpiSig_v1.0.h"
#include "genome_size.h"

void load_bin_file( char* fn_bin, float* bin_data, int size )
{
	FILE* fh;
	const int cnt = size/100+1;
	int r_cnt;

	if( (fh=fopen(fn_bin,"rb")) == NULL )
	{
		printf( "Cannot open file: %s\n", fn_bin );
		exit(1);
	}

	r_cnt = fread( bin_data, sizeof(float), cnt, fh );
	if( r_cnt != cnt )
	{
		printf( "Read number error for file: %s\t expected %d\t actual %d\n", fn_bin, cnt, r_cnt );
		exit(1);
	}
	fclose(fh);

	printf( "Load file %s completed\n", fn_bin );
}

int main( int argc, char **argv )
{
        if( argc < 6 )
        {
                printf( "Usage: %s background_.para filter.idx width wandering_dist genome_name\n", argv[0] );
                exit(1);
        }
        string fn_background(argv[1]);
        string fn_filter(argv[2]);
        int width = atoi(argv[3]);
        int wandering_dist = atoi(argv[4]);

	GS gs;
	GS_init( gs, argv[5] );

	struct _mark *mark;
	struct _filter *filter;

	int mark_num = load_marks(fn_background,&mark);
	int filter_num = load_filter(fn_filter,&filter,&gs);

	if( mark_num<=0 || filter_num<=0 )
	{
		printf( "Fatal errors! Failed to load background_.para or filter.idx\n" );
		exit(1);
	}
	
	int half_wsize = (width+wandering_dist)/2;
	int window_size = width+wandering_dist+1;

	float *loci_data = new float[filter_num*window_size];
	float *bin_data = new float[gs.g_size[0]/100+1]; // temp solution, introduce max_genome_size in future.
	//float *bin_data = new float[hg18[0]/100+1]; // temp solution, introduce max_genome_size in future.
	string dir_name("out_data_ad/");
	string fn;
	
	char fn_bin[1024];
	FILE *fh_dat;
	float *idx_loci;
	int idx_start;
	int idx_end;
	int curr_chr;
	int max_bin = gs.g_size[0]/100+1;
	//int max_bin = hg18[0]/100+1;

	qsort(filter,filter_num,sizeof(_filter),cmp_filter);
	mkdir(dir_name.c_str(),0755);
	for( int i=0; i<mark_num; i++ )
	{
		curr_chr = -1;
		idx_loci = loci_data;
		for( int fi=0; fi<filter_num; fi++, idx_loci+=window_size )
		{
			//printf( "mark %d\t filter %d\tchr %d\tpos %d\n", i, fi, filter[fi].en, filter[fi].pos );
			if( curr_chr != filter[fi].en )	
			{
				sprintf( fn_bin, "out_data/%s.%d.bin", mark[i].name.c_str(), filter[fi].en );
				load_bin_file( fn_bin, bin_data, gs.g_size[filter[fi].en-1] );
				//load_bin_file( fn_bin, bin_data, hg18[filter[fi].en-1] );
				curr_chr = filter[fi].en;
				max_bin = gs.g_size[filter[fi].en-1]/100+1;
				//max_bin = hg18[filter[fi].en-1]/100+1;
			}
			
			idx_start = filter[fi].pos/100 - half_wsize;
			idx_end = idx_start + window_size - 1;
			//printf( "idx_start %d\tidx_end %d\n", idx_start, idx_end );

			if( idx_start>=0 && idx_end<=max_bin )
			{
				memcpy( idx_loci, bin_data+idx_start, window_size*sizeof(float) );
			}
			else if( idx_end>max_bin )
			{
//				printf( "tail bin: idx_start %d\tidx_end%d\n", idx_start, idx_end );
				memcpy( idx_loci, bin_data+idx_start, (max_bin-idx_start+1)*sizeof(float) );
				memset( idx_loci+(max_bin-idx_start+1), 0, (idx_end-max_bin)*sizeof(float) );
			}
			else if( idx_start < 0 )
			{
//				printf( "head bin: idx_start %d\tidx_end%d\n", idx_start, idx_end );
				memset( idx_loci, 0, -idx_start*sizeof(float) );
				memcpy( idx_loci-idx_start, bin_data, (idx_end+1)*sizeof(float) ); 
			}
		}
		
		fn = dir_name + mark[i].name + ".bin";
		if( ( fh_dat=fopen(fn.c_str(),"wb") ) == NULL )
		{
			printf( "Unable to write to: %s\n", fn.c_str() );
			exit(1);
		}
		fwrite( loci_data, sizeof(float), filter_num*window_size, fh_dat );
		fclose(fh_dat);
	}

	if( (fh_dat=fopen("filter_ad.idx","w")) == NULL )
	{
		printf( "Unable to write to: filter_ad.idx\n" );
	}

	for( int i=0; i<filter_num-1; i++ )
	{
		//fprintf( fh_dat, "%s\t%d\n", id2chr(filter[i].en).c_str(), filter[i].pos );
		fprintf( fh_dat, "%s\t%d\n", gs.id_chr_map[filter[i].en-1].c_str(), filter[i].pos );
	}
	fprintf( fh_dat, "%s\t%d", gs.id_chr_map[filter[filter_num-1].en-1].c_str(), filter[filter_num-1].pos );
	//fprintf( fh_dat, "%s\t%d", id2chr(filter[filter_num-1].en).c_str(), filter[filter_num-1].pos );
	fclose(fh_dat);
	
	delete[] bin_data;
	delete[] loci_data;

	return 0;
}
