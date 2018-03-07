#include <iostream>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <math.h>
#include <vector>
#include <string>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "utils.h"
#include "FileUtils.h"
#include "ConfUtils.h"
#include "DataUtils.h"
#include "fast_pmf_cycle.h"
#include "EpiSig_v1.0.h"
#include "filter.h"

using namespace std;

struct timeval tpstart,tpend;
float timeuse;

void print_usage()
{
  string buffer;
  buffer = "Usage: EpiSig_v1.0 [--enriched_region | --global] [--mark_info mark.idx] ...\n";
  buffer += "Options:\n";
  buffer += "-h, --help      display this help and exit\n";
  buffer += "-g, --global     EpiSig_v1.0 will scan the whole genome from input marks\n";
  buffer += "-e, --enriched_region     provide an enriched region file, for example, filter.idx\n";
  buffer += "-m, --mark_info     marks information file. default: mark.idx\n";
  buffer += "-o, --output_name     output file name. default: ChromaSig\n";
  buffer += "-p, --pval_cutoff     p-value cutoff to filter enriched regions. default: 1e-5\n";
  buffer += "-s, --stat_half_window_size   default: 1000\n";
  buffer += "-x, --maxima_half_window_size   default: 2500\n";
buffer += "-r, --overlap_half_window_size   default: 5000\n";
  buffer += "-w, --width       width must between 20 and 200, default: 40\n";
  buffer += "-d, --max_wandering_dist   max_wandering_dist must between 10 and 200, default: 20\n";
  buffer += "-f, --std_factor     default: 1.75\n";
  buffer += "-b, --prior_bg       default: 0.01\n";
  buffer += "-n, --threads      mutiple threads number, default: 4\n";
  buffer += "-i, --genome      genome used in EpiSig_v1.0. Two options: pre-defined genome names: hg18, hg19, hg38, mm9, or mm10, or user-defined tab-delimit file\n";
  printf("%s\n",buffer.c_str());
  fflush(stdout);
}
void print_para(_para para)
{
  printf("EpiSig_v1.0 use following parameters.\n");
  printf("enriched_region: %s\n", para.filter.c_str());
  printf("mark_inf: %s\n", para.mark_info.c_str());
  printf("output_name: %s\n", para.output_name.c_str());
  printf("pval_cutoff: %f\n", para.pval_cutoff);
  printf("stat_half_window_size: %d\n",para.stat_half_window_size);
  printf("maxima_half_window_size: %d\n",para.maxima_half_window_size);
  printf("overlap_half_window_size: %d\n",para.overlap_half_window_size);
  printf("width: %d\n", para.width);
  printf("max_wandering_dist: %d\n", para.wandering_dist);
  printf("std_factor: %f\n", para.std_factor);
  printf("prior_bg: %f\n", para.prior_bg);
  //printf("threads: %d\n", para.max_thread);
  fflush(stdout);
}
void showtime(const char* mesg)
{
  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  timeuse/=1000000;
  int hour = ((int)timeuse) / 3600;
  int minute = ((int)timeuse%3600) / 60;
  float second = timeuse - hour*3600-minute*60;
  printf("%s, elapsed:%dh %dm %.2fs\n",mesg, hour, minute, second);
  fflush(stdout);
}

int compute_motif( _data_size* ds, _gibbs_state* gibbs, int gibbs_size, float* mark_sig, _logr_sums* logr, double* dd, _filter* filter )
{
	int num_elements;
	long d_idx;
	long std_idx;
	long lr_idx;
	double mean;

	for( d_idx=0,std_idx=ds->m_wp; d_idx<ds->m_wp; d_idx++,std_idx++ )
	{
		dd[d_idx] = logr[d_idx].sum;
		dd[std_idx] = logr[d_idx].sum2;
	}

	for( lr_idx=0,std_idx=ds->m_wp; lr_idx<ds->m_wp; lr_idx++,std_idx++ )
	{
		num_elements = logr[lr_idx].count;
		if(num_elements <= 1)
		{
			printf("Zero or negative element motif!\n");
			exit(1);
		}
		
		dd[lr_idx] /= num_elements;

		mean = dd[lr_idx];
		dd[std_idx] -= num_elements*mean*mean;
		dd[std_idx] /= (num_elements - 1);
	}

	return 1;
}

int get_running_stats(_data_size* ds, _gibbs_state* gibbs, int gibbs_size, float* mark_sig, _logr_sums* logr,_filter *filter)
{
	int mi, wi, gi;
	int pol, idx, offset;
	long lr_idx;
	long d_idx;
	double val;

	for( lr_idx=0; lr_idx<ds->m_wp; lr_idx++ )
	{
		logr[lr_idx].sum = 0.5;
		logr[lr_idx].sum2 = 0.25;
		logr[lr_idx].count = 1;
	}

	for( gi=0; gi<gibbs_size; gi++ )
	{
		pol = gibbs[gi].pol;
		idx = gibbs[gi].id;
		offset = (gibbs[gi].loc - filter[idx].pos) / 100; // How offset and pol occurs in gibbs

		d_idx = ds->mark_sig_idx_total[idx]+ds->width_half+ds->wandering_half+offset;
		for( mi=0; mi<ds->mark_num; mi++)
		{
			for( wi=-ds->width_half; wi<=ds->width_half; wi++ )
			{
				lr_idx = ds->logr_idx[mi]+wi;

				val = mark_sig[d_idx+wi*pol];
				logr[lr_idx].sum += val;
				logr[lr_idx].sum2 += val*val;
				logr[lr_idx].count++;
			}
			d_idx += ds->w_size;
		}
	}

	return 1;
}

int compute_motif_full( _data_size* ds, _gibbs_state *gibbs, float* mark_sig, _logr_sums* logr, double* dd, Scratch *ss, _mark *mark,_filter *filter )
{
	get_running_stats( ds, gibbs, ds->gibbs_state_size, mark_sig, logr, filter );
	compute_motif( ds, gibbs, ds->gibbs_state_size, mark_sig, logr, dd, filter );

	int mi, wi;
	int lr_idx;
	int std_idx;

	lr_idx = 0;
	for( mi=0,lr_idx=2*ds->m_wp; mi<ds->mark_num; mi++,lr_idx++ )
	{
		dd[lr_idx] = mark[mi].mean;
	}
	memset( dd+2*ds->m_wp+ds->mark_num, 0, ds->mark_num*sizeof(double) );
	for( mi=0,lr_idx=2*ds->m_wp+ds->mark_num,std_idx=ds->m_wp; mi<ds->mark_num; mi++,lr_idx++ )
	{
		for( wi=0; wi<ds->width_plus; wi++,std_idx++ )
		{
			dd[lr_idx] += sqrt(dd[std_idx]);
		}
		dd[lr_idx] /= ds->width_plus;
	}

	return 1;
}

void fisher_yates_shuffle(_distance* d,int size,THREAD_RAND* thread_rand)
{
	int j;
	int idx=0;
	_distance a;

	for(int i = size-1; i > 0; i--)
	{
		j = thread_rand->rand_num(idx) % (i + 1);
		if(i == j) continue;

		a.id = d[i].id;
		a.distance = d[i].distance;
		d[i].id = d[j].id;
		d[i].distance = d[j].distance;
		d[j].id = a.id;
		d[j].distance = a.distance;
	}
}

void get_euclidean_norm(_data_size* ds, float* mark_sig, double* norm)
{
	long idx;
	double dist_val;
	double val;

	idx = ds->wandering_half;
	for( int fi=0; fi<ds->filter_num; fi++ )
	{
		dist_val = 0;
		for( int mi=0; mi<ds->mark_num; mi++,idx+=para.wandering_dist )
		{
			for( int wi=0; wi<=para.width; wi++,idx++ )
			{
				val = mark_sig[idx];
				dist_val += val*val;
			}
		}
		norm[fi] = dist_val;
	}
}

void tournament_sort(_data_size* ds, int id, int* clusters, _distance* distance, _gibbs_state* gibbs, int& gibbs_size, _filter* filter, THREAD_RAND* thread_rand )
{
	int m_size = ds->norm_size;
	int count;
	int offset;
	int idx1, idx2;
	int *tag = new int[ds->norm_size];
	int *tmp_tag = new int[ds->norm_size];

	fisher_yates_shuffle(distance,ds->norm_size,thread_rand); 
	for(int i = 0; i <ds->norm_size; i++) 
	{
		tag[i] = i;
	}

	while(m_size >= MAX_SEED_SIZE)
	{
		count = 0;
		offset = m_size % 2;
		for(int ii = 0; ii<m_size-offset; ii+=2)
		{
			idx1 = tag[ii];
			idx2 = tag[ii+1];
			if(distance[idx1].distance < distance[idx2].distance) tmp_tag[count]=idx1;
			else tmp_tag[count]=idx2;
			count++;
		}
		if(offset)
		{
			tmp_tag[count] = tag[m_size-1];
			count++;
		}
		m_size = count;
		memcpy(tag,tmp_tag,m_size*sizeof(int));
	}

	for(int i = 0; i < m_size; i++)
	{
		idx1 = tag[i];
		idx2 = distance[idx1].id;  
		gibbs[i].id = idx2;
		gibbs[i].en = filter[idx2].en;
		gibbs[i].loc = filter[idx2].pos;
		gibbs[i].pol = 1;  
	}
	gibbs[m_size].id = id;
	gibbs[m_size].en = filter[id].en;
	gibbs[m_size].loc = filter[id].pos;
	gibbs[m_size].pol = 1;
	gibbs_size = m_size + 1;

	delete[] tag;
	delete[] tmp_tag;
}

double get_seed_score( _data_size* ds, _gibbs_state* gibbs, int gibbs_size, float *mark_sig, _filter* filter )
{
	int mi, wi;
	int d_idx;
	double score;
	double avg_std;
	double differential_factor;
	double* p_mean = new double[ds->width_plus];
	_logr_sums* logr = new _logr_sums[ds->m_wp];
	_motif_dist* mdist = new _motif_dist[ds->m_wp];
	double* dd = new double[2*ds->m_wp];

	get_running_stats( ds, gibbs, gibbs_size, mark_sig, logr, filter );
	compute_motif( ds, gibbs, gibbs_size, mark_sig, logr, dd, filter );

	d_idx = 0;
	score = 0;
	double* idx_std = dd+ds->m_wp;
	for( mi=0; mi<ds->mark_num; mi++ )
	{
		avg_std = 0;
		for( wi=0; wi<=para.width; wi++,d_idx++ )
		{
			avg_std += sqrt(idx_std[d_idx]);
			p_mean[wi] = fabs(dd[d_idx]);
		}
		avg_std /= ds->width_plus;
		qsort(p_mean,ds->width_plus,sizeof(double),cmp_double);

		differential_factor = 0;
		for( int i=0; i<(ds->width_plus)/2; i++)
		{
			differential_factor += p_mean[i];
			differential_factor -= p_mean[para.width-i];
		}
		score += differential_factor/avg_std;
	}

	delete[] p_mean;
	delete[] mdist;
	delete[] logr;

	return score;
}

int get_not_done_ids( _data_size* ds, int* filter_status, int* not_done_list )
{
	int cnt = 0;
	for( int i=0; i<ds->filter_num; i++ )
	{
		if( filter_status[i]==0 )
		{
			not_done_list[cnt] = i;
			cnt++;
		}
	}

	return cnt;
}

void get_node_counts( _data_size* ds, int* filter_status, int* node_name_list, int mpi_size, int* mpi_counts )
{
	ds->not_done_ids_size_total = 0;
	memset( mpi_counts, 0, mpi_size*sizeof(int) );
	
	for( int i=0; i<ds->filter_num_total; i++ )
	{
		if( filter_status[i] == 0 )
		{
			mpi_counts[ node_name_list[i] ]++;
		}
	}

	for( int i=0; i<mpi_size; i++ )
	{
		ds->not_done_ids_size_total += mpi_counts[i];
	}
}

void init_gibbs_seed( _data_size* ds, _gibbs_state* seed, double *all_norm, _filter *filter, float *mark_sig, int* status_list, int* not_done_list, THREAD_RAND* thread_rand, int mpi_rank, int mpi_size, int* mpi_counts, GS* gs )
{
	int num_sites = 100;
	float* mark_sig_seed = new float[num_sites*ds->m_wp+num_sites+1];

	float* idx_seed = NULL;
	_filter* first_sets = NULL;
	_norm* norm = NULL;
 	int* undo_list = NULL;
	int* mpi_recv_cnt = NULL;
	int* mpi_recv_displ = NULL;
	if( mpi_rank == 0 )
	{
		int id;
		int interval;

                int idx = 0;
		float* idx_mark = NULL;
		int size_norm = 0;
		const int percentile = 0;
		const double ratio = 0.25;

		norm = new _norm[ds->filter_num_total];
		undo_list = new int[ds->filter_num_total];
		for( int i=0; i<ds->filter_num_total; i++ )
		{
			if( status_list[i] == 0 )
			{
				norm[size_norm].id = i;
				norm[size_norm].norm = all_norm[i];
				undo_list[size_norm] = i;
				size_norm++;
			}
		}

		ds->sites_size = num_sites;
		ds->norm_size = size_norm;
		qsort(norm,ds->norm_size,sizeof(_norm),cmp_norm);

		interval = int(ratio * ds->norm_size / ds->sites_size);
		if( interval==0 )
		{
			interval = 1; 
		}

		if( ds->sites_size>ds->norm_size )
		{
			ds->sites_size = ds->norm_size;
		}
		mark_sig_seed[0] = ds->sites_size;
		first_sets = new _filter[ds->sites_size]; // DataUtil.h
		for( int i=0; i<ds->sites_size; i++ )
		{
			id = norm[percentile+i*interval].id;
			first_sets[i].en = filter[id].en;
			first_sets[i].pos = id;
		}
		qsort(first_sets,ds->sites_size,sizeof(_filter),cmp_filter);

		idx_seed = mark_sig_seed+1;
		for( int i=0; i<ds->sites_size; i++ )
		{
			*idx_seed = (float)first_sets[i].pos;
			idx_seed++;
			idx_mark = mark_sig + ds->mark_sig_idx_total[first_sets[i].pos] + ds->wandering_half;
			for( int mi=0; mi<ds->mark_num; mi++,idx_seed+=ds->width_plus,idx_mark+=ds->w_size )
			{
				memcpy( idx_seed, idx_mark, ds->width_plus*sizeof(float) );
			}
		}
		idx_seed = NULL;
		idx_mark = NULL;

		mpi_recv_cnt = new int[mpi_size];
                for( int i=0; i<mpi_size; i++ )
                {
                        mpi_recv_cnt[i] = ds->sites_size*mpi_counts[i];
                }

                mpi_recv_displ = new int[mpi_size];
                for( int i=0; i<mpi_size; i++ )
                {
                        mpi_recv_displ[i] = idx;
                        idx += mpi_recv_cnt[i];
                }

		delete[] undo_list;
		delete[] norm;
	}

	MPI_Bcast( mark_sig_seed, num_sites*ds->m_wp+num_sites+1, MPI_FLOAT, 0, MPI_COMM_WORLD );

	ds->sites_size = (int)mark_sig_seed[0];

	_distance* undo_distance = NULL; 
	if( mpi_rank == 0 )
	{
		undo_distance = new _distance[ds->sites_size*ds->not_done_ids_size_total];
	}
	else
	{
		undo_distance = new _distance[ds->sites_size*ds->not_done_ids_size];
	}
	float* idx_mark;
	float* idx_seed_coor;
	double dist;
	double diff;
	int first_id;
	int fi;

	idx_seed = mark_sig_seed+1;
	_distance* idx_undo_d = undo_distance;
	for( int ni=0; ni<ds->sites_size; ni++,idx_seed+=1+ds->m_wp )
	{
		first_id = int(*idx_seed);
		for( int ui=0; ui<ds->not_done_ids_size; ui++,idx_undo_d++ )
		{
			fi = not_done_list[ui];
			if( ds->filter_start+fi == first_id )
			{
				dist = 99999999;
			}
			else
			{
				dist = 0.;
				idx_mark = mark_sig+ds->mark_sig_idx[fi]+ds->wandering_half;
				idx_seed_coor = idx_seed+1;
				for( int mi=0; mi<ds->mark_num; mi++ )
				{
					for( int wi=0; wi<=para.width; wi++,idx_mark++,idx_seed_coor++ )
					{
						diff = ((double)(*idx_mark)) - ((double)(*idx_seed_coor));
						dist += diff*diff;
					}
					idx_mark += para.wandering_dist;
				}
			}
			idx_undo_d->id = fi+ds->filter_start;
			idx_undo_d->distance = dist;
		}
	}
	idx_seed = NULL;
	idx_seed_coor = NULL;
	idx_undo_d = NULL;
	idx_mark = NULL;

	MPI_Datatype MPI_DISTANCE;
	MPI_Datatype oldtype[2] = {MPI_INT,MPI_DOUBLE};
	int blockcounts[2] = {1,1};
	MPI_Aint offsets[2] = {0,8};
	MPI_Type_struct( 2, blockcounts, offsets, oldtype, &MPI_DISTANCE );
	MPI_Type_commit(&MPI_DISTANCE);

	if( mpi_rank!=0 )
        {
                MPI_Send( undo_distance, ds->sites_size*ds->not_done_ids_size, MPI_DISTANCE, 0, mpi_rank, MPI_COMM_WORLD );
        }
        if( mpi_rank == 0 )
        {
		MPI_Status status;
                for( int i=1; i<mpi_size; i++ )
                {
                        MPI_Recv( undo_distance+mpi_recv_displ[i], mpi_recv_cnt[i], MPI_DISTANCE, i, i, MPI_COMM_WORLD, &status );
                }
        }
	MPI_Type_free(&MPI_DISTANCE);

	if( mpi_rank == 0 )
	{
		_distance* idx_data = NULL;
		_distance* data = new _distance[ds->not_done_ids_size_total];
		_gibbs_state gibbs[MAX_SEED_SIZE];
		int gibbs_size;
		double score;
		double best_score = -1;

		for( int i=0; i<ds->sites_size; i++ )
		{
			idx_data = data;
			for( int mpi_i=0; mpi_i<mpi_size; mpi_i++ )
			{
				memcpy( idx_data, undo_distance+mpi_recv_displ[mpi_i]+mpi_counts[mpi_i]*i, mpi_counts[mpi_i]*sizeof(_distance) );
				idx_data += mpi_counts[mpi_i];
			}
			tournament_sort( ds, first_sets[i].pos, not_done_list, data, gibbs, gibbs_size, filter, thread_rand );
			score = get_seed_score( ds, gibbs, gibbs_size, mark_sig, filter );

			if( score>best_score )
			{
				best_score = score;
				ds->seed_size = gibbs_size;
				memcpy(seed,gibbs,sizeof(_gibbs_state)*gibbs_size);
			}
		}
		idx_data = NULL;

		printf("\nBest seed score = %.4f, done:%d, not done:%d, seed set:\n", best_score, ds->done_list_size, ds->not_done_ids_size_total);
		for(int i = 0; i < ds->seed_size; i++)
		{
			printf("%d\t%s\t%d\t%d\n",seed[i].id,gs->id_chr_map[seed[i].en-1].c_str(),seed[i].loc,seed[i].pol);
		}

                delete[] data;
                delete[] mpi_recv_cnt;
                delete[] mpi_recv_displ;
                delete[] first_sets;
	}


	showtime("computing euclidean distance");
	delete[] mark_sig_seed;
        delete[] undo_distance;
}

void load_mark_data( float *mark_sig, string dir_data, _mark *mark, int mark_num, int w_size, long f_start, long filter_num )
{
	unsigned long ival = mark_num*w_size;
	unsigned long cnt = filter_num*w_size;
	unsigned long r_cnt;
	unsigned long sval = w_size*sizeof(float);
	long offset = f_start*w_size*sizeof(float);
	string fn;
	float* mark_idx;
	float* bin_data;
	float* bin_idx;
	FILE* fh;

	bin_data = new float[cnt];
	for( int i=0; i<mark_num; i++ )
	{
		fn = dir_data+"/"+mark[i].name+".bin";
		if( (fh=fopen(fn.c_str(),"rb")) == NULL )
		{
			printf( "Unable to open file: %s\n", fn.c_str() );
			exit(1);
		}
		
		fseek(fh,offset,SEEK_SET);
		r_cnt = fread(bin_data,sizeof(float),cnt,fh);
		if( r_cnt != cnt )
		{
			printf( "Data size %lu not match for %s\n", cnt, fn.c_str() );
			printf( "mark=%d\twindow_size=%d\tfilter_num=%ld\tcnt=%lu\tr_cnt=%lu\n", i, w_size, filter_num, cnt, r_cnt );
			fclose(fh);
			exit(1);
		}
		fclose(fh);

		bin_idx = bin_data;
		mark_idx = mark_sig+i*w_size;
		for( int fi=0; fi<filter_num; fi++,mark_idx+=ival,bin_idx+=w_size )
		{
			memcpy(mark_idx,bin_idx,sval);
		}
	}

	delete[] bin_data;
}

void main_gibbs_cluster()
{
	int mpi_size;
	int mpi_rank;
	MPI_Datatype gtype1;
	MPI_Datatype gtype3;

	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
        MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	GS gs;
	GS_init( gs, para.genome_name );

	THREAD_RAND thread_rand(mpi_size,para.rand_seed);
	_mark* mark;
	_filter* filter;
	Scratch ss;

	int cycle_num = NUM_CYCLES;

	struct _data_size ds;
	ds.mark_num = load_marks("background_.para",&mark);
	ds.filter_num_total = load_filter(para.filter,&filter,&gs);
	if( ds.mark_num<=0 || ds.filter_num_total<=0 )
	{
		printf("Fatal errors! Failed to load background_.para %d or %s %ld on rank %d\n", ds.mark_num, para.filter.c_str(), ds.filter_num_total, mpi_rank );
		return;
	}

	ds_init( &para, &ds, mpi_size, mpi_rank );

	int* filter_status = new int[ds.filter_num];
	int* not_done_ids = new int[ds.filter_num];
	memset( filter_status, 0, ds.filter_num*sizeof(int) );

	int* filter_status_all = NULL;
	int* node_name_all = NULL;
	int* mpi_counts = NULL;
	if( mpi_rank==0  )
	{
		mpi_counts = new int[mpi_size];

		filter_status_all = new int[ds.filter_num_total];
		memset( filter_status_all, 0, ds.filter_num_total*sizeof(int) );

		int idx=0;
		node_name_all = new int[ds.filter_num_total];
		for( int i=0; i<ds.mpi_num_extra; i++ )
		{
			for( int j=0; j<ds.mpi_num_size+1; j++,idx++ )
			{
				node_name_all[idx] = i;
			}
		}
		for( int i=ds.mpi_num_extra; i<mpi_size; i++ )
		{
			for( int j=0; j<ds.mpi_num_size; j++,idx++ )
			{
				node_name_all[idx] = i;
			}
		}
	}

	float* mark_sig;
	if( mpi_rank != 0 )
	{
		mark_sig = new float[ds.mark_num*ds.filter_num*ds.w_size];
		load_mark_data( mark_sig, "out_data_ad", mark, ds.mark_num, ds.w_size, ds.filter_start, ds.filter_num );
	}
	else
	{
		mark_sig = new float[ds.mark_num*ds.filter_num_total*ds.w_size];
		load_mark_data( mark_sig, "out_data_ad", mark, ds.mark_num, ds.w_size, ds.filter_start, ds.filter_num_total );
	}

	init( &ss, &ds );
	double* dd = new double[(2*ds.width_plus+2)*ds.mark_num];
	int* g_size = NULL;
	int* g_displ = NULL;
	double* all_norm = NULL;
	double* norm_node = new double[ds.filter_num];
	
	if( mpi_rank == 0 )
	{
		all_norm = new double[ds.filter_num_total];

		int g_idx = 0;
		g_size = new int[mpi_size];
		g_displ = new int[mpi_size];

		for( int i=0; i<ds.mpi_num_extra; i++,g_idx+=ds.mpi_num_size+1 )
		{
			g_size[i] = ds.mpi_num_size+1;
		}
		for( int i=ds.mpi_num_extra; i<mpi_size; i++,g_idx+=ds.mpi_num_size )
		{
			g_size[i] = ds.mpi_num_size;
		}
		g_displ[0] = 0;
		for( int i=1; i<mpi_size; i++ )
		{
			g_displ[i] = g_displ[i-1]+g_size[i-1];
		}
	}

	get_euclidean_norm( &ds, mark_sig, norm_node );

	MPI_Type_contiguous( ds.filter_num, MPI_DOUBLE, &gtype1 );
	MPI_Type_commit(&gtype1);
	MPI_Gatherv( norm_node, 1, gtype1, all_norm, g_size, g_displ, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Type_free(&gtype1);

	_gibbs_state* gibbs_state = NULL;
	_gibbs_state* new_gibbs_state = NULL;
	_logr_sums *logr_sums = NULL; 
	int* size_gibbs_mpi = NULL;
	int* displ_gibbs_mpi = NULL;
	int* gibbs_data_all = NULL;
	if( mpi_rank == 0 )
	{
		showtime("Computing Euclidean Norm");

		gibbs_state = new _gibbs_state[ds.filter_num_total];
		new_gibbs_state = new _gibbs_state[ds.filter_num_total]; 
		logr_sums = new _logr_sums[ds.mark_num*ds.width_plus];

		size_gibbs_mpi = new int[mpi_size];
		gibbs_data_all = new int[4*ds.filter_num_total];
		displ_gibbs_mpi = new int[mpi_size];
	}

	_gibbs_state copy_gibbs_state[MAX_SEED_SIZE];
	_gibbs_state seed[MAX_SEED_SIZE];
	int* gibbs_data = new int[ds.filter_num*4];
	int round_n = 0;
	int flag_cycle = 1;
	int exclude_id;
	int filter_id;
	int size_gibbs_data;
	int gibbs_size;
	int lr_idx;

	while( flag_cycle )
	{
		round_n++;

		ds.not_done_ids_size = get_not_done_ids( &ds, filter_status, not_done_ids );
		if( mpi_rank==0 )
		{
			get_node_counts( &ds, filter_status_all, node_name_all, mpi_size, mpi_counts );
			printf( "\ncluster = %d\t %d loci clustered, %d loci unfinished\n", ds.cluster_num, ds.done_list_size, ds.not_done_ids_size_total );
		}

		init_gibbs_seed( &ds, seed, all_norm, filter, mark_sig, filter_status_all, not_done_ids, &thread_rand, mpi_rank, mpi_size, mpi_counts, &gs );
		if( mpi_rank==0 )
		{
			printf( "#### END: initial seed round %d size %d ####\n", round_n, ds.seed_size );
			fflush(stdout);
			memcpy(copy_gibbs_state,seed,ds.seed_size*sizeof(_gibbs_state)); 
			ds.copy_gibbs_state_size = ds.seed_size;
			memcpy(gibbs_state,seed,ds.seed_size*sizeof(_gibbs_state));
			ds.gibbs_state_size = ds.seed_size;
		}

		showtime("init_gibbs_seed");
		for(int rd = 1; rd <= cycle_num;rd++)
		{
			if( mpi_rank == 0 )
			{
				compute_motif_full(&ds,gibbs_state, mark_sig, logr_sums, dd, &ss, mark, filter );
			}
			MPI_Bcast( dd, (2*ds.width_plus+2)*ds.mark_num, MPI_DOUBLE, 0, MPI_COMM_WORLD );

			size_gibbs_data = 0;
			for( int ni=0; ni<ds.not_done_ids_size; ni++ ) 
			{	
				exclude_id = not_done_ids[ni];
				lr_idx = ds.mark_sig_idx[exclude_id];
				compute_pmf(dd,mark_sig+lr_idx,&ss,&ds);
				if( !get_votes(&ss,&ds) )
				{
					tally_votes(&ss,&ds);
					filter_id = exclude_id+ds.filter_start;
					gibbs_data[size_gibbs_data] = filter_id;
					gibbs_data[size_gibbs_data+1] = filter[filter_id].en;
					gibbs_data[size_gibbs_data+2] = filter[filter_id].pos+100*get_return_loc(&ss);
					gibbs_data[size_gibbs_data+3] = get_return_pol(&ss);

					size_gibbs_data+=4;
				}
			}

			MPI_Gather( &size_gibbs_data, 1, MPI_INT, size_gibbs_mpi, 1, MPI_INT, 0, MPI_COMM_WORLD );
			if( mpi_rank == 0 )
			{
				displ_gibbs_mpi[0] = 0;
				for( int i=1; i<mpi_size; i++ )
				{
					displ_gibbs_mpi[i] = displ_gibbs_mpi[i-1]+size_gibbs_mpi[i-1];
				}
			}
			if( mpi_rank != 0 )
                        {
                                MPI_Send( gibbs_data, size_gibbs_data, MPI_INT, 0, mpi_rank, MPI_COMM_WORLD );
                        }
                        else
                        {
                                MPI_Status status;
                                for( int i=1; i<mpi_size; i++ )
                                {
                                        MPI_Recv( gibbs_data_all+displ_gibbs_mpi[i], size_gibbs_mpi[i], MPI_INT, i, i, MPI_COMM_WORLD, &status );
                                }
                                memcpy( gibbs_data_all, gibbs_data, size_gibbs_mpi[0]*sizeof(int) );

				ds.gibbs_state_size = 0;
				gibbs_size = displ_gibbs_mpi[mpi_size-1]+size_gibbs_mpi[mpi_size-1];
				for( int j=0; j<gibbs_size; j+=4 )
				{
					gibbs_state[ds.gibbs_state_size].id = gibbs_data_all[j];
					gibbs_state[ds.gibbs_state_size].en = gibbs_data_all[j+1];
					gibbs_state[ds.gibbs_state_size].loc = gibbs_data_all[j+2];
					gibbs_state[ds.gibbs_state_size].pol = gibbs_data_all[j+3];

#if PP_DEBUG
					printf( "Debug gibbs_state: cluster=%d\tcycle=%d\tidx=%d\tid=%d\ten=%d\tloc=%d\tpol=%d\n", ds.cluster_num, rd, j, gibbs_data_all[j], gibbs_data_all[j+1], gibbs_data_all[j+2], gibbs_data_all[j+3] );
#endif
					ds.gibbs_state_size++;
				}
				showtime("This cycle finished");
				printf("Update best alignment set, number:%d\n", ds.gibbs_state_size);

				if(ds.gibbs_state_size <= ds.seed_size)
				{
					printf("Less than seed size, use seed set:%d\n",ds.copy_gibbs_state_size);
					for(int i = 0; i < ds.copy_gibbs_state_size;i++)
					{
						filter_status_all[copy_gibbs_state[i].id] = 1;
						ds.done_list_size++;
					}    
					rd = cycle_num+1;
				}
				printf( "Cycle %d finished, alignment set size %d ...\n", rd, ds.gibbs_state_size );
				fflush(stdout);
			}

			MPI_Bcast( &rd, 1, MPI_INT, 0, MPI_COMM_WORLD );
		}

		if( mpi_rank == 0 )
		{
			if(ds.gibbs_state_size > ds.seed_size)
			{
				printf( "Write data to file\n" );
				fflush(stdout);

				char fname1[128];
				char fname2[128];
				char buf[256];

				sprintf(fname1,"%s",para.output_name.c_str());
				sprintf(fname2,"%s.clusters",para.output_name.c_str());
				if(ds.cluster_num == 0)
				{
					remove(fname1);
					remove(fname2);
				}
				FileUtils file1(20480);
				FileUtils file2(20480);

				strcpy(buf,"");
				if(FileLength(fname1) > 0) file1.WriteLine(fname1,buf,"a+");
				if(FileLength(fname2) > 0) file2.WriteLine(fname2,buf,"a+");

				sprintf(buf,"\ncluster = %d\nmotif at end of %d cycles!",ds.cluster_num, cycle_num);
				file1.WriteLine(fname1, buf, "a+");

				get_running_stats( &ds, gibbs_state, ds.gibbs_state_size, mark_sig, logr_sums, filter );
				compute_motif( &ds, gibbs_state, ds.gibbs_state_size, mark_sig, logr_sums, dd, filter );

				int mi, wi;
				int m_idx;
				long lr_idx;
				for( mi=0; mi<ds.mark_num; mi++)
				{
					sprintf(buf, "\nmod = %s", mark[mi].name.c_str());
					file1.WriteLine(fname1, buf, "a+");
					sprintf(buf, "idx\tmean\tstd");
					file1.WriteLine(fname1, buf, "a+");

					m_idx = mi*ds.width_plus;
					for( wi=0; wi<=para.width; wi++ )
					{
						lr_idx = m_idx+wi;
						sprintf(buf, "%d\t%.15f\t%.15f", wi-ds.width_half, dd[lr_idx], dd[ds.m_wp+lr_idx] );
						file1.WriteLine(fname1, buf, "a+");
					}
				}

				strcpy(buf,"\n");
				file1.WriteLine(fname1, buf, "a+");
				for(int j = 0; j < ds.gibbs_state_size; j++)
				{
					int en = gibbs_state[j].en;
					int loc = gibbs_state[j].loc;
					int pol = gibbs_state[j].pol;
					sprintf(buf,"%s\t%d\t%d\t%d",gs.id_chr_map[en-1].c_str(), loc, pol, ds.cluster_num);
					file1.WriteLine(fname1, buf, "a+");
					file2.WriteLine(fname2, buf, "a+");
				}

				file1.LastWrite();
				file1.Release();
				file2.LastWrite();
				file2.Release();

				for(int j = 0; j < ds.gibbs_state_size;j++)
				{
					filter_status_all[ gibbs_state[j].id ] = 1;
					ds.done_list_size++;
				}    
				ds.cluster_num++;
			}
			if( ds.filter_num_total-ds.done_list_size > MAX_SEED_SIZE )
			{
				flag_cycle = 1;
			}
			else
			{
				flag_cycle = 0;
			}
		}

		MPI_Bcast( &flag_cycle, 1, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Bcast( &ds.cluster_num, 1, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Type_contiguous( ds.filter_num, MPI_INT, &gtype3 );
		MPI_Type_commit(&gtype3);
		MPI_Scatterv( filter_status_all, g_size, g_displ, MPI_INT, filter_status, 1, gtype3, 0, MPI_COMM_WORLD );
		MPI_Type_free(&gtype3);
	}

	delete[] mark;
	delete[] filter;
	delete[] mark_sig;
	delete[] dd;
	delete[] ds.logr_idx;
	delete[] g_size;
	delete[] g_displ;
	delete[] norm_node;
	delete[] gibbs_data;

	if( mpi_rank!=0 )
	{
		delete[] ds.mark_sig_idx;
		delete[] not_done_ids;
		delete[] filter_status;
	}

	if( mpi_rank==0 )
	{
		delete[] ds.mark_sig_idx_total;
		delete[] mpi_counts;
		delete[] node_name_all;
		delete[] all_norm;
		delete[] gibbs_state;
		delete[] new_gibbs_state;
		delete[] logr_sums;
		delete[] size_gibbs_mpi;
		delete[] gibbs_data_all;
		delete[] displ_gibbs_mpi;
	}
}

int main(int argc, char* argv[]) 
{
	int rc;
	int mpi_rank;

	rc = MPI_Init(&argc,&argv);
	if( rc != MPI_SUCCESS )
	{
		printf( "Error starting MPI program. Terminating\n" );
		MPI_Abort(MPI_COMM_WORLD,rc);
	}

	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	para.filter = "filter.idx";
	para.mark_info = "mark.idx";
	para.output_name = "EpiSig";
	para.pval_cutoff = 1e-5; 
	para.stat_half_window_size = 1000;
	para.maxima_half_window_size = 2500;
	para.overlap_half_window_size = 5000;
	para.width = 40; 
	para.wandering_dist = 20; 
	para.std_factor = 1.75; 
	para.prior_bg = 0.01;

	int next_option;
	const char* const short_options = "hgte:m:o:p:s:x:r:w:d:f:b:n:u:i:k:";
	const struct option long_options[]={
		{"help",0,NULL,'h'},
		{"global",0,NULL,'g'},
                {"global_background",0,NULL,'t'},
                {"enriched_region",1,NULL,'e'},
                {"background_region",1,NULL,'k'},
                {"mark_info",1,NULL,'m'},
                {"output_name",1,NULL,'o'},
                {"pval_cutoff",1,NULL,'p'},
                {"stat_half_window_size",1,NULL,'s'},
                {"maxima_half_window_size",1,NULL,'x'},
                {"overlap_half_window_size",1,NULL,'r'},
                {"width",1,NULL,'w'},
                {"max_wandering_dist",1,NULL,'d'},
                {"std_factor",1,NULL,'f'},
                {"prior_bg",1,NULL,'b'},
                {"thread_num",1,NULL,'n'},
                {"random_seed",1,NULL,'u'},
                {"genome_name",1,NULL,'i'},
                {NULL,0,NULL,0}
        };

	para.enriched_region = 0;
	para.global = 0;
	para.flag_genome = 0;
	para.rand_seed = 0;
	do{
		next_option = getopt_long(argc,argv,short_options,long_options,NULL);
		switch(next_option)
		{
			case 'h':
				print_usage();
                                return 0;
			case '?':
				print_usage();
				return 0;
			case ':':
				print_usage();
				return 0;
			case 'g':
				para.global = 1;
				break;
                        case 't':
                                para.global_background = 1;
                                break;
			case 'u':
				para.rand_seed = atoi(optarg);
				break;
                        case 'k':
                                para.region_background = 1;
                                para.background_region = string(optarg);
                                break;
			case 'e':
				para.filter = optarg;
				para.enriched_region = 1;
				break;
			case 'm':
				para.mark_info = optarg;
				break;
			case 'o':
				para.output_name = optarg;
				break;
			case 'p':
				para.pval_cutoff = atof(optarg);
				if(para.pval_cutoff < 0 || para.pval_cutoff > 0.01)
				{
					printf("P-value cutoff should between 0 and 0.01\n");
					return 0;
				}
				break;
			case 's':
				para.stat_half_window_size = atoi(optarg);
				if(para.stat_half_window_size < 0)
				{
					printf("Negative stat_half_window_size!\n");
					return 0;
				}
				if(para.stat_half_window_size % 100 != 0)
				{
					printf("Stat_half_window_size should be divided by 100!\n");
					return 0;
				}
				break;
			case 'x':
				para.maxima_half_window_size = atoi(optarg);
				if(para.maxima_half_window_size < 0)
				{
					printf("Negative maxima_half_window_size!\n");
					return 0;
				}
				if(para.maxima_half_window_size % 100 != 0)
				{
					printf("Maxima_half_window_size should be divided by 100!\n");
					return 0;
				}
				break;
			case 'r':
				para.overlap_half_window_size = atoi(optarg);
				if(para.overlap_half_window_size < 0)
				{
					printf("Negative overlap_half_window_size!\n");
					return 0;
				}
				if(para.overlap_half_window_size % 100 != 0)
				{
					printf("Overlap_half_window_size should be divided by 100!\n");
					return 0;
				}
				break;
			case 'w':
				para.width = atoi(optarg);
				if(para.width < 20 || para.width > 200)
				{
					printf("Width should between 20 and 200!\n");
					return 0;
				}
				if(para.width % 10 != 0)
				{
					printf("Window size should be divided by 10!\n");
					return 0;
				}
				break;
			case 'd':
				para.wandering_dist = atoi(optarg);
				if(para.wandering_dist < 10 || para.wandering_dist > 200)
				{
					printf("Max_wandering_dist should between 10 and 200!\n");
					return 0;
				}
				if(para.wandering_dist % 10 != 0)
				{
					printf("Window size should be divided by 100!\n");
					return 0;
				}
				break;
			case 'f':
				para.std_factor = atof(optarg);
				if(para.std_factor < 0 || para.std_factor > 10)
				{
					printf("std_factor should between 0.0 and 10.0!\n");
					return 0;
				}
				break;
			case 'b':
				para.prior_bg = atof(optarg);
				if(para.prior_bg < 0 || para.prior_bg > 0.5)
				{
					printf("prior_bg should between 0.0 and 0.5!\n");
					return 0;
				}
				break;
			case 'i':
				para.genome_name = string(optarg);
				para.flag_genome = 1;
				break;
			default:
				break;
		}
	}while(next_option !=-1);

	print_para(para);
	if( para.flag_genome == 0 )
	{
		printf( "Genome name not specified, please use -i on running\n" );
                exit(1);
        }

	gettimeofday(&tpstart,NULL);

	if( mpi_rank == 0 )
	{
		GS gs_m;
		GS_init( gs_m, para.genome_name );
		if( para.global_background==0 && para.region_background==0 )
		{
			print_usage();
			printf("You must specify an argument: [--global_background | --background_region]\n");
			return 0;
		}
		if( para.global_background==1 && para.region_background==1 )
		{
			print_usage();
			printf("You must specify only one argument of: [--global_background | --background_region]\n");
			return 0;
		}

		char cmd[256];
		if(para.enriched_region == 0 && para.global == 0)
		{
			print_usage();
			printf("You must specify an argument: [--global | --enriched_region]\n");
			return 0;
		}
		if(para.enriched_region == 1 && para.global == 1)
		{
			print_usage();
			printf("You must specify only one argument of [--global, --enriched_region]\n");
			return 0;
		}
		if(para.enriched_region == 0)
		{
			get_normal_params( para.stat_half_window_size, para.mark_info.c_str(), para.background_region, para.global_background, &gs_m );
			printf("\ncomputing mean and std:%s\n", cmd);
			showtime("");

			get_sig_locs( para.stat_half_window_size, para.maxima_half_window_size, para.overlap_half_window_size, para.pval_cutoff, &gs_m );
			printf("\ncomputing mean and std:%s\n", cmd);
			showtime("");
		}
		printf("\ncomputing background parameters:\n");
		get_background_dist_global( BACKGROUND_BUFFER, para.mark_info.c_str(), para.filter.c_str(), &gs_m );
		showtime("");

		int ret = 0;
		printf("\nextract mark signal for enriched regions:\n");
		sprintf(cmd,"./filter_ad background_.para %s %d %d %s",para.filter.c_str(), para.width, para.wandering_dist, para.genome_name.c_str() );
		printf("%s\n",cmd);
		ret = system(cmd);
		if( ret<0 )
		{
			printf("extract mark signal failed!\n");
			printf("%s\n",cmd);
			exit(1);
		}
		showtime("");
	}
	para.filter = "filter_ad.idx";

	MPI_Barrier( MPI_COMM_WORLD );

	main_gibbs_cluster();

	if( mpi_rank == 0 )
	{
		showtime("Total running time");
	}

	MPI_Finalize();

	return 0;
}
