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

#include "utils.h"
#include "FileUtils.h"
#include "DataUtils.h"
#include "fast_pmf_cycle.h"
#include "EpiSig_v1.0.h"
#include "filter.h"

#define P_DEBUG 1

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
	buffer += "-t, --global_background     EpiSig_v1.0 will use the whole genome as background for enriched regions \n";
	buffer += "-e, --enriched_region     provide an enriched region file, for example, filter.idx\n";
	buffer += "-k, --background_region     provide an background region file for enriched region identification, for example, encode_hg18.bed\n";
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
	buffer += "-u, --random-seed	random seed used in EpiSig_v1.0\n";
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
	printf("threads: %d\n", para.max_thread);
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

int compute_motif( _data_size* ds, int exclude_id, _gibbs_state* gibbs, int gibbs_size, float* mark_sig, _logr_sums* logr, _motif_dist* mdist, _filter* filter )
{
	int mi, wi, gi;
	int is_defined; 
	int id;
	int pol, idx, offset;
	int num_elements;
	long d_idx, lr_idx;
	double val;
	double mean;

	int* num_elements_offset = new int[ds->m_wp];

	memset(num_elements_offset,0,ds->m_wp*sizeof(int));

	for( d_idx=0; d_idx<ds->m_wp; d_idx++ )
	{
		mdist[d_idx].mean = logr[d_idx].sum;
		mdist[d_idx].std = logr[d_idx].sum2;
	}

	is_defined = 0;
	if( exclude_id > -1 )
	{
		for( gi=0; gi<gibbs_size; gi++)
		{
			id = gibbs[gi].id;
			if(id == exclude_id)
			{
				pol = gibbs[gi].pol;
				idx = gibbs[gi].id;
				offset = (gibbs[gi].loc - filter[idx].pos) / 100;

				d_idx = ds->mark_sig_idx[idx]+ds->width_half+ds->wandering_half+offset;
				for( mi=0; mi<ds->mark_num; mi++)
				{
					for( wi=-ds->width_half; wi<=ds->width_half; wi++ )
					{
						lr_idx = ds->logr_idx[mi]+wi;
						val = mark_sig[d_idx+wi*pol];

						mdist[lr_idx].mean -= val;
						mdist[lr_idx].std -= val*val;
						num_elements_offset[lr_idx] = -1;
					}
					d_idx += ds->w_size;
				}
				is_defined = 1;
				break;
			}
		}
	}

	for( lr_idx=0; lr_idx<ds->m_wp; lr_idx++ )
	{
		num_elements = logr[lr_idx].count;
		if(num_elements <= 1)
		{
			printf("Zero or negative element motif!\n");
			exit(1);
		}

		if(is_defined) num_elements += num_elements_offset[lr_idx];
		mdist[lr_idx].mean /= num_elements;

		mean = mdist[lr_idx].mean;
		mdist[lr_idx].std -= num_elements*mean*mean;
		mdist[lr_idx].std /= (num_elements - 1);
	}

	delete[] num_elements_offset;

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
		offset = (gibbs[gi].loc - filter[idx].pos) / 100;

		d_idx = ds->mark_sig_idx[idx]+ds->width_half+ds->wandering_half+offset;
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

int compute_motif_full( _data_size* ds, _gibbs_state *gibbs, int exclude_id, float* mark_sig, _logr_sums* logr, _motif_dist* motif_dist, Data *dd, Scratch *ss, _mark *mark,_filter *filter )
{
	get_running_stats( ds, gibbs, ds->gibbs_state_size, mark_sig, logr, filter );
	compute_motif( ds, exclude_id, gibbs, ds->gibbs_state_size, mark_sig, logr, motif_dist, filter );

	int mi, wi;
	int lr_idx;
	double mean, std;

	lr_idx = 0;
	for( mi=0; mi<ds->mark_num; mi++ )
	{
		for( wi=0; wi<=para.width; wi++,lr_idx++)
		{
			mean = motif_dist[lr_idx].mean;
			std = motif_dist[lr_idx].std;

			for(int j = 0; j < para.max_thread;j++)
			{
				set_motif_mean(&dd[j],&ss[j],mi,wi,mean);
				set_motif_std (&dd[j],&ss[j],mi,wi,std);
			}
		}

		for(int j = 0; j < para.max_thread;j++)
		{
			set_bg_mean(&dd[j],&ss[j],mi,mark[mi].mean);
			set_bg_std(&dd[j],&ss[j],mi,mark[mi].std);
		}
	}

	for(int j = 0; j < para.max_thread;j++)
	{
		set_motif_mean_of_std(&dd[j],&ss[j]);
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

	get_running_stats( ds, gibbs, gibbs_size, mark_sig, logr, filter );
	compute_motif( ds, -1, gibbs, gibbs_size, mark_sig, logr, mdist, filter );

	d_idx = 0;
	score = 0;
	for( mi=0; mi<ds->mark_num; mi++ )
	{
		avg_std = 0;
		for( wi=0; wi<=para.width; wi++,d_idx++ )
		{
			avg_std += sqrt(mdist[d_idx].std);
			p_mean[wi] = fabs(mdist[d_idx].mean);
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

void *distance_thread(void *arg)
{
	_arg_thread* arg_thread = (_arg_thread*) arg;

	THREAD_RAND* thread_rand = arg_thread->thread_rand;
	_data_size* ds = arg_thread->ds;
	int size_sites = arg_thread->ds->sites_size;
	_filter* first_sets = arg_thread->first_sets;
	_filter* filter = arg_thread->filter;
	int size_norm = arg_thread->ds->norm_size;
	int mark_num = arg_thread->ds->mark_num;
	float* mark_sig = arg_thread->mark_sig;
	int* undo_list = arg_thread->undo_list;
	int func_id = arg_thread->func_id;

	int first_id;
	int ni, mi, ui, fi, wi;
	long idx;
	long idx_f;
	double dist;
	double diff;
	double score;
	_distance* undo_distance = new _distance[size_norm];
	_gibbs_state gibbs[MAX_SEED_SIZE];
	int gibbs_size;

	best_score[func_id] = -1;
	for( ni=0; ni<size_sites; ni++ )
	{
		if( ni%para.max_thread != func_id )
			continue;

		first_id = first_sets[ni].pos;
		for( ui=0; ui<size_norm; ui++ )
		{
			fi = undo_list[ui];

			if( fi == first_id )
			{
				dist = 99999999;
			}
			else
			{
				dist = 0;
				idx_f = ds->mark_sig_idx[first_id]+ds->wandering_half;
				idx = ds->mark_sig_idx[fi]+ds->wandering_half;
				for( mi=0; mi<mark_num; mi++ )
				{
					for( wi=0; wi<=para.width; wi++,idx++,idx_f++ )
					{
						diff = ((double)mark_sig[idx])-((double)mark_sig[idx_f]);
						dist += diff*diff;
					}

					idx += para.wandering_dist;
					idx_f += para.wandering_dist;
				}
			}

			(undo_distance+ui)->id = fi;
			(undo_distance+ui)->distance = dist;
		}

		tournament_sort(ds,first_id,undo_list,undo_distance,gibbs,gibbs_size,filter,thread_rand);
		score = get_seed_score(ds, gibbs, gibbs_size, mark_sig, filter);
		if(best_score[func_id] < score)
		{
			best_score[func_id] = score;
			thread_seed_size[func_id] = gibbs_size;
			memcpy(best_gibbs[func_id],gibbs,sizeof(_gibbs_state)*gibbs_size);
		}
	}

	delete[] undo_distance;

	return NULL;
}

int get_not_done_ids( _data_size* ds, int *done_list, int *not_done_ids )
{
	int cnt;
	int* tag = new int[ds->filter_num];

	memset(tag,0,(ds->filter_num)*sizeof(int));

	for(int i=0;i<ds->done_list_size;i++)
	{
		tag[done_list[i]] = 1;
	}

	cnt=0;
	for(int i=0;i<ds->filter_num;i++)
	{
		if(tag[i] == 0)
		{
			not_done_ids[cnt] = i;
			cnt++;
		}
	}
	delete[] tag;

	return cnt;
}

void init_gibbs_seed( _data_size* ds, _gibbs_state* seed, double *all_norm, _filter *filter, float *mark_sig, int* done_list, THREAD_RAND* thread_rand, GS* gs )
{
	_norm* norm = new _norm[ds->filter_num];
	int* tag = new int[ds->filter_num];
	int* undo_list = new int[ds->filter_num];
	int size_norm;
	int id;
	double my_max;
	int my_idx;

	int num_sites = 100;
	const int percentile = 0;
	const double ratio = 0.25;
	int interval;
	_filter* first_sets;

	_arg_thread arg_thread[MAX_THREAD];
	pthread_t id_thread[MAX_THREAD];
	int ret[MAX_THREAD];

	memset(tag,0,ds->filter_num*sizeof(int));

	for( int i=0; i<ds->done_list_size; i++ )
	{
		tag[done_list[i]] = 1;
	}

	size_norm = 0;
	for(int i=0; i<ds->filter_num; i++)
	{
		if(tag[i] == 0)
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
	if (interval == 0)
	{
		interval = 1; 
	}

	if(ds->sites_size > ds->norm_size)
	{
		ds->sites_size = ds->norm_size;
	}
	first_sets = new _filter[ds->sites_size];
	for( int i=0; i<ds->sites_size; i++ )
	{
		id = norm[percentile+i*interval].id;
		first_sets[i].en = filter[id].en;
		first_sets[i].pos = id;
	}

	qsort(first_sets,ds->sites_size,sizeof(_filter),cmp_filter);
	for( int ti=0; ti<para.max_thread; ti++ )
	{
		arg_thread[ti].ds = ds;
		arg_thread[ti].first_sets = first_sets;
		arg_thread[ti].filter = filter;
		arg_thread[ti].mark_sig = mark_sig;
		arg_thread[ti].undo_list = undo_list;
		arg_thread[ti].func_id = ti;
		arg_thread[ti].thread_rand = thread_rand;
	}
	for( int ti=0; ti<para.max_thread; ti++ )
	{
		ret[ti] = pthread_create(&id_thread[ti],NULL,distance_thread,(void*)(&arg_thread[ti]));
	}
	for( int ti=0; ti<para.max_thread; ti++ )
	{
		pthread_join(id_thread[ti],NULL);	
	}

	delete[] norm;
	showtime("computing euclidean distance");

	my_max = 0;
	my_idx = 0;
	for(int i = 0; i < para.max_thread; i++)
	{
		if(my_max < best_score[i])
		{
			my_max = best_score[i];
			my_idx = i;
			ds->seed_size = thread_seed_size[i];
		}
	}

#if P_DEBUG
	printf("\nBest seed score = %.4f, done:%d, not done:%d, seed set:\n", my_max, ds->done_list_size, ds->not_done_ids_size);
	fflush(stdout);
#endif

	memcpy(seed,best_gibbs[my_idx],sizeof(_gibbs_state)*ds->seed_size);

#if P_DEBUG
	for(int i = 0; i < ds->seed_size; i++)
	{
		//printf("%d\t%s\t%d\t%d\n",seed[i].id,id2chr(seed[i].en).c_str(),seed[i].loc,seed[i].pol);
		printf("%d\t%s\t%d\t%d\n",seed[i].id,gs->id_chr_map[seed[i].en].c_str(),seed[i].loc,seed[i].pol);
	}
	fflush(stdout);
#endif

	delete[] first_sets;
	delete[] tag;
	delete[] undo_list;
}

void *func(void* arg)
{
	_arg_thread* arg_thread = (_arg_thread*)arg;
	int mark_num = arg_thread->ds->mark_num;
	_gibbs_state* gibbs = arg_thread->gibbs;
	int* tag = arg_thread->tag;
	int size_not_list = arg_thread->ds->not_done_ids_size;
	int* not_done_ids = arg_thread->undo_list;
	int func_id = arg_thread->func_id;
	float* mark_sig = arg_thread->mark_sig;
	Data* dd = arg_thread->dd;
	Scratch* ss = arg_thread->ss;
	_filter* filter = arg_thread->filter;
	_data_size* ds = arg_thread->ds;

	int ni, mi, wi;
	long lr_idx;
	int cnt;
	int exclude_id;

	cnt = 0;
	for( ni=0; ni<size_not_list; ni++ )
	{
		if( ni%para.max_thread != func_id ) continue;

		exclude_id = not_done_ids[ni];
		lr_idx = ds->mark_sig_idx[exclude_id];
		for( mi=0; mi<mark_num; mi++ )
		{
			for( wi=0; wi<ds->w_size; wi++,lr_idx++ )
			{
				set_region(dd,ss,mi,wi,mark_sig[lr_idx]);
			}
		}
		compute_pmf(dd,ss,ds);
		if( !get_votes(ss,ds) )
		{
			tally_votes(ss,ds);
			gibbs[exclude_id].id = exclude_id;
			gibbs[exclude_id].en =  filter[exclude_id].en;
			gibbs[exclude_id].loc = filter[exclude_id].pos + 100*get_return_loc(ss);
			gibbs[exclude_id].pol = get_return_pol(ss);
			tag[exclude_id] = 1;
		}
		cnt++;
	}

	return NULL;
}

void load_mark_data( float *mark_sig, string dir_data, _mark *mark, _data_size* ds )
{
	unsigned long ival = ds->mark_num*ds->w_size;
	unsigned long cnt = ds->filter_num*ds->w_size;
	unsigned long r_cnt;
	string fn;
	float* mark_idx;
	float* bin_data;
	float* bin_idx;
	FILE* fh;

	bin_data = new float[cnt];
	for( int i=0; i<ds->mark_num; i++ )
	{
		fn = dir_data+"/"+mark[i].name+".bin";
		if( (fh=fopen(fn.c_str(),"rb")) == NULL )
		{
			printf( "Unable to open file: %s\n", fn.c_str() );
			exit(1);
		}

		r_cnt = fread(bin_data,sizeof(float),cnt,fh);
		if( r_cnt != cnt )
		{
			printf( "Data size %lu not match for %s\n", cnt, fn.c_str() );
			printf( "mark=%d\twindow_size=%d\tfilter_num=%ld\tcnt=%lu\tr_cnt=%lu\n", i, ds->w_size, ds->filter_num, cnt, r_cnt );
			fclose(fh);
			exit(1);
		}
		fclose(fh);

		bin_idx = bin_data;
		mark_idx = mark_sig+i*ds->w_size;
		for( int fi=0; fi<ds->filter_num; fi++,mark_idx+=ival,bin_idx+=ds->w_size )
		{
			memcpy(mark_idx,bin_idx,ds->w_size*sizeof(float));
		}
	}

	delete[] bin_data;
}

void main_gibbs_cluster(GS* gs)
{
	THREAD_RAND thread_rand(para.max_thread,para.rand_seed);
	struct _data_size ds; 

	_mark* mark;
	_filter* filter;
	Data* dd = new Data[para.max_thread];
	Scratch* ss = new Scratch[para.max_thread];

	int cycle_num = NUM_CYCLES;
	ds.mark_num = load_marks("background_.para",&mark); 
	ds.filter_num = load_filter(para.filter,&filter,gs);
	if(ds.mark_num <= 0 || ds.filter_num <= 0 )
	{
		printf("Fatal errors! Failed to load background_.para or %s\n",para.filter.c_str());
		return;
	}

	ds.width_half = para.width/2;
	ds.width_plus = para.width+1;
	ds.wandering_half = para.wandering_dist/2;
	ds.wandering_plus = para.wandering_dist+1;
	ds.w_size = para.wandering_dist+para.width+1;
	ds.m_wp = ds.width_plus * ds.mark_num;
	ds.p_ratio = 2*para.prior_bg/(1-para.prior_bg);
	ds.exp_A = exp((para.width+1)*para.std_factor*para.std_factor*-0.5);
	ds.mark_sig_idx = new long[ds.filter_num];
	for( long i=0; i<ds.filter_num; i++ )
	{
		ds.mark_sig_idx[i] = i*ds.w_size*ds.mark_num;
	}
	ds.logr_idx = new long[ds.mark_num];
	for( int i=0; i<ds.mark_num; i++ )
	{
		ds.logr_idx[i] = i*ds.width_plus+ds.width_half;
	}

	float* mark_sig = new float[ds.mark_num*ds.filter_num*ds.w_size];

	load_mark_data( mark_sig, "out_data_ad", mark, &ds );

	for(int i = 0; i < para.max_thread; i++)
	{
		init( &dd[i], &ss[i], &ds );
	}

	qsort(filter, ds.filter_num, sizeof(_filter), cmp_filter);

	double *all_norm = new double[ds.filter_num];
	printf("\nComputing Euclidean norm.................\n");
	fflush(stdout);

	get_euclidean_norm(&ds,mark_sig,all_norm);
	showtime("Computing Euclidean Norm");

	ds.done_list_size = 0;
	int *done_list = new int[ds.filter_num]; 
	ds.gibbs_state_size = 0;
	_gibbs_state *gibbs_state = new _gibbs_state[ds.filter_num];
	ds.not_done_ids_size = ds.filter_num;
	int *not_done_ids = new int[ds.filter_num];
	ds.cluster_num = 0;
	_gibbs_state* new_gibbs_state = new _gibbs_state[ds.filter_num];
	int* tag = new int[ds.filter_num];

	_gibbs_state copy_gibbs_state[MAX_SEED_SIZE];
	_gibbs_state seed[MAX_SEED_SIZE];

	_logr_sums *logr_sums = new _logr_sums[ds.mark_num*ds.width_plus];
	_motif_dist *motif_dist = new _motif_dist[ds.mark_num*ds.width_plus];

	_arg_thread arg_thread[MAX_THREAD];
	pthread_t id[MAX_THREAD];
	int ret[MAX_THREAD];

	int round_n = 0;

	while( ds.filter_num-ds.done_list_size > MAX_SEED_SIZE )
	{
		round_n++;
#if P_DEBUG
		printf("\nDebug: filter_num %ld done_list_size %d MAX_SEED_SIZE %d\n", ds.filter_num, ds.done_list_size, MAX_SEED_SIZE ); 
#endif
		printf("\ncluster = %d\n",ds.cluster_num);
		fflush(stdout);

		ds.not_done_ids_size = get_not_done_ids( &ds, done_list, not_done_ids );

#if P_DEBUG
		printf("\nDebug: init_gibbs_seed round %d\n", round_n );
		fflush(stdout);
#endif

		init_gibbs_seed( &ds, seed, all_norm, filter, mark_sig, done_list, &thread_rand, gs );

		memcpy(copy_gibbs_state,seed,ds.seed_size*sizeof(_gibbs_state));
		ds.copy_gibbs_state_size = ds.seed_size;
		memcpy(gibbs_state,seed,ds.seed_size*sizeof(_gibbs_state)); 
		ds.gibbs_state_size = ds.seed_size;

		showtime("init_gibbs_seed");
		for(int rd = 1; rd <= cycle_num;rd++)
		{
			printf("Cycle %d, alignment set size %d ...\n",rd,ds.gibbs_state_size);
			fflush(stdout);

			compute_motif_full(&ds,gibbs_state,-1, mark_sig, logr_sums, motif_dist, dd, ss, mark, filter );
			memset(tag,0,ds.filter_num*sizeof(int));
			for( int i=0; i<para.max_thread; i++ )
			{
				arg_thread[i].ds = &ds;
				arg_thread[i].dd = dd+i;
				arg_thread[i].ss = ss+i;
				arg_thread[i].filter = filter;
				arg_thread[i].gibbs = new_gibbs_state;
				arg_thread[i].mark_sig = mark_sig;
				arg_thread[i].undo_list = not_done_ids;
				arg_thread[i].tag = tag;
				arg_thread[i].func_id = i;
				arg_thread[i].size_sites = rd;
			}
			for( int i=0; i<para.max_thread; i++ )
			{
				ret[i] = pthread_create(&id[i],NULL,func,(void*)(&arg_thread[i]));
			}
			for( int i=0; i<para.max_thread; i++ )
			{
				pthread_join(id[i], NULL);
			}

			ds.gibbs_state_size = 0;  
			for(int j = 0; j< ds.filter_num; j++)
			{
				if(tag[j] == 1)
				{
					gibbs_state[ds.gibbs_state_size].id = new_gibbs_state[j].id;
					gibbs_state[ds.gibbs_state_size].en = new_gibbs_state[j].en;
					gibbs_state[ds.gibbs_state_size].loc = new_gibbs_state[j].loc;
					gibbs_state[ds.gibbs_state_size].pol = new_gibbs_state[j].pol;
					ds.gibbs_state_size++;
				}    
			}
			showtime("This cycle finished");
			printf("Update best alignment set, number:%d\n", ds.gibbs_state_size);
			fflush(stdout);

			if(ds.gibbs_state_size <= ds.seed_size)
			{
				printf("Less than seed size, use seed set:%d\n",ds.copy_gibbs_state_size);
				fflush(stdout);
				for(int i = 0; i < ds.copy_gibbs_state_size;i++)
				{
					done_list[ds.done_list_size] = copy_gibbs_state[i].id;
					ds.done_list_size++;
				}    
				break;
			}
		}

		if(ds.gibbs_state_size > ds.seed_size)
		{
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
			compute_motif( &ds, -1, gibbs_state, ds.gibbs_state_size, mark_sig, logr_sums, motif_dist, filter );

			int mi, wi;
			int m_idx;
			long lr_idx;
			for( mi=0; mi<ds.mark_num; mi++)
			{
				sprintf(buf, "\nmod = %s", mark[mi].name.c_str());
				file1.WriteLine(fname1, buf, "a+");
				sprintf(buf, "idx\tmean\tstd");
				file1.WriteLine(fname1, buf, "a+");

				m_idx = mi*(para.width+1);
				for( wi=0; wi<=para.width; wi++ )
				{
					lr_idx = m_idx+wi;
					sprintf(buf, "%d\t%.15f\t%.15f", wi-para.width/2, motif_dist[lr_idx].mean, motif_dist[lr_idx].std );
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
				sprintf(buf,"%s\t%d\t%d\t%d",gs->id_chr_map[en-1].c_str(), loc, pol, ds.cluster_num);
				file1.WriteLine(fname1, buf, "a+");
				file2.WriteLine(fname2, buf, "a+");
			}

			file1.LastWrite();
			file1.Release();
			file2.LastWrite();
			file2.Release();

			for(int j = 0; j < ds.gibbs_state_size;j++)
			{
				done_list[ds.done_list_size] = gibbs_state[j].id;
				ds.done_list_size++;
			}    
			ds.cluster_num++;
		}
	}

	delete[] mark;
	delete[] filter;
	delete[] mark_sig;
	delete[] dd;
	delete[] ss;
	delete[] done_list;
	delete[] gibbs_state;
	delete[] not_done_ids;
	delete[] all_norm;
	delete[] tag;
	delete[] new_gibbs_state;
	delete[] logr_sums;
	delete[] motif_dist;
	delete[] ds.mark_sig_idx;
	delete[] ds.logr_idx;
}

int main(int argc, char* argv[]) 
{
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
	para.max_thread = 4;

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

	para.region_background = 0;
	int enriched_region = 0;
	int global = 0;
	para.rand_seed = 0;
	int flag_genome = 0;
	string genome_name;
	// get options
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
				global = 1;
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
				enriched_region = 1;
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
			case 'n':
				para.max_thread = atoi(optarg);
				if(para.max_thread < 1 || para.max_thread > 96)
				{
					printf("max thread number should between 1 and 96!\n");
					return 0;
				}
				break;
			case 'i':
				flag_genome = 1;
				genome_name = string(optarg);
				break;
			default:
				break;
		}
	}while(next_option !=-1);

	print_para(para);
	if( flag_genome == 0 ) 
	{
		printf( "Genome name not specified, please use -i on running\n" );
		exit(1);
	}
	GS gs;
	GS_init( gs, genome_name );

	gettimeofday(&tpstart,NULL);

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
	if(enriched_region == 0 && global == 0)
	{
		print_usage();
		printf("You must specify an argument: [--global | --enriched_region]\n");
		return 0;
	}
	if(enriched_region == 1 && global == 1)
	{
		print_usage();
		printf("You must specify only one argument of [--global, --enriched_region]\n");
		return 0;
	}
	if(enriched_region == 0) 
	{
		get_normal_params( para.stat_half_window_size, para.mark_info.c_str(), para.background_region, para.global_background, &gs );
		printf("\ncomputing mean and std:%s\n", cmd);
		showtime("");

		printf("\nfind enriched regions:\n");
		get_sig_locs( para.stat_half_window_size, para.maxima_half_window_size, para.overlap_half_window_size, para.pval_cutoff, &gs );
		showtime("");
	}
	printf("\ncomputing background parameters:\n");
	get_background_dist_global( BACKGROUND_BUFFER, para.mark_info.c_str(), para.filter.c_str(), &gs );
	showtime("");

	int ret = 0;
	printf("\nextract mark signal for enriched regions:\n");
	sprintf(cmd,"./filter_ad background_.para %s %d %d %s", para.filter.c_str(), para.width, para.wandering_dist, genome_name.c_str() );
	printf("%s\n",cmd);
	ret = system(cmd);
	para.filter = "filter_ad.idx";
	if( ret<0 )
	{
		printf("extract mark signal failed!\n");
		printf("%s\n",cmd);
		exit(1);
	}


	fflush(stdout); 

	showtime("");

	main_gibbs_cluster(&gs);

	showtime("");
	return 0;
}
