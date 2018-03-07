#ifndef DATA_SIZE_H
#define DATA_SIZE_H

#include<string>

using namespace std;

struct _data_size
{
#if CHROMASIG_MPI
	long* mark_sig_idx_total;
	long filter_start;
        long filter_num_total;
	int mpi_num_size;
	int mpi_num_extra;
	int cache_size;
	int not_done_ids_size_total;
#endif
	long* mark_sig_idx;
	long* logr_idx;
	double p_ratio;
	double exp_A;
	long m_wp;
        long filter_num;
        int mark_num;
        int w_size;
        int sites_size;
        int norm_size;
        int done_list_size;
        int gibbs_state_size;
        int not_done_ids_size;
        int cluster_num;
        int copy_gibbs_state_size;
        int seed_size;
	int width_half;
	int width_plus;
	int wandering_plus;
	int wandering_half;
};

struct _para{
        string filter;
        string mark_info;
        string output_name;
	string genome_name;
	string background_region;
        double pval_cutoff;
        int stat_half_window_size;
        int maxima_half_window_size;
        int overlap_half_window_size;
        int width;
        int wandering_dist;
        double std_factor;
        double prior_bg;
	int rand_seed;
	int flag_genome;
        int enriched_region;
        int global;
	int global_background;
	int region_background;
#ifndef CHROMASIG_MPI
        int max_thread; // number of thread in multithreading
#endif
};

#ifdef CHROMASIG_MPI
void ds_init( _para* parameter, _data_size* ds, int mpi_size, int mpi_rank );
#endif

#endif
