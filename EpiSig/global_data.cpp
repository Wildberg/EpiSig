#include<cmath>
#include "global_data.h"

#if CHROMASIG_MPI
void ds_init( _para* parameter, _data_size* ds, int mpi_size, int mpi_rank )
{
	ds->mpi_num_size = ds->filter_num_total/mpi_size;
        ds->mpi_num_extra = ds->filter_num_total%mpi_size;
        if( mpi_rank<ds->mpi_num_extra )
        {
                ds->filter_num = ds->mpi_num_size+1;
                ds->filter_start = mpi_rank * (ds->mpi_num_size+1);
        }
        else
        {
                ds->filter_num = ds->mpi_num_size;
                ds->filter_start = mpi_rank*ds->mpi_num_size + ds->mpi_num_extra;
        }

        ds->width_half = parameter->width/2;
        ds->width_plus = parameter->width+1;
        ds->wandering_half = parameter->wandering_dist/2;
        ds->wandering_plus = parameter->wandering_dist+1;
        ds->w_size = parameter->wandering_dist+parameter->width+1;
        ds->m_wp = long(ds->width_plus) * ds->mark_num;
        ds->sites_size = 100;
        ds->cache_size = 0;
        ds->p_ratio = 2*parameter->prior_bg/(1-parameter->prior_bg);
        ds->exp_A = exp((parameter->width+1)*parameter->std_factor*parameter->std_factor*-0.5);

	ds->done_list_size = 0;
	ds->gibbs_state_size = 0;
	ds->cluster_num = 0;

        ds->mark_sig_idx = new long[ds->filter_num];
        for( long i=0; i<ds->filter_num; i++ )
        {
                ds->mark_sig_idx[i] = i*ds->w_size*ds->mark_num;
        }

        ds->logr_idx = new long[ds->mark_num];
        for( long i=0; i<ds->mark_num; i++ )
        {
                ds->logr_idx[i] = i*ds->width_plus+ds->width_half;
        }

        ds->mark_sig_idx_total = NULL;
        if( mpi_rank==0 )
        {
                ds->mark_sig_idx_total = new long[ds->filter_num_total];
                for( long i=0; i<ds->filter_num_total; i++ )
                {
                        ds->mark_sig_idx_total[i] = i*ds->w_size*ds->mark_num;
                }
        }
}
#endif
