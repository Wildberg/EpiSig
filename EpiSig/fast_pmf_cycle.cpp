////////////////////////////////////////////////////////////////////////////////
//
// fast_pmf_cycle.c
//
////////////////////////////////////////////////////////////////////////////////
#include "fast_pmf_cycle.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern _para para;

////////////////////////////////////////////////////////////////////////////////
//
// GLOBALS
//
////////////////////////////////////////////////////////////////////////////////
int NONE_VAL = -9999;
double ABS_Z_SCORE_MAX = 10;
double ABS_Z_SCORE_MAX_SQ = 100;
double MIN_LOG_RATIO = -100;
double EXP_MAX = -100;

int get_return_pol(Scratch *s) {
  return s->return_pol;
}

int get_return_loc(Scratch *s) {
  return s->return_loc;
}

#ifndef CHROMASIG_MPI
double get_motif_mean(Data *d,Scratch *s,int mod, int pos) {
  return d->motif_mean[mod][pos];
}

void set_motif_mean(Data *d,Scratch *s,int mod, int pos, double val) {
  d->motif_mean[mod][pos] = val;
}

double get_motif_std(Data *d,Scratch *s,int mod, int pos) {
  return d->motif_std[mod][pos];
}

void set_motif_std(Data *d,Scratch *s,int mod, int pos, double val) {
  d->motif_std[mod][pos] = val;
}

double get_bg_mean(Data *d,Scratch *s,int mod) {
  return d->bg_mean[mod];
}

void set_bg_mean(Data *d,Scratch *s,int mod, double val) {
  d->bg_mean[mod] = val;
}

double get_bg_std(Data *d,Scratch *s,int mod) {
  return d->bg_std[mod];
}

void set_bg_std(Data *d,Scratch *s,int mod, double val) {
  d->bg_std[mod] = val;
}

double get_region(Data *d,Scratch *s,int mod, int pos) {
  return d->region[mod][pos];
}

void set_region(Data *d,Scratch *s,int mod, int pos, double val) {
  d->region[mod][pos] = val;
}

void set_motif_mean_of_std(Data *d,Scratch *s)
{
  int mod;
  int motif_index;

  for (mod = 0 ; mod < d->num_mod ; mod++) {
    // init
    d->mean_of_stds[mod] = 0;

    // compute mean
    for(motif_index = 0 ; motif_index <= d->width ; motif_index++) {
      d->mean_of_stds[mod] += sqrt(d->motif_std[mod][motif_index]);
    }
    d->mean_of_stds[mod] /= (d->width + 1);
  }
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
// init
//
////////////////////////////////////////////////////////////////////////////////
#ifndef CHROMASIG_MPI
void init( Data *d,Scratch *s,_data_size* ds )
#else
void init( Scratch *s,_data_size* ds )
#endif
{
	int i;

	////////////////////
	// main struct alloc
	////////////////////
	//d = (Data*) malloc(sizeof(Data));

	// set some values
#ifndef CHROMASIG_MPI
	d->num_mod     = ds->mark_num; // number of marks
	d->width       = para.width; // window size, 4k=40*100
	d->wander_dist = para.wandering_dist; // wander_dist, 20

	// motif_mean alloc
	d->motif_mean = (double**) malloc(sizeof(double*) * ds->mark_num);
	for(i=0 ; i < ds->mark_num ; i++) {
		d->motif_mean[i] = (double*) malloc(sizeof(double) * ds->width_plus);
	}

  // motif_std alloc
  d->motif_std = (double**) malloc(sizeof(double*) * ds->mark_num);
  for(i=0 ; i < ds->mark_num ; i++) {
    d->motif_std[i] = (double*) malloc(sizeof(double) * ds->width_plus);
  }

  // region alloc
  // int total_size = para.width + para.wandering_dist + 1; // 40+20=60
  d->region = (double**) malloc(sizeof(double*) * ds->mark_num);
  for(i=0 ; i < ds->mark_num ; i++) {
    d->region[i] = (double*) malloc(sizeof(double) * ds->w_size);
  }

  // bg_mean alloc
  d->bg_mean      = (double*) malloc(sizeof(double) * ds->mark_num);

  // bg_std alloc
  d->bg_std       = (double*) malloc(sizeof(double) * ds->mark_num);

  // mean_of_stds alloc
  d->mean_of_stds = (double*) malloc(sizeof(double) * ds->mark_num);
  
  // zscore
  d->zscore_bg = (double*) malloc(sizeof(double)*ds->w_size);
  d->val_flag = (int*) malloc(sizeof(int)*ds->w_size);
#else

  s->zscore_bg = (double*) malloc(sizeof(double)*ds->w_size);
  s->val_flag = (int*) malloc(sizeof(int)*ds->w_size);
#endif


  // weights_pos_motif alloc
  s->weights_motif_pos = (double**) malloc(sizeof(double*) * ds->mark_num);
  for(i=0 ; i < ds->mark_num ; i++) {
    s->weights_motif_pos[i] = (double*) malloc(sizeof(double) * ds->wandering_plus);
  }

  // weights_neg_motif alloc
  s->weights_motif_neg = (double**) malloc(sizeof(double*) * ds->mark_num);
  for(i=0 ; i < ds->mark_num ; i++) {
    s->weights_motif_neg[i] = (double*) malloc(sizeof(double) * ds->wandering_plus);
  }

  // weights_pos_bg alloc
  s->weights_bg_pos = (double**) malloc(sizeof(double*) * ds->mark_num);
	  for(i=0 ; i < ds->mark_num ; i++) {
    s->weights_bg_pos[i] = (double*) malloc(sizeof(double) * ds->wandering_plus);
  }

  // weights_neg_bg alloc
  s->weights_bg_neg = (double**) malloc(sizeof(double*) * ds->mark_num);
  for(i=0 ; i < ds->mark_num ; i++) {
    s->weights_bg_neg[i] = (double*) malloc(sizeof(double) * ds->wandering_plus);
  }

  // votes_pos alloc
  s->votes_pos = (double*) malloc(sizeof(double) * ds->wandering_plus);

  // votes_pos alloc
  s->votes_neg = (double*) malloc(sizeof(double) * ds->wandering_plus);
}

#ifndef CHROMASIG_MPI
void compute_pmf(Data *d,Scratch *s,_data_size* ds)
{
	int mod;
	int center_index;
	int motif_index;
	int num_none;
	double temp_sum_log_P;
	int index;
	double logr;
	double mean, std;
	double mean_std_sq;
	double bg_mean;
	double Z_score_P_j;
	double Z_score_B_j;
	double scale_ratio;
	double zscore_max;
	double zscore_max_sq;
	
	// for bg
	for( mod=0; mod<d->num_mod; mod++ )
	{
		bg_mean = d->bg_mean[mod];
		mean_std_sq = d->mean_of_stds[mod]*d->mean_of_stds[mod];
		zscore_max = ABS_Z_SCORE_MAX*d->mean_of_stds[mod];
		zscore_max_sq = -ABS_Z_SCORE_MAX_SQ*mean_std_sq;
		

		for( int wi=0; wi<ds->w_size; wi++ ) 
		{
			logr = d->region[mod][wi];
			if( logr == NONE_VAL )
			{
				d->val_flag[wi] = 0;
				d->zscore_bg[wi] = 0.;
			}
			else
			{
				d->val_flag[wi] = 1;
				if( fabs(Z_score_B_j=(logr-bg_mean)) > zscore_max )
				{
					d->zscore_bg[wi] = zscore_max_sq;
				}
				else
				{
					d->zscore_bg[wi] = -Z_score_B_j*Z_score_B_j;
				}
			}
		}
		/*if( flag )
		{
			for( int i=0; i<ds->w_size; i++ )
			{
				printf( "compute_pmf-zscore: wi=%d\tscore=%lf\tzscore=%lf\tflag=%d\tstd=%lf\n", i, d->zscore_bg[i], d->zscore_bg[i]/mean_std_sq, d->val_flag[i], mean_std_sq );
			}
		}*/

		num_none = 0;
		temp_sum_log_P = 0.;
		for( int mt_i=0; mt_i<=para.width; mt_i++ )
		{
			if( d->val_flag[mt_i] == 0 )
			{
				num_none++;
			}
			else
			{
				temp_sum_log_P += d->zscore_bg[mt_i];
			}
		}
		if( num_none>0 )
		{
			Z_score_B_j = temp_sum_log_P*(para.width+1)/(para.width+1-num_none);
		}
		else
		{
			Z_score_B_j = temp_sum_log_P;
		}
		s->weights_bg_pos[mod][0] = ds->exp_A+ds->p_ratio*exp(Z_score_B_j*0.5/mean_std_sq);
		/*if( flag )
		{
			printf( "compute_pmf-zscore-exp: wi=%d\tzscore_sum=%lf\tZ_score_B_j=%lf\tmean_std_sq=%lf\n", 0, Z_score_B_j*0.5/mean_std_sq, Z_score_B_j, mean_std_sq );
		}*/
		for( int c_i=1; c_i<=para.wandering_dist; c_i++ )
		{
			if( d->val_flag[c_i-1]==0 )
			{
				num_none--;
			}
			else
			{
				temp_sum_log_P -= d->zscore_bg[c_i-1];
			}
			if( d->val_flag[c_i+para.width]==0 )
			{
				num_none++;
			}
			else
			{
				temp_sum_log_P += d->zscore_bg[c_i+para.width];
			}
			if( num_none>0 )
			{
				Z_score_B_j = temp_sum_log_P*(para.width+1)/(para.width+1-num_none);
			}
			else
			{
				Z_score_B_j = temp_sum_log_P;
			}
			s->weights_bg_pos[mod][c_i] = ds->exp_A+ds->p_ratio*exp(Z_score_B_j*0.5/mean_std_sq);
			/*if( flag )
			{
				printf( "compute_pmf-zscore-exp: wi=%d\tzscore_sum=%lf\tZ_score_B_j=%lf\tmean_std_sq=%lf\n", c_i, Z_score_B_j*0.5/mean_std_sq, Z_score_B_j, mean_std_sq );
			}*/
		}
	}

	// for pol=-1
	for (mod = 0 ; mod < d->num_mod ; mod++) // cycle on mark num
	{
		for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) // cycle on offset
		{
			num_none = 0;
			temp_sum_log_P = 0;

			for( motif_index=0; motif_index<=d->width; motif_index++ ) //cycle on position on motif
			{
				index = center_index + d->width - motif_index;
				// motif_index: deal with polarization
				logr = d->region[mod][index];
				if (logr == NONE_VAL)
				{
					num_none++;

				}
				else
				{
					//find P, the probability under the motif model
					mean = d->motif_mean[mod][motif_index];
					std  = d->motif_std[mod][motif_index]; // + d->bg_std[mod] / 2;
					//Z_score_P_j = (logr - mean) / std;
					Z_score_P_j = logr-mean;
					Z_score_P_j *= Z_score_P_j;
					Z_score_P_j /= std;
					//if (fabs(Z_score_P_j) > ABS_Z_SCORE_MAX_SQ)
					if (Z_score_P_j > ABS_Z_SCORE_MAX_SQ)
					{
						temp_sum_log_P += EXP_MAX;
					}
					else
					{
						temp_sum_log_P -= Z_score_P_j;
					}
				} // end else for VAL
			}// end motif_index cycle

			if (num_none > 0)
			{
				scale_ratio = (d->width+1)/(d->width+1-num_none);
				temp_sum_log_P *= scale_ratio;
			}
			s->weights_motif_neg[mod][center_index] = exp(temp_sum_log_P*0.5);
		} // end cycle for offset
	} // end cycle for mark

	// pre-calculation of weights_bg_neg and weights_bg_pos

	// for pol=1
	for (mod = 0 ; mod < d->num_mod ; mod++) // cycle on mark num
	{
		//m_idx = (mod*ds->filter_num+id) * (para.wandering_dist+1);
		for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) // cycle on offset
		{
			num_none = 0;
			temp_sum_log_P = 0;

			for (motif_index = 0 ; motif_index <= d->width ; motif_index++) { //cycle on position on motif
				index = center_index+motif_index; // deal with polarization
				logr = d->region[mod][index];
				if (logr == NONE_VAL)
				{
					num_none++;
				}
				else
				{
					//find P, the probability under the motif model
					mean = d->motif_mean[mod][motif_index];
					std  = d->motif_std[mod][motif_index]; // + d->bg_std[mod] / 2;
					//Z_score_P_j = (logr - mean) / std;
					Z_score_P_j = logr-mean;
					Z_score_P_j *= Z_score_P_j;
					Z_score_P_j /= std;
					//if (fabs(Z_score_P_j) > ABS_Z_SCORE_MAX) {
					if (Z_score_P_j > ABS_Z_SCORE_MAX_SQ)
					{
						temp_sum_log_P += EXP_MAX;
					}
					else
					{
						temp_sum_log_P -= Z_score_P_j;
					}
				}
			}

			if (num_none > 0) {
				scale_ratio = (d->width+1)/(d->width+1-num_none);
				temp_sum_log_P *= scale_ratio;
			}
			s->weights_motif_pos[mod][center_index] = exp(temp_sum_log_P*0.5);
		} // end cycle for offset
	}
}
#else
void compute_pmf(double* d, float* mark_sig, Scratch* s,_data_size* ds)
{
        int mod;
        int center_index;
        int motif_index;
        int num_none;
        double temp_sum_log_P;
        int index;
        double logr;
        double mean, std;
        double mean_std_sq;
        double bg_mean;
        double Z_score_P_j;
        double Z_score_B_j;
        double scale_ratio;
        double zscore_max;
        double zscore_max_sq;

        double* d_motif_mean = d;
        double* d_motif_std = d+ds->m_wp;
        double* d_bg_mean = d+2*ds->m_wp;
        double* d_bg_std = d+2*ds->m_wp+ds->mark_num;
        double* idx_motif_mean;
        double* idx_motif_std;
        float* idx_mark;

        // for bg
        for( mod=0,idx_mark=mark_sig; mod<ds->mark_num; mod++ )
        {
                bg_mean = d_bg_mean[mod];
                mean_std_sq = d_bg_std[mod]*d_bg_std[mod];
                zscore_max = ABS_Z_SCORE_MAX*d_bg_std[mod];
                zscore_max_sq = -ABS_Z_SCORE_MAX_SQ*mean_std_sq;

                for( int wi=0; wi<ds->w_size; wi++,idx_mark++ )
                {
                        logr = *idx_mark;
                        if( logr == NONE_VAL )
                        {
                                s->val_flag[wi] = 0;
                                s->zscore_bg[wi] = 0.;
                        }
                        else
                        {
                                s->val_flag[wi] = 1;
                                if( fabs(Z_score_B_j=(logr-bg_mean)) > zscore_max )
                                {
                                        s->zscore_bg[wi] = zscore_max_sq;
                                }
                                else
                                {
                                        s->zscore_bg[wi] = -Z_score_B_j*Z_score_B_j;
                                }
                        }
                }

                num_none = 0;
                temp_sum_log_P = 0.;
                for( int mt_i=0; mt_i<=para.width; mt_i++ )
                {
                        if( s->val_flag[mt_i] == 0 )
                        {
                                num_none++;
                        }
                        else
                        {
                                temp_sum_log_P += s->zscore_bg[mt_i];
                        }
                }
                if( num_none>0 )
                {
                        Z_score_B_j = temp_sum_log_P*ds->width_plus/(ds->width_plus-num_none);
                }
                else
                {
                        Z_score_B_j = temp_sum_log_P;
                }
                s->weights_bg_pos[mod][0] = ds->exp_A+ds->p_ratio*exp(Z_score_B_j*0.5/mean_std_sq);
                for( int c_i=1; c_i<=para.wandering_dist; c_i++ )
                {
                        if( s->val_flag[c_i-1]==0 )
                        {
                                num_none--;
                        }
                        else
                        {
                                temp_sum_log_P -= s->zscore_bg[c_i-1];
                        }
                        if( s->val_flag[c_i+para.width]==0 )
                        {
                                num_none++;
                        }
                        else
                        {
                                temp_sum_log_P += s->zscore_bg[c_i+para.width];
                        }
                        if( num_none>0 )
                        {
                                Z_score_B_j = temp_sum_log_P*(para.width+1)/(para.width+1-num_none);
                        }
                        else
                        {
                                Z_score_B_j = temp_sum_log_P;
                        }
                        s->weights_bg_pos[mod][c_i] = ds->exp_A+ds->p_ratio*exp(Z_score_B_j*0.5/mean_std_sq);
                }
        }

        // for pol=-1
        idx_mark = mark_sig;
        idx_motif_mean = d_motif_mean;
        idx_motif_std = d_motif_std;
        for (mod = 0 ; mod < ds->mark_num ; mod++,idx_mark+=ds->w_size,idx_motif_mean+=ds->width_plus,idx_motif_std+=ds->width_plus) // cycle on mark num
        {
                for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) // cycle on offset
                {
                        num_none = 0;
                        temp_sum_log_P = 0;

                        for( motif_index=0; motif_index<=para.width; motif_index++ ) //cycle on position on motif
                        {
                                index = center_index + para.width - motif_index;
                                // motif_index: deal with polarization
                                //logr = d->region[mod][index];
                                logr = *(idx_mark+index);
                                if (logr == NONE_VAL)
                                {
                                        num_none++;

                                }
                                else
                                {
                                        //find P, the probability under the motif model
                                        mean = *(idx_motif_mean+motif_index);
                                        std  = *(idx_motif_std+motif_index); // + d->bg_std[mod] / 2;
                                        //Z_score_P_j = (logr - mean) / std;
                                        Z_score_P_j = logr-mean;
                                        Z_score_P_j *= Z_score_P_j;
                                        Z_score_P_j /= std;
                                        //if (fabs(Z_score_P_j) > ABS_Z_SCORE_MAX_SQ)
                                        if (Z_score_P_j > ABS_Z_SCORE_MAX_SQ)
                                        {
                                                temp_sum_log_P += EXP_MAX;
                                        }
                                        else
                                        {
                                                temp_sum_log_P -= Z_score_P_j;
                                        }
                                } // end else for VAL
                        }// end motif_index cycle

                        if (num_none > 0)
                        {
                                scale_ratio = ds->width_plus/(ds->width_plus-num_none);
                                temp_sum_log_P *= scale_ratio;
                        }
                        s->weights_motif_neg[mod][center_index] = exp(temp_sum_log_P*0.5);
                } // end cycle for offset
        } // end cycle for mark

        // pre-calculation of weights_bg_neg and weights_bg_pos

        // for pol=1
        idx_mark = mark_sig;
        idx_motif_mean = d_motif_mean;
        idx_motif_std = d_motif_std;
        for (mod = 0 ; mod < ds->mark_num ; mod++,idx_mark+=ds->w_size,idx_motif_mean+=ds->width_plus,idx_motif_std+=ds->width_plus) // cycle on mark num
        {
                //m_idx = (mod*ds->filter_num+id) * (para.wandering_dist+1);
                for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) // cycle on offset
                {
                        num_none = 0;
                        temp_sum_log_P = 0;

                        for (motif_index = 0 ; motif_index <= para.width ; motif_index++) { //cycle on position on motif
                                index = center_index+motif_index; // deal with polarization
                                //logr = d->region[mod][index];
                                logr = *(idx_mark+index);
                                if (logr == NONE_VAL)
                                {
                                        num_none++;
                                }
                                else
                                {
                                        //find P, the probability under the motif model
                                        mean = *(idx_motif_mean+motif_index);
                                        std  = *(idx_motif_std+motif_index); // + d->bg_std[mod] / 2;
                                        //Z_score_P_j = (logr - mean) / std;
                                        Z_score_P_j = logr-mean;
                                        Z_score_P_j *= Z_score_P_j;
                                        Z_score_P_j /= std;
                                        //if (fabs(Z_score_P_j) > ABS_Z_SCORE_MAX) {
                                        if (Z_score_P_j > ABS_Z_SCORE_MAX_SQ)
                                        {
                                                temp_sum_log_P += EXP_MAX;
                                        }
                                        else
                                        {
                                                temp_sum_log_P -= Z_score_P_j;
                                        }
                                }
                        }

                        if (num_none > 0) {
                                scale_ratio = ds->width_plus/(ds->width_plus-num_none);
                                temp_sum_log_P *= scale_ratio;
                        }
                        s->weights_motif_pos[mod][center_index] = exp(temp_sum_log_P*0.5);
                } // end cycle for offset
        }
}
#endif

////////////////////////////////////////////////////////////////////////////////
//
// get_votes
//
////////////////////////////////////////////////////////////////////////////////
int get_votes(Scratch *s,_data_size* ds)
{
        int mod;
        //int pol;
        int center_index;

        // if any single mark's best alignment is not as good as the background,
        // then reject
        for (mod = 0 ; mod < ds->mark_num ; mod++) {
                int reject = 1;
                for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) {

                        // pol == 1
                        if (s->weights_motif_pos[mod][center_index] > s->weights_bg_pos[mod][center_index]) {
                                reject = 0;
                                break;
                        }

                        // pol == -1
                        if (s->weights_motif_neg[mod][center_index] > s->weights_bg_pos[mod][center_index]) {
                                reject = 0;
                                break;
                        }
                }
                if( reject )
                {
                        return 1;
                }
        }

        double norm;
        for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) {
                s->votes_pos[center_index] = 0;
                s->votes_neg[center_index] = 0;

                for (mod = 0 ; mod <ds->mark_num  ; mod++) {

                        // pol = 1
                        if (s->weights_motif_pos[mod][center_index] != 0) {
                                norm = log( s->weights_motif_pos[mod][center_index] /
                                                s->weights_bg_pos[mod][center_index] );
                        } else {
                                norm = MIN_LOG_RATIO;
                        }
                        s->votes_pos[center_index] += norm;

                        // pol = -1
                        if (s->weights_motif_neg[mod][center_index] != 0) {
                                norm = log( s->weights_motif_neg[mod][center_index] /
                                                s->weights_bg_pos[mod][center_index] );
                        } else {
                                norm = MIN_LOG_RATIO;
                        }
                        s->votes_neg[center_index] += norm;
                }
        }

        return 0;
}
////////////////////////////////////////////////////////////////////////////////
//
// tally_votes
//
////////////////////////////////////////////////////////////////////////////////
//void tally_votes(int has_exclude_max) {
void tally_votes(Scratch *s,_data_size* ds) {
	int center_index;

	/*
	// one of the modifications had the most probable state being exclude
	if (has_exclude_max) {
	s->return_pol = 1;
	s->return_loc = NONE_VAL;
	return;
	}
	*/
	int max_index = NONE_VAL;
	int max_pol   = NONE_VAL;
	double max_val   = NONE_VAL;

	// pol = -1
	for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++)
	{
		if ((max_val == NONE_VAL) || (max_val < s->votes_neg[center_index]))
		{
			max_index = center_index;
			max_pol = -1;
			max_val = s->votes_neg[center_index];
		}
	}

	// pol = 1
	for (center_index = 0 ; center_index <= para.wandering_dist ; center_index++) {
    if ((max_val == NONE_VAL) ||
	(max_val < s->votes_pos[center_index])) {
      max_index = center_index;
      max_pol = 1;
      max_val = s->votes_pos[center_index];
    }
  }

  s->return_pol = max_pol;
  s->return_loc = max_index - ds->wandering_half;
}
