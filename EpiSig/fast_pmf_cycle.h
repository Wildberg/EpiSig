#ifndef fast_pmf_cycle_header
#define fast_pmf_cycle_header

#include "global_data.h"

////////////////////////////////////////////////////////////////////////////////
//
// fast_pmf_cycle.h
//
////////////////////////////////////////////////////////////////////////////////

////////////////////
// structs
////////////////////
#ifndef CHROMASIG_MPI 
typedef struct {
  int num_mod;
  int width;
  int wander_dist;

  double **motif_mean;
  double **motif_std;
  double **region;
  double *bg_mean;
  double *bg_std;
  double *mean_of_stds;
  double *zscore_bg;
  int *val_flag;
} Data;
#endif

typedef struct {
  int return_pol;
  int return_loc;

#ifdef CHROMASIG_MPI
  double *zscore_bg;
  int *val_flag;
#endif

  double **weights_motif_pos;
  double **weights_motif_neg;
  double **weights_bg_pos;
  double **weights_bg_neg;

  double *votes_pos;
  double *votes_neg;
} Scratch;

extern int NONE_VAL;
extern double ABS_Z_SCORE_MAX;
extern double MIN_LOG_RATIO;

////////////////////
// accessors
////////////////////
#ifndef CHROMASIG_MPI
double get_motif_mean(Data *d,Scratch *s,int mod, int pos);
void set_motif_mean(Data *d,Scratch *s,int mod, int pos, double val);

double get_motif_std(Data *d,Scratch *s,int mod, int pos);
void set_motif_std(Data *d,Scratch *s,int mod, int pos, double val);

double get_bg_mean(Data *d,Scratch *s,int mod);
void set_bg_mean(Data *d,Scratch *s,int mod, double val);

double get_bg_std(Data *d,Scratch *s,int mod);
void set_bg_std(Data *d,Scratch *s,int mod, double val);

double get_region(Data *d,Scratch *s,int mod, int pos);
void set_region(Data *d,Scratch *s,int mod, int pos, double val);
#endif

int get_return_pol(Scratch *s);
int get_return_loc(Scratch *s);
int get_votes(Scratch *s,_data_size* ds);
void tally_votes(Scratch *s,_data_size* ds);

////////////////////
// init
////////////////////
#ifndef CHROMASIG_MPI
void init(Data *d,Scratch *s, _data_size* ds );

void set_motif_mean_of_std(Data *d,Scratch *s);
void compute_pmf(Data *d,Scratch *s, _data_size* ds);
#else
void init(Scratch *s, _data_size* ds );

void compute_pmf(double *d, float* mark_sig, Scratch *s, _data_size* ds);
#endif

#endif
