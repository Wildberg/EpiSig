#ifndef __UPDATE_GIBBS_STATE_H
#define __UPDATE_GIBBS_STATE_H

#include "global_data.h"
#include "thread_random.h"

//using namespace std;

const int NUM_CYCLES = 5;
const int MAX_SEED_SIZE = 20;
const int MAX_THREAD = 20;
const int BACKGROUND_BUFFER = 20000;
const int BUFFER_LIMIT = 6000;

struct _para para;

#ifndef CHROMASIG_MPI
struct _arg_thread
{
	_data_size* ds;
	THREAD_RAND* thread_rand;
	Data *dd;
	Scratch *ss;
	_filter* filter;
	_filter* first_sets;
	_gibbs_state* gibbs;
	double* bg_prob;
	float* mark_sig;
	int* undo_list;
	int* tag;
	int size_sites;
	int func_id;
};

struct _func_arg{
        int mark_num;
        int width;
        _gibbs_state* gibbs;
        int* tag;
        int not_list;
        int* not_done_ids;
        float* undo_buffer;
        int page_num;
        int func_id;
	Data *dd;
	Scratch *ss;
	_filter* filter;
};

double best_score[MAX_THREAD];
int thread_seed_size[MAX_THREAD];
_gibbs_state best_gibbs[MAX_THREAD][MAX_SEED_SIZE+1];
#endif

#endif
