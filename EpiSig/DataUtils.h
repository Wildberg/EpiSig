#ifndef __DATA_UTILS_H
#define __DATA_UTILS_H

#include <string>
#include "genome_size.h"
using namespace std;
#pragma once
#define HG_CHR_NUM 24
// marker information structure
struct _mark{
	int id; // id, in mark.idx file
	string name; // file name/marker name in mark.idx file
	double mean;
	double std;
	double z_max;
	double exp_max;
};
struct _gibbs_state{
// state for each region/motif
	int id;
	int en; // chr id
	int loc; // location
	int pol; // direction +/-
};
struct _logr_sums{
	double sum;
	double sum2;
	int count;
};
struct _encode{
	int chrom;
	int min;
	int max;
	string name;
};
//struct _background{
//	int id;
//	string mark;
//	double mean;
//	double std;
//};
struct _motif_dist{
	double mean;
	double std;
};
struct _sig{
	int pos;
	double sum;
	int marks;
};

struct _filter{
	int en; // chr id
	int pos; // position, chromosome loci
};
struct _distance{
	int id; // filtered region id
	double distance; // distance to the seed
};
struct _norm{
        int id; // id to filter_num
        double norm; // norm of the region
};
struct _mark_sig{
        int pos;
        double z;
};

extern int hg18[];

//extern struct _mark *mark;
extern struct _encode *encode;
//extern struct _filter *filter;

#define MAX_MARK 1024
extern float* mark_buf[MAX_MARK];
extern float* mark_buf_ref[MAX_MARK];

extern int cmp_filter(const void *a, const void *b);
extern int cmp_norm(const void *a, const void *b);
extern int cmp_double(const void *a, const void *b);
extern int get_marks(string fname,_mark **mark);
extern int load_marks(string fname,_mark **mark);
extern int load_filter(string fname,_filter **filter,GS* gs);
extern int get_encode(string fname,_encode **encode,GS* gs);
extern int alloc_mem(float** mark_buf, int mark_num,int chr_idx,int old_chr,_mark *mark,GS* gs);
#endif
