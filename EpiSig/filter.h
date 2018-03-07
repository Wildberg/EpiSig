#ifndef __READ_TEST_DATA_H
#define __READ_TEST_DATA_H

#include <string>
using namespace std;
#pragma once

void get_local_maxima(_mark_sig* mark_sig,int num,int half_size,int* list, int& list_num);
void print_non_overlapping(_sig* sig,int sig_num,int chr,int half_size,GS* gs);
void get_sig_locs(int stat_half_size, int maxima_half_size, int overlap_half_size, double pval_cutoff, GS* gs);
void get_normal_params(int stat_half_window_size, string mark_file, string encode_file, int flag_global_background, GS* gs);
void get_background_dist_global(int bg_buf_size, string mark_file, string filter_file, GS* gs);

#endif
