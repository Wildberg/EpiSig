#ifndef __PROCESSDATA_H
#define __PROCESSDATA_H

#include <string>
using namespace std;
#pragma once

const int MAX_THREAD = 96;

//int hg18[] = {247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737,
//		 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964,
//		 46944323, 49691432, 154913754, 57772954};

struct _para{
        int max_thread; // multi-threading threads
	int format;
	string data_dir;
	int skip;
	string genome_name;
};
struct _para para;

void ExtractBed(string fname,const char* mark,char delimiter);
#endif
