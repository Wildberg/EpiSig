#ifndef __GENOME_SIZE_H
#define __GENOME_SIZE_H

#include<vector>
#include<map>
#include<string>

using namespace std;

struct GS
{
	map<string,int> chr_id_map;
	vector<string> id_chr_map; 
	vector<long> g_size;
	int g_cnt;
};

void GS_init( GS& gs, string& name );
void GS_init( GS& gs, const char* name );
void GS_init_preset( GS& gs, const char* name );
#endif
