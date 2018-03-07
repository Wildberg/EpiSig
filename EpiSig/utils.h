#ifndef __UTILS_H
#define __UTILS_H

//#include <string>
using namespace std;
#pragma once

#define ISSPACE(x) ((x)==' '||(x)=='\r'||(x)=='\n'||(x)=='\f'||(x)=='\b'||(x)=='\t')

struct _meta
{
	char chrom[16];
	int start;
	int step;
	int span;
	int steptype;
};

extern int cmp_int(const void *a, const void *b);
int FileLength(string fname);
char* MakeLower(char* str);
char* TrimLeft(char* str);
char* TrimRight(char* str);
char* Trim(char* str);
string basename_(char* path);
void spliter(const char* line,char delimiter,string* rt,int& cnt);
#endif
