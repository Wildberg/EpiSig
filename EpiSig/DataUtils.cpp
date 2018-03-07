#include <cstdio>
#include <iostream>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "utils.h"
#include "FileUtils.h"
#include "DataUtils.h"
#include <string.h>
#include <stdlib.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
extern struct _mark *mark;
extern struct _encode *encode;
extern struct _filter *filter;

extern float* mark_buf[MAX_MARK];
extern float* mark_buf_ref[MAX_MARK];

double F_ABS_Z_SCORE_MAX = 10;

// sort on chromosome id
int cmp_filter(const void *a, const void *b)
{
  int id1 = ((_filter*)a)->en;
  int id2 = ((_filter*)b)->en;
  if(id1 > id2) return 1;
  else if(id1 < id2) return -1;
  else return 0;
}

// sort on the value of norm
int cmp_norm(const void *a, const void *b)
{
  double aa = ((_norm*)a)->norm;
  double bb = ((_norm*)b)->norm;
  if(aa > bb) return -1;
  else if(aa < bb) return 1;
  else return 0;
}
int cmp_double(const void *a, const void *b)
{
  double aa = *((double*)a);
  double bb = *((double*)b);
  if(aa > bb) return -1;
  else if(aa < bb) return 1;
  else return 0;
}

int get_marks(string fname,_mark **mark)
{
  FileUtils file;
  FILE* fp;
  
  int len = file.CountLines(fname); // CountLines, FileUtils.cpp
  if(len < 0) return -1;

  *mark = new  _mark[len+1]; // mark-DataUtil.cpp global variable for marker information
  if(!mark) 
  {
    cout<<"Allocate memory failed!"<<endl;
    return -1;
  }

  fp = fopen(fname.c_str(),"r");
  if(!fp)
  {
    cout<<"Can not open mark file!"<<endl;
    return -1;
  }
  char s[512];
  string fields[3];
  int cnt = 0;
  while(!feof(fp))
  {
    fgets(s,512,fp);
    int num = 2;
    spliter(s,'\t',fields,num); // string split, utils.cpp
    (*mark)[cnt].id = atoi(fields[0].c_str()); // id, integer number
    (*mark)[cnt].name = fields[1]; // file name
    (*mark)[cnt].mean = 0;
    (*mark)[cnt].std = 0;
    cnt++;
  }
  fclose(fp);
  return len;
}
int load_marks(string fname,_mark **mark)
{
	FileUtils file;
	FILE* fp;
  
	int len = file.CountLines(fname);
	if(len < 0) return -1;

	//mark = new  _mark[len+1];
	*mark = new _mark[len];
	if(!mark) 
	{
		cout<<"Allocate memory failed!"<<endl;
		return -1;
	}

	fp = fopen(fname.c_str(),"r");
	if(!fp)
	{
		cout<<"Can not open mark file!"<<endl;
		return -1;
	}
	char s[512];
	string fields[5];
	int cnt = 0;
	while(!feof(fp))
	{
		fgets(s,512,fp);
		int num = 4;
		spliter(s,'\t',fields,num);
		(*mark)[cnt].id = atoi(fields[0].c_str());
		(*mark)[cnt].name = fields[1];
		(*mark)[cnt].mean = atof(fields[2].c_str());
		(*mark)[cnt].std = atof(fields[3].c_str());
		(*mark)[cnt].z_max = F_ABS_Z_SCORE_MAX*(*mark)[cnt].std;
		(*mark)[cnt].exp_max = -(*mark)[cnt].z_max * (*mark)[cnt].z_max;
		cnt++;
	}
	fclose(fp);

	return len;
}

int load_filter(string fname,_filter **filter, GS* gs)
{
  FileUtils file;
  FILE* fp;
  
  int len = file.CountLines(fname);
  if(len < 0) return -1;

//  filter = new  _filter[len+1];
  *filter = new  _filter[len+1];
  if(!filter) 
  {
    cout<<"Allocate memory failed!"<<endl;
    return -1;
  }

  fp = fopen(fname.c_str(),"r");
  if(!fp)
  {
    cout<<"Can not open mark file!"<<endl;
    return -1;
  }
  char s[512];
  string fields[5];
  int cnt = 0;
  while(!feof(fp))
  {
    fgets(s,512,fp);
    int num = 2;
    spliter(s,'\t',fields,num);
    //(*filter)[cnt].en = chr2id(fields[0]);
    (*filter)[cnt].en = gs->chr_id_map[fields[0]]+1;
    (*filter)[cnt].pos = atoi(fields[1].c_str());
    cnt++;
  }
  fclose(fp);
  return len;
}
int get_encode(string fname,_encode **encode,GS* gs)
{
  FileUtils file;
  FILE* fp;
  
  int len = file.CountLines(fname);
  if(len < 0) return -1;

  *encode = new  _encode[len+1]; // encode global variable, DataUtil.cpp
  if(!encode) 
  {
    cout<<"Allocate memory failed!"<<endl;
    return -1;
  }

  fp = fopen(fname.c_str(),"r");
  if(!fp)
  {
    cout<<"Can not open mark file!"<<endl;
    return -1;
  }
  char s[512];
  string fields[5];
  int cnt = 0;
  while(!feof(fp))
  {
    fgets(s,512,fp);
    int num = 4;
    spliter(s,'\t',fields,num);
    //(*encode)[cnt].chrom = chr2id(fields[0]); //chr2id
    (*encode)[cnt].chrom = gs->chr_id_map[fields[0]]+1; //chr2id
    (*encode)[cnt].min = atoi(fields[1].c_str());
    (*encode)[cnt].max = atoi(fields[2].c_str());
    (*encode)[cnt].name = fields[3];
    cnt++;
  }
  fclose(fp);

  return len;
}


int alloc_mem(float** mark_buf, int mark_num,int chr_idx,int old_chr,_mark *mark, GS* gs)
// mark_buf[MAX_MARK] see filter.cpp
// mark_num, number of mark (number of files)
// chr_idx, 1-24 index
// old_chr, the last chr processed in encode list, see filter.cpp
{
  char currfile[512];
  int buf_size = gs->g_size[chr_idx-1]/100; // number of bins used throughout chromosome
  //int buf_size = hg18[chr_idx-1]/100; // number of bins used throughout chromosome

  // if on same chromosome, do not assign
  if(old_chr == chr_idx && old_chr >= 0) return 1;

  for(int j = 0; j < mark_num; j++) // loop over all mark number
  {
     if(old_chr >=0 ) free(mark_buf[j]); // if region assigned previously
     mark_buf[j] = (float*)malloc(sizeof(float)*(buf_size+1)); // assign memory for each mark on every chromosome
     if(mark_buf[j] == NULL)
     {
       cout<<"Allocate memory failed when opening files!"<<endl;
     }

     // get file name for mark[j] and chromosome chr_idx
     sprintf(currfile,"out_data/%s.%d.bin",mark[j].name.c_str(),chr_idx);
     FILE* fp = fopen(currfile,"rb");
     if(!fp)
     {
        cout<<"Cannot open "<<currfile<<endl;
        exit(-1);
     }
     // read all data from pre-EpiSig binary file
     int sizeread = fread(mark_buf[j], sizeof(float), buf_size, fp);
     if(sizeread != buf_size)
     {
        cout<<"Fatal error occurs when reading data files!"<<endl;
        exit(-1);
     }
     fclose(fp);
     fp = NULL;
  }
  return 1;
}
