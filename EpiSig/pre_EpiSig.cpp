#include <iostream>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <math.h>

#include "utils.h"
#include "FileUtils.h"
#include "ConfUtils.h"
#include "pre_EpiSig.h"
#include "genome_size.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// for multi-threading purpose
struct _arg{
  char directory[512];
        int func_id;
	GS* gs;
};

void print_usage()
{
        string buffer;
        buffer = "Usage: pre_EpiSig [--data-dir] [--format] ...\n";
        buffer += "Options:\n";
        buffer += "-h, --help      display this help and exit\n";
        buffer += "-d, --data-dir      directory of input data marks\n";
        buffer += "-f, --format      [BED | WIG | BIN] format\n";
        buffer += "-s, --skip-norm      skip normalization(normalized, or subtract control data already)\n";
        buffer += "-n, --threads      mutiple threads number, default: 6\n";
        buffer += "-i, --genome      genome used in pre_EpiSig. Two options: pre-defined genome names: hg18, hg19, hg38, mm9, or mm10, or user-defined tab-delimit file\n";
        printf("%s\n",buffer.c_str());
}
void print_para(_para para)
{
        printf("pre_EpiSig use following parameters.\n");
        printf("Data directory: %s\n", para.data_dir.c_str());
  if(para.format == 0) printf("format: BED\n");
  if(para.format == 1) printf("format: WIG\n");
  if(para.format == 2) printf("format: BIN\n");
  if(para.skip == 0) printf("Normalization: sigmoid function\n");
  if(para.skip == 1) printf("Normalization: skipped\n");
        printf("threads: %d\n", para.max_thread);
        printf("genome: %s\n", para.genome_name.c_str());
}

void GetMeta(char*line,_meta& meta)
{
  Trim(line);
  MakeLower(line);

  char* p = strstr(line,"fixedstep");
  if(p) meta.steptype = 1;
  p = strstr(line,"variablestep");
  if(p) meta.steptype = 0;

  char value[128];
  ConfUtils conf;
  int ret;

  ret = conf.GetPara(line,"chrom",value);
  if(ret) strcpy(meta.chrom,value);
  ret = conf.GetPara(line,"start",value);
  if(ret) meta.start = atoi(value);
  ret = conf.GetPara(line,"step",value);
  if(ret) meta.step = atoi(value);
  ret = conf.GetPara(line,"span",value);
  if(ret) meta.span = atoi(value);
  if(meta.span <= 1) meta.span = meta.step;
}

void normalize(int** bucket, const char* mark, GS* gs) // mark, mark name/file name of the bed file
{
  //normalize data
  int read_count[32768*10]; // histogram by count
  double cdf[32768*10];
  int num_bins = 0; // total number of bins of all chromosomes
  int total_reads = 0; // total number of reads of all chromosomes
  cdf[0] = cdf[1] = 0;
  memset(read_count,0,32768*sizeof(int)*10); // initialization of the array
  for(int i=0; i<gs->g_cnt; i++) // loop over all chromosomes
  {
    //int size = hg18[i] / 100 + 1;
    int size = gs->g_size[i] / 100 + 1;
    for(int j=0;j<size;j++)
    {
      int count = bucket[i][j]; // bucket[i][j], i-chromosome, j-bin number, bucket[i][j]-signal of the bin
      if(count >= 2) // only bins with read counts more than 1 is used for the statistics
      {
        num_bins++;
        total_reads += count;
        read_count[count] ++;
      }
    }
  }

  double mean = total_reads*1.0/num_bins; // mean of signal among all bins
  double std = 0;
  double median = 0;
  int flag = 1;
  // calculate mean, std, and median
  for(int i=2; i<32768*10; i++) // loop over all bins more than 2
  {
    std += read_count[i]*(i - mean)*(i - mean);
    cdf[i] = (read_count[i]*i*1.0)/total_reads + cdf[i-1];
    if(flag && cdf[i] >= 0.5)
    {
      median = i;
      flag = 0;
    }
  }
  std = std/(num_bins-1);
  std = sqrt(std);
  // finish calculation

  double prob;
  for(int i=0; i<gs->g_cnt; i++)
  {
    //int size = hg18[i] / 100 + 1;
    int size = gs->g_size[i] / 100 + 1;

    float* out = new float[size];
    if(!out)
    {
      printf("Allocate output buffer failed!\n");
      break;
    }
    for(int j=0;j<size;j++)
    {
      int count = bucket[i][j];
      if(count == 0) { prob = 0; }
      else{
        prob = 1.0/(1+exp(-1*(count-median)/std)); // sigmoid function normalization on all signal of all chromosomes
      }
      char bbb[32];
      sprintf(bbb,"%.4f",prob);
      if(para.skip == 1) *(out+j) = count;
      else *(out+j) = (float)atof(bbb);
    }

    //write whole chromosome signal to binary files
    char fn1[1024];
    sprintf(fn1,"out_data/%s.%d.bin",mark,i+1);
  
    cout<<"Write to "<<fn1<<"......";
    FILE* ff = fopen(fn1, "wb+");
    if(!ff)
    {
      cout<<"Failed to open file for write!"<<endl;
      delete[] out;
      break;
      
    }
    fwrite(out, sizeof(float), size, ff);
    fclose(ff);
    cout<<"OK!"<<endl;
    delete[] out;
  }
}
void ExtractWig(string fname,const char* mark,char delimiter,GS* gs)
{
  //int* bucket[24];
  int** bucket;
        char line[1024];
        int ret;
  _meta meta;

  //allocate memory space
	bucket = new int*[gs->g_cnt];
  for(int i=0; i<gs->g_cnt; i++)
  {
    //int size = hg18[i] / 100 + 1;
    int size = gs->g_size[i] / 100 + 1;
    bucket[i] = new int[size]; 
    memset(bucket[i],0,size*sizeof(int));
    if(bucket[i] == NULL) return;
  }

        FileUtils file;
        ret = file.OpenFileForRead(fname);
        if(!ret)
        {
                cout<<"Can't open file "<<fname.c_str()<<" for read!"<<endl;
                exit(-1);
        }
        ret = file.AllocBuffer();
        if(!ret)
        {
                cout<<"can't allocate memory!"<<endl;
                exit(-1);
        }
        ret = file.FirstRead();
        if(!ret)
        {
                file.Release();
                cout<<"Read file failed!"<<endl;
                exit(-1);
        }

  int count = 0;
        int total = 0;
        string chrom;

        strcpy(meta.chrom,"");
        meta.start = 1;
        meta.step = 1;
        meta.span = 1;
        meta.steptype = 1;
  string fields[5];
        while(1)
        {
                int not_end = file.GetLine(line,1024);
    if(!not_end)
    {
      Trim(line);
      if(line[0] == '\0') break;
    }
    
    if(strstr(line,"="))
                {
                        GetMeta(line,meta);
      count = 0;
                        continue;
                }

    int num = 2;
    int start,end;
    int value = 0;

    spliter(line,delimiter,fields,num);
    //int index = chr2id(meta.chrom);
    int index = gs->chr_id_map[meta.chrom]+1;

    //get start information  
    if(meta.steptype == 1) start = meta.start + count*meta.step;
    else start = atoi(fields[0].c_str());

    //get end information
    end = start + meta.span - 1 ;
    int pos = (start + end)/200;
  
    //get value information
    if(meta.steptype == 1) value = atoi(fields[0].c_str());
    if(meta.steptype == 0) value = atoi(fields[1].c_str());

    total++; //total 

    bucket[index-1][pos] += value;
  
                if(total%2000000 == 0)
                {
      printf("%s: %d records were scanned!\n", mark, total);
                }
  
                if(!not_end) break;
        }

        file.Release();
  
  normalize(bucket, mark, gs);

  for(int i=0; i<gs->g_cnt; i++) delete[] bucket[i];
	delete[] bucket;
}
void ExtractBed(string fname,const char* mark,char delimiter,GS* gs)
{
  //int* bucket[24]; // 24-number of chromosomes
	int** bucket;
        char line[1024];
        int ret;

  //allocate memory space
	bucket = new int*[gs->g_cnt];
  for(int i=0; i<gs->g_cnt; i++) // loop over all chromosomes
  {
    //int size = hg18[i] / 100 + 1; // 100bp bin
    int size = gs->g_size[i] / 100 + 1; // 100bp bin
    bucket[i] = new int[size]; 
    memset(bucket[i],0,size*sizeof(int));
    if(bucket[i] == NULL) return;
  }

        FileUtils file;
        ret = file.OpenFileForRead(fname); // open file for read in FileUtils
        if(!ret)
        {
                cout<<"Can't open file "<<fname.c_str()<<" for read!"<<endl;
                exit(-1);
        }
        ret = file.AllocBuffer(); // allocate buffer array in FileUtils
        if(!ret)
        {
                cout<<"can't allocate memory!"<<endl;
                exit(-1);
        }
        ret = file.FirstRead(); // read in contents of the file, FileUtils
        if(!ret)
        {
                file.Release(); // de-construct FileUtils
                cout<<"Read file failed!"<<endl;
                exit(-1);
        }

        int total = 0;
	string fields[5];
	while(1)
	{
		int not_end = file.GetLine(line,1024);
		if(!not_end)
		{
			Trim(line);
			if(line[0] == '\0') break;
		}

		int num = 3;
		spliter(line,delimiter,fields,num);

		// get the chromsome index
		//int index = chr2id(fields[0]);
		int index = gs->chr_id_map[fields[0]]+1;
		int start = atoi(fields[1].c_str());  
		int end = atoi(fields[2].c_str());
		int pos = (start + end)/200; // bin number

		// calculate the signal for each bin
		//if(pos <= hg18[index-1]/100){
		if(pos <= gs->g_size[index-1]/100){
			bucket[index-1][pos]++; 
		}

                total++; //total count
  
                if(total%2000000 == 0)
                {
			printf("%s: %d records were scanned!\n", mark, total);
                }
  
                if(!not_end) break;
        }

        file.Release();

  normalize(bucket, mark, gs);

  for(int i=0; i<gs->g_cnt; i++) delete[] bucket[i];
	delete[] bucket;
}
void ExtractBin(string fname,const char* mark,char delimiter,GS* gs)
{
  //int* bucket[24];
  int** bucket;
        char line[1024];
        int ret;

  //allocate memory space
	bucket = new int*[gs->g_cnt];
  for(int i=0; i<gs->g_cnt; i++)
  {
    //int size = hg18[i] / 100 + 1;
    int size = gs->g_size[i] / 100 + 1;
    bucket[i] = new int[size]; 
    memset(bucket[i],0,size*sizeof(int));
    if(bucket[i] == NULL) return;
  }

        FileUtils file;
        ret = file.OpenFileForRead(fname);
        if(!ret)
        {
                cout<<"Can't open file "<<fname.c_str()<<" for read!"<<endl;
                exit(-1);
        }
        ret = file.AllocBuffer();
        if(!ret)
        {
                cout<<"can't allocate memory!"<<endl;
                exit(-1);
        }
        ret = file.FirstRead();
        if(!ret)
        {
                file.Release();
                cout<<"Read file failed!"<<endl;
                exit(-1);
        }

        int total = 0;
  string fields[5];
        while(1)
        {
                int not_end = file.GetLine(line,1024);
    if(!not_end)
    {
      Trim(line);
      if(line[0] == '\0') break;
    }

    int num = 3;
    spliter(line,delimiter,fields,num);

    //get the chromsome index
    //int index = chr2id(fields[0]);
    int index = gs->chr_id_map[fields[0]]+1;
    int pos = atoi(fields[1].c_str());  
    int val = atoi(fields[2].c_str());
    if(pos % 100 == 0) pos = pos / 100;

    //if(pos <= hg18[index-1]/100)
    if(pos <= gs->g_size[index-1]/100)
    bucket[index-1][pos] += val;

                total++; //total count
  
                if(total%2000000 == 0)
                {
      printf("%s: %d records were scanned!\n", mark, total);
                }
  
                if(!not_end) break;
        }

        file.Release();

  normalize(bucket, mark, gs);

  for(int i=0; i<gs->g_cnt; i++) delete[] bucket[i];
	delete[] bucket;
}
int make_mark_idx(const char* directory)
{
  FileUtils file;
  DIR *dir;
        char currfile[512];
        struct dirent *s_dir;
        struct stat file_stat;
  char buf[1024];

        stat(directory,&file_stat);
        if(!S_ISDIR(file_stat.st_mode))
        {
    strcpy(currfile,directory);
    sprintf(buf,"1\t%s",basename_(currfile).c_str());
    file.WriteLine("mark.idx",buf);
    cout<<"Building mark.idx.............[OK!]"<<endl;
        }
  else{
          dir = opendir(directory);
          if(dir == NULL) return 0;
    int cnt = 0;
          while((s_dir=readdir(dir))!=NULL)
          {
                  if((strcmp(s_dir->d_name,".")==0)|| (strcmp(s_dir->d_name,"..")==0))
                  continue;

                  sprintf(currfile,"%s/%s",directory,s_dir->d_name);
                  stat(currfile,&file_stat);

                  if(!S_ISDIR(file_stat.st_mode))
                  {
        sprintf(buf,"%d\t%s",cnt+1,s_dir->d_name);
        file.WriteLine("mark.idx",buf);
        cnt++;
                  }
          }
          closedir(dir);
    cout<<"Building mark.idx.............[OK!]"<<endl;
  }
  file.LastWrite();
  file.Release();
  return 1;
}
void *func(void* arg)
{
        char delimiter='\t';
        DIR *dir;
        char currfile[512];
        struct dirent *s_dir;
        struct stat file_stat;

        _arg arg1 = *((_arg*)arg);
        char* directory = arg1.directory;
        int func_id = arg1.func_id;
	GS* gs = arg1.gs;

        stat(directory,&file_stat);
	// common file
        if(!S_ISDIR(file_stat.st_mode) && func_id == 0)
        {
                strcpy(currfile,directory);
                cout<<"processingF "<<currfile<<endl;
    if(para.format == 0)
                  ExtractBed(currfile,basename_(currfile).c_str(),delimiter,gs);
    if(para.format == 1)
                  ExtractWig(currfile,basename_(currfile).c_str(),delimiter,gs);
    if(para.format == 2)
                  ExtractBin(currfile,basename_(currfile).c_str(),delimiter,gs);
        }
        // directory
        else{
                dir = opendir(directory);
                if(dir == NULL) return 0;
                int cnt = 0;
                while((s_dir=readdir(dir))!=NULL)
                {
                        if((strcmp(s_dir->d_name,".")==0)|| (strcmp(s_dir->d_name,"..")==0))
                        continue;

                        sprintf(currfile,"%s/%s",directory,s_dir->d_name);
                        stat(currfile,&file_stat);

                        if(!S_ISDIR(file_stat.st_mode))
                        {
                                if(cnt % para.max_thread == func_id) // assign job according to mod
                                {
                                        cout<<"processingD "<<currfile<<endl;
          if(para.format == 0)
                                          ExtractBed(currfile, s_dir->d_name, delimiter,gs);
          if(para.format == 1)
                        ExtractWig(currfile, s_dir->d_name, delimiter,gs);
          if(para.format == 2)
                        ExtractBin(currfile, s_dir->d_name, delimiter,gs);
                                }
                                cnt++;
                        }
                }
                closedir(dir);
        }
  
        return NULL;
}

// program pre_EpiSig starts here
int main(int argc, char* argv[]) 
{
        char directory[512]; // directory name for processing

  // default para definition, see .h file
  para.max_thread = 6;
  para.format = 0;
  para.skip = 0;

  // get options
        int next_option;
        const char* const short_options = "hd:f:sn:i:";
        const struct option long_options[]={
                {"help",0,NULL,'h'},
                {"data-dir",1,NULL,'d'},
                {"format",1,NULL,'f'},
                {"skip-norm",0,NULL,'s'},
                {"threads",1,NULL,'n'},
                {"genome",1,NULL,'i'},
                {NULL,0,NULL,0}
        };
  int data_dir = 0;
  int data_format = 0;
	int flag_genome = 0;
  char buf[256];
         do{
                next_option = getopt_long(argc,argv,short_options,long_options,NULL);
                switch(next_option)
                {
                        case 'h':
                                print_usage();
                                return 0;
                        case '?':
                                print_usage();
                                return 0;
                        case ':':
                                print_usage();
                                return 0;
                        case 'd':
                                para.data_dir = optarg;
        data_dir = 1;
                                break;
                        case 'f':
        strcpy(buf, optarg);
        Trim(buf);
        MakeLower(buf);
        if(strncmp(buf,"bed",3) == 0) para.format = 0;
        else if(strncmp(buf,"wig",3) == 0) para.format = 1;
        else if(strncmp(buf,"bin",3) == 0) para.format = 2;
        else{
                      print_usage();
          printf("You should specify [BED | WIG | BIN] format\n");
          return 0;
        }
        data_format = 1;
                                break;
                        case 's':
        para.skip = 1;
                                break;
                        case 'n':
                                para.max_thread = atoi(optarg);
                                if(para.max_thread < 1 || para.max_thread > 16)
                                {
                                        printf("max thread number should between 1 and 16!\n");
                                        return 0;
                                }
                               break;
			case 'i':
				para.genome_name = string(optarg);
				flag_genome = 1;
				break;
                        default:
                                break;
                }
        }while(next_option !=-1);

  if(data_dir == 0)
        {
                print_usage();
                printf("You must specify input data directory: [--data-dir]\n");
                return 0;
        }
        if(data_format == 0)
        {
                print_usage();
                printf("You must specify argument: [--format BED|WIG|BIN]\n");
                return 0;
        }
	if( flag_genome == 0 )
	{
		print_usage();
		printf("You must specify genome name: -g\n" );
		return 0;
	}

  print_para(para);
  // end get options

	GS gs;
	GS_init(gs,para.genome_name);

  mkdir("out_data",0755); // make directory
  // re-format directory name of data source directory
  strncpy(directory,para.data_dir.c_str(),512); // copy directory name
  Trim(directory); // remove space, like chomp in perl
  int len = strlen(directory);
  if(len < 1) return 0;
  if(directory[len-1] == '/') directory[len-1] = '\0';

  // running time setup
  struct timeval tpstart,tpend;
  float timeuse;
  gettimeofday(&tpstart,NULL); 
  
  struct _arg arg[MAX_THREAD];
  for(int i = 0; i < MAX_THREAD; i++)
  {
    arg[i].func_id = i;
    strcpy(arg[i].directory,directory);
	arg[i].gs = &gs;
  }
  pthread_t id[MAX_THREAD];
  int ret[MAX_THREAD];
  //int *ret = new int[para.max_thread];
  for(int i = 0; i < para.max_thread; i++)
  ret[i] = pthread_create(&id[i],NULL,func,(void*)(&arg[i]));
  for(int i = 0; i < para.max_thread; i++)
  pthread_join(id[i], NULL);

  make_mark_idx(directory);

  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  timeuse/=1000000;

  cout<<"......................finished! Total time:"<<timeuse<<" seconds!"<<endl;
  //delete[] ret;
  return 0;
} 
