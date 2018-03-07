#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>

#include "utils.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
int cmp_int(const void *a, const void *b)
{
	int a_int = *((int*)a);
	int b_int = *((int*)b);
	if(a_int > b_int) return 1;
	if(a_int < b_int) return -1;
	return 0;
}

void spliter(const char* line,char delimiter,string* rt,int& cnt)
{
        char field[128];
        char ch;

        int i = 0;
        int pos = 0;
        int count = 0;
        field[0] ='\0';

        while((ch=*(line+i)) != '\0')
        {
                //find the column
                if(ch == delimiter)
                {
                        field[pos] = 0;
                        pos = 0;
                        i++;

                        Trim(field);
			if(count < cnt)
			{
				*(rt+count) = field;
                        	count++;
			}
			if(count == cnt) break;
                }
                else
                {
			if(pos < 127)
			{
                        	field[pos] = ch;
                        	pos++;
			}
                        i++;
                }
        }
        field[pos] = 0;
        Trim(field);
	if(count < cnt && field[0] != '\0') 
	{
		*(rt+count) = field;
		count++;
	}
	cnt = count;
}
int FileLength(string fname)
{
  int fd;
  struct stat stat;
  
  fd = open(fname.c_str(), O_RDONLY);
  if(fd < 0)
  {
	return -1;
  }	

  if( fstat(fd,&stat) < 0 )
  {
	return -1;
  }
  else
  {
	return stat.st_size;
  }
}

char* MakeLower(char* str)
{
	if(!str) return NULL;
	int i = 0;
	char ch;
        while((ch=*(str+i)) != '\0')
	{
		*(str+i) = tolower(ch);
		i++;
	}
	return str;
}
char* TrimLeft(char* str)
{
	if(!str) return NULL;

	int i = 0;
	char ch;
        while((ch=*(str+i)) != '\0')
	{
		if(!ISSPACE(ch)) break;
		i++;
	}

	if(i > 0)
	{
		char *p = str+i;
		i = 0;
		while(*(p+i) != '\0') 
		{
			*(str+i) = *(p+i);
			i++;
		}
		*(str+i) = '\0';
	}

	return str;
} 
char* TrimRight(char* str)
{
	if(!str) return NULL;

	int len = strlen(str);
	int i = len-1;
	while(i>=0)
	{
		char ch = *(str+i);
		if(!ISSPACE(ch)) break;
		else
		{
			*(str+i) = '\0';
			i--;
		}
	}

	return str;
}
char* Trim(char* str)
{
	TrimRight(str);
	TrimLeft(str);
	return str;
}
string basename_(char* path)
{
	if(!path) return "";
	string basename_ = "";

	char* p = path;
	char* p1 = NULL;
	while(p != NULL)
	{
		if(p1 != NULL) p++;
		p1 = p;
		p = strstr(p,"/");
	}
	basename_ = p1;
	return basename_;
}
