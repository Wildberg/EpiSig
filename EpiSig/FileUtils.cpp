#include <cstdio>
#include <iostream>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "utils.h"
#include "FileUtils.h"
#include <stdlib.h>
#include <string.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
FileUtils::FileUtils()
{
	fp = NULL;
	out_fp = NULL;		
	g_buf_size = 4*1024*1024;
	real_buf_size = 0;
	g_fp_pos = 0;
	g_buf_pos = 0;
	out_buf_pos = 0;
	g_buffer = NULL;
	out_buffer = NULL;
}
FileUtils::FileUtils(int buf_size)
{
	fp = NULL;
	out_fp = NULL;
	if(buf_size < 4096) g_buf_size = 4096;
	else g_buf_size = buf_size;

	real_buf_size = 0;
	g_fp_pos = 0;
	g_buf_pos = 0;
	out_buf_pos = 0;
	g_buffer = NULL;
	out_buffer = NULL;
}
FileUtils::~FileUtils()
{
	if(g_buffer) delete[] g_buffer;
	if(out_buffer) delete[] out_buffer;
}
bool FileUtils::OpenFileForRead(string fname)
{
	fp = fopen(fname.c_str(), "r");
	if(!fp) return 0;

	return 1;
}
bool FileUtils::FirstRead()
{
	real_buf_size = fread(g_buffer, sizeof(char), g_buf_size, fp);
	if(real_buf_size < 0) return 0;

	return 1;
}
bool FileUtils::AllocBuffer()
{
	g_buffer = new char[g_buf_size];
	out_buffer = new char[g_buf_size];

	if(g_buffer && out_buffer) return 1;
	else return 0;
}
void FileUtils::Release()
{
	if(fp) fclose(fp);
	if(out_fp) fclose(out_fp);
	if(g_buffer)delete[] g_buffer;
	if(out_buffer) delete[] out_buffer;
	fp = out_fp = NULL;
	g_buffer = out_buffer = NULL;
	real_buf_size = 0;
	g_fp_pos = 0;
	g_buf_pos = 0;
	out_buf_pos = 0;
}
int FileUtils::CountLines(string fname)
{
	FILE* fp = fopen(fname.c_str(), "r");
	if(!fp) return -1;

	char* buffer=NULL;
	int buf_size=1024*1024;

	buffer = new char[buf_size];
	if(!buffer) return -1;

	fseek(fp,0,SEEK_SET);
	int sizeread = -1;

	int count=0;
	while(1)
	{
		sizeread = fread(buffer, sizeof(char), buf_size, fp);

		int old_pos = 0;
		for(int i=0;i<sizeread;i++)
		{
			if(*(buffer+i) == '\n') 
			{
				old_pos = i+1;
				count++;
			}
		}
		//EOF or error occurs
		if(feof(fp)) 
		{
			if(old_pos < sizeread) count++; // the rest of lines until error/eof
			break;
		}
	}
	fclose(fp);
	delete[] buffer;

	return count;
}
int FileUtils::GetLine(char* line, int line_size)
{
	int i;
	int count = 0;

	if(line == NULL) return 0;
	while(1)
	{
		for(i=g_buf_pos;i<real_buf_size;i++)
		{
			char ch = *(g_buffer+i);
			if(ch == '\n')
			{
				if(count >= line_size) count = line_size - 1;
				*(line+count) = '\0';
				g_buf_pos = i+1;
				return 1;
			}
			else
			{
				if(count < line_size-1) *(line+count) = ch;
				count++;
			}
		}
		//eof of the file
		if(feof(fp))
		{
			if(count < line_size-1) *(line+count) = '\0';
			else  *(line+line_size-1) = '\0';
			return 0;
		}

		//not eof of file and did not get a line,refill the buffer and continue fill line buffer
		g_buf_pos = 0;
		real_buf_size = fread(g_buffer, sizeof(char), g_buf_size, fp);
	}

	return 1;
}
int FileUtils::WriteLine(char* line)
{
	if(!line) return -1;
	int len = strlen(line);

	//buffer is full
	if(out_buf_pos+len+1 > g_buf_size)
	{
		//write to file
		fwrite(out_buffer, sizeof(char), out_buf_pos, out_fp);
		out_buf_pos = 0;
	}

	memcpy(out_buffer+out_buf_pos,line,len);
	out_buf_pos += len;
	*(out_buffer + out_buf_pos) = '\n';
	out_buf_pos++;

	return 1;
}
bool FileUtils::WriteLine(string fname,char* line)
{
	if(out_fp == NULL)
	{
		out_fp = fopen(fname.c_str(), "w+");
		if(!out_fp) return 0;
        	out_buffer = new char[g_buf_size];
        	if(!out_buffer) 
		{
			fclose(out_fp);
			out_fp = NULL;
			return 0;
		}
	}

	WriteLine(line);
	return 1;
}
bool FileUtils::WriteLine(string fname,char* line,const char* tag)
{
	if(out_fp == NULL)
	{
		out_fp = fopen(fname.c_str(), tag);
		if(!out_fp) return 0;
        	out_buffer = new char[g_buf_size];
        	if(!out_buffer) 
		{
			fclose(out_fp);
			out_fp = NULL;
			return 0;
		}
	}

	WriteLine(line);
	return 1;
}
bool FileUtils::LastWrite()
{
	if(out_buf_pos > 1) fwrite(out_buffer, sizeof(char), out_buf_pos-1, out_fp);
	out_buf_pos = 0;

	return 1;
}
