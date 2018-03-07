#include <cstdio>
#include <iostream>
#include "ConfUtils.h"

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include "utils.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
ConfUtils::ConfUtils()
{
}
ConfUtils::~ConfUtils()
{
}
int ConfUtils::GetPara(char* line,const char* para,char* value)
{
	if(line == NULL) return 0;

	Trim(line);
	MakeLower(line);

	char* p =strstr(line,para);
	if(!p) return 0;

	char* l = strstr(p,"=");
	if(!l) return 0;

	p = l + 1;

	int i = 0;
	while(isspace(*(p+i)) && (*(p+i) != '\0')) i++;

	int j = 0;
	while((!isspace(*(p+i+j))) && (*(p+i+j) != '\0')) j++;
	if(j > 0) memcpy(value,p+i,j);
	value[j] = '\0';

	return 1;
}
