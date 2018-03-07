#include<cstring>
#include<cstdlib>
#include "thread_random.h"

THREAD_RAND::THREAD_RAND( int thread_num, int seed )
{
	srandom(seed);

        size = 1000000;
        rand = new long[size];
        nb = new int[size];

        for( int i=0; i<size; i++ )
        {
                *(rand+i) = random();
                *(nb+i) = i+1;
        }
        *(nb+size-1) = 0;
}

THREAD_RAND::~THREAD_RAND()
{
	delete[] rand;
	delete[] nb;
}

long THREAD_RAND::rand_num( int& idx )
{
	long r = *(rand+idx);
	idx = *(nb+idx);

	return r;
}
