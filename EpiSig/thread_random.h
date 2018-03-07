#ifndef __THREAD_RAND_H
#define __THREAD_RAND_H

class THREAD_RAND
{
public:
	THREAD_RAND( int thread_num, int seed );
	~THREAD_RAND();
	long rand_num(int& idx);
private:
	long* rand;
	int* nb;
	int size;
};
#endif
