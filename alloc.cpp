// Onur G. Guleryuz 1995, 1996, 1997,
// University of Illinois at Urbana-Champaign,
// Princeton University,
// Polytechnic University.

#include <stdlib.h>
#include "macros.h"

float **allocate_2d_float(int N,int M,char zero)

{
	int i;
	float **mymat;
	
	mymat=(float **)malloc(N*sizeof(float *));
	check_ptr(mymat,"allocate_2d_float");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(float *)malloc(M*sizeof(float));
			check_ptr(mymat[i],"allocate_2d_float");
	}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(float *)calloc(M,sizeof(float));
			check_ptr(mymat[i],"allocate_2d_float");
	}
	return(mymat);
}

void free_2d_float(float **a,int N)

{
	int i;
	
	for(i=0;i<N;i++)
		free((void *)a[i]); 
	free((void *)a);
}

double **allocate_2d_double(int N,int M,char zero)

{
	int i;
	double **mymat;
	
	mymat=(double **)malloc(N*sizeof(double *));
	check_ptr(mymat,"allocate_2d_double");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(double *)malloc(M*sizeof(double));
			check_ptr(mymat[i],"allocate_2d_double");
		}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(double *)calloc(M,sizeof(double));
			check_ptr(mymat[i],"allocate_2d_double");
		}
	return(mymat);
}

void free_2d_double(double **a,int N)

{
	int i;
	
	for(i=0;i<N;i++)
		free((void *)a[i]); 
	free((void *)a);
}


int **allocate_2d_int(int N,int M,char zero)

{
	int i;
	int **mymat;
	
	mymat=(int **)malloc(N*sizeof(int *));
	check_ptr(mymat,"allocate_2d_int");
	if(!zero)
		for(i=0;i<N;i++) {
			mymat[i]=(int *)malloc(M*sizeof(int));
			check_ptr(mymat[i],"allocate_2d_int");
	}
	else
		for(i=0;i<N;i++) {
			mymat[i]=(int *)calloc(M,sizeof(int));
			check_ptr(mymat[i],"allocate_2d_int");
	}
	return(mymat);
}

void free_2d_int(int **a,int N)

{
	int i;
	
	for(i=0;i<N;i++)
		free((void *)a[i]); 
	free((void *)a);
}

float *allocate_1d_float(int N,char zero)

{
	float *arr;
	
	if(!zero)
		arr=(float *)malloc(N*sizeof(float));
	else
		arr=(float *)calloc(N,sizeof(float));
	check_ptr(arr,"allocate_1d_float");
	return(arr);
}

double *allocate_1d_double(int N,char zero)

{
	double *arr;
	
	if(!zero)
		arr=(double *)malloc(N*sizeof(double));
	else
		arr=(double *)calloc(N,sizeof(double));
	check_ptr(arr,"allocate_1d_double");
	return(arr);
}

int *allocate_1d_int(int N,char zero)

{
	int *arr;
	
	if(!zero)
		arr=(int *)malloc(N*sizeof(int));
	else
		arr=(int *)calloc(N,sizeof(int));
	check_ptr(arr,"allocate_1d_int");
	return(arr);
}

unsigned char *allocate_1d_uchar(int N,char zero)

{
	unsigned char *arr;
	
	if(!zero)
		arr=(unsigned char *)malloc(N*sizeof(unsigned char));
	else
		arr=(unsigned char *)calloc(N,sizeof(unsigned char));
	check_ptr(arr,"allocate_1d_uchar");
	return(arr);
}

