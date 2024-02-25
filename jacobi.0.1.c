/************************************************************************/
/*  jacobi.0.1.c .					                */
/*  opera con una matrice di esempio in formato non compatto (non CSR)  */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 10

double A[N][N]={
		{10,1,1,1,1,1,1,1,1,1},
	       	{1,10,1,1,1,1,1,1,1,1},
	       	{1,1,10,1,1,1,1,1,1,1},
		{1,1,1,10,1,1,1,1,1,1},
		{1,1,1,1,10,1,1,1,1,1},
		{1,1,1,1,1,10,1,1,1,1},
	       	{1,1,1,1,1,1,10,1,1,1},
	       	{1,1,1,1,1,1,1,10,1,1},
		{1,1,1,1,1,1,1,1,10,1},
		{1,1,1,1,1,1,1,1,1,10},
	      };

double b[N]={2,3,4,1,5,7,8,1,2,5};
double x0[N]={0,0,0,0,0,0,0,0,0,0};
double x1[N];

int main()
{	
	/***************************/
	/*  vairabili di supporto  */
	/***************************/	
	const unsigned long int NMAX=pow(N, 4);
	const double TOL_QUAD=1e-20;
	double e_num=0, e_den=0;
	unsigned long int i,j,k=0;
	double* x_k=x0; 
	double* x_k_1=x1; 
	double* x_tmp;
	double sum;
	
	/********************************************/
	/*	       		                    */
	/*                  JACOBI                  */
	/*		                            */
	/********************************************/		
	while(k < NMAX)
	{
		/***********************************calcolo di x(k+1)**********************************/
		/*										      */
		/*   xi(k+1) = (1/aij)*[ bi - ( somm( aij xj(k) ) ) {j=1..n and j!=i} ]  { i=1..n }   */
		/*										      */
		/**************************************************************************************/ 
		
		for (i=0; i<N; i++)
		{
			sum=0;
			for(j=0; j<N; j++)
			{
				if(i!=j)
				{
					sum = sum + A[i][j]*x_k[j];					
				}
			}
			x_k_1[i]=(1/A[i][i])*(b[i] - sum);		
		}


		/***********condizione d'uscita***********/
		/* 					 */
		/*   ||x(k+1)-x(k)|| / ||x(k+1)|| < TOL  */
		/*					 */
		/*****************************************/
		e_num=0;
		for(i=0; i<N; i++)
		{
			e_num = e_num + pow(x_k_1[i] - x_k[i], 2);
		}
	
		e_den=0;
		for(i=0; i<N; i++)
		{
			e_den = e_den + pow(x_k_1[i], 2);
		}

		if(e_num/e_den < TOL_QUAD)
		{
			break;
		}

		/*************nuova iterazione************/
		/*					 */
		/*	       x(k) = x(k+1)             */
		/* 					 */
		/*****************************************/
		x_tmp=x_k;
		x_k=x_k_1;
		x_k_1=x_tmp;	
		k++;	
	}

	//printf("double: %lu\n", sizeof(double));
	//printf("long double: %lu\n", sizeof(long double));


	/********************************/
	/*				*/
	/*    stampa della soluzione	*/
	/*				*/
	/********************************/
	for(i=0; i<N; i++)
	{
		printf("%le\n", x_k_1[i]);
	}
	printf("iterazioni:%ld\n", k+1);
	printf("errore: %le\n", e_num/e_den);

	exit(0);
}

