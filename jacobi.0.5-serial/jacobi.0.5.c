/*****************************************************************************/
/* jacobi.0.5.c .					                     */
/* - a differenza della versione 0.4 l'output avviene su file		     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct
{
	unsigned long int n_row, n_col, n_elem; //dimensioni matrice
	
	unsigned int simm; //simmetria matrice	
		
	double* le;	        //vettori CSR			
	unsigned long int* lr;	//		
	unsigned long int* lc;	//	
}
CSR_matrix;

int main(int argc, char *argv[])
{
	/********************************/
	/*  allocazione struttura CSR   */
	/********************************/
	CSR_matrix A;
	
	/***************************/
	/*  vairabili di supporto  */
	/***************************/
	double e_num, e_den;
	unsigned long int i=0, j, k=0, row, iter;
	double* x_tmp;
	double sum, x_k_1_quad, aii;

	/*****************************************************************/
	/*  apertura file e raccolta info matrice A del sistema lineare  */
	/*****************************************************************/	
	FILE *f=fopen(argv[1], "r");
	fscanf(f,"%lu,%lu,%lu,%d\n", &(A.n_row), &(A.n_col), &(A.n_elem), &(A.simm)); //dimensioni matrice di input
	const unsigned NMAX=pow(A.n_row,4); //numero massimo di iterazioni dell'algoritmo
	const double TOL_QUAD=1e-10; //tolleranza

	/******************************/
	/*  allocazione vettori CSR   */
	/******************************/
	A.le=(double*)malloc(A.n_elem*sizeof(double));			     //le = vettore degli elementi no nulli di A[N][N] (ordinati per riga)
	A.lr=(unsigned long int*)malloc(A.n_row*sizeof(unsigned long int));  //lr = vettore degli indici, in le, dei primi elementi, non nulli, di riga, della matrice A[N][N]
	A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int)); //lc = vettore degli indici colonna, in A[N][N], degli elementi di "le"

	/***************************************************/
	/*  - allocazione vettore "b" del sistema lineare  */ 
	/*  - allocazione del vettore "x(k)" e "x(k+1)"    */
	/***************************************************/
	double* b=(double*)malloc(A.n_col*sizeof(double));
	double* x_k=(double*)malloc(A.n_row*sizeof(double));
	double* x_k_1=(double*)malloc(A.n_row*sizeof(double));

	/*********************************/
	/*  definizione dei vettori CSR  */
	/*********************************/	
	for(k=0;k<A.n_elem; k++)
	{
		fscanf(f, "%lu,%lu,%le\n", &row, &(A.lc[k]), &(A.le[k]));
		if(row == i)
		{
			A.lr[row]=k;
			i++;
		}
	}
	A.lr[i]=A.n_elem+1;

	/******************************************/
	/*  definizione del vettore "b" e "x(0)"  */
	/******************************************/
	for(i=0; i<A.n_row; i++)
	{
		fscanf(f, "%le\n", &b[i]);
		x_k[i]=0;
	}
	fclose(f);	
	
	/********************************************/
	/*	       		                    */
	/*                  JACOBI                  */
	/*		                            */
	/********************************************/	
	for(iter=0; iter<NMAX; iter++)
	{
		/***********************************calcolo di x(k+1)**********************************/
		/*										      */
		/*   xi(k+1) = (1/aij)*[ bi - ( somm( aij xj(k) ) ) {j=1..n and j!=i} ]  { i=1..n }   */
		/*										      */
		/**************************************************************************************/ 
		e_num=0; e_den=0;
		
		for(i=0; i<A.n_row; i++)
		{
			x_k_1[i] = b[i];
			
			for(k=A.lr[i]; k<A.lr[i+1]; k++)
			{
				j=A.lc[k];
				if(i!=j)
				{
					x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[j] );
				}
				else
				{
					aii=A.le[k];
				}		
			}
			x_k_1[i] = x_k_1[i]/aii;

			/****rapporto per condizione d'uscita*****/
			/* 					 */
			/*      ||x(k+1)-x(k)|| / ||x(k+1)||	 */
			/*					 */
			/*****************************************/
			x_k_1_quad = x_k_1[i] * x_k_1[i];
			e_den=e_den + x_k_1_quad;
			e_num=e_num + x_k_1_quad + x_k[i]*x_k[i] - 2*x_k_1[i]*x_k[i];		
		}

		/*****valutazione condizione d'uscita*****/
		/* 					 */
		/*   ||x(k+1)-x(k)|| / ||x(k+1)|| < TOL  */
		/*					 */
		/*****************************************/
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
	}

	/********************************/
	/*				*/
	/*    stampa della soluzione	*/
	/*				*/
	/********************************/
	f=fopen(strcat(argv[1],"_out.dat"), "w");

	for(i=0; i<A.n_row; i++)
	{
		fprintf(f,"%le\n", x_k_1[i]);
	}
	fclose(f);
	printf("\niterazioni:%ld\n", iter+1);
	printf("tolleranza: %le\n", e_num/e_den);
	
	exit(0);
}




