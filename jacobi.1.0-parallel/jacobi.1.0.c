/*****************************************************************************/
/* jacobi.1.0.c parallelo						     */
/*    	il coordinatore ora è il processo con rank più alto		     */
/*	in modo da sfruttare meglio la funzione "gather" per la		     */
/*	sincronizzazione						     */
/*****************************************************************************/

#include <mpi.h>
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

void stampa_CSR(unsigned long int n_row, unsigned long int n_col, unsigned long int n_elem, double* le, unsigned long int* lr, unsigned long int* lc, double* b)
{
	int k;

	printf("n_row: %lu\n", n_row);
	printf("n_col: %lu\n", n_col);
	printf("n_elem: %lu\n", n_elem);
	printf("le: ");
	for(k=0; k<n_elem; k++)
	{
		printf("%le ", le[k]);
	}
	printf("\nlr: ");
	for(k=0; k<n_row+1; k++)
	{
		printf("%lu ", lr[k]);
	}
	printf("\nlc: ");
	for(k=0; k<n_elem; k++)
	{
		printf("%lu ", lc[k]);
	}
	printf("\nb: ");
	for(k=0; k<n_row; k++)
	{
		printf("%le ", b[k]);
	}
	printf("\n\n");
}

int main(int argc, char *argv[])
{

	//==================================================workers===================================================//

	CSR_matrix A;//matrice CSR
	int size;//numero processi
	int rank;//id processo
	unsigned long int k, i, j;//contatori
	unsigned long int dim;//dimensione buffer
	int position;//posizione elemento corrente nel buffer
	char* buffer1;//buffer 
	char* buffer2;//buffer

	double e_num, e_den;
	double* x_tmp;
	double sum, x_k_1_quad, aii;
	unsigned long int iter;

	double* b;
	double* x_k;
	double* x_k_1;


	unsigned NMAX;//numero massimo di iterazioni dell'algoritmo
	double TOL_QUAD; //tolleranza

	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD);
	
	//=================================================coordinatore===============================================//
	if(rank==(size-1))
	{	
	
		//
		//___________________creazione struttura CSR generale
		//

		//apertura file e raccolta info matrice A del sistema lineare
		FILE *f=fopen(argv[1], "r");
		fscanf(f,"%lu,%lu,%lu,%d\n", &(A.n_row), &(A.n_col), &(A.n_elem), &(A.simm));
		const unsigned NMAX=pow(A.n_row,4);	

		//le = vettore degli elementi non nulli di A[N][N] (ordinati per riga)
		A.le=(double*)malloc(A.n_elem*sizeof(double));

		//lr = vettore degli indici, in le, dei primi elementi, non nulli, di riga, della matrice A[N][N]
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));

		//lc = vettore degli indici colonna, in A[N][N], degli elementi di "le"   
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));			

		//definizione dei vettori CSR
		unsigned long int row;
		i=0; 
		for(k=0; k<A.n_elem; k++)
		{
			fscanf(f, "%lu,%lu,%le\n", &row, &(A.lc[k]), &(A.le[k]));
			if(row == i)
			{
				A.lr[row]=k;
				i++;
			}
		}
		A.lr[i]=A.n_elem;
		
		//
		//___________________allocazione e definizione vettore b
		//

		//allocazione vettore "b" del sistema lineare
		b=(double*)malloc(A.n_col*sizeof(double));

		//definizione del vettore b
		for(k=0; k<A.n_row; k++)
		{
			fscanf(f, "%le\n", &b[k]);
		}	
		fclose(f);

		printf("STRUTTURA CSR GENERALE:\n\n");
		stampa_CSR(A.n_row, A.n_col, A.n_elem, A.le, A.lr, A.lc, b);

				
		//
		//_________________________creazione e invio delle strutture CSR e delle parti di vettore b, per i worker 
		//
		//

		//definizione del job dei worker
		unsigned long int job_worker=A.n_row/size;	
		unsigned long int job_max;
		unsigned long int job_min;
		int worker; int j_worker; int i_worker;
		unsigned long int worker_nelem;

		//vettori da inviare ai worker afinchè inizializzino
		//la propria struttura CSR
		unsigned long int* v_nrow;
		unsigned long int* v_ncol;
		unsigned long int* v_nelem;
		unsigned int* v_simm;
		unsigned long int** v_lr;
		unsigned long int** v_lc;
		double** v_le;

		//allocazione vettori che conterranno i vettori CSR per i worker
		v_nrow = (unsigned long int*)malloc(size*sizeof(unsigned long int));
		v_ncol = (unsigned long int*)malloc(size*sizeof(unsigned long int));
		v_nelem = (unsigned long int*)malloc(size*sizeof(unsigned long int));
		v_simm = (unsigned int*)malloc(size*sizeof(unsigned int));
		v_le = (double**)malloc(size*sizeof(double*));
		v_lr = (unsigned long int**)malloc(size*sizeof(unsigned long int*));
		v_lc = (unsigned long int**)malloc(size*sizeof(unsigned long int*));

	
		//ciclo su tutti i worker
		for(worker=0; worker<(size-1); worker++)
		{
	
			//
			//___________________________________________definizione del job per il worker corrente
			//
			
			//limiti del job per il worker corrente
			job_min = worker*job_worker /*- job_worker*/;
			job_max = job_min + job_worker;
			worker_nelem=A.lr[job_max]-A.lr[job_min];			
								
			//
			//__________________________________definizione e invio valori n_row, n_col, n_elem e simm per il worker corrente
			//
		
			//definizione "n_row", "n_col", "n_elem" e "simm" per il worker corrente
			v_nrow[worker]=job_worker;
			v_ncol[worker]=A.n_col;
			v_nelem[worker]=worker_nelem;
			v_simm[worker]=0;
	
			//allocazione buffer info
			dim=
			3*sizeof(unsigned long int)+//n_row, n_col, n_elem
			sizeof(unsigned int);//simm
			buffer1=(char*)malloc(dim);

			//packaging info struttura CSR
			position=0;
			MPI_Pack(&(v_nrow[worker]), 1, MPI_UNSIGNED_LONG, buffer1, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(v_ncol[worker]), 1, MPI_UNSIGNED_LONG, buffer1, dim, &position, MPI_COMM_WORLD);		
			MPI_Pack(&(v_nelem[worker]), 1, MPI_UNSIGNED_LONG, buffer1, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(v_simm[worker]), 1, MPI_UNSIGNED, buffer1, dim, &position, MPI_COMM_WORLD);
			
			//invio al worker le info per allocare la sua struttura
			MPI_Send(buffer1, dim, MPI_PACKED, worker, 0, MPI_COMM_WORLD);				
						
			//deallocazione del buffer per le info sulla struttura CSR	
			free(buffer1);	
	
			//
			//___________________________________allocazione, definizione ed invio vettori le, lr, lc e b per il worker corrente
			//
			
			//allocazione dei vettori "le", "lc" e "lr" per il worker corrente
			v_le[worker]=(double*)malloc(worker_nelem*sizeof(double) );
			v_lc[worker]=(unsigned long int*)malloc(worker_nelem*sizeof(unsigned long int) );	
			v_lr[worker]=(unsigned long int*)malloc( (job_worker+1)*sizeof(unsigned long int) );

			//definizione dei vettori "le", "lr" e "lc" per il worker corrente
			for(i=job_min, i_worker=0, j_worker=0; i<job_max; i++, i_worker++)
			{	
				//lr
				v_lr[worker][i_worker]=j_worker;		

				//le, lc
				for(j=A.lr[i]; j<A.lr[i+1]; j++)
				{		
					v_lc[worker][j_worker]=A.lc[j];
					v_le[worker][j_worker]=A.le[j];
					j_worker++;
				}		
			}
			//[lr]=n+1
			v_lr[worker][i_worker]=j_worker;

			//printf("STRUTTURA CSR WORKER %d:\n\n", worker);
			//stampa_CSR(v_nrow[worker], v_ncol[worker], v_nelem[worker], v_le[worker], v_lr[worker], v_lc[worker]);
									
			//allocazione buffer vettori CSR e parte del vettore b
			dim=
			worker_nelem*sizeof(double)+//le
			(job_worker+1)*sizeof(unsigned long int)+//lr	
			worker_nelem*sizeof(unsigned long int)+//lc
			job_worker*sizeof(double);//parte del vettore b
			buffer2=(char*)malloc(dim);
	
			//packaging vettori CSR e vettore b
			position=0;
			MPI_Pack(v_le[worker], worker_nelem, MPI_DOUBLE, buffer2, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(v_lr[worker], job_worker+1, MPI_UNSIGNED_LONG, buffer2, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(v_lc[worker], worker_nelem, MPI_UNSIGNED_LONG, buffer2, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(b[job_min]), job_worker, MPI_DOUBLE, buffer2, dim, &position, MPI_COMM_WORLD);
			
			//invio al worker i vettori CSR e il vettore b
			MPI_Send(buffer2, dim, MPI_PACKED, worker, 1, MPI_COMM_WORLD);		

			//deallocazione buffer vettori CSR
			free(buffer2);
			free(v_le[worker]);
			free(v_lr[worker]);
			free(v_lc[worker]);			
		}
		
		//deallocazione info per i workers
		free(v_nrow);
		free(v_ncol);
		free(v_nelem);
		free(v_simm);

		//printf("\n-processo %d: inviate le strutture CSR ai worker\n", rank);	
	}
	
	//==================================================workers===================================================//
	if(rank!=(size-1)){
		
		//
		//__________________________ricezione struttura CSR 
		//
		
		//ricezione info dimensioni struttura
		dim=
		3*sizeof(unsigned long int)+//n_row, n_col, n_elem
		sizeof(unsigned int);//simm
		buffer1=(char*)malloc(dim);	

		MPI_Recv(buffer1, dim, MPI_PACKED, size-1, 0, MPI_COMM_WORLD, &status);
		
		//definizione "n_row", "n_col" e "n_elem" della struttura CSR del worker
		position=0;
		MPI_Unpack(buffer1, dim, &position, &(A.n_row), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer1, dim, &position, &(A.n_col), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer1, dim, &position, &(A.n_elem), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer1, dim, &position, &(A.simm), 1, MPI_UNSIGNED, MPI_COMM_WORLD);

		//definizione del massimo numero di iterazioni e 
		//della tolleranza dell'algoritmo di jacobi
		NMAX=pow(A.n_row,4);
		TOL_QUAD=1e-8;

		//allocazione vettori "le", "lr" "lc"  e "b" del worker
		A.le=(double*)malloc(A.n_elem*sizeof(double));
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));
		b=(double*)malloc(A.n_col*sizeof(double));

		//deallocazione buffer
		free(buffer1);

		//allocazione buffer vettori CSR e del vettore b
		dim=
		(A.n_elem)*sizeof(double)+//le
		(A.n_row+1)*sizeof(unsigned long int)+//lr	
		(A.n_elem)*sizeof(unsigned long int)+//lc
		(A.n_row)*sizeof(double);//parte del vettore b
		buffer2=(char*)malloc(dim);

		//allocazione vettori x(k) e x(k+1)
		x_k_1=(double*)malloc(A.n_row*sizeof(double));
		x_k=(double*)malloc(A.n_col*sizeof(double));
		
		//x(k)=0;
		for(k=0; k<A.n_row; k++){
			x_k[k]=0;
		}

		//ricezione vettori CSR e vettore b
		position=0;
		MPI_Recv(buffer2, dim, MPI_PACKED, size-1, 1, MPI_COMM_WORLD, &status);
		MPI_Unpack(buffer2, dim, &position, A.le, A.n_elem, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(buffer2, dim, &position, A.lr, A.n_row+1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer2, dim, &position, A.lc, A.n_elem, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer2, dim, &position, b, A.n_row, MPI_DOUBLE, MPI_COMM_WORLD);
	
		//deallocazione buffer
		free(buffer2);

		//printf("\nprocesso %d: ricevuta la struttura CSR\n", rank);
		printf("STRUTTURA CSR WORKER %d:\n\n", rank);
		stampa_CSR(A.n_row, A.n_col, A.n_elem, A.le, A.lr, A.lc, b);
		
		/*	
		//
		//_______________________________________JACOBI                  
		//
		for(iter=0; iter<NMAX; iter++)
		{
			//_______________________________________________________________calcolo di x(k+1)
			//
			//   xi(k+1) = (1/aij)*[ bi - ( somm( aij xj(k) ) ) {j=1..n and j!=i} ]  { i=1..n }   
			//
			//________________________________________________________________________________
			e_num=0; e_den=0;
			
			for(i=0; i<A.n_row; i++)
			{
				x_k_1[i] = b[i];
	
				for(k=A.lr[i]; A.lc[k]<i; k++)
				{
					x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[A.lc[k]] );
				}
						
				for(aii=A.le[k++]; k<A.lr[i+1]; k++)
				{
					x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[A.lc[k]] );
				}
	
				x_k_1[i] = x_k_1[i]/aii;
	
				//__________rapporto per condizione d'uscita
				//
				//      ||x(k+1)-x(k)|| / ||x(k+1)||
				//
				//__________________________________________
				x_k_1_quad = x_k_1[i] * x_k_1[i];
				e_den=e_den + x_k_1_quad;
				e_num=e_num + x_k_1_quad + x_k[i]*x_k[i] - 2*x_k_1[i]*x_k[i];		
			}

			//__________valutazione condizione d'uscita
			//
			//   ||x(k+1)-x(k)|| / ||x(k+1)|| < TOL
			//
			//_________________________________________
			if(e_num/e_den < TOL_QUAD)
			{
				break;
			}
			
						
	
			//_________________________nuova iterazione
			//
			//	       x(k) = x(k+1)
			//
			//_________________________________________
			x_tmp=x_k;
			x_k=x_k_1;
			x_k_1=x_tmp;		
		}
		*/	
	}	
	
	//=================================================coordinatore===============================================//
	if(rank==(size-1))
	{	
		//
		//________________________________creazione struttura CSR per il coordinatore
		//

		CSR_matrix A_coord;
		double* b_coord;
	
		//definizione del job del coordinatore
		//NB: il coordinatore lavora anche ciò che resta
		int job_coord=(A.n_row/size)+(A.n_row%size);
		int j_coord; int i_coord;

		//limiti del job per il coordinatore
		unsigned long int job_worker=A.n_row/size;
		unsigned long int job_min = (size-1)*job_worker;
		unsigned long int job_max = A.n_row;

		//definizione "n_row", "n_col", "n_elem" e "simm" per il coordinatore
		A_coord.n_row=job_coord;
		A_coord.n_col=A.n_col;
		A_coord.n_elem=A.lr[job_max]-A.lr[job_min];
		A_coord.simm=0;

		//definizione del massimo numero di iterazioni e 
		//della tolleranza dell'algoritmo di jacobi
		NMAX=pow(A_coord.n_row,4);
		TOL_QUAD=1e-8;

		//allocazione dei vettori "le", "lc" "lr" e "b" per il coordinatore
		A_coord.le=(double*)malloc(A_coord.n_elem*sizeof(double) );
		A_coord.lc=(unsigned long int*)malloc(A_coord.n_elem*sizeof(unsigned long int));
		A_coord.lr=(unsigned long int*)malloc((A_coord.n_row+1)*sizeof(unsigned long int));
		b_coord=(double*)malloc(A_coord.n_col*sizeof(double));

		//allocazione vettori x(k) e x(k+1)
		x_k_1=(double*)malloc(A_coord.n_col*sizeof(double));
		x_k=(double*)malloc(A_coord.n_col*sizeof(double));
		
		//x(k)=0;
		for(k=0; k<A.n_row; k++){
			x_k[k]=0;
		}
	
		//definizione dei vettori "le", "lr" e "lc" per il coordinatore
		for(i=job_min, j_coord=0, i_coord=0; i<job_max; i++, i_coord++)
		{
			//lr
			A_coord.lr[i_coord]=j_coord;	
			
			//b
			b_coord[i_coord]=b[i];
			
			//le, lc
			for(j=A.lr[i]; j<A.lr[i+1]; j++)
			{		
				A_coord.lc[j_coord]=A.lc[j];
				A_coord.le[j_coord]=A.le[j];
				j_coord++;
			}			
		}
		//[lr]=n+1
		A_coord.lr[i_coord]=j_coord;	

		printf("\nSTRUTTURA CSR COORDINATORE:\n\n");
		stampa_CSR(A_coord.n_row, A_coord.n_col, A_coord.n_elem, A_coord.le, A_coord.lr, A_coord.lc, b_coord);
	}

	MPI_Finalize();/*_____________________________________fine programma distribuito*/	

	exit(0);
}

