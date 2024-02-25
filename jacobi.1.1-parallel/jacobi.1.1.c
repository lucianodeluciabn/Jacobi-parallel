/*****************************************************************************/
/* jacobi.1.1.c parallelo						     */
/*    	
	prima versione funzionante	     			             */
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

	unsigned long int NMAX=2000000;//numero massimo di iterazioni dell'algoritmo
	double TOL_QUAD=1e-8; //tolleranza

	
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
		
		CSR_matrix A_gen;
		double* b_gen;

		//apertura file e raccolta info matrice A del sistema lineare
		FILE *f=fopen(argv[1], "r");
		fscanf(f,"%lu,%lu,%lu,%d\n", &(A_gen.n_row), &(A_gen.n_col), &(A_gen.n_elem), &(A_gen.simm));
		NMAX=pow(A_gen.n_row,4);	

		//le = vettore degli elementi non nulli di A[N][N] (ordinati per riga)
		A_gen.le=(double*)malloc(A_gen.n_elem*sizeof(double));

		//lr = vettore degli indici, in le, dei primi elementi, non nulli, di riga, della matrice A[N][N]
		A_gen.lr=(unsigned long int*)malloc((A_gen.n_row+1)*sizeof(unsigned long int));

		//lc = vettore degli indici colonna, in A[N][N], degli elementi di "le"   
		A_gen.lc=(unsigned long int*)malloc(A_gen.n_elem*sizeof(unsigned long int));			

		//definizione dei vettori CSR
		unsigned long int row;
		i=0; 
		for(k=0; k<A_gen.n_elem; k++)
		{
			fscanf(f, "%lu,%lu,%le\n", &row, &(A_gen.lc[k]), &(A_gen.le[k]));
			if(row == i)
			{
				A_gen.lr[row]=k;
				i++;
			}
		}
		A_gen.lr[i]=A_gen.n_elem;
		
		//
		//___________________allocazione e definizione vettore b
		//

		//allocazione vettore "b" del sistema lineare
		b_gen=(double*)malloc(A_gen.n_col*sizeof(double));

		//definizione del vettore b
		for(k=0; k<A_gen.n_row; k++)
		{
			fscanf(f, "%le\n", &b_gen[k]);
		}	
		fclose(f);

		printf("STRUTTURA CSR GENERALE:\n\n");
		stampa_CSR(A_gen.n_row, A_gen.n_col, A_gen.n_elem, A_gen.le, A_gen.lr, A_gen.lc, b_gen);

				
		//
		//_________________________creazione e invio delle strutture CSR e delle parti di vettore b, per i worker 
		//
		//

		//definizione del job dei worker
		unsigned long int job_worker=A_gen.n_row/size;	
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
			worker_nelem=A_gen.lr[job_max]-A_gen.lr[job_min];			
								
			//
			//__________________________________definizione e invio valori n_row, n_col, n_elem e simm per il worker corrente
			//
		
			//definizione "n_row", "n_col", "n_elem" e "simm" per il worker corrente
			v_nrow[worker]=job_worker;
			v_ncol[worker]=A_gen.n_col;
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
			
			//invio al worker corrente le info per allocare la sua struttura
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
				for(j=A_gen.lr[i]; j<A_gen.lr[i+1]; j++)
				{		
					v_lc[worker][j_worker]=A_gen.lc[j];
					v_le[worker][j_worker]=A_gen.le[j];
					j_worker++;
				}		
			}
			//[lr]=n+1
			v_lr[worker][i_worker]=j_worker;
									
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
			MPI_Pack(&(b_gen[job_min]), job_worker, MPI_DOUBLE, buffer2, dim, &position, MPI_COMM_WORLD);
			
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
	
		//
		//________________________________creazione struttura CSR per il coordinatore
		//
	
		//definizione del job del coordinatore
		//NB: il coordinatore lavora anche ciò che resta
		int job_coord=(A_gen.n_row/size)+(A_gen.n_row%size);
		int j_coord; int i_coord;

		//limiti del job per il coordinatore
		job_min = (size-1)*job_worker;
		job_max = A_gen.n_row;

		//definizione "n_row", "n_col", "n_elem" e "simm" per il coordinatore
		A.n_row=job_coord;
		A.n_col=A_gen.n_col;
		A.n_elem=A_gen.lr[job_max]-A_gen.lr[job_min];
		A.simm=A_gen.simm;

		//allocazione dei vettori "le", "lc" "lr" e "b" per il coordinatore
		A.le=(double*)malloc(A.n_elem*sizeof(double) );
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));
		b=(double*)malloc(A.n_row*sizeof(double));

		//allocazione vettori x(k) e x(k+1)
		x_k_1=(double*)malloc(A.n_col*sizeof(double));
		x_k=(double*)malloc(A.n_col*sizeof(double));	
	
		//x(k)=0;
		for(k=0; k<A.n_col; k++){
			x_k[k]=0;
		}
	
		//definizione dei vettori "le", "lr" e "lc" per il coordinatore
		for(i=job_min, j_coord=0, i_coord=0; i<job_max; i++, i_coord++)
		{
			//lr
			A.lr[i_coord]=j_coord;	
			
			//b
			b[i_coord]=b_gen[i];
			
			//le, lc
			for(j=A_gen.lr[i]; j<A_gen.lr[i+1]; j++)
			{		
				A.lc[j_coord]=A_gen.lc[j];
				A.le[j_coord]=A_gen.le[j];
				j_coord++;
			}			
		}
		//[lr]=n+1
		A.lr[i_coord]=j_coord;	

		printf("\nSTRUTTURA CSR COORDINATORE:\n\n");
		stampa_CSR(A.n_row, A.n_col, A.n_elem, A.le, A.lr, A.lc, b);

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

		//allocazione vettori "le", "lr" "lc"  e "b" del worker
		A.le=(double*)malloc(A.n_elem*sizeof(double));
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));
		b=(double*)malloc(A.n_row*sizeof(double));

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
		x_k_1=(double*)malloc(A.n_col*sizeof(double));
		x_k=(double*)malloc(A.n_col*sizeof(double));
		
		//x(k)=0;
		for(k=0; k<A.n_col; k++){
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

		printf("STRUTTURA CSR WORKER %d:\n\n", rank);
		stampa_CSR(A.n_row, A.n_col, A.n_elem, A.le, A.lr, A.lc, b);
			
	}	
	

	//==================================================all=======================================================//

	
	//
	//_______________________________________________________________________________allocazione vettori per Allgaterv
	//
	//per far si che i processi, ad ogni iterazione aggiornino il valore della
	//soluzione correntemente calcolata, bisogna definire le posizioni in cui
	//tali vettori andranno ad esser scritti nei buffer di ricezione di ogni
	//worker
	int* recvcounts=(int*)malloc(size*sizeof(int));//vettore delle dimensioni di ogni x(k+1) di ogni worker
	int* dspls=(int*)malloc(size*sizeof(int));//vettore posizionamenti degli x(k+1) dei worker nel buffer di ricezione
	int dspl=0;
	//il coordinatore determina le posizioni
	//in modo differente rispetto ai worker
	if(rank==(size-1))
	{
		for(k=0; k<(size-1); k++)
		{
			recvcounts[k]=A.n_col/size;
			dspls[k]=dspl;
			dspl=dspl+recvcounts[k];
		}
		recvcounts[size-1]=A.n_row;
		dspls[size-1]=dspl;
	}
	else
	{
		for(k=0; k<(size-1); k++)
		{
			recvcounts[k]=A.n_row;
			dspls[k]=dspl;
			dspl=dspl+recvcounts[k];
		}
		recvcounts[size-1]=(A.n_row)+(A.n_col%size);
		dspls[size-1]=dspl;
	}
	//____________________________________________________________________________________________________________________	

	
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
		
		for(i=0; i<A.n_row; i++)
		{
			x_k_1[i] = b[i];

			for(k=A.lr[i]; A.lc[k]<i+((A.n_col/size)*rank); k++)//<------aggiunto: +(A.ncol/size)*rank
			{
				x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[A.lc[k]] );
			}
					
			for(aii=A.le[k++]; k<A.lr[i+1]; k++)
			{
				x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[A.lc[k]] );
			}

			x_k_1[i] = x_k_1[i]/aii;
		
		}

		//________________________________________________________________<<sincronizzazione>> 
		//tutti i processi aggiornano la x(k+1) scrivendola in x(k)
		//per la prossima iterazione
		if(iter==0){
			printf("%d, x(k+1): ", rank);
			for(i=0;i<A.n_row;i++)
			{
				printf("%le ", x_k_1[i]);
			}
			printf("\n");
		}

		MPI_Allgatherv(x_k_1, A.n_row, MPI_DOUBLE, x_k_1, recvcounts, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
		if(iter==0){
			printf("%d, x(k+1): ", rank);
			for(i=0;i<A.n_col;i++)
			{
				printf("%le ", x_k_1[i]);
			}
			printf("\n");
		}
		//____________________________________________________________________________________
		/*
		printf("%d, x(k+1): ", rank);
		for(i=0;i<A.n_col;i++)
		{
			printf("%le ", x_k_1[i]);
		}
		printf("\n");
		*/
		//___________________________________calcolo rapporto per condizione d'uscita
		//
		//      ||x(k+1)-x(k)|| / ||x(k+1)||
		//
		e_num=0; e_den=0;
		for(i=0; i<A.n_col; i++)
		{
			x_k_1_quad = x_k_1[i] * x_k_1[i];
			e_den=e_den + x_k_1_quad;
			e_num=e_num + x_k_1_quad + x_k[i]*x_k[i] - 2*x_k_1[i]*x_k[i];
		}
		//____________________________________________________________________________


		//__________valutazione condizione d'uscita
		//
		//   ||x(k+1)-x(k)|| / ||x(k+1)|| < TOL
		//
		if(e_num/e_den < TOL_QUAD)
		{
			break;
		}				
		//_________________________________________

		
		//
		//____________________nuova iterazione: x(k)=x(k+1)
		//
		x_tmp=x_k;
		x_k=x_k_1;
		x_k_1=x_tmp;	

	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	//=================================================coordinatore===============================================//

	if(rank==(size-1))
	{
					
		//________________stampa della soluzione	
		//				
		//______________________________________

		//f=fopen(strcat(argv[1],"_out.dat"), "w");
	
		for(i=0; i<A.n_col; i++)
		{
			/*f*/printf(/*f,*/"%le\n", x_k_1[i]);
		}
		//fclose(f);
		//printf("\niterazioni:%ld\n", iter+1);
		//printf("tolleranza: %le\n", e_num/e_den);
	}

	MPI_Finalize();/*_____________________________________fine programma distribuito*/	

	exit(0);
}

