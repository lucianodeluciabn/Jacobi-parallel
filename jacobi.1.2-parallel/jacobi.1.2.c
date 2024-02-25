/*****************************************************************************/
/* jacobi.1.2.c parallelo						     */
/*  					     			             */
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
	//==================================================all======================================================//

	CSR_matrix A;//matrice CSR
	int size;//numero processi
	int rank;//id processo
	unsigned long int k, i, j;//contatori
	unsigned long int dim;//dimensione buffer
	int position;//posizione elemento corrente nel buffer
	char* buffer;//buffer 

	double e_num, e_den;//numeratore e denominatore per la condizione d'uscita
	double* x_tmp;//puntatore per swap dei vettori risultato x(k) e x(k+1)
	double x_k_1_quad;//variabile destinata a contenere il valore x(k+1)*x(k+1) durante l'algoritmo
	double aii;//variabile destinata a contenere i valori sulla diagonale della matrice A
	unsigned long int iter;//contatore iterazioni

	double* b;//puntatore al vettore dei termini noti
	double* x_k;//risultato al passo "k": x(k)
	double* x_k_1;//risultato al passo "k+1": x(k+1)

	unsigned long int job_worker;//job identico per ciascun worker (fetta uguale di lavoro)

	unsigned long int NMAX;//numero massimo di iterazioni dell'algoritmo
	double TOL_QUAD=1e-8; //tolleranza

	FILE *f;

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
		f=fopen(argv[1], "r");
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
		A_gen.lr[i]=A_gen.n_elem;//[lr]=n+1
	
	
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
				
		//
		//______________________________creazione e invio delle strutture CSR e delle parti di vettore b, per i worker 
		//
		//

		//job dei worker e del coordinatore
		job_worker=A_gen.n_row/size;//job = numero di righe di A per ciascun worker (fetta del lavoro)
		unsigned long int job_coord=(A_gen.n_row/size)+(A_gen.n_row%size);//il coordinatore lavora anche ciò che resta	
		unsigned long int job_max;//indice massimo fetta di job
		unsigned long int job_min;//indice minimo fetta di job
		int worker;//rank del worker corrente nel ciclo di distribuzione dei job
		unsigned long int j_worker, i_worker, j_coord, i_coord;//indici di scorrimento
		unsigned long int worker_nelem;//elementi di ciascun worker
		
		//vettori le lr e lc per i worker
		double* le;
		unsigned long int* lr;
		unsigned long int* lc;
	
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
			//__________________________________invio valori n_row, n_col, n_elem e simm per il worker corrente
			//
			//allocazione buffer info
			dim=
			3*sizeof(unsigned long int);//n_row, n_col, n_elem
			buffer=(char*)malloc(dim);

			//packaging info struttura CSR
			position=0;
			MPI_Pack(&job_worker, 1, MPI_UNSIGNED_LONG, buffer, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(A_gen.n_col), 1, MPI_UNSIGNED_LONG, buffer, dim, &position, MPI_COMM_WORLD);		
			MPI_Pack(&worker_nelem, 1, MPI_UNSIGNED_LONG, buffer, dim, &position, MPI_COMM_WORLD);
			
			//___________________________________________________________________________<<sincronizzazione>>	
			//invio al worker corrente le info per allocare la sua struttura
			MPI_Send(buffer, dim, MPI_PACKED, worker, 0, MPI_COMM_WORLD);				
			//_______________________________________________________________________________________________
						
			//deallocazione del buffer per le info sulla struttura CSR	
			free(buffer);	
	
			//
			//______________________________________________________vettori le, lr, lc e b per il worker corrente
			//
			//allocazione dei vettori "le", "lc" e "lr" per il worker corrente
			le=(double*)malloc(worker_nelem*sizeof(double) );
			lc=(unsigned long int*)malloc(worker_nelem*sizeof(unsigned long int) );	
			lr=(unsigned long int*)malloc( (job_worker+1)*sizeof(unsigned long int) );

			//definizione dei vettori "le", "lr" e "lc" per il worker corrente
			for(i=job_min, i_worker=0, j_worker=0; i<job_max; i++, i_worker++)
			{	
				//lr
				lr[i_worker]=j_worker;		

				//le, lc
				for(j=A_gen.lr[i]; j<A_gen.lr[i+1]; j++)
				{	
					le[j_worker]=A_gen.le[j];	
					lc[j_worker]=A_gen.lc[j];
					j_worker++;
				}		
			}
			//[lr]=n+1
			lr[i_worker]=j_worker;
									
			//allocazione buffer vettori CSR e parte del vettore b
			dim=
			worker_nelem*sizeof(double)+//le
			(job_worker+1)*sizeof(unsigned long int)+//lr	
			worker_nelem*sizeof(unsigned long int)+//lc
			job_worker*sizeof(double);//parte del vettore b
			buffer=(char*)malloc(dim);
	
			//packaging vettori CSR e vettore b
			position=0;
			MPI_Pack(le, worker_nelem, MPI_DOUBLE, buffer, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(lr, job_worker+1, MPI_UNSIGNED_LONG, buffer, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(lc, worker_nelem, MPI_UNSIGNED_LONG, buffer, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(b_gen[job_min]), job_worker, MPI_DOUBLE, buffer, dim, &position, MPI_COMM_WORLD);
		
			//___________________________________________________________________________<<sincronizzazione>>	
			//invio al worker i vettori le, lr, lc e b
			MPI_Send(buffer, dim, MPI_PACKED, worker, 1, MPI_COMM_WORLD);		
			//_______________________________________________________________________________________________

			//deallocazione buffer vettori CSR
			free(buffer);
			free(le);
			free(lr);
			free(lc);			
		}
			
		//
		//________________________________creazione struttura CSR per il coordinatore
		//
		//limiti del job per il coordinatore
		job_min = (size-1)*job_worker;
		job_max = A_gen.n_row;

		//definizione "n_row", "n_col" e "n_elem" per il coordinatore
		A.n_row=job_coord;
		A.n_col=A_gen.n_col;
		A.n_elem=A_gen.lr[job_max]-A_gen.lr[job_min];

		//allocazione dei vettori "le", "lc" "lr" e "b" per il coordinatore
		A.le=(double*)malloc(A.n_elem*sizeof(double) );
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));
		b=(double*)malloc(A.n_row*sizeof(double));

		//allocazione vettori x(k) e x(k+1)
		x_k_1=(double*)malloc(A.n_col*sizeof(double));
		x_k=(double*)malloc(A.n_col*sizeof(double));	
	
		//x(0)=0;
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
	}
	
	//==================================================workers===================================================//
	if(rank!=(size-1)){
		
		//
		//__________________________ricezione struttura CSR 
		//
		
		//allocazione buffer per ricevere n_row, n_col e n_elem
		dim=
		3*sizeof(unsigned long int);//n_row, n_col, n_elem
		buffer=(char*)malloc(dim);	

		//__________________________________________________________________________________<<sincronizzazione>>
		//ricezione n_row, n_col e n_elem
		MPI_Recv(buffer, dim, MPI_PACKED, size-1, 0, MPI_COMM_WORLD, &status);
		//______________________________________________________________________________________________________
		
		//definizione di n_row, n_col e n_elem
		position=0;
		MPI_Unpack(buffer, dim, &position, &(A.n_row), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer, dim, &position, &(A.n_col), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer, dim, &position, &(A.n_elem), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

		//definizione della fetta di lavoro uguale per tutti i worker.
		//serve per l'esecuzione dell'algoritmo di jacobi
		job_worker=A.n_row;	

		//allocazione vettori le, lr, lc e b
		A.le=(double*)malloc(A.n_elem*sizeof(double));
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));
		b=(double*)malloc(A.n_row*sizeof(double));

		//deallocazione buffer
		free(buffer);

		//allocazione buffer per la ricezione dei vettori le, lr, lc e b
		dim=
		(A.n_elem)*sizeof(double)+//le
		(A.n_row+1)*sizeof(unsigned long int)+//lr	
		(A.n_elem)*sizeof(unsigned long int)+//lc
		(A.n_row)*sizeof(double);//parte del vettore b
		buffer=(char*)malloc(dim);

		//allocazione vettori x(k) e x(k+1)
		x_k_1=(double*)malloc(A.n_col*sizeof(double));
		x_k=(double*)malloc(A.n_col*sizeof(double));
		
		//x(k)=0;
		for(k=0; k<A.n_col; k++){
			x_k[k]=0;
		}

		//ricezione vettori CSR e vettore b
		position=0;
		MPI_Recv(buffer, dim, MPI_PACKED, size-1, 1, MPI_COMM_WORLD, &status);
		MPI_Unpack(buffer, dim, &position, A.le, A.n_elem, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(buffer, dim, &position, A.lr, A.n_row+1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer, dim, &position, A.lc, A.n_elem, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer, dim, &position, b, A.n_row, MPI_DOUBLE, MPI_COMM_WORLD);
	
		//deallocazione buffer
		free(buffer);
			
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
	NMAX=pow(A.n_col, 4);
	for(iter=0; iter<NMAX; iter++)
	{
		//_______________________________________________________________calcolo di x(k+1)
		//
		//   xi(k+1) = (1/aij)*[ bi - ( somm( aij xj(k) ) ) {j=1..n and j!=i} ]  { i=1..n }   
		//		
		for(i=0; i<A.n_row; i++)
		{
			x_k_1[i] = b[i];

			for(k=A.lr[i]; A.lc[k]<i+(job_worker*rank); k++)//<------aggiunto: "+(job_worker)*rank"
									//NB: l'elemento aii si trova in posizioni
									//    diverse a seconda del worker (essendo
									//    esso sulla diagonale della matrice A)
			{
				x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[A.lc[k]] );
			}
					
			for(aii=A.le[k++]; k<A.lr[i+1]; k++)
			{
				x_k_1[i] = x_k_1[i] - ( A.le[k]*x_k[A.lc[k]] );
			}

			x_k_1[i] = x_k_1[i]/aii;
		
		}
		//________________________________________________________________________________


		//________________________________________________________________<<sincronizzazione>> 
		//tutti i processi aggiornano la x(k+1) scrivendola in x(k)
		//per la prossima iterazione

		MPI_Allgatherv(x_k_1, A.n_row, MPI_DOUBLE, x_k_1, recvcounts, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
		//____________________________________________________________________________________


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

		
		
		//____________________nuova iterazione: x(k)=x(k+1)
		//
		x_tmp=x_k;
		x_k=x_k_1;
		x_k_1=x_tmp;
		//_________________________________________________	

	
	}
		
	//=================================================coordinatore===============================================//

	if(rank==(size-1))
	{				
		//________________stampa della soluzione	
		//				
		//______________________________________

		f=fopen(strcat(argv[1],"_out_parallel.dat"), "w");
	
		for(i=0; i<A.n_col; i++)
		{
			fprintf(f, "%le\n", x_k_1[i]);
		}
		fclose(f);
	}

	MPI_Finalize();/*_____________________________________fine programma distribuito*/	

	exit(0);
}

