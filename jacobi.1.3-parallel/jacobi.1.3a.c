/*****************************************************************************/
/* jacobi.1.3a.c-parallelo-con send asincrona			             */
/*	versione con send asincrona rispetto alla precedente
	NB: con la send asincrona non si può fare free sui buffer            
	altrimenti il supporto alla messaggistica causa una seg. fault.
        in quanto legge un buffer liberato dal sender che ha fatto free      */
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

int main(int argc, char *argv[])
{
	//==================================================all======================================================//

	CSR_matrix A;//matrice CSR
	int size;//numero processi
	int rank;//id processo
	unsigned long int k, i, j;//contatori
	unsigned long int dim;//dimensione buffer
	int position;//posizione elemento corrente nel buffer
	char *buffer1, *buffer2;//buffer 

	double e_num, e_den;//numeratore e denominatore per la condizione d'uscita
	double* x_tmp;//puntatore per swap dei vettori risultato x(k) e x(k+1)
	double x_k_1_quad;//variabile destinata a contenere il valore x(k+1)*x(k+1) durante l'algoritmo
	double aii;//variabile destinata a contenere i valori sulla diagonale della matrice A
	unsigned long int iter;//contatore iterazioni

	double* b;//puntatore al vettore dei termini noti
	double* x_k;//risultato al passo "k": x(k)
	double* x_k_1;//risultato al passo "k+1": x(k+1)

	unsigned long int job_worker;//job identico per ciascun worker (fetta uguale di lavoro)
	unsigned long int job_coord;//job del coordinatore. pari al giob di un worker più ciò che resta

	unsigned long int NMAX;//numero massimo di iterazioni dell'algoritmo
	double TOL_QUAD=1e-8; //tolleranza

	int* recvcounts;//vettore delle dimensioni di ogni x(k+1) di ogni worker
	int* dspls;//vettore posizionamenti degli x(k+1) dei worker nel buffer di ricezione
	int dspl;//spiazzamento corrente


	FILE *f;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Request req;
	
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
		job_coord=(A_gen.n_row/size)+(A_gen.n_row%size);//il coordinatore lavora anche ciò che resta	
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
			buffer1=(char*)malloc(dim);

			//packaging info struttura CSR
			position=0;
			MPI_Pack(&job_worker, 1, MPI_UNSIGNED_LONG, buffer1, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(A_gen.n_col), 1, MPI_UNSIGNED_LONG, buffer1, dim, &position, MPI_COMM_WORLD);		
			MPI_Pack(&worker_nelem, 1, MPI_UNSIGNED_LONG, buffer1, dim, &position, MPI_COMM_WORLD);
			
			//___________________________________________________________________________<<sincronizzazione>>	
			//invio al worker corrente le info per allocare la sua struttura
			MPI_Isend(buffer1, dim, MPI_PACKED, worker, 0, MPI_COMM_WORLD, &req);				
			//_______________________________________________________________________________________________
			

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
			buffer2=(char*)malloc(dim);
	
			//packaging vettori CSR e vettore b
			position=0;
			MPI_Pack(le, worker_nelem, MPI_DOUBLE, buffer2, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(lr, job_worker+1, MPI_UNSIGNED_LONG, buffer2, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(lc, worker_nelem, MPI_UNSIGNED_LONG, buffer2, dim, &position, MPI_COMM_WORLD);
			MPI_Pack(&(b_gen[job_min]), job_worker, MPI_DOUBLE, buffer2, dim, &position, MPI_COMM_WORLD);
		
			//___________________________________________________________________________<<sincronizzazione>>	
			//invio al worker i vettori le, lr, lc e b
			MPI_Isend(buffer2, dim, MPI_PACKED, worker, 1, MPI_COMM_WORLD, &req);		
			//_______________________________________________________________________________________________

			//deallocazione vettori CSR
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
		buffer1=(char*)malloc(dim);	

		//__________________________________________________________________________________<<sincronizzazione>>
		//ricezione n_row, n_col e n_elem
		MPI_Recv(buffer1, dim, MPI_PACKED, size-1, 0, MPI_COMM_WORLD, &status);
		//______________________________________________________________________________________________________
		
		//definizione di n_row, n_col e n_elem
		position=0;
		MPI_Unpack(buffer1, dim, &position, &(A.n_row), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer1, dim, &position, &(A.n_col), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer1, dim, &position, &(A.n_elem), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

		//definizione della fetta di lavoro del worker e del coordinatore
		job_worker=A.n_row;
		job_coord=(A.n_row)+(A.n_col%size);	

		//allocazione vettori le, lr, lc e b
		A.le=(double*)malloc(A.n_elem*sizeof(double));
		A.lr=(unsigned long int*)malloc((A.n_row+1)*sizeof(unsigned long int));
		A.lc=(unsigned long int*)malloc(A.n_elem*sizeof(unsigned long int));
		b=(double*)malloc(A.n_row*sizeof(double));

		//deallocazione buffer
		free(buffer1);

		//allocazione buffer per la ricezione dei vettori le, lr, lc e b
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
	
		//__________________________________________________________________________________<<sincronizzazione>>
		//ricezione di le lr lc e b
		MPI_Recv(buffer2, dim, MPI_PACKED, size-1, 1, MPI_COMM_WORLD, &status);
		//______________________________________________________________________________________________________
		
		//definizione di le lr lc e b
		position=0;
		MPI_Unpack(buffer2, dim, &position, A.le, A.n_elem, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(buffer2, dim, &position, A.lr, A.n_row+1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer2, dim, &position, A.lc, A.n_elem, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		MPI_Unpack(buffer2, dim, &position, b, A.n_row, MPI_DOUBLE, MPI_COMM_WORLD);
	
		//deallocazione buffer
		free(buffer2);
			
	}	
	

	//==================================================all=======================================================//

	
	//____________________________________________allocazione e definizione vettori recvcount e dspls
	//
	//allocazione recvcount e dspls 
	recvcounts=(int*)malloc(size*sizeof(int));
	dspls=(int*)malloc(size*sizeof(int));
	dspl=0;

	//definizione recvcount e dspls
	for(k=0; k<(size-1); k++)
	{
		recvcounts[k]=job_worker;
		dspls[k]=dspl;
		dspl=dspl+recvcounts[k];
	}
	recvcounts[size-1]=job_coord;
	dspls[size-1]=dspl;
	//_________________________________________________________________________________________________
	

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
		f=fopen(strcat(argv[1],"_out_parallel.dat"), "w");
	
		for(i=0; i<A.n_col; i++)
		{
			fprintf(f, "%le\n", x_k_1[i]);
		}
		fclose(f);
		//______________________________________
	}

	MPI_Finalize();	

	exit(0);
}

