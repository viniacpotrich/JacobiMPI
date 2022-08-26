#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

/*----------------------------------------------------------*/
#define ERRO 1e-10

double calcErro(double *a, double *b, int tamA){
    double somaA = 0, somaB = 0, resultA = 0, resultB = 0;
    for(int i = 0; i < tamA; i++){
        somaA += pow(a[i], 2);
        somaB += pow(b[i], 2);
    }
    resultA = sqrt(somaA);
    resultB = sqrt(somaB);
    return resultA - resultB;
}

void escreveVetorResultado(double *b, int tamA);
void escreveVetorMatriz(double *b, int tamLinha, int tamColuna);

/*----------------------------------------------------------*/
void jacobi(double *A, double *b, double *x, int tamA, int nIteracoes, int nProcessos, int id, int *tamanhoPosicaoVetor, int *posicaoVetor, int *tamanhoPosicaoMatriz, int *posicaoMatriz){
    int i, j, marcador=0;
    int tamanhoLinha = tamA;
    int inicioBloco = (tamanhoLinha/nProcessos)*(id);
    int diagonalPrincipal;
    int fimBloco = tamanhoPosicaoVetor[id];
    int countIteracoes = 0;
    double soma;
    double *xAtual, *xPart, *xNovo;
    xAtual = (double *)malloc(tamanhoLinha * sizeof(double));
    xPart = (double *)malloc(tamanhoPosicaoVetor[id] * sizeof(double));
    xNovo = (double *)malloc(tamanhoLinha * sizeof(double));
    for (int i = 0; i<tamanhoLinha; i++){
        xAtual[i] = 0;
        xNovo[i] = 0;
    }
    for (int i = 0; i<tamanhoPosicaoVetor[id]; i++){
        xPart[i] = 0;
    }
    while (countIteracoes < nIteracoes){
        diagonalPrincipal =0 ;
        for (i = 0; i < fimBloco; i++){
            diagonalPrincipal = posicaoVetor[id] + (tamanhoLinha * i) + i;
            soma = 0;
            
            for (j = 0; j < tamanhoLinha; j++){
                if (diagonalPrincipal != (tamanhoLinha * i)+j){
                    soma += A[i*tamanhoLinha+j] * xAtual[j];
                }
            }
            xPart[i] = (b[i] - soma)/A[diagonalPrincipal];
        }

	    MPI_Allgatherv(xPart, tamanhoPosicaoVetor[id], MPI_DOUBLE, xNovo, tamanhoPosicaoVetor, posicaoVetor, MPI_DOUBLE,MPI_COMM_WORLD);
        
        if ( fabs( calcErro(xNovo, xAtual, tamanhoLinha) ) < ERRO){
            for (i = 0; i < tamanhoLinha; i++){
                x[i] = xNovo[i];
            }
            marcador = 1;
        } else {
            for (i = 0; i < tamanhoLinha; i++){
                x[i] = xAtual[i];
                xAtual[i] = xNovo[i];
            }
        }
        countIteracoes++;
        if (marcador == 1){
            break;
        }
    }
    printf("iteracoes: %d\n", countIteracoes);
}


/*----------------------------------------------------------*/
void geraVetorMatriz(double *A, int tamA){
    for (int i = 0; i<tamA; i++){ // linha
        for (int j = 0; j<tamA; j++){ // coluna
            if ((i == 0 && j == 0) || (i == tamA-1 && j == tamA-1)){
                A[i*tamA+j] = 6;
            } else if (i == j){
                A[i*tamA+j] = 4;
            } else if (((j == i-1) || (j == i+1)) && i != 0 && i != tamA-1){
                A[i*tamA+j] = 1;
            } else {
                A[i*tamA+j] = 0;
            }
        }
    }
    return;
}

void escreveVetorMatriz(double *A, int tamLinha, int tamColuna){
    for (int i = 0; i<tamColuna; i++){ // linha
        for (int j = 0; j<tamLinha; j++){ // coluna
            printf("%lf ", A[i*tamLinha+j]);
        }
        printf("\n");
    }
    return;
}

void geraVetorResultado(double *b, int tamA){
    for (int i = 0; i < tamA; i++){
        //printf("i = %d\n",i);
        if ((i == 0) || (i == tamA-1)){
            b[i] = 0;
        } else if ((i == 1) || (i == tamA-2)){
            b[i] = 1;
        } else if((i == 2) || (i == tamA-3)){
            b[i] = 2;
        } else{
            b[i] = -6;
        }
    }
}

void escreveVetorResultado(double *b, int tamA){
    for (int i = 0; i < tamA; i++){
        printf("%lf ",b[i]);
    }
    printf("\n");
    return;
}


void zeraX(double *x, int tamA){
    for (int i = 0; i < tamA; i++){
        x[i]= 0;
    }
}



/*----------------------------------------------------------*/
int main(int argc, char **argv ){
    int tamanhoMatriz,tamanhoRecorte,nIteracoes,np,id;
    double *A = NULL;
    double *PartA = NULL;
    double *PartB = NULL;
    double *b = NULL;
    double *x = NULL;
    double ti,tf=0;
    int *posicaoMatriz = NULL;
    int *posicaoVetor = NULL;
    int *tamanhoPosicaoMatriz = NULL;
    int *tamanhoPosicaoVetor = NULL;


    //parametros do programa ------------------------------------------------------/
	if ( argc != 3 ){
		printf("%s < Ordem da Matriz > < Max Iteracoes >\n", argv[0]);
		exit(0);
	}
    tamanhoMatriz = atoi(argv[1]);
	nIteracoes = atoi(argv[2]);

    
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);


    tamanhoPosicaoMatriz = (int *)malloc(np*sizeof(int));
    tamanhoPosicaoVetor = (int *)malloc(np*sizeof(int));
    posicaoMatriz = (int *)malloc(np*sizeof(int));
    posicaoVetor = (int *)malloc(np*sizeof(int));


    int posicaoLivre = 0;
	for(int i=0; i<np; i++){
		if ( i == 0 ){
			tamanhoPosicaoMatriz[i] = (tamanhoMatriz/np + tamanhoMatriz%np)*tamanhoMatriz;      
            tamanhoPosicaoVetor[i] = tamanhoMatriz/np + tamanhoMatriz%np;
			posicaoMatriz[i] = posicaoLivre * tamanhoMatriz;   
            posicaoVetor[i] = posicaoLivre;                               
			posicaoLivre += tamanhoMatriz/np + tamanhoMatriz%np;            
		}
		else{
            tamanhoPosicaoMatriz[i] = (tamanhoMatriz/np)*tamanhoMatriz;
			tamanhoPosicaoVetor[i] = tamanhoMatriz/np;
            posicaoMatriz[i] = posicaoLivre * tamanhoMatriz;   
			posicaoVetor[i] = posicaoLivre;
			posicaoLivre += tamanhoMatriz/np;
		}
	}

    if ( id == 0 ){
		tamanhoRecorte = tamanhoMatriz/np + tamanhoMatriz%np;
	}
	else{
		tamanhoRecorte = tamanhoMatriz/np;
	}

    if(id==0){
        A = (double *)malloc(tamanhoMatriz * tamanhoMatriz * sizeof(double));
        b = (double *)malloc(tamanhoMatriz * sizeof(double));
        geraVetorMatriz(A,tamanhoMatriz);
        geraVetorResultado(b,tamanhoMatriz);
    }
    x = (double *)malloc(tamanhoMatriz * sizeof(double));
    PartA = (double *)malloc((tamanhoMatriz*tamanhoRecorte)*sizeof(double));
    PartB = (double *)malloc(tamanhoRecorte*sizeof(double));
    zeraX(x,tamanhoMatriz);
    
    MPI_Scatterv(A,tamanhoPosicaoMatriz,posicaoMatriz,MPI_DOUBLE,PartA,(tamanhoMatriz*tamanhoRecorte),MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatterv(b,tamanhoPosicaoVetor,posicaoVetor,MPI_DOUBLE,PartB,tamanhoRecorte,MPI_DOUBLE,0,MPI_COMM_WORLD); 

    ti = MPI_Wtime();
    jacobi(PartA, PartB, x, tamanhoMatriz, nIteracoes, np, id, tamanhoPosicaoVetor, posicaoVetor, tamanhoPosicaoMatriz, posicaoMatriz);
    tf = MPI_Wtime();
    
    if(id==0){
        printf("resultado, id= %d\n", id);
        for (int i = 0; i<tamanhoMatriz; i++){
            printf("%f\n", x[i]);
        }
		printf("Tempo: %f\n", tf - ti);
        free(A);
        free(b);
    }

    MPI_Finalize();
    
    free(x);
    free(PartA);
    free(PartB);
    free(tamanhoPosicaoMatriz);
    free(tamanhoPosicaoVetor);
    free(posicaoVetor);
}
/*----------------------------------------------------------*/

