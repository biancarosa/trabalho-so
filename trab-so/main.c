#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<semaphore.h>
#include<string.h>

#define nLerCatalogo         1
#define nLerArquivos         3
#define nMultiplicarMatrizes 4
#define nDeterminanteMatriz  5
#define nEscreverArquivo     3
#define OUT_FILE "saida.out"



int debug, contLeitura, contEscrita;

typedef double **mat;

typedef struct {
    char Nome[100];
    int Ordem;
    mat A, B, C;
    double E;
} S;

typedef struct{
    S *buffer[5];
    sem_t full, empty, mutex;
    int bufferIn, bufferOut;
} sharedStruct;

sharedStruct shared[4];

void mat_show(mat matrix, int size) {
    int i, j;
    printf("\n\n");
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            printf(" %f", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void mat_print_file(mat matrix, int size, FILE * f) {
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            fprintf(f," %f", matrix[i][j]);
        }
        fprintf(f, "\n");
    }
}


void wolfranFormat(mat matrix, int size) {
    int i, j;
    printf("{");
    for (i = 0; i < size; i++) {
        printf("{");
        for (j = 0; j < size; j++) {
            printf("%d", (int) matrix[i][j]);
            if (j < size - 1) {
                printf(",");
            }

        }
        printf("}");
        if (i < size - 1) {
            printf(",");
        }
    }
    printf("}");
}

mat newMatrix(int size) {
    int i;
    mat matrix = malloc(size * sizeof ( double *));
    for (i = 0; i < size; i++) {
        matrix[i] = malloc(size * sizeof ( double));
    }
    return matrix;
}

mat cloneMatrix( mat matrix, int size ){
    mat clone = newMatrix( size );
    int i, j;
    for( i=0; i<size; i++ ){
        for( j=0; j<size; j++ ){
            clone[i][j] = matrix[i][j];
        }
    }
    return clone;
}

void zeroMatrix(mat matrix, int size) {
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            matrix[i][j] = 0;
        }
    }
}

void randomFill(mat matrix, int size) {
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            matrix[i][j] = (double) (rand() % 10);
        }
    }
}


int multiplyMatrices(mat matA, mat matB, mat resolt, int size){
    int i,j,k;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(2);
    #pragma omp parallel shared(matA,matB,resolt) private(i,j,k) 
    {
        #pragma omp for  schedule(static)
        for (i=0; i<size; i=i+1){
            for (j=0; j<size; j=j+1){
                resolt[i][j]=0.;
                for (k=0; k<size; k=k+1){
                    resolt[i][j] += (matA[i][k])*(matB[k][j]);
                }
            }
        }
    }
    return 0;
}

void mat_pivot(mat a, mat p, int size) {
    int i, j, k;
    double aux;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            p[i][j] = (i == j);
        }
    }
    for (i = 0; i < size; i++) {
        int max_j = i;
        for (j = i; j < size; j++) {
            if (fabs(a[j][i]) > fabs(a[max_j][i])) max_j = j;
        }
        if (max_j != i)
            for (k = 0; k < size; k++) {
                aux = p[i][k];
                p[i][k] = p[max_j][k];
                p[max_j][k] = aux;
            }
    }
}

void mat_LU(mat matrix, mat matL, mat matU, mat matP, int size) {
    int i, j, k;
    mat Aprime = newMatrix(size);
    double s;
    zeroMatrix(matL, size);
    zeroMatrix(matU, size);
    mat_pivot(matrix, matP, size);
    multiplyMatrices(matP, matrix, Aprime, size);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (j <= i) {
                s = 0;
                for (k = 0; k < size; k++) {
                    s += matL[j][k] * matU[k][i];
                }
                matU[j][i] = Aprime[j][i] - s;
            }
            if (j >= i) {
                s = 0;
                for (k = 0; k < size; k++) {
                    s += matL[j][k] * matU[k][i];
                }
                matL[j][i] = (double) ((Aprime[j][i] - s) / matU[i][i]);
            }
        }
    }
    for (i = 0; i < size; i++) {
        matL[i][i] = 1;
    }

}

double triangularMatrixDeterminant(mat matrix, int size) {
    int i;
    double det = 1;
    for (i = 0; i < size; i++) {
        det *= matrix[i][i];
    }
    return det;
}

double determinant(mat matrix, int size) {
    double det;
    int i, pivotCount = 0;
    mat matL = newMatrix(size), matU = newMatrix(size), matP = newMatrix(size);
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(2);
    #pragma omp parallel shared(matrix, matL, matU, matP, size) private(i)
    {
        mat_LU(matrix, matL, matU, matP, size);
        det = triangularMatrixDeterminant(matU, size);
    
        for (i = 0; i < size; i++) {
            if (matP[i][i] == 0) {
                pivotCount++;
            }
        }
        if (pivotCount % 2 == 0 && pivotCount != 0) {
            det *= -1;
        }
    }
    free(matL);
    free(matU);
    free(matP);
    return det;
}

/*LC - Thread Leitora de Catálogo - Lê um arquivo (entrada.in) contendo uma lista de arquivos de entrada (um nome de arquivo por linha). A cada linha lida, a thread LC  cria dinamicamente uma estrutura S, preenche o nome do arquivo de entrada e coloca o ponteiro para a estrutura S em shared[0].buffer[in] para serem processada pela etapa seguinte. Só teremos 1 instância desta thread.*/

void lerCatalogo(){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int i;
    fp = fopen("entrada.in", "r");
    if (fp == NULL){
        printf("\n\nErro na leitura do catalogo, entrada.in nao encontrado.\n\n");
        exit(EXIT_FAILURE);        
    }
    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
        if (read > 100 ){
            printf("\n\nErro na leitura do catalogo, maximo de 100 caracteres no nome do arquivo.\n\n");
            exit(EXIT_FAILURE);        
        }
        contLeitura++;
        if (read > 1 ){
            sem_wait(&shared[0].empty);
            sem_wait(&shared[0].mutex);
            
            

            strtok( line, "\n" );
            for( i=0; i<read; i++ ){
                shared[0].buffer[shared[0].bufferIn]->Nome[i] = line[i];
            }
            
            

            //printf("%s", shared[0].buffer[shared[0].bufferIn]->Nome);
           // sem_getvalue( &shared[0].empty, &debug);
           // printf("\n\nempty:%d\n\n", debug);
            shared[0].bufferIn = (shared[0].bufferIn+1)%5;
            
            sem_post(&shared[0].full);
            sem_post(&shared[0].mutex);
        }
    }
    fclose(fp);
}

/*LA - Thread Leitora de Arquivos - Recebe a estrutura preenchida por LC (recebe de shared[0].buffer[out]) e abre o arquivo indicado pela variável nome e lê o arquivo que contém um par matriz quadradas de doubles com ordem variando entre 500 e 1000. As duas matizes em um arquivo possuem a mesma ordem, que é indicada na primeira linha do arquivo. Elementos de uma linha são separados por espaços e as linhas são separadas por quebras de linha. A arquivo lido, a thread LA completa a estrutura recebida de LC, preenchendo as matrizes A e B e a variável ordem. Depois, a thread LA coloca o ponteiro para a estrutura S em shared[1].buffer[in] para serem processada pela etapa seguinte. Teremos 3 instâncias desta thread.*/

void lerArquivos(){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int index = 0, i, j, ordem;
    S *temp = ( S* )malloc( sizeof( S ) );
    while( 1 ){
        line = NULL;
        len = 0;
        index = 0;
        //printf("\n\nEsperando full\n\n");
        //sem_getvalue( &shared[0].full, &debug);
        //printf("\n\nfull:%d\n\n", debug);
        sem_wait(&shared[0].full);
        //printf("\n\nEsperando mutex\n\n");
        sem_wait(&shared[0].mutex);
        //printf("\n\nEsperando mais nada\n\n");
        
        fp = fopen(shared[0].buffer[shared[0].bufferOut]->Nome, "r");
        printf("\n\n*%s*\n\n", shared[0].buffer[shared[0].bufferOut]->Nome);
        if (fp == NULL){
            printf("\n\nErro na leitura do arquivo %s.\n\n", shared[0].buffer[shared[0].bufferOut]->Nome);
            exit(EXIT_FAILURE);        
        }
        //temp->Nome = shared[0].buffer[shared[0].bufferOut]->Nome;
        for( i=0; i<100; i++ ){
            temp->Nome[i] = shared[0].buffer[shared[0].bufferOut]->Nome[i];
        }
        shared[0].bufferOut = ( shared[0].bufferOut +1 )%5;

        sem_post(&shared[0].empty);
        sem_post(&shared[0].mutex);


        fscanf( fp, "%d", &ordem );
        temp->Ordem = ordem;
        
        temp->A = newMatrix(ordem);
        temp->B = newMatrix(ordem);
        temp->C = newMatrix(ordem);
        for( i=0; i<ordem; i++ ){
            for( j=0; j<ordem; j++ ){
                fscanf( fp, "%lf", &temp->A[i][j] );
            }
        }
        for( i=0; i<ordem; i++ ){
            for( j=0; j<ordem; j++ ){
                fscanf( fp, "%lf", &temp->B[i][j] );
            }
        }
        fclose(fp);

        sem_wait(&shared[1].empty);
        sem_wait(&shared[1].mutex);
        shared[1].buffer[shared[1].bufferIn] = temp;        
        shared[1].bufferIn  = ( shared[1].bufferIn + 1 )%5;
        sem_post(&shared[1].full);
        sem_post(&shared[1].mutex);
        
    }
}

/*MM - Thread Multiplicadora de Matrizes - Move shared[1]->buffer[out] para um ponteiro temporário, calcula C=A*B no elemento temporário e move o ponteiro temporário para shared[2]->buffer[in]. Teremos 4 instâncias desta thread.*/

void multiplicarMatrizes(){
    S *temp = ( S* )malloc( sizeof( S ) );
    int i;
    while(1){
        printf("\n\nEsperando full 1\n\n");
        sem_wait(&shared[1].full);
        printf("\n\nEsperando mutex 1\n\n");
        sem_wait(&shared[1].mutex);
        printf("\n\nEsperando mais nada 1\n\n");

        int ordem = shared[1].buffer[shared[1].bufferOut]->Ordem;
        //temp->Nome = shared[1].buffer[shared[1].bufferOut]->Nome;
        
        for( i=0; i<100; i++ ){
            temp->Nome[i] = shared[1].buffer[shared[1].bufferOut]->Nome[i];
        }
        

        temp->Ordem = ordem;
        temp->A = cloneMatrix ( shared[1].buffer[shared[1].bufferOut]->A, ordem );
        temp->B = cloneMatrix ( shared[1].buffer[shared[1].bufferOut]->B, ordem );
        temp->C = cloneMatrix ( shared[1].buffer[shared[1].bufferOut]->C, ordem );
        multiplyMatrices(temp->A, temp->B, temp->C, ordem);

        sem_post(&shared[1].empty);
        sem_post(&shared[1].mutex);

        sem_wait(&shared[2].empty);
        sem_wait(&shared[2].mutex);
        shared[2].buffer[shared[2].bufferIn] = temp;
        shared[1].bufferOut = ( shared[1].bufferOut +1 )%5;
        shared[2].bufferIn =  ( shared[2].bufferIn + 1 )%5;
        sem_post(&shared[2].full);
        sem_post(&shared[2].mutex);
    }
}


/*DM - Thread de Cálculo de determinante - Move shared[2]->buffer[out] para um ponteiro temporário, calcula E = det(C)  e move o ponteiro temporário para shared[3]->buffer[in]. Teremos 5 instâncias desta thread.*/



void determinanteMatriz(){
    S *temp = ( S* )malloc( sizeof( S ) );
    int i;
    while(1){
        printf("\n\nEsperando full 2\n\n");
        sem_wait(&shared[2].full);
        printf("\n\nEsperando mutex 2\n\n");
        sem_wait(&shared[2].mutex);
        printf("\n\nEsperando mais nada 2\n\n");

        int ordem = shared[2].buffer[shared[2].bufferOut]->Ordem;

        for( i=0; i<100; i++ ){
            temp->Nome[i] = shared[2].buffer[shared[2].bufferOut]->Nome[i];
        }
        printf("\n\ndetermiant %s\n\n", temp->Nome);
        temp->Ordem = ordem;
        temp->A = cloneMatrix ( shared[2].buffer[shared[2].bufferOut]->A, ordem );
        temp->B = cloneMatrix ( shared[2].buffer[shared[2].bufferOut]->B, ordem );
        temp->C = cloneMatrix ( shared[2].buffer[shared[2].bufferOut]->C, ordem );
        temp->E = determinant(temp->C, ordem);

        

        sem_post(&shared[2].empty);
        sem_post(&shared[2].mutex);

        sem_wait(&shared[3].empty);
        sem_wait(&shared[3].mutex);
        shared[3].buffer[shared[3].bufferIn] = temp;
        shared[2].bufferOut = ( shared[2].bufferOut +1 )%5;
        shared[3].bufferIn =  ( shared[3].bufferIn + 1 )%5;
        sem_post(&shared[3].full);
        sem_post(&shared[3].mutex);
    }
}



/*EA - Thread Escritora de Arquivo - Escreve um arquivo nome.out (nome = nome de entrada). O arquivo contém o nome do arquivo de entrada, a ordem das matrizes e os valores de A, B, C e E. Depois de imprimir, faz free na estrutura S.  Teremos 3 instâncias desta thread.*/

void escreverArquivo(){
    FILE * f = fopen(OUT_FILE, "w");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    while(1){
        printf("\n\nEsperando full 3\n\n");
        sem_wait(&shared[3].full);
        printf("\n\nEsperando mutex 3\n\n");
        sem_wait(&shared[3].mutex);
        printf("\n\nEsperando mais nada 3\n\n");

         S *temp =  shared[3].buffer[shared[3].bufferOut];
        shared[3].bufferOut = ( shared[3].bufferOut +1 )%5;


        printf("\n\n Arquivo %s\n",temp->Nome);
        printf("\n\n Ordem %d\n",temp->Ordem);


        fprintf(f,"================================\n");
        fprintf(f,"%s - %d\n", temp->Nome, temp->Ordem);
        fprintf(f, "--------------------------------\n");
        mat_print_file(temp->A, temp->Ordem, f);
        fprintf(f, "--------------------------------\n");
        mat_print_file(temp->B, temp->Ordem, f);
        fprintf(f, "--------------------------------\n");
        mat_print_file(temp->C, temp->Ordem, f);
        fprintf(f, "--------------------------------\n");
        fprintf(f, "%f\n",temp->E);
        fprintf(f,"================================\n");

        contEscrita++;

        printf("Escrita %d  Leitura %d\n", contEscrita, contLeitura);

        if (contEscrita == contLeitura) break;

        sem_post(&shared[3].empty);
        sem_post(&shared[3].mutex);

    }

    fclose(f);

}

void initShared(){
    int i, j;
    for( i=0; i<4; i++ ){
        shared[i].bufferIn = 0;
        shared[i].bufferOut = 0;
        //shared[i].full, shared[i].empty, shared[i].mutex;
        sem_init(&shared[i].mutex, 1, 1);
        sem_init(&shared[i].empty, 1, 5);
        sem_init(&shared[i].full, 1, 0);
        for( j=0;j<5;j++ ){
            shared[i].buffer[j] = ( S* )malloc( sizeof( S ) );
        }
    }
}

int main() {
    int tLerCatalogo[nLerCatalogo],
        tLerArquivos[nLerArquivos],
        tMultiplicarMatrizes[nMultiplicarMatrizes],
        tDeterminanteMatriz[nDeterminanteMatriz],
        tEscreverArquivo[nEscreverArquivo];
    pthread_t idLerCatalogo, 
              idLerArquivos,
              idMultiplicarMatrizes,
              idDeterminanteMatriz,
              idEscreverArquivo;
    int index;

    initShared();

    index = 0;
    contLeitura = 0;
    contEscrita = 0;
    
    while ((index<nLerCatalogo) || (index<nLerArquivos) || ( index < nMultiplicarMatrizes ) ) {
        if ( index < nLerCatalogo){
            tLerCatalogo[index]=index;
       		pthread_create(&idLerCatalogo, NULL, lerCatalogo, &tLerCatalogo[index]);
        }
        if ( index < nLerArquivos ){
            tLerArquivos[index]=index;
       		pthread_create(&idLerArquivos, NULL, lerArquivos, &tLerArquivos[index]);
        }
        if ( index < nMultiplicarMatrizes ){
            tMultiplicarMatrizes[index]=index;
            pthread_create(&idMultiplicarMatrizes, NULL, multiplicarMatrizes, &tMultiplicarMatrizes[index]);
        }
        if ( index < nDeterminanteMatriz ){
            tDeterminanteMatriz[index]=index;
            pthread_create(&idDeterminanteMatriz, NULL, determinanteMatriz, &tDeterminanteMatriz[index]);
        }
        if ( index < 1){//nEscreverArquivo ){
            tEscreverArquivo[index]=index;
            pthread_create(&idEscreverArquivo, NULL, escreverArquivo, &tEscreverArquivo[index]);
        }
        //if (index<numThreads){
        //    threadArray[index]=index;
        //    pthread_create(&idthread, NULL, threadFunc, &threadArray[index]);
        //}
        index++;
    }

    

    pthread_exit(NULL);
    
}

