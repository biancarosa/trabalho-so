#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <semaphore.h>
#include <string.h>

#define nLerCatalogo         1
#define nLerArquivos         3
#define nMultiplicarMatrizes 4
#define nDeterminanteMatriz  5
#define nEscreverArquivo     3
#define bufferSize           5



int debug, contLeitura, contEscrita;

sem_t mutexEscrita;

typedef double **mat;

typedef struct {
    char Nome[100];
    int Ordem;
    mat A, B, C;
    double E;
} S;

typedef struct{
    S *buffer[bufferSize];
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
    #pragma omp parallel shared( matrix, clone, size ) private(i,j)
    {
        #pragma omp for schedule(static)
        for( i=0; i<size; i++ ){
            for( j=0; j<size; j++ ){
                clone[i][j] = matrix[i][j];
            }
        }
    }
    return clone;
}

void zeroMatrix(mat matrix, int size) {
    int i, j;
    #pragma omp parallel shared( matrix, size ) private(i,j)
    {
        #pragma omp for schedule(static)
        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                matrix[i][j] = 0;
            }
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
        #pragma omp for schedule(static)
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

double determinant(mat matrix, int size) {
    double det = 1;
    int i, j, k;
    mat matL = newMatrix(size), matU = newMatrix(size);
    zeroMatrix(matL, size);
    zeroMatrix(matU, size);
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(2);
    #pragma omp parallel shared( matrix, matL, matU, size) private( i, j, k)
    {
        #pragma omp for schedule( static, 10 )
        for( k=0; k<size; k++ ){
            matL[k][k] = 1;
            for( i=k+1; i<size; i++ ){
                matL[i][k] = matrix[i][k] / matrix[k][k];
                for( j=k+1; j<size; j++ ){
                    matrix[i][j] = matrix[i][j] - matL[i][k] * matrix[k][j];
                }
                #pragma omp flush( matrix )
            }
            #pragma omp nowait
            for( j=k; j<size; j++ ){
                matU[k][j] = matrix[k][j];
            }
        }
        for (i = 0; i < size; i++) {
            det *= matrix[i][i];
        }
    }
    
    //printf( "\n\ndet: %lf\n\n", det );
    //mat_show( matU, size );
    //mat_show( matL, size );

    free(matL);
    free(matU);
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
        
        if( contLeitura < 0 ){
            contLeitura = 0;
        }
        contLeitura++;

        if (read > 1 ){
            //printf("\nesperando empty lerCatalogo");
            sem_wait(&shared[0].empty);
            //printf("\nesperando mutex lerCatalogo");
            sem_wait(&shared[0].mutex);
            //printf("\nescrita lerCatalogo");
            
            

            strtok( line, "\n" );
            strcpy( shared[0].buffer[shared[0].bufferIn]->Nome, line );
            //for( i=0; i<bufferSize; i++ ){
               // printf("\n arquivo %s",shared[0].buffer[i]->Nome);
           // }
            

            //printf("%s", shared[0].buffer[shared[0].bufferIn]->Nome);
           // sem_getvalue( &shared[0].empty, &debug);
           // printf("\n\nempty:%d\n\n", debug);
            shared[0].bufferIn = (shared[0].bufferIn+1)%bufferSize;
            
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
    S *temp;
    while( 1 ){
        line = NULL;
        len = 0;
        index = 0;
        //printf("\nesperando full lerArquivos");
        sem_wait(&shared[0].full);
        //printf("\nesperando mutex lerArquivos");
        sem_wait(&shared[0].mutex);
        //printf("\nleitura lerArquivos");
        temp = ( S* )malloc( sizeof( S ) );        
        fp = fopen(shared[0].buffer[shared[0].bufferOut]->Nome, "r");
        if (fp == NULL){
            printf("\n\nErro na leitura do arquivo %s.\n\n", shared[0].buffer[shared[0].bufferOut]->Nome);
            exit(EXIT_FAILURE);        
        }

        strcpy( temp->Nome, shared[0].buffer[shared[0].bufferOut]->Nome );        
        
        shared[0].bufferOut = ( shared[0].bufferOut +1 )%bufferSize;

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
        //printf("\nesperando empty lerArquivos");
        sem_wait(&shared[1].empty);
        //printf("\nesperando mutex lerArquivos");
        sem_wait(&shared[1].mutex);
        //printf("\nescrita lerArquivos");
        shared[1].buffer[shared[1].bufferIn] = temp;


        //printf("\n\n%s   %d", temp->Nome, shared[1].bufferIn);
        //printf("\nbufferDump lerArquivos\n\n");
        //for( i=0; i<bufferSize; i++ ){
            //rintf("\n%d arquivo %s",i ,shared[1].buffer[i]->Nome);
        //}
        temp = NULL;

        shared[1].bufferIn = (shared[1].bufferIn+1)%bufferSize;
        sem_post(&shared[1].full);
        sem_post(&shared[1].mutex);
        
    }
}

/*MM - Thread Multiplicadora de Matrizes - Move shared[1]->buffer[out] para um ponteiro temporário, calcula C=A*B no elemento temporário e move o ponteiro temporário para shared[2]->buffer[in]. Teremos 4 instâncias desta thread.*/

void multiplicarMatrizes(){
    S *temp;
    int i;
    while(1){
        //printf("\nesperando full multiplicarMatrizes");
        sem_wait(&shared[1].full);
        //printf("\nesperando mutex multiplicarMatrizes");
        sem_wait(&shared[1].mutex);
        //printf("\nleitura multiplicarMatrizes");

        temp = ( S* )malloc( sizeof( S ) );
        int ordem = shared[1].buffer[shared[1].bufferOut]->Ordem;        
        strcpy( temp->Nome, shared[1].buffer[shared[1].bufferOut]->Nome );
        temp->Ordem = ordem;
        temp->A = cloneMatrix ( shared[1].buffer[shared[1].bufferOut]->A, ordem );
        temp->B = cloneMatrix ( shared[1].buffer[shared[1].bufferOut]->B, ordem );
        temp->C = cloneMatrix ( shared[1].buffer[shared[1].bufferOut]->C, ordem );
        multiplyMatrices(temp->A, temp->B, temp->C, ordem);
        shared[1].bufferOut = ( shared[1].bufferOut +1 )%bufferSize;


        
        sem_post(&shared[1].empty);
        sem_post(&shared[1].mutex);

        
        //printf("\nesperando empty multiplicarMatrizes");
        sem_wait(&shared[2].empty);
        //printf("\nesperando mutex multiplicarMatrizes");
        sem_wait(&shared[2].mutex);
        //printf("\nescrita multiplicarMatrizes");

        shared[2].buffer[shared[2].bufferIn] = temp;            
        shared[2].bufferIn =  ( shared[2].bufferIn + 1 )%bufferSize;

        sem_post(&shared[2].full);
        sem_post(&shared[2].mutex);
    }
}


/*DM - Thread de Cálculo de determinante - Move shared[2]->buffer[out] para um ponteiro temporário, calcula E = det(C)  e move o ponteiro temporário para shared[3]->buffer[in]. Teremos 5 instâncias desta thread.*/



void determinanteMatriz(){
    S *temp;
    int i;
    while(1){
        //printf("\nesperando full determinanteMatriz");
        sem_wait(&shared[2].full);
        //printf("\nesperando mutex determinanteMatriz");
        sem_wait(&shared[2].mutex);
        //printf("\nleitura determinanteMatriz");

        temp = ( S* )malloc( sizeof( S ) );
        int ordem = shared[2].buffer[shared[2].bufferOut]->Ordem;
        strcpy( temp->Nome, shared[2].buffer[shared[2].bufferOut]->Nome );
        temp->Ordem = ordem;
        temp->A = cloneMatrix ( shared[2].buffer[shared[2].bufferOut]->A, ordem );
        temp->B = cloneMatrix ( shared[2].buffer[shared[2].bufferOut]->B, ordem );
        temp->C = cloneMatrix ( shared[2].buffer[shared[2].bufferOut]->C, ordem );
        temp->E = determinant(temp->C, ordem);

        //printf("\n\n bufferDump determinanteMatriz\n\n");
        //for( i=0; i<bufferSize; i++ ){
            //printf("\n arquivo %s",shared[2].buffer[i]->Nome);
        //}
        shared[2].bufferOut = ( shared[2].bufferOut +1 )%bufferSize;
        sem_post(&shared[2].empty);
        sem_post(&shared[2].mutex);

        //printf("\nesperando empty determinanteMatriz");
        sem_wait(&shared[3].empty);
        //printf("\nesperando mutex determinanteMatriz");
        sem_wait(&shared[3].mutex);
        //printf("\nescrita determinanteMatriz");
        shared[3].buffer[shared[3].bufferIn] = temp;
        temp = NULL;
        shared[3].bufferIn =  ( shared[3].bufferIn + 1 )%bufferSize;
        sem_post(&shared[3].full);
        sem_post(&shared[3].mutex);
    }
}



/*EA - Thread Escritora de Arquivo - Escreve um arquivo nome.out (nome = nome de entrada). O arquivo contém o nome do arquivo de entrada, a ordem das matrizes e os valores de A, B, C e E. Depois de imprimir, faz free na estrutura S.  Teremos 3 instâncias desta thread.*/

void escreverArquivo(){
    char arquivoSaida[100];
    int i;
    while(1){

        //printf("\nesperando full escreverArquivo");
        sem_wait(&shared[3].full);
        //printf("\nesperando mutex escreverArquivo");
        sem_wait(&shared[3].mutex);
        //printf("\nescrita escreverArquivo");


        S *temp =  shared[3].buffer[shared[3].bufferOut];
        shared[3].bufferOut = ( shared[3].bufferOut +1 )%bufferSize;

        
        strcpy( arquivoSaida,"" );
        strcat( arquivoSaida, temp->Nome);
        strcat( arquivoSaida, ".out");
        FILE * f = fopen(arquivoSaida, "wb");
        if (f == NULL) {
            printf("Error opening file!\n");
            exit(1);
        }


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

        fclose(f);
    
        sem_post(&shared[3].empty);
        sem_post(&shared[3].mutex);
        
        //printf("\nesperando mutexEscrita escreverArquivo");
        sem_wait(&mutexEscrita);
        //printf("\ncontEscrita++");
        contEscrita++;
        sem_post(&mutexEscrita);
    }

}

void initShared(){
    int i, j;
    for( i=0; i<4; i++ ){
        shared[i].bufferIn = 0;
        shared[i].bufferOut = 0;
        //shared[i].full, shared[i].empty, shared[i].mutex;
        sem_init(&shared[i].mutex, 1, 1);
        sem_init(&shared[i].empty, 1, bufferSize);
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

    contLeitura = -1;
    contEscrita = 0;
    
    while( ( index < nLerCatalogo) || ( index < nLerArquivos ) || ( index < nMultiplicarMatrizes )
                                   || ( index < nDeterminanteMatriz ) || ( index < nEscreverArquivo )   ) {
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
        if ( index < nEscreverArquivo ){
            tEscreverArquivo[index]=index;
            pthread_create(&idEscreverArquivo, NULL, escreverArquivo, &tEscreverArquivo[index]);
        }
        index++;
        //if (index<numThreads){
        //    threadArray[index]=index;
        //    pthread_create(&idthread, NULL, threadFunc, &threadArray[index]);
        //}
    }
    
    sem_init(&mutexEscrita, 1, 1);
    while(1) {
        //printf("\nesperando mutexEscrita main");
        sem_wait(&mutexEscrita);
        sem_wait( &shared[0].mutex );
        //printf("\n\n%d..%d\n\n", contEscrita, contLeitura);
        if (contEscrita == contLeitura && contLeitura >=0 ) {
            //printf("\nmatando tudo\n\n");
            sem_post(&mutexEscrita);
            //printf("\n\nAcabou!!! %d..%d\n\n", contEscrita, contLeitura);

            int retvalue = pthread_cancel(idEscreverArquivo);
            //printf("Retorno / thread_cancel EscreverArquivo:: %d \n",retvalue);
            retvalue = pthread_cancel(idDeterminanteMatriz);
            //printf("Retorno / thread_cancel DeterminanteMatriz:: %d \n",retvalue);
            retvalue = pthread_cancel(idMultiplicarMatrizes);
            //printf("Retorno / thread_cancel MultiplicarMatrizes:: %d \n",retvalue);
            retvalue = pthread_cancel(idLerArquivos);
            //printf("Retorno / thread_cancel LerArquivos:: %d \n",retvalue);
            retvalue = pthread_cancel(idLerCatalogo);
            //printf("Retorno / thread_cancel LerCatalogo:: %d \n",retvalue);

            return 0;
        };
        sem_post(&mutexEscrita);
        sem_post( &shared[0].mutex );
    }


    return 1;
    
}
