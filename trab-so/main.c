#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double **mat;

void mat_show( mat matrix, int size ){
    int i, j;
    printf("\n\n");
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            printf(" %f", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void wolfranFormat( mat matrix, int size ){
    int i, j;
    printf("{");
    for( i=0;i<size;i++ ){
        printf("{");
        for( j=0; j<size;j++){
            printf("%d",(int)matrix[i][j]);
            if( j<size-1 ){
                printf( "," );
            }
        
        }
        printf("}");
        if( i<size-1 ){
            printf( "," );
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
            matrix[i][j] = rand()%10;
        }
    }
}

void multplyMatrices(mat matA, mat matB, mat resolt, int size) {
    int i, j, k;
    zeroMatrix(resolt, size);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            for (k = 0; k < size; k++) {
                resolt[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }
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
    multplyMatrices(matP, matrix, Aprime, size);
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
                matL[j][i] = (Aprime[j][i] - s) / matU[i][i];
            }
        }
    }
    for (i = 0; i < size; i++) {
        matL[i][i] = 1;
    }
    
}

double triangularMatrixDeterminant(mat matrix, int size) {
    int i, det = 1;
    for (i = 0; i < size; i++) {
        det *= matrix[i][i];
    }
    return det;
}

double determinant(mat matrix, int size) {
    double det;
    mat matL = newMatrix(size), matU = newMatrix(size),matP = newMatrix(size);
    mat_LU(matrix, matL, matU, matP, size);
    det = triangularMatrixDeterminant(matU, size);
    
    mat debug = newMatrix(size);
    multplyMatrices(matL, matU, debug, size);
    
    printf( "\n\nL" );
    mat_show(matL,size);
    printf( "u" );
    mat_show(matU,size);
    printf( "LU" );
    mat_show(debug,size);
    wolfranFormat(debug, size);
    
    
    free(matL);
    free(matU);
    free(matP);
    return det;
}

void mat_show( mat matrix, int size ){
    int i, j;
    printf("\n\n");
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            printf(" %f", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void wolfranFormat( mat matrix, int size ){
    int i, j;
    printf("{");
    for( i=0;i<size;i++ ){
        printf("{");
        for( j=0; j<size;j++){
            printf("%d",(int)matrix[i][j]);
            if( j<size-1 ){
                printf( "," );
            }
        
        }
        printf("}");
        if( i<size-1 ){
            printf( "," );
        }
    }
    printf("}");
}

int main() {
    int size = 4;
    srand( time(NULL) );
    mat matrix = newMatrix(size);
    randomFill(matrix, size);
    printf( "original" );
    mat_show(matrix, size);
    wolfranFormat(matrix, size);
    printf("Determinant: %f", determinant(matrix, size) );
}