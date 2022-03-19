//
// Created by Aviv Avraham on 2/27/2022.
//

#include "nsclustering.h"
#include "spkmeans.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
    1.1.1 - calculate_weighted_matrix
    1.1.2 - calculate_diagonal_degree_matrix
    1.1.3 -
    1.2 -
    1.3 - calculate_eigengap_heuristic
*/

/*
 this function takes the data points matrix and creates the weighted matrix W[nXn]
 that represent the weight between every two data points
 while n represent the amount of data points received
 the function returns the W matrix (1.1.1)
*/

void calculate_weighted_matrix(double ** weighted_matrix, double** data_points_array, int num_rows, int d){
    int i,j;
    for (i = 0; i < num_rows; ++i) {
        for (j = 0; j < num_rows; ++j) {
            if(i==j){
                weighted_matrix[i][j] = 0;
            }
            if(j < i){
                weighted_matrix[i][j] = weighted_matrix[j][i];
            }
            if(i < j){
                weighted_matrix[i][j] = calculate_weight(data_points_array[i],data_points_array[j],d);
            }
        }
    }
}

/*
 this is the weight function as followed-
  wij = exp(-(l_norm||dt1,dt2|| / 2))
  while l_norm is as followed-
  l_norm||dt1,dt2|| = sqrt_root(sum_i((dt1[i] - dt2[i])^2))
*/

double calculate_weight(const double* data_point_1, const double* data_point_2,int d){
    double sum = 0;
    int i;
    for(i=0; i < d; i++){
        double sub = data_point_1[i] - data_point_2[i];
        sum += pow(sub,2);
    }
    sum = sqrt(sum);
    sum = -sum / 2;
    return exp(sum);
}

/*
 * The diagonal degree matrix D ∈ Rn×n
is defined as D = (dij)(i,j) = 1,...,n, such that:
dij =
        sum(Wiz)z=1 until z=n, if i = j
        0, otherwise
*/
void **calculate_degree_matrix(double **degree_matrix, double ** weighted_matrix, int N){
    for (int i = 0; i < N; ++i) {
        degree_matrix[i][i] = sum_of_row(weighted_matrix[i],N);
    }
}
/*
 * returns the sum of values in a given row vector
*/
double sum_of_row(const double *row_vector, int len){
    double sum=0;
    for (int i = 0; i < len; ++i) {
        sum += row_vector[i];
    }
    return sum;
}
/*
 * prints the values of a given matrix
*/
void print_matrix(double **matrix, int size_row, int size_column){
    int i,j;
    printf("[");
    for (i = 0; i < size_row ; ++i) {
        printf("[");
        for (j = 0; j < size_column; ++j) {
            printf("%f", matrix[i][j]);
            if (j != size_column-1){printf(",");}
        }
        if(i == size_row-1){printf("]");}
        printf("]\n");
    }
    printf("\n");
}

/*
 * this function raising every element of a given matrix by the power of a given integer.
 * optional: if diag_pow_optional = 1(true), only the element on the diagonal
 * of the matrix are raised by the given power.
*/
void **pow_matrix(double ** matrix,int size_row, int size_column , double pow_, int diag_pow_optional){
    if(diag_pow_optional){
        for (int i = 0; i < size_row; ++i) {
            matrix[i][i] = pow(matrix[i][i],pow_);
        }
    } else{
        for (int i = 0; i < size_row; ++i) {
            for (int j = 0; j < size_column; ++j) {
                matrix[i][j] = pow(matrix[i][j],pow_);
            }
        }
    }
}

/*
 * this function return the D matrix normalized by the power of -1/2.
*/
double **degree_matrix_normalized(double ** degree_matrix, int num_rows){
    pow_matrix(degree_matrix,num_rows,num_rows,-0.5,1);
    return degree_matrix;
}

/*
 * this function multiply matrix1 by matrix2 and return the product matrix.
 * assumption: the dimensions for the multiplication are valid,there for,column1_size = row2_size.
*/
void matrix_multiplication(double ** matrix1, int row1_size, int column1_size,
                               double ** matrix2, int row2_size, int column2_size, double **result){
    for (int i = 0; i < row1_size; ++i) {
        for (int j = 0; j < column2_size; ++j) {
            double sum = 0;
            for (int k = 0; k < column1_size; ++k) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = sum;
        }
    }
}

/*
 * this function subtract matrix1 by matrix2 and return the product matrix.
 * assumption: the dimensions for the subtraction are valid,
 * there for,row_size for matrix 1 = row_size for matrix 2. and so for columns sizes.
*/
void matrix_subtraction(double ** matrix1, double ** matrix2, int row_size, int column_size, double **result){
    for (int i = 0; i < row_size; ++i) {
        for (int j = 0; j < column_size; ++j) {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
}

/*
 * this function calculates the Diagonal Degree Matrix (1.1.2)
 *
 */

void calculate_diagonal_degree_matrix(double **matrix, double** weight_matrix,int len){
    calculate_degree_matrix(matrix,weight_matrix,len);
    degree_matrix_normalized(matrix,len); //returns D^-(0.5)
}

/*
 * TODO: finish,
*/
void calculate_lnorm_matrix(double **lnorm_matrix,double** weight_matrix, double **diagonal_matrix, int n){
    double **Identity_matrix = allocate_array_2d(n,n);
    double **result1 = allocate_array_2d(n,n);
    double **result2 = allocate_array_2d(n,n);
    error_occurred(Identity_matrix == NULL);
    error_occurred(result1 == NULL);
    for (int i = 0; i < n; ++i) {
        Identity_matrix[i][i] = 1;
    }

    // Calculate W x D(^-0.5)
    matrix_multiplication(weight_matrix,n,n,diagonal_matrix,n,n,result1);

    // Calculate (D^-0.5)*W*D(^-0.5)
    matrix_multiplication(diagonal_matrix,n,n,result1,n,n,result2);

    // Calculate lnorm_matrix = I - (D^-0.5)*W*D(^-0.5)
    matrix_subtraction(Identity_matrix,result2,n,n,lnorm_matrix);

    free_array_2d(Identity_matrix,n);
    free_array_2d(result1,n);
    free_array_2d(result2,n);
}
/*
 * the function gets *sorted* eigenvalues array (1.3)
 */

int calculate_eigengap_heuristic(double **eigenvalues,int n){
    int i,j;
    int max=0,argmax=0;
    double **eigengap = allocate_array_2d(n,n);
    error_occurred(eigengap == NULL);

    for(i=0;i<n-1;i++){ //TODO: verify this calculation
        for(j=0;j<n;j++ ){
            eigengap[i] = eigenvalues[i+1] - eigenvalues[i];
        }
    }

    for(i=0;i < floor(n/2);i++){
        if(eigengap[i] > max){
            max = eigengap[i];
            argmax = i;
        }
    }

    free_array_2d(eigengap,n);
    return argmax;
}


