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
 NEW CODE
*/

/*
 this function takes the data points matrix and creates the weighted matrix W[nXn]
 that represent the weight between every two data points
 while n represent the amount of data points received
 the function returns the W matrix
*/

double** weighted_matrix(double** data_points_array){
    double** weighted_matrix = allocate_array_2d(num_rows, num_rows);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_rows; ++j) {
            if(i==j){
                weighted_matrix[i][j] = 0;
            }
            if(j < i){
                weighted_matrix[i][j] = weighted_matrix[j][i];
            }
            if(i < j){
                weighted_matrix[i][j] = weight(data_points_array[i],data_points_array[j]);
            }
        }
    }
    return weighted_matrix;
}

/*
 this is the weight function as followed-
  wij = exp(-(l_norm||dt1,dt2|| / 2))
  while l_norm is as followed-
  l_norm||dt1,dt2|| = sqrt_root(sum_i((dt1[i] - dt2[i])^2))
*/

double weight(double* data_point_1, double* data_point_2){
    double sum = 0;
    for(int i=0; i < d; i++){
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
double **degree_matrix(double ** weighted_matrix){
    double** degree_matrix = allocate_array_2d(num_rows, num_rows);
    for (int i = 0; i < num_rows; ++i) {
        degree_matrix[i][i] = sum_of_row(weighted_matrix[i]);
    }
    return degree_matrix;
}
/*
 * returns the sum of values in a given row vector
*/
double sum_of_row(double *row_vector){
    double sum=0;
    for (int i = 0; i < num_rows; ++i) {
        sum += row_vector[i];
    }
    return sum;
}
/*
 * prints the values of a given matrix
*/
void print_matrix(double **matrix, int size_row, int size_column){
    printf("[");
    for (int i = 0; i < size_row ; ++i) {
        printf("[");
        for (int j = 0; j < size_column; ++j) {
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
double **degree_matrix_normalized(double ** degree_matrix){
    pow_matrix(degree_matrix,num_rows,num_rows,-0.5,1);
    return degree_matrix;
}

/*
 * this function multiply matrix1 by matrix2 and return the product matrix.
 * assumption: the dimensions for the multiplication are valid,there for,column1_size = row2_size.
*/
double **matrix_multiplication(double ** matrix1, int row1_size, int column1_size,
                               double ** matrix2, int row2_size, int column2_size){
    double ** product = allocate_array_2d(row1_size,column2_size);
    for (int i = 0; i < row1_size; ++i) {
        for (int j = 0; j < column2_size; ++j) {
            double sum = 0;
            for (int k = 0; k < column1_size; ++k) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            product[i][j] = sum;
        }
    }
    return product;
}

/*
 * this function subtract matrix1 by matrix2 and return the product matrix.
 * assumption: the dimensions for the subtraction are valid,
 * there for,row_size for matrix 1 = row_size for matrix 2. and so for columns sizes.
*/
double **matrix_subtraction(double ** matrix1, double ** matrix2, int row_size, int column_size){
    double ** product = allocate_array_2d(row_size,column_size);
    for (int i = 0; i < row_size; ++i) {
        for (int j = 0; j < column_size; ++j) {
            product[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
    return product;
}

/*
 * TODO: finish
*/
double **lnorm_matrix(double ** weight_matrix, double ** diagonal_matrix, int row_size, int column_size){
    double ** product = allocate_array_2d(row_size,column_size);
    for (int i = 0; i < row_size; ++i) {
        product[i][i] = 1; //turning the product matrix to the unit matrix,I.
    }
    degree_matrix_normalized(diagonal_matrix);
    product = matrix_subtraction(product,matrix_multiplication(matrix_multiplication(
                                                                       diagonal_matrix,row_size,column_size,
                                                                       weight_matrix,row_size,column_size),
                                                               row_size,column_size,
                                                               diagonal_matrix,row_size,
                                                               column_size),row_size,column_size);
    // product = I - (D^-0.5)*W*D(^-0.5)
    return product;
}
