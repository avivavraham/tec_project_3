//
// Created by aviva on 2/27/2022.
//

#ifndef TEC_PROJECT_3_NSCLUSTERING_H
#define TEC_PROJECT_3_NSCLUSTERING_H

typedef struct
{
    int i;
    int j;
} Point;

#include "libraries.h"
#include "spkmeans.h"

double calculate_weight(const double* data_point_1, const double* data_point_2, int);
void calculate_weighted_matrix(double ** weighted_matrix, double** data_points_array, int num_rows, int d);
double **allocate_array_2d(int r, int c);
void **calculate_degree_matrix(double **degree_matrix, double ** weighted_matrix, int N);
double sum_of_row(const double *row_vector, int);
void print_matrix(double **matrix, int size_row, int size_column);
void **pow_matrix(double ** matrix,int size_row, int size_column , double pow_, int diag_pow_optional);
double **degree_matrix_normalized(double ** degree_matrix, int);
double **matrix_multiplication(double ** matrix1, int row1_size, int column1_size,
                               double ** matrix2, int row2_size, int column2_size);
double **matrix_subtraction(double ** matrix1, double ** matrix2, int row_size, int column_size);
double **calculate_lnorm_matrix(double **lnorm_matrix,double** weight_matrix, double ** diagonal_matrix, int row_size, int column_size);
Point find_Jacobi_ij(double ** matrix,int n);
int get_sign(double num);
double Jacobi_find_diff_off(double **A,double **A_new, int n);
void set_matrix_A_new(double **A_new, double ** A,int n,int i,int j,double c,double s);


#endif //TEC_PROJECT_3_NSCLUSTERING_H
