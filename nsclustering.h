//
// Created by aviva on 2/27/2022.
//

#ifndef TEC_PROJECT_3_NSCLUSTERING_H
#define TEC_PROJECT_3_NSCLUSTERING_H

double weight(double* data_point_1, double* data_point_2, int);
double** weighted_matrix(double** data_points_array, int);
double **allocate_array_2d(int r, int c);
double **degree_matrix(double ** weighted_matrix, int);
double sum_of_row(double *row_vector, int);
void print_matrix(double **matrix, int size_row, int size_column);
void **pow_matrix(double ** matrix,int size_row, int size_column , double pow_, int diag_pow_optional);
double **degree_matrix_normalized(double ** degree_matrix, int);
double **matrix_multiplication(double ** matrix1, int row1_size, int column1_size,
                               double ** matrix2, int row2_size, int column2_size);
double **matrix_subtraction(double ** matrix1, double ** matrix2, int row_size, int column_size);
double **lnorm_matrix(double ** weight_matrix, double ** diagonal_matrix, int row_size, int column_size);

#endif //TEC_PROJECT_3_NSCLUSTERING_H
