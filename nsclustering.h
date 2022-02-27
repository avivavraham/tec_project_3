//
// Created by aviva on 2/27/2022.
//

#ifndef TEC_PROJECT_3_NSCLUSTERING_H
#define TEC_PROJECT_3_NSCLUSTERING_H

double weight(double* data_point_1, double* data_point_2);
double** weighted_matrix(double** data_points_array);
double **allocate_array_2d(int r, int c);
double **degree_matrix(double ** weighted_matrix);
double sum_of_row(double *row_vector);
void print_matrix(double **matrix, int size_row, int size_column);
void **pow_matrix(double ** matrix,int size_row, int size_column , double pow_, int diag_pow_optional);
double **degree_matrix_normalized(double ** degree_matrix);

#endif //TEC_PROJECT_3_NSCLUSTERING_H
