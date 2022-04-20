#ifndef TEC_PROJECT_3_NSCLUSTERING_H
#define TEC_PROJECT_3_NSCLUSTERING_H

#include "helper.h"

typedef struct
{
    int i;
    int j;
} Point;

typedef struct
{
    double value;
    double *vector;
} Eigen;

double calculate_weight(const double* data_point_1, const double* data_point_2, int);
void calculate_weighted_matrix(double ** weighted_matrix, double** data_points_array, int num_rows, int d);
void calculate_degree_matrix(double **degree_matrix, double ** weighted_matrix, int N);
double sum_of_row(const double *row_vector, int);
void print_matrix(double **matrix, int size_row, int size_column);
void set_matrix_to_Identity(double ** matrix, int n);
void pow_matrix(double ** matrix,int size_row, int size_column , double pow_, int diag_pow_optional);
void **degree_matrix_normalized(double ** degree_matrix, int num_rows);
void calculate_diagonal_degree_matrix(double **matrix, double** weight_matrix,int len);
void matrix_multiplication(double ** matrix1, int row1_size, int column1_size,
                           double ** matrix2, int column2_size, double **result);
void matrix_subtraction(double ** matrix1, double ** matrix2, int row_size, int column_size, double **result);
void calculate_lnorm_matrix(double **lnorm_matrix,double** weight_matrix, double **diagonal_matrix, int n);
Point Jacobi_find_ij(double ** matrix,int n);
int Jacobi_get_sign(double num);
double Jacobi_find_diff_off(double **A,double **A_new, int n);
void transpose_matrix(double **transposed,double ** matrix, int rows, int columns);
void Jacobi_set_matrix_A_new(double **A_new, double ** A,int n,int i,int j,double c,double s);
void Jacobi_set_matrix_P(double **matrix, int n, int i,int j, double c, double s);
void Jacobi_algorithm(double **laplacian, int n, Eigen *eigen);
int calculate_eigengap_heuristic(Eigen *eigens ,int n);
void init_Eigen_struct(Eigen *eigen, int n);
void free_Eigen_struct(Eigen *eigen, int n);
void Jacobi_set_Eigen(Eigen *eigen, int n, double **vectors, double **values);
void sort_eigen(Eigen *eigen, int n);
void set_U_matrix(double ** U,Eigen *eigen, int n, int k);
double get_squared_sum_of_row(double ** matrix,int n, int c);
void calculate_T_matrix(double ** T,double ** U,int n,int k);
void print_eigenvalues(Eigen *eigen, int num_rows);
void get_first_k_columns(double **to,double **from, int n, int k);


#endif
