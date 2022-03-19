//
// Created by Aviv Avraham on 2/27/2022.
//

#include "nsclustering.h"

/*
 NEW CODE
*/

/*
 this function takes the data points matrix and creates the weighted matrix W[nXn]
 that represent the weight between every two data points
 while n represent the amount of data points received
 the function returns the W matrix
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
 * TODO: finish, where is the diagonal degree Matrix? need to create function for calculating it
*/
double **calculate_lnorm_matrix(double **lnorm_matrix,double** weight_matrix, double ** diagonal_matrix, int row_size, int column_size){
    double ** product = allocate_array_2d(row_size,column_size);
    for (int i = 0; i < row_size; ++i) {
        product[i][i] = 1; //turning the product matrix to the unit matrix,I.
    }
    degree_matrix_normalized(diagonal_matrix, row_size);
    product = matrix_subtraction(product,matrix_multiplication(matrix_multiplication(
                                                                       diagonal_matrix,row_size,column_size,
                                                                       weight_matrix,row_size,column_size),
                                                               row_size,column_size,
                                                               diagonal_matrix,row_size,
                                                               column_size),row_size,column_size);
    // product = I - (D^-0.5)*W*D(^-0.5)
    return product;
}

/*
 * this function sets the i,j for them:
 * max{|A_ij|} for any j>i
 */

Point find_Jacobi_ij(double ** matrix,int n) {
    int i,j;
    double max=0, num;
    Point p;

    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            num = matrix[i][j];
            if (num > max){
                max = num;
                p.i = i;
                p.j = j;
            }
        }
    }
    return p;
}

double Jacobi_find_diff_off(double **A,double **A_new, int n){
    double sum_A=0,sum_A_new=0;
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i == j){
                continue;
            }
            sum_A += pow(A[i][j],2);
            sum_A_new += pow(A_new[i][j],2);
        }
    }

    return sum_A - sum_A_new;
}

int get_sign(double num){
    if(num<0){
        return -1;
    } else{
        return 1;
    }
}

void set_matrix_A_new(double **A_new, double ** A,int n,int i,int j,double c,double s){
    int row,column;

    for(row=0;row<n;row++){
        for(column=0;column<n;column++){
            if(row != i && row != j){
                A_new[row][i] = c * A[row][i] - s * A[row][j];
                A_new[row][j] = c * A[row][j] + s * A[row][i];
            }
        }
    }

    A_new[i][i] = pow(c,2) * A[i][i] + pow(s,2) * A[j][j] - 2 * s * c * A[i][j];
    A_new[j][j] = pow(s,2) * A[i][i] + pow(c,2) * A[j][j] + 2 * s * c * A[i][j];
    A_new[i][j] = 0;
}

void Jacobi_algorithm(double **laplacian, int n){
    double diff,theta,t,c,s;
    double epsilon = pow(10,-5);
    int i,j;
    int num_iter=100;
    Point p;
    double ** V = allocate_array_2d(n,n);
    double ** A = allocate_array_2d(n,n);
    double ** A_new = allocate_array_2d(n,n);
    double ** P_matrix = allocate_array_2d(n,n);
    set_equal_array_2d(A,laplacian,n,n);
    set_equal_array_2d(A_new,laplacian,n,n);



    while( num_iter > 0 && diff < epsilon){
        p = find_Jacobi_ij(A,n);
        i = p.i;
        j = p.j;
        theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
        t = get_sign(theta) / (fabs(theta) + sqrt(pow(theta,2) + 1 ));
        c = 1 / (sqrt(pow(t,2)) + 1 );
        s = t * c;



        set_matrix_A_new(A_new,A,n,i,j,c,s);


        diff = Jacobi_find_diff_off(A,A_new,n);
        num_iter--;
    }
}
