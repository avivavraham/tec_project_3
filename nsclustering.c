#include "nsclustering.h"

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
void calculate_degree_matrix(double **degree_matrix, double ** weighted_matrix, int N){
    int i;
    for (i = 0; i < N; ++i) {
        degree_matrix[i][i] = sum_of_row(weighted_matrix[i],N);
    }
}
/*
 * returns the sum of values in a given row vector
*/
double sum_of_row(const double *row_vector, int len){
    double sum=0;
    int i;
    for (i = 0; i < len; ++i) {
        sum += row_vector[i];
    }
    return sum;
}
/*
 * prints the values of a given matrix
*/
void print_matrix(double **matrix, int size_row, int size_column){
    int i,j;
    for (i = 0; i < size_row ; ++i) {
        for (j = 0; j < size_column; ++j) {
            printf("%.4f", matrix[i][j]);
            if (j != size_column-1){printf(",");}
        }
        printf("\n");
    }
}

/*
 * this function raising every element of a given matrix by the power of a given integer.
 * optional: if diag_pow_optional = 1(true), only the element on the diagonal
 * of the matrix are raised by the given power.
*/
void pow_matrix(double ** matrix,int size_row, int size_column , double pow_, int diag_pow_optional){
    int i,j;
    if(diag_pow_optional){
        for (i = 0; i < size_row; ++i) {
            matrix[i][i] = pow(matrix[i][i],pow_);
        }
    } else{
        for (i = 0; i < size_row; ++i) {
            for (j = 0; j < size_column; ++j) {
                matrix[i][j] = pow(matrix[i][j],pow_);
            }
        }
    }
}

/*
 * this function return the D matrix normalized by the power of -1/2.
*/
void **degree_matrix_normalized(double ** degree_matrix, int num_rows){
    pow_matrix(degree_matrix,num_rows,num_rows,-0.5,1);
    return NULL;
}

/*
 * this function multiply matrix1 by matrix2 and return the product matrix.
 * assumption: the dimensions for the multiplication are valid,there for,column1_size = row2_size.
*/
void matrix_multiplication(double ** matrix1, int row1_size, int column1_size,
                               double ** matrix2, int column2_size, double **result){
    int i,j,k;
    for (i = 0; i < row1_size; ++i) {
        for (j = 0; j < column2_size; ++j) {
            double sum = 0;
            for (k = 0; k < column1_size; ++k) {
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
    int i,j;
    for (i = 0; i < row_size; ++i) {
        for (j = 0; j < column_size; ++j) {
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
}

void set_matrix_to_Identity(double ** matrix, int n){
    int i;
    for (i = 0; i < n; ++i) {
        matrix[i][i] = 1;
    }
}

/*
 *
*/
void calculate_lnorm_matrix(double **lnorm_matrix,double** weight_matrix, double **diagonal_matrix, int n){
    double **Identity_matrix = allocate_array_2d(n,n);
    double **result1 = allocate_array_2d(n,n);
    double **result2 = allocate_array_2d(n,n);
    set_matrix_to_Identity(Identity_matrix,n);

    degree_matrix_normalized(diagonal_matrix,n); /* returns D^-(0.5) */

    /* Calculate W x D(^-0.5) */
    matrix_multiplication(weight_matrix,n,n,diagonal_matrix,n,result1);

    /* Calculate (D^-0.5)*W*D(^-0.5) */
    matrix_multiplication(diagonal_matrix,n,n,result1,n,result2);

    /* Calculate lnorm_matrix = I - (D^-0.5)*W*D(^-0.5) */
    matrix_subtraction(Identity_matrix,result2,n,n,lnorm_matrix);

    free_array_2d(Identity_matrix,n);
    free_array_2d(result1,n);
    free_array_2d(result2,n);
}

/*
 * this function sets the i,j for them:
 * max{|A_ij|} for any j>i (1.2.1.3)
 */

Point Jacobi_find_ij(double ** matrix,int n) {
    int i,j;
    double max=0, num;
    Point p;

    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            num = fabs(matrix[i][j]);
            if (num > max){
                max = num;
                p.i = i;
                p.j = j;
            }
        }
    }
    return p;
}

/*
 * this function gets the difference between A and A' (1.2.1.5)
 */

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

/*
 * retrieves the sign of the number
 * returns 1 if num >= 0 or -1 if negative
 */

int Jacobi_get_sign(double num){
    if(num<0){
        return -1;
    } else{
        return 1;
    }
}

/*
 * sets the matrix A' according to the rules specified in section 1.2.1.6
 */

void Jacobi_set_matrix_A_new(double **A_new, double ** A,int n,int i,int j,double c,double s){
    int row;

    for(row=0;row<n;row++){
        if(row != i && row != j){
            A_new[row][i] = c * A[row][i] - s * A[row][j];
            A_new[row][j] = c * A[row][j] + s * A[row][i];
        }
    }

    A_new[i][i] = pow(c,2) * A[i][i] + pow(s,2) * A[j][j] - 2 * s * c * A[i][j];
    A_new[j][j] = pow(s,2) * A[i][i] + pow(c,2) * A[j][j] + 2 * s * c * A[i][j];
    A_new[i][j] = 0;
    A_new[j][i] = 0;
}

void Jacobi_set_matrix_P(double **matrix, int n, int i,int j, double c, double s){
    set_matrix_to_Identity(matrix,n);
    matrix[i][i] = c;
    matrix[j][j] = c;
    matrix[i][j] = s;
    matrix[j][i] = -s;
}

void init_Eigen_struct(Eigen *eigen, int n){
    int i;
    for(i=0;i<n;i++){
        eigen[i].vector = calloc(n, sizeof(double));
        error_occurred(eigen[i].vector == NULL);
    }
}

void free_Eigen_struct(Eigen *eigen, int n){
    int i;
    for(i=0;i<n;i++){
        free(eigen[i].vector);
    }
    free(eigen);
}

/*
 * this function sets the eigenvalues and eigenvectors
 * the vectors matrix is the V matrix when each column is eigenvector
 * the values matrix is the A matrix when each A[i][i] is eigenvalue
 */

void Jacobi_set_Eigen(Eigen *eigen, int n, double **vectors, double **values){
    int i,j;
    for(j=0;j<n;j++){
        eigen[j].value = values[j][j];
        for(i=0;i<n;i++){
            eigen[j].vector[i] = vectors[i][j]; /* transfer the column in vectors to row */
        }
    }
}

void transpose_matrix(double **transposed,double ** matrix, int rows, int columns){
    int i,j;
    for(i=0; i < rows; i ++){
        for(j=0; j< columns; j++){
            transposed[j][i] = matrix[i][j];
        }
    }
}

void Jacobi_algorithm(double **laplacian, int n, Eigen *eigen){
    double diff=1,theta,t,c,s;
    double epsilon = pow(10,-5);
    int num_iter=100;
    double ** V = allocate_array_2d(n,n);
    double ** A = allocate_array_2d(n,n);
    double ** A_new = allocate_array_2d(n,n);
    double ** P_matrix = allocate_array_2d(n,n);
    double ** P_matrix_transposed = allocate_array_2d(n,n);
    double ** result = allocate_array_2d(n,n);
    set_equal_array_2d(A,laplacian,n,n);
    set_equal_array_2d(A_new,laplacian,n,n);

    set_matrix_to_Identity(V,n);

    while(num_iter > 0 && diff >= epsilon){
        Point p = Jacobi_find_ij(A,n); /* 1.2.1.3 */
        /* 1.2.1.4 */
        int i = p.i;
        int j = p.j;
        theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
        t = Jacobi_get_sign(theta) / (fabs(theta) + sqrt(pow(theta,2) + 1 ));
        c = 1 / sqrt(pow(t,2) + 1);
        s = t * c;

        zero_array_2d(P_matrix,n,n);
        zero_array_2d(result,n,n);

        Jacobi_set_matrix_P(P_matrix,n,i,j,c,s); /* 1.2.1.2 */
        matrix_multiplication(V,n,n,P_matrix,n,result); /* 1.2.1.1.e */
        set_equal_array_2d(V,result,n,n);

        transpose_matrix(P_matrix_transposed,P_matrix,n,n);
        matrix_multiplication(P_matrix_transposed,n,n,A,n,result); /* 1.2.1.1.e */
        matrix_multiplication(result,n,n,P_matrix,n,A_new); /* 1.2.1.1.e */

        diff = Jacobi_find_diff_off(A,A_new,n); /* 1.2.1.5 */
        num_iter--;

        set_equal_array_2d(A,A_new,n,n);
    }

    Jacobi_set_Eigen(eigen,n,V,A);
    set_equal_array_2d(laplacian,V,n,n);

    free_array_2d(V,n);
    free_array_2d(A,n);
    free_array_2d(A_new,n);
    free_array_2d(P_matrix,n);
    free_array_2d(result,n);
}

void sort_eigen(Eigen *eigen, int n){
    int i,j;
    Eigen temp;

    for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (eigen[i].value > eigen[j].value)
            {
                temp =  eigen[i];
                eigen[i] = eigen[j];
                eigen[j] = temp;
            }
        }
    }
}
/*
 * the function gets Eigen array (1.3)
 * it calls the sort function and performs the heuristic
 */

int calculate_eigengap_heuristic(Eigen *eigens ,int n){
    int i,argmax;
    double max=0;
    double *eigengap = calloc(n, sizeof(double));;
    error_occurred(eigengap == NULL);

    sort_eigen(eigens,n);

    for(i=0;i<n-1;i++){
        eigengap[i] = eigens[i+1].value - eigens[i].value;
    }

    for(i=0;i < n/2;i++){
        if(eigengap[i] > max){
            max = eigengap[i];
            argmax = i;
        }
    }

    free(eigengap);
    /* TODO: verify this */
    return argmax + 1; /* because according to the presentation we calculate the eigenvalues from 1 to n */
}

/*
 * U is a matrix nxk, and we set each of the first k eigenvectors to be a column in U matrix
 */

void set_U_matrix(double ** U,Eigen *eigen, int n, int k){
    int i,j;

    for (i = 0; i < k; ++i){
        for (j = 0; j < n; ++j){
            U[j][i] = eigen[i].vector[j];
        }
    }
}

void get_first_k_columns(double **to,double **from, int n, int k){
    int i,j;

    for (i = 0; i < k; ++i){
        for (j = 0; j < n; ++j){
            to[j][i] = from[j][i];
        }
    }
}

/*
 * the function return the squared sum of the column
 */

double get_squared_sum_of_row(double ** matrix,int n, int r){
    int c;
    double sum=0;
    for(c=0;c<n;c++){
        sum += pow(matrix[r][c],2);
    }
    return sum;
}

/*
 * the function takes the U matrix and creates from it T matrix
 */

void calculate_T_matrix(double ** T,double ** U,int n,int k){
    int i,j;

    for (i = 0; i < k; ++i){
        for (j = 0; j < n; ++j){
            T[i][j] = U[i][j] / sqrt(get_squared_sum_of_row(U,n,i));
        }
    }
}

void print_eigenvalues(Eigen *eigen, int num_rows){
    int i;
    for(i=0;i<num_rows-1;i++){
        printf("%.4f", eigen[i].value);
        printf(",");
    }
    printf("%.4f\n", eigen[i].value);
}
