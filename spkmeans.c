
#include "spkmeans.h"

int main(int argc, char **argv) {

    int d,num_rows;
    char *input_file, *goal;
    double **data_points;
    Eigen *eigen;
    double **diagonal_degree_matrix, **weighted_matrix, **lnorm_matrix, **eigenvectors;
    FILE *ifp;

    input_valid(argc == 3);
    goal = argv[1];
    input_file = argv[2];

    ifp = fopen(input_file, "r");
    error_occurred(ifp == NULL);

    num_rows = find_vectors_len(ifp);
    d = find_d_of_vector(ifp);

    data_points = allocate_array_2d(num_rows, d);

    init_data_frame(ifp, d,data_points);

    if (strcmp(goal, "wam") == 0){
        weighted_matrix = allocate_array_2d(num_rows, num_rows);
        calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);
        print_matrix(weighted_matrix,num_rows,num_rows);
        free_array_2d(weighted_matrix,num_rows);
    }
    else if (strcmp(goal, "ddg") == 0){
        weighted_matrix = allocate_array_2d(num_rows, num_rows);
        calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);

        diagonal_degree_matrix = allocate_array_2d(num_rows, num_rows);
        calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,num_rows);
        print_matrix(diagonal_degree_matrix, num_rows, num_rows);
        free_array_2d(diagonal_degree_matrix,num_rows);

        free_array_2d(weighted_matrix, num_rows);
    }
    else if (strcmp(goal, "lnorm") == 0){
        weighted_matrix = allocate_array_2d(num_rows, num_rows);
        calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);
        diagonal_degree_matrix = allocate_array_2d(num_rows, num_rows);
        calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,num_rows);

        lnorm_matrix = allocate_array_2d(num_rows, num_rows);
        calculate_lnorm_matrix(lnorm_matrix,weighted_matrix,diagonal_degree_matrix,num_rows);
        print_matrix(lnorm_matrix, num_rows, num_rows);
        free_array_2d(lnorm_matrix, num_rows);

        free_array_2d(diagonal_degree_matrix,num_rows);
        free_array_2d(weighted_matrix, num_rows);
    }
    else if (strcmp(goal, "jacobi") == 0){
        input_valid(num_rows == d);
        eigen = calloc(num_rows, sizeof (Eigen));
        error_occurred(eigen == NULL);
        init_Eigen_struct(eigen,num_rows);
        Jacobi_algorithm(data_points,num_rows,eigen);
        print_eigenvalues(eigen,num_rows);

        eigenvectors = allocate_array_2d(num_rows, num_rows);
        set_U_matrix(eigenvectors,eigen,num_rows,num_rows);
        print_matrix(eigenvectors,num_rows,num_rows);

        free_Eigen_struct(eigen,num_rows);

    } else{
        free_array_2d(data_points,num_rows);
        input_valid(0);
    }

    free_array_2d(data_points,num_rows);

    return 0;
}
