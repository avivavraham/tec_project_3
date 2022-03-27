
#include "main.h"

int main(int argc, char **argv) {

    int *d,*num_rows;
    char *input_file, *goal;
    double **data_points;
    Eigen *eigen;

    input_valid(argc == 3);
    goal = argv[1];
    input_file = argv[2];

    init_data_frame(d,num_rows,input_file,data_points);

    if (strcmp(goal, "wam") == 0){ //input datapoins
        double** weighted_matrix = allocate_array_2d(*num_rows, *num_rows);
        calculate_weighted_matrix(weighted_matrix,data_points,*num_rows,*d);
        print_matrix(weighted_matrix,*num_rows,*num_rows);
        free_array_2d(weighted_matrix,*num_rows);
    }
    else if (strcmp(goal, "ddg") == 0){ //input datapoins
        double** weighted_matrix = allocate_array_2d(*num_rows, *num_rows);
        calculate_weighted_matrix(weighted_matrix,data_points,*num_rows,*d);

        double** diagonal_degree_matrix = allocate_array_2d(*num_rows, *num_rows);
        calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,*num_rows);
        print_matrix(diagonal_degree_matrix, *num_rows, *num_rows);
        free_array_2d(diagonal_degree_matrix,*num_rows);

        free_array_2d(weighted_matrix, *num_rows);
    }
    else if (strcmp(goal, "lnorm") == 0){ //input datapoins
        double** weighted_matrix = allocate_array_2d(*num_rows, *num_rows);
        calculate_weighted_matrix(weighted_matrix,data_points,*num_rows,*d);
        double** diagonal_degree_matrix = allocate_array_2d(*num_rows, *num_rows);
        calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,*num_rows);

        double** lnorm_matrix = allocate_array_2d(*num_rows, *num_rows);
        calculate_lnorm_matrix(lnorm_matrix,weighted_matrix,diagonal_degree_matrix,*num_rows);
        print_matrix(lnorm_matrix, *num_rows, *num_rows);
        free_array_2d(lnorm_matrix, *num_rows);

        free_array_2d(diagonal_degree_matrix,*num_rows);
        free_array_2d(weighted_matrix, *num_rows);

    }
    else if (strcmp(goal, "jacobi") == 0){ //input symmetric matrix, what is the output?

        input_valid(*num_rows == *d); //validates that this is squared matrix
        eigen = calloc(*num_rows, sizeof (Eigen *));
        error_occurred(eigen == NULL);
        init_Eigen_struct(eigen,*num_rows);
        Jacobi_algorithm(data_points,*num_rows,eigen);

        print_eigenvalues(eigen,*num_rows);

        double** eigenvectors = allocate_array_2d(*num_rows, *num_rows);
        // U matrix takes the k first vectors as columns,
        // so if we insert k=n we get matrix that each column of the n columns is eigenvector
        set_U_matrix(eigenvectors,eigen,*num_rows,*num_rows);
        print_matrix(eigenvectors,*num_rows,*num_rows);

    } else{
        free_array_2d(data_points,*num_rows);
        input_valid(0);
    }

    free_array_2d(data_points,*num_rows);

    return 0;
}

void input_valid(int condition) {
    if (condition == 0) {
        printf("Invalid Input!\n");
        exit(1);
    }
}

int find_vectors_len(FILE *fp) {
    char line[1024];
    int num_rows=0;
    while (fgets(line, sizeof line, fp) != NULL) {
        num_rows++;
    }
    fseek(fp, 0, SEEK_SET);
    return num_rows;
}

int find_d_of_vector(FILE *fp) {
    char line[1024];
    int i = 0;
    int d=0;
    fseek(fp, 0, SEEK_SET);
    if (fgets(line, sizeof line, fp) != NULL) {
        line[i] = line[0];
        while (line[i] != '\n') {
            if (line[i] == ',') {
                d++;
            }
            i++;
        }
    }
    fseek(fp, 0, SEEK_SET);
    return d;
}

void init_data_frame(const int *d,const int *num_rows, char *input_file,double **data_points) {
    FILE *ifp = NULL;
    int i = 0, j;
    char line[1024];
    int local_num_rows, local_d;

    ifp = fopen(input_file, "r");
    error_occurred(ifp == NULL);

    local_num_rows = find_vectors_len(ifp);
    num_rows = &local_num_rows;
    local_d = find_d_of_vector(ifp);
    d = &local_d;

    data_points = allocate_array_2d(*num_rows, *d);

    while (fgets(line, sizeof line, ifp) != NULL) {
        line[strlen(line) - 1] = '\0';
        for (j = 0; j < *d; j++) {
            if (j == 0) {
                data_points[i][j] = atof(strtok(line, ",")); //TODO: check the warning about atof
            } else {
                data_points[i][j] = atof(strtok(NULL, ","));
            }
        }
        i++;
    }
    fclose(ifp);
}
