#include <stdio.h>
#include <spkmeans.c>

int main(int argc, char **argv) {

    int k, max_iter, num_rows, d = 1;
    int *num_elements_in_cluster;
    char *input_file, *output_file;
    double **data_points, **centroids, **new_centroids, **clusters;

    input_valid(argc == 2);
    input_file = argv[1];

    init_data_frame(input_file,data_points);

    return 0;
}

void init_data_frame(char *input_file,double **data_points) {
    FILE *ifp = NULL;
    int i = 0, j;
    char line[1024];

    ifp = fopen(input_file, "r");
    error_occured(ifp == NULL);

    find_vectors_len(ifp);
    find_d_of_vector(ifp);

    data_points = allocate_array_2d(num_rows, d);
    error_occured(data_points == NULL);

    while (fgets(line, sizeof line, ifp) != NULL) {
        line[strlen(line) - 1] = '\0';
        for (j = 0; j < d; j++) {
            if (j == 0) {
                data_points[i][j] = atof(strtok(line, ","));
            } else {
                data_points[i][j] = atof(strtok(NULL, ","));
            }
        }
        i++;
    }
    fclose(ifp);
}
