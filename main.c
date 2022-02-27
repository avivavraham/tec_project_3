#include <spkmeans.h>
#include "nsclustering.h"
#include <main.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {

    int *k,*d,*num_rows;
    char *input_file, *goal;
    double **data_points;

    input_valid(argc == 3);
    goal = argv[1];
    input_file = argv[2];



    init_data_frame(k,d,num_rows,input_file,data_points);

    if (strcmp(goal, "wam") == 0){

    }
    else if (strcmp(goal, "ddg") == 0){

    }
    else if (strcmp(goal, "lnorm") == 0){

    }
    else if (strcmp(goal, "jacobi") == 0){

    } else{
        input_valid(0);
    }

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

void init_data_frame(const int *k,int *d,int *num_rows, char *input_file,double **data_points) {
    FILE *ifp = NULL;
    int i = 0, j;
    char line[1024];

    ifp = fopen(input_file, "r");
    error_occurred(ifp == NULL);

    *num_rows = find_vectors_len(ifp);
    input_valid(*k < *num_rows);
    *d = find_d_of_vector(ifp);

    data_points = allocate_array_2d(*num_rows, *d);
    error_occurred(data_points == NULL);

    while (fgets(line, sizeof line, ifp) != NULL) {
        line[strlen(line) - 1] = '\0';
        for (j = 0; j < *d; j++) {
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
