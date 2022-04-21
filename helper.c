#include "helper.h"

/*
this function handles error,
it checks the condition if it is true it will print "An Error Has Occurred" and exit
 */

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
    int d=1;
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

void init_data_frame(FILE *ifp, int d,double **data_points) {
    int i = 0, j;
    char line[1024];

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

void error_occurred(int condition) {
    if (condition == 1) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
}

/*
the function gets data points and calculates the closest cluster
by the smallest squared distance from a centroid
it returns the index of the cluster in the array
*/

int find_closets_cluster(int k,int d,double *data_point, double **centroids) {
    double *difference, *current_mu, sum, min_sum;
    int index = 0, i;

    difference = calloc(k, sizeof(int *));
    error_occurred(difference == NULL);
    for (i = 0; i < k; i++) {
        current_mu = centroids[i];
        sum = calculate_distance_squared(current_mu, data_point,d);
        difference[i] = sum;
    }
    min_sum = difference[0];
    for (i = 0; i < k; i++) {
        if (difference[i] < min_sum) {
            min_sum = difference[i];
            index = i;
        }
    }
    free(difference);
    return index;
}

/*
the function iterate through all the data points and assigns each to the closest cluster
it calculates the num elements in each cluster
*/

void set_clusters(int k,int d, int num_rows, double *num_elements_in_cluster,double **centroids, double **clusters, double **data_points) {
    int index, i, a;

    for (i = 0; i < k; i++) {
        num_elements_in_cluster[i] = 0;
    }
    for (i = 0; i < num_rows; i++) {
        index = find_closets_cluster(k,d,data_points[i],centroids);
        for (a = 0; a < d; a++) {
            clusters[index][a] += data_points[i][a];
        }
        num_elements_in_cluster[index]++;
    }
}

/*
the function calculate the centroids by the new clusters that was already calculated
*/
void calculate_new_centroids(int k,int d,const double *num_elements_in_cluster,double** new_centroids, double **clusters) {
    double *cluster;
    int len_cluster, i, j;

    for (i = 0; i < k; i++) {
        cluster = clusters[i];
        len_cluster = num_elements_in_cluster[i];
        for (j = 0; j < d; j++) {
            new_centroids[i][j] = cluster[j] / len_cluster;
        }
    }
}

/*
the function calculate the squared distance of two points
*/

double calculate_distance_squared(double *a, double *b, int d) {
    double distance = 0;
    int i = d;
    while (i--) {
        distance += pow(a[i] - b[i], 2);
    }
    return distance;
}

/*
pre-condition: centroids and data points arrays already initialized
each centroid represent cluster (the average of the cluster)
every iteration it calculates the clusters and then updates the centroids
when the difference squared distance of all the centroids is smaller than epsilon we stop and return the centroids
*/

void algorithm(int k, int d, int num_rows, int max_iter, double **data_points, double **centroids) {
    int i, j;
    double diff, sq_diff, max, epsilon = 0.001;
    int *num_elements_in_cluster;
    double *sum_diff_centroids;
    double **new_centroids,**clusters;

    num_elements_in_cluster = calloc(k, sizeof(int));
    error_occurred(num_elements_in_cluster == NULL);
    clusters = allocate_array_2d(k, d);
    new_centroids = allocate_array_2d(k, d);
    sum_diff_centroids = calloc(k, sizeof(double));
    error_occurred(sum_diff_centroids == NULL);
    while (max_iter > 0) {
        for (j = 0; j < k; j++) {
            sum_diff_centroids[j] = 0;
        }

        zero_array_2d(new_centroids, k, d);
        zero_array_2d(clusters, k, d);

        error_occurred(num_elements_in_cluster == NULL);

        set_clusters(k,d,num_rows,num_elements_in_cluster,centroids,clusters,data_points);
        calculate_new_centroids(k,d,num_elements_in_cluster,new_centroids,clusters);

        max_iter--;
        error_occurred(sum_diff_centroids == NULL);

        for (i = 0; i < k; i++) {
            diff = calculate_distance_squared(centroids[i], new_centroids[i],d);
            sq_diff = pow(diff, 0.5);
            sum_diff_centroids[i] = sq_diff;
        }
        max = sum_diff_centroids[0];
        for (i = 0; i < k; i++) {
            if (sum_diff_centroids[i] > max) {
                max = sum_diff_centroids[i];
            }
        }
        set_equal_array_2d(centroids, new_centroids, k, d);
        if (max <= epsilon) {
            max_iter = 0;
        }
    }
    free(sum_diff_centroids);
    free_array_2d(new_centroids, k);
    free_array_2d(clusters, k);
    free(num_elements_in_cluster);
}

/*
the function allocate memory for 2 dimensional array of type double
*/

double **allocate_array_2d(int r, int c) {
    double **arr;
    arr = calloc(r, sizeof(double *));
    error_occurred(arr == NULL);
    while (r--) {
        arr[r] = calloc(c, sizeof(double));
        error_occurred(arr[r] == NULL);
    }
    return arr;
}

/*
the function free memory for 2 dimensional array
*/

void free_array_2d(double **arr, int r) {
    while (r--) {
        free(arr[r]);
    }
    free(arr);
}

/*
the function gets two 2 dimensional arrays and num of rows and columns
it sets the first array to be equal to the second array
*/


void set_equal_array_2d(double **new, double **current, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            new[i][j] = current[i][j];
        }
    }
}

/*
the function gets 2 dimensional array and num of rows and columns
it sets all the values of the array to zero
*/

void zero_array_2d(double **arr, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            arr[i][j] = 0;
        }
    }
}
