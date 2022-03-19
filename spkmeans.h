#ifndef CAP_H_
#define CAP_H_

#include "libraries.h"

int find_closets_cluster(int k,int d,double *data_point, double **centroids);
void algorithm(int k, int d, int num_rows, int max_iter, double **data_points);
void free_array_2d(double **arr, int r);
void set_equal_array_2d(double **new, double **current, int r, int c);
void zero_array_2d(double **arr, int r, int c);
void calculate_new_centroids(int k,int d,const double *num_elements_in_cluster,double** new_centroids, double **clusters);
double **allocate_array_2d(int r, int c);
double calculate_distance_squared(double *a, double *b, int d);
void error_occurred(int condition);
void set_clusters(int k,int d, int num_rows, double *num_elements_in_cluster,double **centroids, double **clusters, double **data_points);
void init_centroids(int k,int d,double **centroids, double **data_points);

#endif
