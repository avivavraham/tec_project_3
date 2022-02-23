#ifndef CAP_H_
#define CAP_H_

int find_closets_cluster(double *data_point);
void algorithm();
void free_array_2d(double **arr, int r);
void set_equal_array_2d(double **new, double **current, int r, int c);
void zero_array_2d(double **arr, int r, int c);
void calculate_new_centroids();
double **allocate_array_2d(int r, int c);
double calculate_distance_squared(double *a, double *b);
void error_occurred(int condition);
void set_clusters();

#endif
