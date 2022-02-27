//
// Created by Yonatan Romm on 26/02/2022.
//

#ifndef TEC_PROJECT_3_MAIN_H
#define TEC_PROJECT_3_MAIN_H

void input_valid(int condition);
int find_vectors_len(FILE *fp);
int find_d_of_vector(FILE *fp);
void init_data_frame(const int *k,int *d,int *num_rows, char *input_file,double **data_points);

#endif //TEC_PROJECT_3_MAIN_H
