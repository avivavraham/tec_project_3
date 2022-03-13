//
// Created by Yonatan Romm on 26/02/2022.
//

#ifndef TEC_PROJECT_3_MAIN_H
#define TEC_PROJECT_3_MAIN_H

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void input_valid(int condition);
int find_vectors_len(FILE *fp);
int find_d_of_vector(FILE *fp);
void init_data_frame(const int *k,const int *d,const int *num_rows, char *input_file,double **data_points);

#endif //TEC_PROJECT_3_MAIN_H
