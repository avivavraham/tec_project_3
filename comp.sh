#!/bin/bash

gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c nsclustering.c helper.c -lm -o spkmeans