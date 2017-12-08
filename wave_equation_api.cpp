/*
Copyright (C) 2017 the team of Jim00000, ActKz and pityYo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/**
 * @file        wave_equation_api.c
 * @author      Jim00000
 * @date        12.6.2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "wave_equation_api.h"

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

void c_sequential_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)
{
    // Check that data is not a null pointer
    if(data == NULL) {
        fprintf(stderr, "[ERROR] array 'data' is null \n");
        exit(EXIT_FAILURE);
    }

    if(olddata == NULL) {
        fprintf(stderr, "[ERROR] array 'olddata' is null \n");
        exit(EXIT_FAILURE);
    }

    if(newdata == NULL) {
        fprintf(stderr, "[ERROR] array 'newdata' is null \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that row_size is a positive number
    if(row_size <= 0) {
        fprintf(stderr, "[ERROR] row_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that col_size is a positive number
    if(col_size <= 0) {
        fprintf(stderr, "[ERROR] col_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    // Central part
    for(int i = 1; i < row_size - 1; i++) {
        for(int j = 1; j < col_size - 1; j++) {
            double potential = data[(i + 1) * col_size + j] + data[(i - 1) * col_size + j] + data[i * col_size + j + 1] +
                               data[i * col_size + j - 1] - 4 * data[i * col_size + j];
            newdata[i * col_size + j] = ( pow(C * dt, 2) * potential * 2 + 4 * data[i * col_size + j] - olddata[i * col_size + j] *
                                          (2 - K * dt) ) / (2 + K * dt);
        }
    }

    // Four edges
    for(int i = 1; i < row_size - 1; i++) {
        double P1 = data[(i + 1) * col_size] + data[(i - 1) * col_size] + data[i * col_size + 1] - 3 * data[i * col_size];
        double P2 = data[(i + 1) * col_size + col_size - 1] + data[(i - 1) * col_size + col_size - 1] +
                    data[i * col_size + col_size - 2] - 3 * data[i * col_size + col_size - 1];
        double P3 = data[col_size + i] + data[i + 1] + data[i - 1] - 3 * data[i];
        double P4 = data[(row_size - 2) * col_size + i] + data[(row_size - 1) * col_size + i + 1] +
                    data[(row_size - 1) * col_size + i - 1] - 3 * data[(row_size - 1) * col_size + i];
        newdata[i * col_size] = ( pow(C * dt, 2) * P1 * 2 + 4 * data[i * col_size] - olddata[i * col_size] *
                                  (2 - K * dt) ) / (2 + K * dt);
        newdata[i * col_size + col_size - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * data[i * col_size + col_size - 1] -
                                               olddata[i * col_size + col_size - 1] * (2 - K * dt) ) / (2 + K * dt);
        newdata[i] = ( pow(C * dt, 2) * P3 * 2 + 4 * data[i] - olddata[i] * (2 - K * dt) ) / (2 + K * dt);
        newdata[(row_size - 1) * col_size + i] = ( pow(C * dt, 2) * P4 * 2 + 4 * data[(row_size - 1) * col_size + i] -
                olddata[(row_size - 1) * col_size + i] * (2 - K * dt) ) / (2 + K * dt);
    }

    // Four corners
    double P1 = data[col_size] + data[1] - 2 * data[0];
    double P2 = data[col_size + col_size - 1] + data[col_size - 2] - 2 * data[col_size - 1];
    double P3 = data[(row_size - 2) * col_size] + data[(row_size - 1) * col_size + 1] - 2 * data[(row_size - 1) * col_size];
    double P4 = data[(row_size - 2) * col_size + col_size - 1] + data[(row_size - 1) * col_size + col_size - 2] - 2 *
                data[(row_size - 1) * col_size + col_size - 1];
    newdata[0] = ( pow(C * dt, 2) * P1 * 2 + 4 * data[0] - olddata[0] * (2 - K * dt) ) / (2 + K * dt);
    newdata[col_size - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * data[col_size - 1] - olddata[col_size - 1] * (2 - K * dt) ) /
                            (2 + K * dt);
    newdata[(row_size - 1) * col_size] = ( pow(C * dt, 2) * P3 * 2 + 4 * data[(row_size - 1) * col_size] - olddata[(row_size - 1)
                                           * col_size] * (2 - K * dt) ) / (2 + K * dt);
    newdata[(row_size - 1) * col_size + col_size - 1] = ( pow(C * dt, 2) * P4 * 2 +
            4 * data[(row_size - 1) * col_size + col_size - 1] - olddata[(row_size - 1) * col_size + col_size - 1] * (2 - K * dt) )
            / (2 + K * dt);
}

void c_openmp_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)
{
    // Check that data is not a null pointer
    if(data == NULL) {
        fprintf(stderr, "[ERROR] array 'data' is null \n");
        exit(EXIT_FAILURE);
    }

    if(olddata == NULL) {
        fprintf(stderr, "[ERROR] array 'olddata' is null \n");
        exit(EXIT_FAILURE);
    }

    if(newdata == NULL) {
        fprintf(stderr, "[ERROR] array 'newdata' is null \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that row_size is a positive number
    if(row_size <= 0) {
        fprintf(stderr, "[ERROR] row_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that col_size is a positive number
    if(col_size <= 0) {
        fprintf(stderr, "[ERROR] col_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    // Central part
    #pragma omp parallel for shared(row_size, col_size, newdata, data, olddata, C, K, dt) collapse(2)
    for(int i = 1; i < row_size - 1; i++) {
        for(int j = 1; j < col_size - 1; j++) {
            double potential = data[(i + 1) * col_size + j] + data[(i - 1) * col_size + j] + data[i * col_size + j + 1] +
                               data[i * col_size + j - 1] - 4 * data[i * col_size + j];
            newdata[i * col_size + j] = ( pow(C * dt, 2) * potential * 2 + 4 * data[i * col_size + j] - olddata[i * col_size + j] *
                                          (2 - K * dt) ) / (2 + K * dt);
        }
    }

    // Four edges
    #pragma omp parallel for shared(row_size, col_size, newdata, data, olddata, C, K, dt)
    for(int i = 1; i < row_size - 1; i++) {
        double P1 = data[(i + 1) * col_size] + data[(i - 1) * col_size] + data[i * col_size + 1] - 3 * data[i * col_size];
        double P2 = data[(i + 1) * col_size + col_size - 1] + data[(i - 1) * col_size + col_size - 1] +
                    data[i * col_size + col_size - 2] - 3 * data[i * col_size + col_size - 1];
        double P3 = data[col_size + i] + data[i + 1] + data[i - 1] - 3 * data[i];
        double P4 = data[(row_size - 2) * col_size + i] + data[(row_size - 1) * col_size + i + 1] +
                    data[(row_size - 1) * col_size + i - 1] - 3 * data[(row_size - 1) * col_size + i];
        newdata[i * col_size] = ( pow(C * dt, 2) * P1 * 2 + 4 * data[i * col_size] - olddata[i * col_size] *
                                  (2 - K * dt) ) / (2 + K * dt);
        newdata[i * col_size + col_size - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * data[i * col_size + col_size - 1] -
                                               olddata[i * col_size + col_size - 1] * (2 - K * dt) ) / (2 + K * dt);
        newdata[i] = ( pow(C * dt, 2) * P3 * 2 + 4 * data[i] - olddata[i] * (2 - K * dt) ) / (2 + K * dt);
        newdata[(row_size - 1) * col_size + i] = ( pow(C * dt, 2) * P4 * 2 + 4 * data[(row_size - 1) * col_size + i] -
                olddata[(row_size - 1) * col_size + i] * (2 - K * dt) ) / (2 + K * dt);
    }

    // Four corners
    double P1 = data[col_size] + data[1] - 2 * data[0];
    double P2 = data[col_size + col_size - 1] + data[col_size - 2] - 2 * data[col_size - 1];
    double P3 = data[(row_size - 2) * col_size] + data[(row_size - 1) * col_size + 1] - 2 * data[(row_size - 1) * col_size];
    double P4 = data[(row_size - 2) * col_size + col_size - 1] + data[(row_size - 1) * col_size + col_size - 2] - 2 *
                data[(row_size - 1) * col_size + col_size - 1];
    newdata[0] = ( pow(C * dt, 2) * P1 * 2 + 4 * data[0] - olddata[0] * (2 - K * dt) ) / (2 + K * dt);
    newdata[col_size - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * data[col_size - 1] - olddata[col_size - 1] * (2 - K * dt) ) /
                            (2 + K * dt);
    newdata[(row_size - 1) * col_size] = ( pow(C * dt, 2) * P3 * 2 + 4 * data[(row_size - 1) * col_size] - olddata[(row_size - 1)
                                           * col_size] * (2 - K * dt) ) / (2 + K * dt);
    newdata[(row_size - 1) * col_size + col_size - 1] = ( pow(C * dt, 2) * P4 * 2 +
            4 * data[(row_size - 1) * col_size + col_size - 1] - olddata[(row_size - 1) * col_size + col_size - 1] * (2 - K * dt) )
            / (2 + K * dt);
}
