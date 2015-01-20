/******************************************************************************
* This file is part of CHAOTICB.
* Copyright (C) <2012-2014> <Xiaocan Li> <xl0009@uah.edu>
* 
* CHAOTICB is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* CHAOTICB is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with CHAOTICB.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "quick_sort.h"

/******************************************************************************
 * Main function of quick sorting and keep tracking the indices of the array.
 *
 * Input:
 *  data_array: input array.
 *  data_size: the size of input array.
 * Output:
 *  Sorted data_array.
******************************************************************************/
void quicksort(double *data_array, int *index_array, int data_size)
{
    partition(data_array, index_array, 0, data_size-1);
}

/******************************************************************************
 * Swap the values of two int or double variables.
******************************************************************************/
void swap_int(int *a, int *b) 
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
void swap_double(double *a, double *b) 
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

/******************************************************************************
 * Get the median of three numbers. And the index of the middle of the those
 * three numbers will be updated.
 *
 * Input:
 *  a, b, c: 3 numbers.
 *  index_a, index_b, index_c: indices of an array.
 * Output:
 *  index_b will be updated.
 * Return:
 *  b: the median of the 3 numbers.
******************************************************************************/
double median(int index_a, int *index_b, int index_c, 
        double a, double b, double c)
{
    if(b < a) {
        swap_double(&a, &b);
        swap_int(&index_a, index_b);
    }
    if(c < a) {
        swap_double(&a, &c);
        swap_int(&index_a, &index_c);
    }
    if(c < b) {
        swap_double(&b, &c);
        swap_int(index_b, &index_c);
    }
    return b;
}

/******************************************************************************
 * Divide and conquer algorithm of Quicksort. No sub-array is needed to copy
 * part of the input array.
 *
 * Input:
 *  input: the array to sorted.
 *  index: the change of the indices of the array will be tracked.
 *  left: the index of the leftmost element of the "sub-array".
 *  right: the index of the rightmost element of the "sub-array".
 * Output:
 *  Updated input and index.
******************************************************************************/
void partition(double *input, int *index, int left, int right)
{
    double pivot;
    int i, j, size;
    int pivotIndex;

    size = right - left + 1;
    i = left + 1;
    if(size % 2 == 0) {
        pivotIndex = left + size / 2 - 1;
    }
    else {
        pivotIndex = left + size / 2;
    }
    pivot = median(left, &pivotIndex, right, 
            input[left],input[pivotIndex],input[right]);

    swap_double(&input[pivotIndex], &input[left]);
    swap_int(&index[pivotIndex], &index[left]);

    for(j = left+1; j <= right; j++) {
        if(input[j] < pivot) {
            swap_double(&input[i],&input[j]);
            swap_int(&index[i],&index[j]);
            ++i;
        }
    }
    swap_double(&input[left],&input[i-1]);
    swap_int(&index[left],&index[i-1]);
    if(left < i-1) {
        partition(input, index, left, i-2);
    }
    if(right > i-1) {
        partition(input, index, i, right);
    }
}
