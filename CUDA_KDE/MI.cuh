/*
 * Author: Malin Prematilake
 * Date: 04.04.2016 14:43:52
 * Version: 1.0
 * Description: This code contains all the definitions related with MI.c
*/
#ifndef MI_H
#define MI_H

__device__ int normaliseArray(double *doubleArray, int *intArray, int length);

__device__ double MI(  double *MIPointer, double *dataArray, double *targetArray, int start, int Length,
					 int *normalisedArray1, int *normalisedArray2,
					 int *firstStateProbs, int *secondStateProbs, int *jointStateProbs);
#endif
