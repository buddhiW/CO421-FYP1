/*
 * Author: Malin Prematilake
 * Date: 23.03.2016 00:28:11
 * Version: 1.0
 * Description: This code contains all the definitions related with arrays
*/
#ifndef ARRAYOPS_H
#define ARRAYOPS_H

void printArrayDouble(double *array, int length);
void printArrayInt(int *array, int length);
//void printArrayByte(byte *array, int length);
int isInArray(double *array, double key, int length);
int totalUniqueVals(double *array, int length);
int normaliseArray(double *doubleArray, int *intArray, int length);
void riseUp(double *array, double factor, int length);
#endif
