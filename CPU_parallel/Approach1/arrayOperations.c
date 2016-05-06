/*
 * Author: Malin Prematilake
 * Date: 22.03.2016 17:23:22
 * Version: 1.0
 * Description: This code contains all the operations related with arrays
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "arrayOperations.h"

//prints array of doubles
void printArrayDouble(double *array, int length){
	int i;
	
	for (i=0;i<length;i++){
		double val = array[i];
		printf("%lf  ",val);
	}
	printf("\n");
}

//prints array of ints
void printArrayInt(int *array, int length){
	int i;
	
	for (i=0;i<length;i++){
		int val = array[i];
		printf("%d  ",val);
	}
	printf("\n");
}
/*
//prints array of bytes
void printArrayByte(byte *array, int length){
	int i;
	
	for (i=0;i<length;i++){
		int val = (int)array[i];
		printf("%d  ",val);
	}
	printf("\n");
}*/

int isInArray(double *array, double key, int length){
	
	int i;
	
	for(i=0;i<length;i++){
		if (array[i] == key) 
			return 1;
	}
	return 0;
}

int totalUniqueVals(double *array, int length){
	
	double *vals = (double *)malloc(sizeof(double)*length);
	
	int i,j=0;
	
	for(i=0; i<length; i++){
		if(isInArray(vals, array[i], i))
			continue;		
		else{
			vals[j] = array[i];
			j++;
		}
	}
	free(vals);
	return j;
}

int normaliseArray(double *doubleArray, int *intArray, int length){

	int minVal = 0;
	int maxVal = 0;
	
	int currentValue;
	int i;

	if (length > 0){
		minVal = (int) floor(doubleArray[0]);
		maxVal = (int) floor(doubleArray[0]);

		for (i = 0; i < length; i++){
			
			currentValue = (int) floor(doubleArray[i]);
			intArray[i] = currentValue;

			if (currentValue < minVal)
			minVal = currentValue;

			else if (currentValue > maxVal)
			maxVal = currentValue;

			}/*for loop over vector*/

		for (i = 0; i < length; i++)
			intArray[i] = intArray[i] - minVal;

		maxVal = (maxVal - minVal) + 1;
	}
	return maxVal;
}

void riseUp(double *array, double factor, int length){
	int i;
	
	for(i=0;i<length;i++){
		array[i] = array[i]*factor;
	}
}
/*
int main(){
	
	double fd[] = {2.3,34.2,56.7,45.0,45.0};
	int *te = malloc(sizeof(int)*5);
	
	int y = totalUniqueVals(fd, 5);
	printf("--%d \n",y);
	printArrayDouble(fd, 5);
	
	int yy = normaliseArray(fd, te, 5);
	
	printf("----%d \n",yy);
	printArrayInt(te, 5);
	
	riseUp(fd,10.0,5);
	printArrayDouble(fd,5);
	
	//--4 
	//2.300000  34.200000  56.700000  45.000000  45.000000  
	//----55 
	//0  32  54  43  43
	
	
	free(te);
	return 0;
}
*/
