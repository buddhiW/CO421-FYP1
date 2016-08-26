/*
 * Author: Malin Prematilake
 * Date: 22.03.2016 17:23:22
 * Version: 1.0
 * Description: This code contains all the operations related with arrays
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "helpers.cuh"
#include "arrayOperations.cuh"

//number of sample points in one encryption
/*#define SAMPLEPOINTS 100000
#define SAMPLES 200

//defined for 128 bit AES
#define KEYBYTES 16
#define KEYS 256*/

//#define WIDTH 5.0 //No. of elements in a bin

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
void printArrayInt(unsigned int *array, int length){
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

__global__ void normaliseWaveDataKernel(double *doubleArray, unsigned int *intArray, unsigned int* firstNumStates){

	/* BlockSize = (1000,1,1);
	 * GridSize = (100, 1, 1);
	 * 
	 * doubleArray = input
	 * intArray = result
	 * */

	int samplePoint = blockIdx.x*blockDim.x + threadIdx.x; //Directly gives sample point
	
	//int y = blockDim.y; //Directly gives sample point
	
	int minVal = 0;
	int maxVal = 0;
	
	int currentValue;
	int i;
	int length = SAMPLES;
	
	int startPosition = samplePoint*length;
	
	int bins;
	//int binNumber;
	int temp;

	if (samplePoint < SAMPLEPOINTS){
		
		//minVal = (int) floor(doubleArray[0]); //??
		//maxVal = (int) floor(doubleArray[0]); //??

		minVal = (int) floor(doubleArray[startPosition]);
		maxVal = (int) floor(doubleArray[startPosition]);

		for (i = startPosition; i < startPosition+length; i++){
			
			currentValue = (int) floor(doubleArray[i]);
			intArray[i] = currentValue;		//?????? index
			
			if (currentValue < minVal)
			minVal = currentValue;

			else if (currentValue > maxVal)
			maxVal = currentValue;

		}/*for loop over vector*/

		bins = ceil(((maxVal-minVal)+1)/(double)WIDTH); //No. of bins
		maxVal = (maxVal - minVal)+1;
		//for (i = startPosition; i < startPosition+length; i++)
		for (i = startPosition; i < startPosition+length; i++){
			//temp2 = intArray[i];
			temp = intArray[i] - minVal;
			//intArray[i] = temp/(maxVal/BINS);
			intArray[i] = temp/WIDTH;
			//if(samplePoint == 99999)
				//printf("actual: %d  min: %d temp: %d bin: %d\n", temp2, minVal, temp, intArray[i]);	
		}
				
		firstNumStates[samplePoint] = bins;
				
		//__syncthreads();
	}
	
}

/* needs error checking, badly */
/*int *buildHist(int bins, double min, double max, int n, double *data){
   double *hist=malloc(bins*sizeof(int));
   if (hist == NULL) return hist;
   for (int i=0; i<n; ++i){
      int bin=int( (data[i]-min)/((max-min)/(bins)) );
      if ( (bin>=0) && (bin<n) ) hist[bin]++;
   }
   return hist;
}*/




__global__ void normaliseHammingKernel(double *doubleArray, unsigned int *intArray, unsigned int* secondNumStates){
	
	/*	GridSize = (256, 1, 1)
	 *	BlockSize = (16, 1, 1)
	 * 
	 * */
	 
	//int key =  blockIdx.x*blockDim.x + threadIdx.x;
	int key =  blockIdx.x; // There are 256 blocks 
	int keyByte = threadIdx.x; 
	
	int minVal = 0;
	int maxVal = 0;
	
	int currentValue;
	int i;

	int startPosition = key*SAMPLES*KEYBYTES + keyByte*SAMPLES;
	int length = SAMPLES;

	if ((key < KEYS) && (keyByte < KEYBYTES)){
		
		minVal = (int) floor(doubleArray[startPosition]);
		maxVal = (int) floor(doubleArray[startPosition]);

		for (i = startPosition; i < startPosition+length; i++){
			
			currentValue = (int) floor(doubleArray[i]);
			intArray[i] = currentValue;
			
			//intArray[i] = key;
			
			if (currentValue < minVal)
			minVal = currentValue;

			else if (currentValue > maxVal)
			maxVal = currentValue;

		}//for loop over vector

		for (i = startPosition; i < startPosition+length; i++)
			intArray[i] = intArray[i] - minVal;

		maxVal = (maxVal - minVal) + 1;
		
		secondNumStates[key*KEYBYTES+keyByte] = maxVal;
		
		//__syncthreads();
	}
	
	
}



void riseUp(double *array, double factor, int length){
	int i;
	
	for(i=0;i<length;i++){
		array[i] = array[i]*factor;
	}
}

void printMatInt(unsigned int *mat, int row, int col){
	int i,j;
	
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			int post = i*col+j;
			//printf("%02X  ",mat[post]);
			printf("%d  ",mat[post]);
		}
		printf("\n");
	}
	return;
}

void printMatDouble(double *mat, int row, int col){
	int i,j;
	
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			int post = i*col+j;
			printf("%lf  ",mat[post]);
		}
		printf("\n");
	}
	return;
}
/*
int main(){
	
	//double fd[] = {2.3,34.2,56.7,45.0,45.0};
	unsigned int *te = malloc(sizeof(unsigned int)*40);
	
	int i;
	for (i=0;i<40;i++){
		te[i] = i;
	}
	printMatInt(te,8,5);
	free(te);
	return 0;
}
*/
