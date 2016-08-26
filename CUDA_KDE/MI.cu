#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "arrayOperations.cuh"
#include "MI.cuh"
#include "helpers.cuh"

/*__global__ void MIKernel(unsigned int *firstNumStates, unsigned int *secondNumStates, 
							double *waveProbs, double *hammingProbs, double *jointProbs){
	
	int fNumState = firstNumStates[xx];
	int sNumState = secondNumStates[key*KEYBYTES + yy];
		
		int jNumState = fNumState*sNumState;
		double MI = 0.0;
		
		for (i=0; i<jNumState; i++){
			firstIndx = i%fNumState;
			secondIndx = i/fNumState;
			
			double c1 = waveProbs[xx*maxOfFirst + firstIndx];
			double c2 = hammingProbs[key*KEYBYTES*maxOfSecond + yy*maxOfSecond + secondIndx];
			double c3 = jointProbs[yy];
			if ((c3 > 0.0) && (c1 > 0.0) && (c2 > 0.0)){
					double temp = c3 * log10f(c3 / c1 / c2);
					MI = MI + temp;
			}
		}
		MI = MI/log10f(2.0);
		//MIvals[key*KEYBYTES*SAMPLEPOINTS + yy*SAMPLEPOINTS + xx] = MI;
		MIvals[yy*SAMPLEPOINTS + (xx+startPosition)] = MI; //MIvals is for one key
	
	
	
}*/


__device__ int normaliseArray(double *doubleArray, int *intArray, int length){

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

__device__ double MI(double *dataArray, double *targetArray, int startD, int startT, int Length,
					 int *normalisedArray1, int *normalisedArray2,
					 int *firstStateProbs, int *secondStateProbs, int *jointStateProbs){
	
	double MI = 0.0;
	int firstIndex=0,secondIndex=0;

	int firstNumStates = normaliseArray(dataArray,normalisedArray1,Length);
	int secondNumStates = normaliseArray(targetArray,normalisedArray2,Length);
	int jointNumStates = firstNumStates * secondNumStates;
		
	double length = Length;
	
	int i;
	
	for (i = 0; i < Length; i++){
		firstStateProbs[normalisedArray1[i]] += 1.0;
		secondStateProbs[normalisedArray2[i]] += 1.0;
		jointStateProbs[normalisedArray2[i] * firstNumStates + normalisedArray1[i]] += 1.0;
	}

	for (i = 0; i < firstNumStates; i++){
		firstStateProbs[i] = firstStateProbs[i] / length;
	}
	
	for (i = 0; i < secondNumStates; i++){
		secondStateProbs[i] = secondStateProbs[i] / length;
	}
	
	for (i = 0; i < jointNumStates; i++){
		jointStateProbs[i] = jointStateProbs[i] / length;
	}
	
	/*normalisedArray1 = NULL;
	normalisedArray2 = NULL;
	
	firstStateCounts = NULL;
	secondStateCounts = NULL;
	jointStateCounts = NULL;*/
	
	/*Probability calculatoins ends*/
	
	/* MI calculation begins
	** I(X;Y) = sum sum p(xy) * log (p(xy)/p(x)p(y))
	*/
	
	for (i = 0; i < jointNumStates; i++){
		
		firstIndex = i % firstNumStates;
		secondIndex = i / firstNumStates;
		
		if ((jointStateProbs[i] > 0.0) && (firstStateProbs[firstIndex] > 0.0) && (secondStateProbs[secondIndex] > 0.0)){
			/*
			 * double division is probably more stable than multiplying two small numbers together
			 * mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / (state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex]));
			*/
			double temp = jointStateProbs[i] * log10f(jointStateProbs[i] / firstStateProbs[firstIndex] / secondStateProbs[secondIndex]);
			MI = MI + temp;
		}
	}

	MI /= log(2.0);
	
	//cudaFree(firstStateProbs);
	firstStateProbs = NULL;
	
	//cudaFree(secondStateProbs);
	secondStateProbs = NULL;
	
	//cudaFree(jointStateProbs);
	jointStateProbs = NULL;
		
	return MI;
}

int main(){
	
	printf("Working\n");
	return 0;
}
