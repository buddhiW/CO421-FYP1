/*
 * Author: Buddhi Wickramasinghe
 * Date: 08.07.2016 
 * Version: 1.0
 * Description: This file contains functions for calculating PDF using KDE 
 * 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "helpers.h"
#include "probOps.h"
#include "arrayOperations.h"

double find_max(double * data, int length){
	
	double max = 0;
	int i;
	
	for(i=0; i<length; i++){	
		if(data[i] > max)
			max = data[i];
	}	
	
	return max;
}


double find_min(double * data, int length){
	
	double min = 100000.0;
	int i;
	
	for(i=0; i<length; i++){	
		if(data[i] < min)
			min = data[i];
	}	
	
	return min;
}



double gauss(double x, double data, double h){
	
	//double h = bandWidth;
	double two_h_square=2*h*h;
	double pi=3.14159265358979;
	
	double temp = x - data;
	double norm = temp*temp; 
	
	return (exp(-norm/two_h_square)/(h*sqrt(2*pi)));
}

double kde(double sourceArray[2][SAMPLES], double *target, int sourceLength, int vars, double *bandWidth){

	//double h = bandWidth;
	//double two_h_square=2*h*h;
	//double pi=3.14159265358979;
	//double q=(pow(-1,r))/(sqrt(2*pi)*N*(pow(h,(r+1))));

	int i,j;
	double prob=0.0;

	for(i=0; i<sourceLength; i++){
		
		double mid = 1.0;
		
		for(j=0; j<vars; j++){
		
			double temp=target[j]-sourceArray[j][i];
			
			double temp1 = gauss(target[j], sourceArray[j][i], bandWidth[j]);

			mid = mid * temp1;	
			
		}
		
		prob = prob + mid;
		//printf("%lf ", prob);	
	}
	
	//printf("prob: %lf\n", prob);
	
	return (prob/SAMPLES);//?????
}

void gauss_kde_marginal(double *sourceArray, double bandWidth, int targetLength, double * probArray){
	
	double x_increment = 0.0;
	
	double x_max = find_max(sourceArray, SAMPLES) + (3*bandWidth); 
	double x_min = find_min(sourceArray, SAMPLES) - (3*bandWidth); 

	int points = 200;
	
	x_increment =  (x_max - x_min)/points;
	
	int i=0,j;
	
	double x = x_min;

	for(x = x_min; x<x_max; i++,x += x_increment){

		double prob=0.0;

		printf("%d\n", i);
		for(j=0; j<SAMPLES; j++){

			double temp = gauss(x, sourceArray[j], bandWidth);	
			prob = prob + temp;	
			
		}	
	
		probArray[i] = prob/SAMPLES;
	}
	
}



void gauss_kde_joint(double sourceMatrix[2][SAMPLES], double *targetArray1, double bandWidth1, int targetLength1, 
					double *targetArray2, double bandWidth2, int targetLength2, double *jointProbArray){
						
	int i,j;
	int jointTargetLength = targetLength1 * targetLength2;
	double target[2] = {0};
	double bandWidth[2] = {0};
	
	bandWidth[0] = bandWidth1;
	bandWidth[1] = bandWidth2;
	int vars = 2;
	
	double x_increment = 0.0;
	double y_increment = 0.0;
	
	double x_max = find_max(sourceMatrix[0], SAMPLES) + (3*(bandWidth[0])); 
	double x_min = find_min(sourceMatrix[0], SAMPLES) - (3*(bandWidth[0])); 
	double y_max = find_max(sourceMatrix[1], SAMPLES) + (3*(bandWidth[1])); 
	double y_min = find_min(sourceMatrix[1], SAMPLES) - (3*(bandWidth[1])); 
	int points = 200;
	x_increment =  (x_max - x_min)/points;
	y_increment =  (y_max - y_min)/points;
	
	double x = x_min;
	double y = y_min; 
		
	printf("incx: %lf, incy: %lf\n", x_increment, y_increment);
	printf("minX: %lf, minY: %lf\n", x_min, y_min);
	
	//for(i = 0; i < 200; i++, x += x_increment){
		//	printf("%lf ", x);			
	//}
	
	for(j=0; j<targetLength1; j++, x += x_increment){
		
		double y = y_min; 

		for(i=0; i<targetLength2; i++, y += y_increment){
			
			target[0] = x;
			target[1] = y;

			double probVal = kde(sourceMatrix, target, 200, vars, bandWidth);
			
			//if(probVal >= 1)
				//printf("damn! %d\n", j);
			jointProbArray[j*SAMPLES+i] = probVal;
			//printf("%lf ", jointProbArray[j]);	
		}	
	}						
}



JointProbabilityState* jointProbability(double *array1, double *array2, int Length){
	
	JointProbabilityState *state = malloc(sizeof(JointProbabilityState));

	int *normalisedArray1 = (int *) malloc(Length*sizeof(int));
	checkMalloc(normalisedArray1,0);
	int *normalisedArray2 = (int *) malloc(Length*sizeof(int));
	checkMalloc(normalisedArray2,1);
	
	int firstNumStates = normaliseArray(array1,normalisedArray1,Length);
	int secondNumStates = normaliseArray(array2,normalisedArray2,Length);
	int jointNumStates = firstNumStates * secondNumStates;
	
	//printf("This is the size of first num states: %d\n",firstNumStates);
	//printf("This is the size of joint num states: %d\n",jointNumStates);
	
	int *firstStateCounts = (int *) malloc(firstNumStates*sizeof(int));
	checkMalloc(firstStateCounts,2);
	int *secondStateCounts = (int *) malloc(secondNumStates*sizeof(int));
	checkMalloc(secondStateCounts,3);
	int *jointStateCounts = (int *) malloc(jointNumStates*sizeof(int));
	checkMalloc(jointStateCounts,4);
	
	int ii;
	
	for (ii=0;ii<firstNumStates;ii++)
		firstStateCounts[ii] = 0;
	for (ii=0;ii<secondNumStates;ii++)
		secondStateCounts[ii] = 0;
	for (ii=0;ii<jointNumStates;ii++)
		jointStateCounts[ii] = 0;
	
	double *firstStateProbs = (double *) malloc(firstNumStates*sizeof(double));
	checkMalloc(firstStateProbs,5);
	
	double *secondStateProbs = (double *) malloc(secondNumStates*sizeof(double));
	checkMalloc(secondStateProbs,6);
	
	double *jointStateProbs = (double *) malloc(jointNumStates*sizeof(double));
	checkMalloc(jointStateProbs,7);//??????????????????
	
	double length = Length;
	int i;
	
	for (i = 0; i < Length; i++){
		firstStateCounts[normalisedArray1[i]] += 1;
		secondStateCounts[normalisedArray2[i]] += 1;
		jointStateCounts[normalisedArray2[i] * firstNumStates + normalisedArray1[i]] += 1;
	}

	for (i = 0; i < firstNumStates; i++){
		//printf("First state count: %d\n",firstStateCounts[i]);
		firstStateProbs[i] = firstStateCounts[i] / length;
		//assert(firstStateProbs[i] >= 0.00000000);
	}
	
	for (i = 0; i < secondNumStates; i++){
		secondStateProbs[i] = secondStateCounts[i] / length;
		//assert(secondStateProbs[i] >= 0.0000000);
	}
	
	for (i = 0; i < jointNumStates; i++){
		jointStateProbs[i] = jointStateCounts[i] / length;
		//assert(jointStateProbs[i] >= 0.00000000);
	}
	
	free(normalisedArray1);
	free(normalisedArray2);
	
	normalisedArray1 = NULL;
	normalisedArray2 = NULL;
	
	free(firstStateCounts);
	free(secondStateCounts);
	free(jointStateCounts);

	
	firstStateCounts = NULL;
	secondStateCounts = NULL;
	jointStateCounts = NULL;
	

	
	state->jointProbabilityVector = jointStateProbs;
	state->numJointStates = jointNumStates;
	state->firstProbabilityVector = firstStateProbs;
	state->numFirstStates = firstNumStates;
	state->secondProbabilityVector = secondStateProbs;
	state->numSecondStates = secondNumStates;

	return state;
}

ProbabilityState probability(double *array, int Length){

	ProbabilityState state;

	double length = (double)Length;

	int *normalisedVector = (int *) malloc(Length*sizeof(int));

	int numStates = normaliseArray(array,normalisedVector,Length);//??????????
	
	//printf("Normalised array: ");
	
	//printf("Num states= %d\n=",numStates);
	
	int *stateCounts = (int *) malloc(numStates*sizeof(int));// Uninitialised value was created by a heap allocation
	double *stateProbs = (double *) malloc(numStates*sizeof(double));

	int i;
	for (i = 0; i < Length; i++){
		stateCounts[normalisedVector[i]] += 1;
	}
	
	for (i = 0; i < numStates; i++){
		stateProbs[i] = stateCounts[i] / length;
	}

	free(stateCounts);
	free(normalisedVector);
	
	state.probabilityVector = stateProbs;
	state.numStates = numStates;

	return state;
}

/*
int main(){
	
	double fe[] = {5.280000,5.320000,5.320000,5.360000,5.360000,5.360000,5.360000,5.320000,5.320000,5.320000,5.320000,5.280000,5.320000,5.320000,5.360000,5.360000,5.360000,5.360000,5.360000,5.400000,5.400000,5.360000,5.400000,5.400000,5.440000,5.440000,5.440000,5.440000,5.440000,5.480000,5.480000,5.480000,5.480000,5.480000,5.480000,5.480000,5.480000,5.480000,5.440000,5.440000,5.400000,5.360000,5.360000,5.320000,5.320000,5.360000,5.360000,5.400000,5.360000,5.320000,5.280000,5.320000,5.320000,5.320000,5.360000,5.360000,5.400000,5.400000,5.400000,5.400000,5.400000,5.400000,5.400000,5.400000,5.400000,5.440000,5.440000,5.440000,5.440000,5.480000,5.480000,5.480000,5.480000,5.480000,5.480000,5.520000,5.520000,5.520000,5.520000,5.520000,5.480000,5.400000,5.320000,5.320000,5.320000,5.360000,5.400000,5.400000,5.400000,5.440000,5.440000,5.400000,5.360000,5.320000,5.320000,5.360000,5.360000,5.360000,5.360000,5.360000,5.360000,5.400000,5.400000,5.400000,5.400000,5.440000,5.480000,5.480000,5.480000,5.440000,5.440000,5.440000};
	
	double df[] =  {0, 2, 4, 8, 8, 1, 8, 8, 4, 7, 1, 3, 8, 7, 8, 5, 0, 7, 8, 6, 
					6, 6, 3, 5, 1, 6, 0, 2, 0, 0, 7, 6, 2, 8, 0, 3, 3, 6, 7, 1, 
					4, 4, 5, 6, 6, 2, 6, 5, 1, 1, 4, 8, 3, 5, 2, 6, 2, 4, 6, 8, 
					8, 4, 1, 1, 2, 7, 2, 7, 2, 8, 3, 1, 2, 5, 4, 3, 7, 5, 4, 8, 
					2, 6, 6, 3, 5, 0, 0, 4, 7, 8, 1, 5, 4, 0, 3, 1, 7, 2, 4, 1, 
					5, 2, 5, 6, 6, 4, 0, 2, 8, 1, 7, 4};
					
	riseUp(fe,100,112);
	//printArrayDouble(fe,112);
	
	//ProbabilityState fds = probability(fe,112);
	
	//printf("Testing probability\n");
	//printf("P(fe)= ");
	//printArrayDouble(fds.probabilityVector,112);//????????????????????????????????
	//printf("States= %d\n",fds.numStates);
	
	JointProbabilityState *sre = jointProbability(fe, df, 112);*/
	  
	  /*double *jointProbabilityVector;
		int numJointStates;
		double *firstProbabilityVector;
		int numFirstStates;
		double *secondProbabilityVector;
		int numSecondStates;
	  */
	 /*
	printf("Testing joint probability\n");
	printArrayDouble(sre->jointProbabilityVector,112);
	printf("States= %d\n",sre->numJointStates);
	return 0;
}
*/
