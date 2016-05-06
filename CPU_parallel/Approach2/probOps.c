/*
 * Author: Malin Prematilake
 * Date: 23.03.2016 00:28:11
 * Version: 1.0
 * Description: This code contains all the functions related with calculating probability
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "helpers.h"
#include "probOps.h"
#include "arrayOperations.h"

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
