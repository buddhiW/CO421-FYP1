#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "probOps.h"
#include "arrayOperations.h"
#include "MI.h"

double MI_KDE(double *dataArray, double *targetArray, int length){
	//printf("Inside MI \n");
	double MI = 0.0;
	//printf("MI at the begining = %lf\n",MI);
	int firstIndex=0,secondIndex=0;
	int i;
	
	//JointProbabilityState *state = malloc(sizeof(JointProbabilityState));
	JointProbabilityState *state = NULL;
	
	state = jointProbability(dataArray,targetArray,length);
	
	//JointProbabilityState state = jointProbability(dataArray,targetArray,length);
	
	/*
	** I(X;Y) = sum sum p(xy) * log (p(xy)/p(x)p(y))
	*/
	
	for (i = 0; i < state->numJointStates; i++){
		
		firstIndex = i % state->numFirstStates;
		secondIndex = i / state->numFirstStates;
		
		//printf("Condition0= %d\n",i);
		//printf("Condition1= %.10lf\n",state->jointProbabilityVector[i]);
		//printf("[%d]Condition2= %.10lf\n",firstIndex,state->firstProbabilityVector[firstIndex] );
		//printf("[%d]Condition3= %.10lf\n",secondIndex,state->secondProbabilityVector[secondIndex]);
		
		if ((state->jointProbabilityVector[i] > 0.0) && (state->firstProbabilityVector[firstIndex] > 0.0) && (state->secondProbabilityVector[secondIndex] > 0.0)){
			/*
			 * double division is probably more stable than multiplying two small numbers together
			 * mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / (state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex]));
			*/
			double temp = state->jointProbabilityVector[i] * log(state->jointProbabilityVector[i] / state->firstProbabilityVector[firstIndex] / state->secondProbabilityVector[secondIndex]);
			//printf("Temp(%d)= %lf\n",i,temp);
			MI = MI + temp;
		}
	}

	MI /= log(2.0);

	
	free(state->firstProbabilityVector);
	state->firstProbabilityVector = NULL;
	free(state->secondProbabilityVector);
	state->secondProbabilityVector = NULL;
	free(state->jointProbabilityVector);
	state->jointProbabilityVector = NULL;
	
	free(state);
	
	return MI;
}


double MI(double *dataArray, double *targetArray, int length){
	//printf("Inside MI \n");
	double MI = 0.0;
	//printf("MI at the begining = %lf\n",MI);
	int firstIndex=0,secondIndex=0;
	int i;
	
	//JointProbabilityState *state = malloc(sizeof(JointProbabilityState));
	JointProbabilityState *state = NULL;
	
	state = jointProbability(dataArray,targetArray,length);
	
	//JointProbabilityState state = jointProbability(dataArray,targetArray,length);
	
	/*
	** I(X;Y) = sum sum p(xy) * log (p(xy)/p(x)p(y))
	*/
	
	for (i = 0; i < state->numJointStates; i++){
		
		firstIndex = i % state->numFirstStates;
		secondIndex = i / state->numFirstStates;
		
		//printf("Condition0= %d\n",i);
		//printf("Condition1= %.10lf\n",state->jointProbabilityVector[i]);
		//printf("[%d]Condition2= %.10lf\n",firstIndex,state->firstProbabilityVector[firstIndex] );
		//printf("[%d]Condition3= %.10lf\n",secondIndex,state->secondProbabilityVector[secondIndex]);
		
		if ((state->jointProbabilityVector[i] > 0.0) && (state->firstProbabilityVector[firstIndex] > 0.0) && (state->secondProbabilityVector[secondIndex] > 0.0)){
			/*
			 * double division is probably more stable than multiplying two small numbers together
			 * mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / (state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex]));
			*/
			double temp = state->jointProbabilityVector[i] * log(state->jointProbabilityVector[i] / state->firstProbabilityVector[firstIndex] / state->secondProbabilityVector[secondIndex]);
			//printf("Temp(%d)= %lf\n",i,temp);
			MI = MI + temp;
		}
	}

	MI /= log(2.0);

	
	free(state->firstProbabilityVector);
	state->firstProbabilityVector = NULL;
	free(state->secondProbabilityVector);
	state->secondProbabilityVector = NULL;
	free(state->jointProbabilityVector);
	state->jointProbabilityVector = NULL;
	
	free(state);
	
	return MI;
}
/*
int main(){
	double dataA[] = {	5.280000,5.320000,5.320000,5.360000,5.360000,5.360000,5.360000,
					5.320000,5.320000,5.320000,5.320000,5.280000,5.320000,5.320000,
					5.360000,5.360000,5.360000,5.360000,5.360000,5.400000,5.400000,
					5.360000,5.400000,5.400000,5.440000,5.440000,5.440000,5.440000,
					5.440000,5.480000,5.480000,5.480000,5.480000,5.480000,5.480000,
					5.480000,5.480000,5.480000,5.440000,5.440000,5.400000,5.360000,
					5.360000,5.320000,5.320000,5.360000,5.360000,5.400000,5.360000,
					5.320000,5.280000,5.320000,5.320000,5.320000,5.360000,5.360000,
					5.400000,5.400000,5.400000,5.400000,5.400000,5.400000,5.400000,
					5.400000,5.400000,5.440000,5.440000,5.440000,5.440000,5.480000,
					5.480000,5.480000,5.480000,5.480000,5.480000,5.520000,5.520000,
					5.520000,5.520000,5.520000,5.480000,5.400000,5.320000,5.320000,
					5.320000,5.360000,5.400000,5.400000,5.400000,5.440000,5.440000,
					5.400000,5.360000,5.320000,5.320000,5.360000,5.360000,5.360000,
					5.360000,5.360000,5.360000,5.400000,5.400000,5.400000,5.400000,
					5.440000,5.480000,5.480000,5.480000,5.440000,5.440000,5.440000};
	
	double trgtA[] =  {0.0, 2.0, 4.0, 8.0, 8.0, 1.0, 8.0, 8.0, 4.0, 7.0, 1.0, 3.0, 8.0, 7.0, 8.0, 5.0, 0.0, 7.0, 8.0, 6.0, 
					6.0, 6.0, 3.0, 5.0, 1.0, 6.0, 0.0, 2.0, 0.0, 0.0, 7.0, 6.0, 2.0, 8.0, 0.0, 3.0, 3.0, 6.0, 7.0, 1.0, 
					4.0, 4.0, 5.0, 6.0, 6.0, 2.0, 6.0, 5.0, 1.0, 1.0, 4.0, 8.0, 3.0, 5.0, 2.0, 6.0, 2.0, 4.0, 6.0, 8.0, 
					8.0, 4.0, 1.0, 1.0, 2.0, 7.0, 2.0, 7.0, 2.0, 8.0, 3.0, 1.0, 2.0, 5.0, 4.0, 3.0, 7.0, 5.0, 4.0, 8.0, 
					2.0, 6.0, 6.0, 3.0, 5.0, 0.0, 0.0, 4.0, 7.0, 8.0, 1.0, 5.0, 4.0, 0.0, 3.0, 1.0, 7.0, 2.0, 4.0, 1.0, 
					5.0, 2.0};
	
	double sde = MI(trgtA,dataA,112);
	printf("MI= %lf\n",sde);
	
	return 0;
}
*/
