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

#include "helpers.cuh"
#include "probOps.cuh"
#include "arrayOperations.cuh"

/*#define SAMPLEPOINTS 100000
#define SAMPLES 200

//defined for 128 bit AES
#define KEYBYTES 16
#define KEYS 256*/

__global__ void findProbsWave(unsigned int *dev_normaliseArray1, double * dev_firstStateProbs, unsigned int *dev_firstNumState, int length){
	
	//int firstStateCounts[SAMPLES];
	int i,j,postt;
	double add = 1.0/(double)SAMPLES;
	//int samplePoint = (blockIdx.x*blockDim.x+blockDim.y) + threadIdx.x; //Directly gives sample point
	int samplePoint = (gridDim.x*blockDim.x*blockIdx.y) + threadIdx.x + blockIdx.x*blockDim.x;
	
	int startPositionProb = samplePoint*length;
	int startPositionNorm = samplePoint*SAMPLES; //(sizeof(dev_normaliseArray1) != dev_firstStateProbs)
	
	if (samplePoint < SAMPLEPOINTS){
	
		//int end = startPositionProb+dev_firstNumState[samplePoint];
		
		//if(samplePoint>8192 && samplePoint<8192)
			//printf("Index: (%d,%d)\n",startPositionNorm,end);
		//	printf("sample: %d\n", samplePoint);
		
		
		
		for (i = startPositionProb; i < startPositionProb+length; i++){ //changed end to length
			dev_firstStateProbs[i] = 0.0;
			//dev_firstStateProbs[i] = samplePoint;
		}
		
		for (j = startPositionNorm; j < startPositionNorm+SAMPLES ; j++){
			postt = startPositionProb + dev_normaliseArray1[j];
			dev_firstStateProbs[postt] = dev_firstStateProbs[postt] + add;	//????????	
		}
		
	}
	
} 

__global__ void findProbsHamming(unsigned int *dev_normaliseArray2, double * dev_secondStateProbs, unsigned int *dev_secondNumState, int length){
	
	int key =  blockIdx.x; // There are 256 blocks 
	int keyByte = threadIdx.x; 

	int i;

	int startPositionNorm = key*SAMPLES*KEYBYTES + keyByte*SAMPLES;
	int startPositionProb = key*length*KEYBYTES + keyByte*length;
	//int length = SAMPLES;

	if ((key < KEYS) && (keyByte < KEYBYTES)){
	
		int end = dev_secondNumState[key*KEYBYTES+keyByte];
		
		//if(key==0)
		//	printf("(%d, %d)\n", startPositionProb,startPositionNorm);
		
		for (i = startPositionProb; i < startPositionProb+end; i++){
			dev_secondStateProbs[i] = 0; //??
		}
		
		for (i = startPositionNorm; i < startPositionNorm+SAMPLES ; i++){
			dev_secondStateProbs[startPositionProb + dev_normaliseArray2[i]] += 1;	//????????
		}
	
	
		for (i = startPositionProb; i < startPositionProb+end; i++){
			dev_secondStateProbs[i] = dev_secondStateProbs[i] / SAMPLES;
		}
	
	}	
}

__global__ void findJointProbs(double *MIvals, unsigned int *normaliseWave, unsigned int *normaliseHamming, unsigned int *firstNumStates, 
								unsigned int *secondNumStates, double *waveProbs, double *hammingProbs, double *jointProbs, 
								int maxOfFirst, int maxOfSecond, int key, int startPosition){
	
	double add = 1.0/(double)SAMPLES;
	
	int xx = blockIdx.x*blockDim.x + threadIdx.x;
	int yy = blockIdx.y;
	
	//int index = yy*gridDim.x*blockDim.x + xx;
	int maxOfJoint = maxOfFirst*maxOfSecond;
	int startJoint = xx*KEYBYTES*maxOfJoint + yy*maxOfJoint;
	
	//int startJoint = (xx+startPosition)*KEYBYTES*maxOfJoint + yy*maxOfJoint;
	//int endPosition = startPosition + SAMPLEPOINTS/REPEAT;
	int firstIndx, secondIndx;
	
	if (xx < SAMPLEPOINTS/REPEAT){
		int i;
		
		for (i=0; i<SAMPLES; i++){
			int postt = normaliseHamming[key*KEYBYTES*SAMPLES + yy*SAMPLES +i] * firstNumStates[xx] + normaliseWave[(xx+startPosition)*SAMPLES + i];
			jointProbs[startJoint + postt] = jointProbs[startJoint + postt] + add;
			
			//if(xx==0 && yy==0)
			//	printf("postt: %d\n", (startJoint + postt));
		}
		
	/*	int fNumState = firstNumStates[xx];
		int sNumState = secondNumStates[key*KEYBYTES + yy];
		
		int jNumState = fNumState*sNumState;
		double MI = 0.0;
		
		//if(xx==0 && yy==0)
		//	printf("sNum %d\n", sNumState);
		
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
		//if(xx==0 && yy==0)
		//	printf("MI: %d\n", MI);
		//MIvals[key*KEYBYTES*SAMPLEPOINTS + yy*SAMPLEPOINTS + xx] = MI;
		MIvals[yy*SAMPLEPOINTS + (xx+startPosition)] = MI; //MIvals is for one key*/
	}
	
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
