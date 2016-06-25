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

//#define WIDTH 5

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
		
		//if(samplePoint>99000)
			//printf("Index: (%d,%d)\n",startPositionNorm,end);
			//printf("binNum: %d\n", dev_firstNumState[samplePoint]);
		
		
		
		for (i = startPositionProb; i < startPositionProb+length; i++){ //changed end to length
			dev_firstStateProbs[i] = 0.0;
			//dev_firstStateProbs[i] = samplePoint;
		}
		
		for (j = startPositionNorm; j < startPositionNorm+SAMPLES ; j++){
			postt = startPositionProb + dev_normaliseArray1[j];
			//postt = startPositionProb + (dev_normaliseArray1[j]/WIDTH); //(dev_normaliseArray1[j]/WIDTH) = index of the bin
			//if(samplePoint==0)
				//printf("val: %d %d\n", dev_normaliseArray1[j], (dev_normaliseArray1[j]/WIDTH));
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
			int postt = normaliseHamming[key*KEYBYTES*SAMPLES + yy*SAMPLES +i] * firstNumStates[xx+startPosition] + normaliseWave[(xx+startPosition)*SAMPLES + i];
			jointProbs[startJoint + postt] = jointProbs[startJoint + postt] + add;
			
			//if(xx==0 && yy==0)
			//	printf("postt: %d\n", (startJoint + postt));
		}
		int fNumState = firstNumStates[xx+startPosition];
		int jNumState = fNumState*secondNumStates[key*KEYBYTES + yy];
		double MI = 0.0;
		//printf("MI at 1:- %lf\n", MI);
		
		for (i=0; i<jNumState; i++){
			firstIndx = i%fNumState;
			secondIndx = i/fNumState;
			
			double c1 = waveProbs[(xx+startPosition)*maxOfFirst + firstIndx];
			double c2 = hammingProbs[key*KEYBYTES*maxOfSecond + yy*maxOfSecond + secondIndx];
			double c3 = jointProbs[xx*KEYBYTES*maxOfJoint + yy*maxOfJoint + i];
			
			if ((c3 > 0.0) && (c1 > 0.0) && (c2 > 0.0)){
					
					double temp = c3 * log(c3 / (c1 * c2));
					MI = MI + temp;
					//printf("MI: %lf, temp: %lf\n",MI,temp);
			}			
		}
		
		MI = MI/log(2.0);
		//if (MI>0.5) printf("MI at 2:- %lf\n", MI);
		
		//MIvals[xx + startPosition] = MI;
		//if ((KEYBYTES*SAMPLEPOINTS + yy*SAMPLEPOINTS + xx + startPosition)<10)  printf("The MI in CUDA = %lf\n",MI);
		MIvals[key*KEYBYTES*SAMPLEPOINTS + yy*SAMPLEPOINTS + (xx+startPosition)] = MI;
		//MIvals[yy*SAMPLEPOINTS + (xx+startPosition)] = MI; //MIvals is for one key
	}	
}


__global__ void findJointProbs2(double *MIvals, unsigned int *normaliseWave, unsigned int *normaliseHamming, unsigned int *firstNumStates, 
								unsigned int *secondNumStates, double *waveProbs, double *hammingProbs, double *jointProbs, 
								unsigned long *jointSizes, int maxOfFirst, int maxOfSecond, int turn){
	
	double add = 1.0/(double)SAMPLES;
	
	int xx = blockIdx.x*blockDim.x + threadIdx.x;//16 KEYBYTES
	int yy = blockIdx.y*blockDim.y + threadIdx.y;//8 KEYS/32
	int zz = blockIdx.z*blockDim.z + threadIdx.z;//100,000 SAMPLEPOINTS
	
	int posjSizes = yy*KEYBYTES*SAMPLEPOINTS + xx*SAMPLEPOINTS + zz; //OK
	unsigned long startPos;
	int i;
	
	if(zz<SAMPLEPOINTS && yy<8 && xx<16){
	
		if(posjSizes == 0)
			startPos = 0;
		else 
			startPos = jointSizes[posjSizes-1];//starting position for jointProbs
		
		int jNumState = (int)(jointSizes[posjSizes] - startPos);
		
		int normHammingPost = (turn+yy)*KEYBYTES*SAMPLES + xx*SAMPLES; //Added SAMPLES
		int fNumStatePost = zz;
		int normWavePost = zz * SAMPLES; //Added SAMPLES
		
		for (i=0; i<SAMPLES; i++){
				
				int postt = normaliseHamming[normHammingPost +i] * firstNumStates[fNumStatePost] + normaliseWave[normWavePost + i];
				jointProbs[startPos + postt] = jointProbs[startPos + postt] + add;
		}
		
		int fNumState = firstNumStates[fNumStatePost];
		
		double MI = 0.0;
		int firstIndx, secondIndx;
		
		for (i=0; i<jNumState; i++){
			firstIndx = i%fNumState;
			secondIndx = i/fNumState;
			
			double c1 = waveProbs[zz*maxOfFirst + firstIndx];
			double c2 = hammingProbs[(turn+yy)*KEYBYTES*maxOfSecond + xx*maxOfSecond + secondIndx];
			double c3 = jointProbs[startPos + i];
			
			if ((c3 > 0.0) && (c1 > 0.0) && (c2 > 0.0)){
					
					double temp = c3 * log(c3 / (c1 * c2));
					MI = MI + temp;
			}		
		}
		MI = MI/log(2.0);
		MIvals[(turn+yy)*KEYBYTES*SAMPLEPOINTS + xx*SAMPLEPOINTS + zz] = MI;
	}
}
