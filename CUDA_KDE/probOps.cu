/*
 * Version info
 * MP		Initial version 23.03.2016
 * BW		Changed methods to do KDe based probability calculation 
 *
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

/**********************************************************************************************************************************/

__device__ double gauss(double x, double data, double h){
	
	//double h = bandWidth;
	double two_h_square=2*h*h;
	double pi=3.14159265358979;
	
	double temp = x - data;
	double norm = temp*temp; 
	
	double temp1 = exp(-norm/two_h_square)/(h*sqrt(2*pi));
	
	
	if(isnan(temp1)!=0){
		__threadfence();
		asm("trap;");
	}
				
	return temp1;
}


__global__ void KDE_findProbsWave(double *dev_waveData, double * dev_waveDataProbs, 
									double *bandWidthArray, int targetLength, double * maxMin){
	
	//int firstStateCounts[SAMPLES];
	int i=0,j=0;
	
	//int samplePoint = (blockIdx.x*blockDim.x+blockDim.y) + threadIdx.x; //Directly gives sample point
	int samplePoint = (gridDim.x*blockDim.x*blockIdx.y) + threadIdx.x + blockIdx.x*blockDim.x; //directly gives samplepoint
	
	//int startPositionProb = samplePoint*length;
	//int startPositionNorm = samplePoint*SAMPLES; //(sizeof(dev_normaliseArray1) != dev_firstStateProbs)
	
	
	double x_increment = 0.0;
	
	if (samplePoint < SAMPLEPOINTS){

		double min = maxMin[samplePoint*2];
		double max = maxMin[samplePoint*2+1];
		double bandWidth = bandWidthArray[samplePoint];
		double x_max = max + (3*bandWidth); 
		double x_min = min - (3*bandWidth); 
	
		x_increment =  (x_max - x_min)/targetLength;
	
		double x = x_min;
		double prob = 0.0;
		
		//if(isnan(bandWidth)!=0){
			//printf("errorNan %d %d\n",samplePoint, __LINE__);
			//__threadfence();
			//asm("trap;");
		//}
		
		//if(samplePoint == 0)
			//printf("%lf %lf %lf %lf\n", x_increment, bandWidth, x_max, x_min);
		
		for(x = x_min; x<x_max; i++,x += x_increment){

			if(i<targetLength){
				
				prob=0.0;

				//printf("%d\n", i);
				for(j=0; j<SAMPLES; j++){

					double temp = gauss(x, dev_waveData[samplePoint*SAMPLES+j], bandWidth);	
					prob = prob + temp;	
					
				}	
				
				double probFinal = prob/SAMPLES; 
				dev_waveDataProbs[samplePoint*targetLength+i] = probFinal;
				//if(probFinal < 0)
					//printf("%d ", samplePoint);
			}
			
			//if(i==targetLength)
				//printf("%lf\n", prob/SAMPLES);
			
		}
	}	
} 

__global__ void KDE_findProbsHamming(double *dev_hammingData, double * dev_hammingDataProbs,
										double *bandWidthArray, int targetLength, double * maxMin){
	
	int key =  blockIdx.x; // There are 256 blocks 
	int keyByte = threadIdx.x; 

	int i=0,j=0;
	
	int startPositionData = key*SAMPLES*KEYBYTES + keyByte*SAMPLES;
	int startPositionProb = key*targetLength*KEYBYTES + keyByte*targetLength;
	//int length = SAMPLES

	if ((key < KEYS) && (keyByte < KEYBYTES)){
	
		double bandWidth = bandWidthArray[key*KEYBYTES + keyByte];

		double min = maxMin[key*2*KEYBYTES+2*keyByte];
		double max = maxMin[key*2*KEYBYTES+2*keyByte+1];
		
		double x_max = max + (3*bandWidth); 
		double x_min = min - (3*bandWidth); 
		
		//if(key<100)
			//printf("%lf, %lf \n", x_max,x_min);
			
		double x_increment =  (x_max - x_min)/targetLength;
	
		double x = x_min;
		
		for(x = x_min; x<x_max; i++,x += x_increment){

			if(i<targetLength){
				double prob=0.0;

				//printf("%d\n", i);
				for(j=0; j<SAMPLES; j++){

					double temp = gauss(x, dev_hammingData[startPositionData+j], bandWidth);	
					prob = prob + temp;	
					
				}	
				dev_hammingDataProbs[startPositionProb+i] = prob/SAMPLES;
			}
		}
	
	}	
}

//Without MI calculation
/*	dim3 block3d(16,2,4);
	dim3 grid3d(1,4,25000);*/
__global__ void KDE_findJointProbs(double *dev_hammingData, double *dev_waveData, unsigned int targetLengthW, 
									unsigned int targetLengthH, double *jointProbs, double * bandwidthW, double * bandwidthH, 
									double * maxMinH, double * maxMinW, int turn){
	
	int xx = blockIdx.x*blockDim.x + threadIdx.x;//16 KEYBYTES
	int yy = blockIdx.y*blockDim.y + threadIdx.y;//8 KEYS/32
	int zz = blockIdx.z*blockDim.z + threadIdx.z;//100,000 SAMPLEPOINTS
	
	int key = yy + turn;
	int keyByte = xx;
	int samplePoint = zz;
	
	//int startJoint = (xx+startPosition)*KEYBYTES*maxOfJoint + yy*maxOfJoint;
	//int endPosition = startPosition + SAMPLEPOINTS/REPEAT;
	
	int i,j;
	
	//int jointTargetLength = targetLengthH * targetLengthW;
	
	double bandWidthW = bandwidthW[zz];
	double bandWidthH = bandwidthH[key*KEYBYTES + keyByte];
	
	int vars = 2;
	
	double x_increment = 0.0;
	double y_increment = 0.0;
	
	double x_min = maxMinW[samplePoint*2] - (3*bandWidthW);
	double x_max = maxMinW[samplePoint*2+1] + (3*bandWidthW);
	double y_min = maxMinH[key*2*KEYBYTES+2*keyByte] - (3*bandWidthH); 
	double y_max = maxMinH[key*2*KEYBYTES+2*keyByte+1] + (3*bandWidthH);
	
	x_increment =  (x_max - x_min)/targetLengthW;
	y_increment =  (y_max - y_min)/targetLengthH;

	if (zz<SAMPLEPOINTS && yy<KEYS/32 && xx<KEYBYTES){

		double x = x_min;
		
		for(j=0; j<targetLengthW; j++, x += x_increment){
		
			double y = y_min; 
			double prob;
			for(i=0; i<targetLengthH; i++, y += y_increment){
				
				for(i=0; i<SAMPLES; i++){
		
					double mid = 1.0;
		
					for(j=0; j<vars; j++){
		
						double tempW = gauss(x, dev_waveData[samplePoint+i], bandWidthW);
						double tempH = gauss(y, dev_hammingData[key*KEYBYTES + keyByte], bandWidthH);

						mid = tempW * tempW;				
					}		
					prob = prob + mid;
				}
				
				prob = (prob/SAMPLES);

				jointProbs[key*SAMPLEPOINTS*KEYBYTES + samplePoint*KEYBYTES+keyByte+targetLengthW*targetLengthH] = prob;

			}	
		}					
		
	}	
}

__global__ void KDE_findJointProbs2(double *dev_hammingData, double *dev_waveData, unsigned int targetLengthW, 
									unsigned int targetLengthH, double * bandwidthW, double * bandwidthH, 
									double *dev_hammingTarget, double *dev_waveTarget, double *jointProbs,
									int turn, int key, int keyByte){
	
	int xx = blockIdx.x*blockDim.x + threadIdx.x;//SAMPLEPOINTS
	int yy = blockIdx.y*blockDim.y + threadIdx.y;//Index of wave target array
	int zz = blockIdx.z*blockDim.z + threadIdx.z;//Index of Hamming target array
	
	int jointIndex = zz + blockDim.z*gridDim.z * yy; // Index in the joint prob array for the given sample point
	int i;
	
	if (xx<SAMPLEPOINTS/DIVIDE && yy<targetLengthW && zz<targetLengthH){

			int samplePoint = xx+turn; //turn goes from 0 to 90,000, in steps of 10,000
		
			double bandWidthW = bandwidthW[samplePoint];
			double bandWidthH = bandwidthH[key*KEYBYTES + keyByte];
		
			double prob=0.0;	
			unsigned int i1 = key*KEYBYTES*targetLengthH+keyByte*targetLengthH+zz;	
			unsigned int i2 = samplePoint*targetLengthW+yy;
			double targetHamming = dev_hammingTarget[i1]; 
			double targetWave = dev_waveTarget[i2];  	
			int startPositionHamming = key*SAMPLES*KEYBYTES + keyByte*SAMPLES;	
		
			//if(isnan(bandWidthW)!=0){
					//printf("errorNan %d %d\n",samplePoint, __LINE__);
				//	__threadfence();
				//	asm("trap;");
			//}
		
				for(i=0; i<SAMPLES; i++){
		
					double mid = 0.0;
		
					double tempW = gauss(targetWave, dev_waveData[samplePoint*SAMPLES+i], bandWidthW);
					double tempH = gauss(targetHamming, dev_hammingData[startPositionHamming+i], bandWidthH);

					mid = tempW * tempH;				
					
									
					prob = prob + mid;
				}
				
				prob = (prob/SAMPLES);

				int indexGlobal = jointIndex + targetLengthW*targetLengthH*xx;
				
				//prob = prob * 1000000;
				jointProbs[indexGlobal] = prob;
				//jointProbs[indexGlobal] = xx;
				
				//if(samplePoint < 2000){
					//if(isnan(jointProbs[indexGlobal])!=0){
						//__threadfence();
						//asm("trap;");
					//}
				//}
	}	
}

/* ***Different grid*** */
__global__ void KDE_MI(double *MIvals, double *waveProbs, double *hammingProbs, double *jointProbs, 
								int targetLengthH, int targetLengthW, int key, int keyByte, int turn){

	int xx = blockIdx.x*blockDim.x + threadIdx.x; //Samplepoint value
	
	int i;
	int firstIndx, secondIndx;
	double MI = 0.0;
	
	int samplePoint = xx+turn;
	int targetLengthJoint = targetLengthW*targetLengthH;
	
	if(xx<SAMPLEPOINTS/DIVIDE){
	
		for (i=0; i<targetLengthJoint; i++){
				//firstIndx = i%targetLengthW;
				//secondIndx = i/targetLengthW;
				firstIndx = i/targetLengthH;
				secondIndx = i%targetLengthH;
				
				
				double c1 = waveProbs[samplePoint*targetLengthW + firstIndx];
				double c2 = hammingProbs[key*KEYBYTES*targetLengthH + keyByte*targetLengthH + secondIndx];
				double c3 = jointProbs[xx*targetLengthJoint + i];
				
				//if(samplePoint==0)
					//printf("%lf %lf %lf\n", c1, c2, c3);
					
				if ((c3 > 0.0000001) && (c1 > 0.0000001) && (c2 > 0.0000001)){
						
						double temp = c3 * log(c3 / (c1 * c2));
						//if(samplePoint==0)
							//printf("%lf\n",(c3 / (c1 * c2)));
						MI = MI + temp;
						
						//if(samplePoint == 0)
							//printf("%lf\n",MI);
				}
				
				
		}
			
			MI = MI/log(2.0);
			
			MIvals[key*KEYBYTES*SAMPLEPOINTS + keyByte*SAMPLEPOINTS + samplePoint] = MI;
									
	}
}
/**********************************************************************************************************************************/

/*__global__ void findProbsWave(unsigned int *dev_normaliseArray1, double * dev_firstStateProbs, unsigned int *dev_firstNumState, int length){
	
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
}*/

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
