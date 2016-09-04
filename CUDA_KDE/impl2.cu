//
//	Version information:
//
//	Contains normalising functions for hammingMat and waveData
//
//	MP			Initial version
//	BW			Added makeAllHamming() method
//				Added printMatDouble() to arrayOperations.cu file
//				Changed hammingMat size
//				Changed argument order in cudaMemcpyDeviceToHost
//				Added sampleNum argument to hamming() method
//				
//	Grid for hamming matrix calculation
//	------------------------------------------------------------------  
//	|	256	|	256	|	256	|	256	|	256	|	256	| .....	|	256	|		
//	------------------------------------------------------------------
//	<-------------------------- 16 blocks ---------------------------->
//
//	Hamming result vector
//	--------------------------------------------------------------------
//	|<-200->			|   				|					|		
//	|   0  |...|   15   |		|	|		|		|	|		|	
//	|<-----16*200----->	|					|					|	
//	--------------------------------------------------------------------
// <-------------------------- 16*200*256 ---------------------------->
//	index of an element = key*KEYBYTES*SAMPLES + keyByte*SAMPLES + sampleNumber

#include <stdlib.h>
#include <stdio.h>
#include <syslog.h>
#include <limits.h>
#include "data.cuh"
#include "arrayOperations.cuh"
#include "probOps.cuh"
#include "helpers.cuh"

#define PLAINFILE "/storage/buddhi/data/plain.txt"
//#define WAVEFILE "wave.txt"
#define WAVEFILE "/storage/buddhi/data/wave.txt"

//number of sample points in one encryption
/*#define SAMPLEPOINTS 100000
#define SAMPLES 200

//defined for 128 bit AES
#define KEYBYTES 16
#define KEYS 256*/

//calculates hamming weight of a 8 bit number
__device__ byte hammingweight(byte H){

	// Count the number of set bits
	byte dist=0;
	while(H){
		dist++; 
		H &= H - 1;
	}
	return dist;
}

//find hamming weight for the selection function
/*__device__ byte hamming(unsigned int *sample, unsigned int i,unsigned int n,unsigned int key){ //n is byteno  i is the sample
    //byte inter = (byte)sbox[sample[i*KEYBYTES+n] ^ key];
    byte inter = (byte)sbox[sample[i] ^ key];
    //byte inter = 0x00;
    byte dist = hammingweight(inter);	  
	return dist;
}*/

__device__ byte hamming(unsigned int *sample, unsigned int index, unsigned int key, unsigned int sampleNum){ //n is byteno  i is the sample
	byte inter = sbox[sample[sampleNum*KEYBYTES+index] ^ key]; 
    //byte inter = 0xAB;
    byte dist = hammingweight(inter);	  
	return dist;
}

int arrayMax(unsigned int * data, int length){
	
	int max = data[0];
	for(int i=0; i<length; i++){
		if(data[i]>max)
			max = data[i];
	}
	return max;
}


/******************************************************************/
/**                   Kernels          							 **/
/******************************************************************/

//-------------------------------------------------------------------------------------------------------
//	Calculate all bandwidths for all keys and keybytes for hamming data
//  Calculate maximum and minimum of each data set
//-------------------------------------------------------------------------------------------------------
__global__ void KDE_bandwidthHammingKernel(double * dataArray, double * bandwidthArray, double * maxMinArray,
												int targetLengthH, double * targetArray){          
	
	/* Block size = KEYS/16, KEYS/16,1
	 * Grid size = KEYBYTE, 1,1
	 * */
	
	
	int j;
	int max = INT_MIN;
	int min = INT_MAX; 
	int keyByte = threadIdx.y;
	int key = blockIdx.x*blockDim.x + threadIdx.x;
	int targetLength = targetLengthH;
	double sum_array=0;
	double sum_2_array=0;

	
	
	if ((keyByte < KEYBYTES) && (key < KEYS)){
		for(j=0; j<SAMPLES; j++){

			int index = key*KEYBYTES*SAMPLES + keyByte*SAMPLES + j;
			sum_array += dataArray[index];
			sum_2_array += dataArray[index] * dataArray[index];
		
			if(dataArray[index] > max){
				max = dataArray[index];
			}
			
			if(dataArray[index] < min){
				min = dataArray[index];
			}
		
		}
		
		double x  = sum_array/SAMPLES;
		double x2 = sum_2_array/SAMPLES;
		//printf("%lf %lf ", x, x2);
		double sigma = sqrt(x2 - (x*x));
		double b = sigma*(pow((3.0*SAMPLES/4.0),(-1.0/5.0)));
		bandwidthArray[key*KEYBYTES + keyByte] = b;
		
		maxMinArray[key*2*KEYBYTES+2*keyByte] = min;
		maxMinArray[key*2*KEYBYTES+2*keyByte+1] = max;
		
		/* **************** calculate target arrays ***************** */
		
		double x_max = max + (3*b); 
		double x_min = min - (3*b); 	
		double x_increment =  (x_max - x_min)/targetLength;		
		double xVal=x_min;
		j=0;
		for(xVal = x_min; xVal<x_max; j++,xVal += x_increment){
			if(j<targetLength){
				
				int index = key*KEYBYTES*targetLength + keyByte*targetLength+j;
				targetArray[index] = xVal;
				if((keyByte == 10) && (key == 255))
					printf("%d\n", j);
					
			}
			
		}
		
		
	}
}

//-------------------------------------------------------------------------------------------------------
//	Calculate all bandwidths for all samplepoints for wave data
//-------------------------------------------------------------------------------------------------------
__global__ void KDE_bandwidthWaveKernel(double * dataArray, double * bandwidthArray, double * maxMinArray, 
											int targetLengthW, double * targetArray){
	
	/* block Size = SAMPLEPOINTS/1000,1,1
	 * grid Size = SAMPLEPOINTS/100,1,1
	 * */
	
	
	int j;
	int max = INT_MIN;
	int min = INT_MAX; 
	int samplePoint = blockIdx.x*blockDim.x + threadIdx.x;
	int targetLength = targetLengthW;
	double sum_array=0;
	double sum_2_array=0;

	if (samplePoint < SAMPLEPOINTS){
		
	
		
		for(j=0; j<SAMPLES; j++){
			int index = samplePoint*SAMPLES + j;
			sum_array += dataArray[index];
			sum_2_array += dataArray[index] * dataArray[index];
			
			//if(isnan(sum_2_array)!=0){ //x2 ERROR!!
			//printf("error at: %d\n", samplePoint);
				//__threadfence();
				//asm("trap;");
			//}
			
			if(dataArray[index] > max){
				max = dataArray[index];
			}
			
			if(dataArray[index] < min){
				min = dataArray[index];
			}
		
		}
		
		double x  = sum_array/SAMPLES;
		double x2 = sum_2_array/SAMPLES;
		
		double sigma = sqrt(x2 - (x*x));
		
		double b = sigma*(pow((3.0*SAMPLES/4.0),(-1.0/5.0)));
		bandwidthArray[samplePoint] = b;
		
		//if(isnan(x2)!=0){ //x2 ERROR!!
			//printf("error at: %d\n", samplePoint);
			//__threadfence();
			//asm("trap;");
		//}
		
		maxMinArray[2*samplePoint] = min;
		maxMinArray[2*samplePoint+1] = max;
		
		/* **************** calculate target arrays ***************** */
		
		double x_max = max + (3*b); 
		double x_min = min - (3*b); 	
		double x_increment =  (x_max - x_min)/targetLength;		
		double xVal=x_min;
		j=0;
		for(xVal = x_min; xVal<x_max; j++,xVal += x_increment){
			if(j<targetLength){
				int index = samplePoint*targetLength+j;
				targetArray[index] = xVal;
			}
		}
		
		
	}
}

//----------------------------------------------------------------------------------------------------------
//	Method to calculate hamming weights of all keyBytes of all samples w.r.t. all key guesses
//----------------------------------------------------------------------------------------------------------
__global__ void makeAllHamming(unsigned int *plaintxt, double *hammingRes, int n){//n is the size of samples 
	
	int keyByte = threadIdx.y;
	int key = blockIdx.x*KEYBYTES + threadIdx.x;
	int i;

	if ((keyByte < KEYBYTES) && (key < KEYS)){
		for (i=0; i<SAMPLES; i++){
			double temp = (double)hamming(plaintxt,keyByte, key, i); 			
			hammingRes[key*KEYBYTES*SAMPLES + keyByte*SAMPLES + i] = temp;
		}
	}
} 

//--------------------------------------------------------------------------------------------------------
//	Sort each MI set
//--------------------------------------------------------------------------------------------------------
__global__ void simpleSort(double *origMat, int cols, double *maxOfEach, int key){
	
	int index = threadIdx.x;
	int i;
	double max = 0.0;
	for (i=0;i<cols;i++){
		if (max<origMat[key*KEYBYTES*cols + index*cols + i])
			max = origMat[key*KEYBYTES*cols + index*cols + i];
	}
	maxOfEach[key*KEYBYTES + index] = max;
}

int main(){

	//cudaSetDevice(1);
	//int nDevices;
	
	//cudaGetDeviceCount(&nDevices);
	//for (int i = 0; i < nDevices; i++) {
    //cudaDeviceProp prop;
    //cudaGetDeviceProperties(&prop, i);
    //printf("Device Number: %d\n", i);
    //printf("  Device name: %s\n", prop.name);
    //printf("  Max Blocks: %d\n",
    //       *prop.maxGridSize);
  //}


	//plainText at host
	int width = 256*16;
	int i,j;
	unsigned int *dev_plainTexts;
	unsigned int *plainTexts = (unsigned int*)malloc(KEYBYTES*SAMPLES*sizeof(unsigned int));
	//checkMalloc(plainTexts);

	/*** Reading plaintext file ***/
	FILE *plainT = fopen(PLAINFILE,"r");
	char temp1[32];	

	for(i=0; i<SAMPLES ;i++){
		fscanf(plainT,"%s",temp1);
		for(j=0; j<KEYBYTES; j++){
			int post = i*KEYBYTES + j;
			sscanf(&temp1[2*j],"%02X",&plainTexts[post]);		
		}
	}
	
	/*** reading wave file ***/
	FILE *waveStuff = fopen(WAVEFILE,"r");
	
	double *waveData= (double *)malloc(sizeof(double) * SAMPLEPOINTS * SAMPLES);

	float dat;

	
	for(i=0; i<SAMPLES ;i++){
		for(j=0; j<SAMPLEPOINTS; j++){
			fscanf(waveStuff,"%f",&dat); //?? 
			int pos = j*SAMPLES + i;
			//waveData[j][i]=(double)dat*1000.0;
			waveData[pos]=(double)dat*1000.0;
		}
	}
	
	fclose(waveStuff);
	
	//openlog("CUDAmia", NULL, 0);
	printf("WIDTH: %d\n", WIDTH);
	/****************************Calling makeAllHamming kernel*************************/
	
	//hammingMat at host
	//double *hammingMat = (double *)malloc(KEYBYTES*SAMPLES*sizeof(double));
	double *hammingMat = (double *)malloc(KEYS*KEYBYTES*SAMPLES*sizeof(double));
	double *dev_hammingMat;
	cudaMalloc((void **)&dev_plainTexts, KEYBYTES*SAMPLES*sizeof(unsigned int));
	cudaMemcpy(dev_plainTexts, plainTexts, KEYBYTES*SAMPLES*sizeof(unsigned int), cudaMemcpyHostToDevice);
	//cudaMalloc((void **)&dev_hammingMat, KEYBYTES*SAMPLES*sizeof(double));
	cudaMalloc((void **)&dev_hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double));
	
	dim3 numBlocks(16,1); //Blocks in the grid
	dim3 numThreads(16,16); // Threads per block
	
	makeAllHamming<<<numBlocks,numThreads>>>(dev_plainTexts, dev_hammingMat, width);
	cudaDeviceSynchronize();
	
	//cudaMemcpy(hammingMat, dev_hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double), cudaMemcpyDeviceToHost);
	
	//cudaFree(dev_hammingMat);
	cudaFree(dev_plainTexts);
	printf("make hamming done\n");
	
	/*****************************************************************************************
	 ******************************************************************************************/
	

	 /* **************************** KDE_bandwidthHammingKernel ***************************** */
	int targetLengthH = 12; 
	double *dev_banwWidthHamming;
	cudaMalloc((void **)&dev_banwWidthHamming, KEYBYTES*KEYS*sizeof(double)); checkCudaError();
	
	double *dev_targetHamming;
	cudaMalloc((void **)&dev_targetHamming, KEYBYTES*KEYS*targetLengthH*sizeof(double)); checkCudaError();
	
	double *dev_maxMinHamming;
	cudaMalloc((void **)&dev_maxMinHamming, KEYBYTES*KEYS*2*sizeof(double)); checkCudaError();
	
	
	
	double *banwWidthHamming = (double *)malloc(KEYS*KEYBYTES*sizeof(double));
	double *targetHamming = (double *)malloc(KEYBYTES*KEYS*targetLengthH*sizeof(double));
	
	dim3 numBlocksBWH(16,1); //Blocks in the grid
	dim3 numThreadsBWH(16,16); // Threads per block
	
	KDE_bandwidthHammingKernel<<<numBlocksBWH, numThreadsBWH>>>(dev_hammingMat, dev_banwWidthHamming, dev_maxMinHamming,
																targetLengthH, dev_targetHamming); checkCudaError();
	//SEEMS FINE!!															
	cudaMemcpy(banwWidthHamming, dev_banwWidthHamming, KEYBYTES*KEYS*sizeof(double), cudaMemcpyDeviceToHost);checkCudaError();
	
	/*FILE * fpP = freopen("results.txt", "w", stdout);
	int ii, jj;
	for(ii=0;ii<KEYS;ii++){
		for(jj=0;jj<KEYBYTES;jj++){
			int post = ii*KEYBYTES+jj;
			printf("%lf  ",banwWidthHamming[post]);
		}
			printf("\n");
	}
	fclose(fpP);*/
	
	printf("bandwidth Hamming done\n");
	free(banwWidthHamming);
	free(targetHamming);
	/* **************************** KDE_bandwidthWaveKernel ***************************** */
	int targetLengthW = 256; 
	double * dev_waveData;
	cudaMalloc((void **)&dev_waveData, SAMPLEPOINTS * SAMPLES*sizeof(double)); checkCudaError();
	cudaMemcpy(dev_waveData, waveData, SAMPLEPOINTS * SAMPLES*sizeof(double), cudaMemcpyHostToDevice);checkCudaError();
	
	double *dev_banwWidthWave;
	cudaMalloc((void **)&dev_banwWidthWave, SAMPLEPOINTS*sizeof(double)); checkCudaError();
	
	double *dev_targetWave;
	cudaMalloc((void **)&dev_targetWave, SAMPLEPOINTS*targetLengthW*sizeof(double)); checkCudaError();
	
	double *dev_maxMinWave;
	cudaMalloc((void **)&dev_maxMinWave, SAMPLEPOINTS*2*sizeof(double)); checkCudaError();

	double *banwWidthWave = (double *)malloc(SAMPLEPOINTS*sizeof(double));
	double *targetWave = (double *)malloc(SAMPLEPOINTS*targetLengthW*sizeof(double));
	
	dim3 numBlocksBWW(100, 1, 1);
	dim3 numThreadsBWW(1000, 1, 1);
	
	KDE_bandwidthWaveKernel<<<numBlocksBWW, numThreadsBWW>>>(dev_waveData, dev_banwWidthWave, dev_maxMinWave,
																targetLengthW, dev_targetWave);checkCudaError();
	cudaMemcpy(banwWidthWave, dev_banwWidthWave, SAMPLEPOINTS*sizeof(double), cudaMemcpyDeviceToHost);checkCudaError();
	
	/*FILE * fpP = freopen("results.txt", "w", stdout);
	int ii, jj;
	for(ii=0;ii<SAMPLEPOINTS;ii++){
		//for(jj=0;jj<SAMPLES;jj++){
			//int post = ii*SAMPLES+jj;
			printf("%lf  ",banwWidthWave[ii]);
		//}
			printf("\n");
	}
	fclose(fpP);*/
	printf("bandwidth Wave done\n");
	free(banwWidthWave);
	/* **************************** KDE_findProbsHamming ******************************** */
	
	double *dev_hammingDataProbs;
	cudaMalloc((void **)&dev_hammingDataProbs, KEYBYTES*KEYS*targetLengthH*sizeof(double)); checkCudaError();
	
	double *hammingDataProbs = (double *)malloc(KEYBYTES*KEYS*targetLengthH*sizeof(double));
	
	dim3 numBlocksProbsHKDE(256, 1, 1); 
	dim3 threadsPerBlocksProbsHKDE(16,1,1);
	
	KDE_findProbsHamming<<<numBlocksProbsHKDE, threadsPerBlocksProbsHKDE>>>(dev_hammingMat, dev_hammingDataProbs, 
							dev_banwWidthHamming, targetLengthH, dev_maxMinHamming);
	
	//cudaMemcpy(hammingDataProbs, dev_hammingDataProbs, KEYBYTES*KEYS*targetLengthH*sizeof(double), cudaMemcpyDeviceToHost);
	
	//FILE * fpP = freopen("results.txt", "w", stdout);
	//printMatDouble(hammingDataProbs, KEYS, KEYBYTES*targetLengthH);
	printf("prob hamming done\n");
	//fclose(fpP);
	free(hammingDataProbs);
	
	/* **************************** KDE_findProbsWave ******************************** */
	
	double *dev_waveDataProbs;
	cudaMalloc((void **)&dev_waveDataProbs, SAMPLEPOINTS*targetLengthW*sizeof(double)); checkCudaError();
	
	double *waveDataProbs = (double *)malloc(SAMPLEPOINTS*targetLengthW*sizeof(double));
	
	dim3 numBlocksProbsWKDE(16, 16, 1); 
	dim3 threadsPerBlocksProbsWKDE(512, 1, 1);
	
	KDE_findProbsWave<<<numBlocksProbsWKDE, threadsPerBlocksProbsWKDE>>>(dev_waveData, dev_waveDataProbs, 
							dev_banwWidthWave, targetLengthW, dev_maxMinWave);
	checkCudaError();						
	//cudaMemcpy(waveDataProbs, dev_waveDataProbs, SAMPLEPOINTS*targetLengthW*sizeof(double), cudaMemcpyDeviceToHost);
	
	//FILE * fpP = freopen("results.txt", "w", stdout);
	//printMatDouble(waveDataProbs, SAMPLEPOINTS, targetLengthW);
	printf("prob wave done\n");
	//fclose(fpP);
	free(waveDataProbs);
	
	cudaFree(dev_maxMinWave);
	cudaFree(dev_maxMinHamming);
	
	/* **************************** KDE_findProbsJoint2 ******************************** */
	int targetLengthJoint = targetLengthW*targetLengthH;
	double *dev_jointProbs;
	cudaMalloc((void **)&dev_jointProbs, (SAMPLEPOINTS/DIVIDE)*targetLengthJoint*sizeof(double)); checkCudaError();
	
	double *jointProbs = (double *)malloc((SAMPLEPOINTS/DIVIDE)*targetLengthJoint*sizeof(double));checkCudaError();
	
	int turn = 0;
	int key = 0;
	int keyByte = 0;
	
	dim3 blocksProbsJointKDE(ceil(targetLengthH/4.0),ceil(targetLengthW/16.0),SAMPLEPOINTS/(DIVIDE*4));
	//dim3 blocksProbsJointKDE(3,16,12500);
	dim3 threadsProbsJointKDE(4,16,4);
	
	/* KDE_MI kernel */
	dim3 blocksMIJointKDE(200,0,0);
	dim3 threadsProbsJointKDE(256,0,0);
	
	cudaEvent_t start,stop;
	float elapsedtime;
	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	
	for(key=0; key<1; key++){
		for(i=0; i<1; i++){
			
			KDE_findJointProbs2<<<blocksProbsJointKDE, threadsProbsJointKDE>>>(dev_hammingMat, dev_waveData, targetLengthW, 
									targetLengthH, dev_banwWidthWave, dev_banwWidthHamming, dev_targetHamming, dev_targetWave, 
									dev_jointProbs,turn, key, keyByte);
			cudaThreadSynchronize();						
			checkCudaError();
			
			turn = turn + 50000;
			//cudaMemset(dev_jointProbs, 0, (SAMPLEPOINTS/DIVIDE)*targetLengthJoint*sizeof(double)); checkCudaError();
			
			cudaMemcpy(jointProbs, dev_jointProbs, (SAMPLEPOINTS/DIVIDE)*targetLengthJoint*sizeof(double), cudaMemcpyDeviceToHost);
			checkCudaError();
			
		}
	}
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedtime,start,stop);
	fprintf(stdout,"Time spent for operation : %.10f seconds\n",elapsedtime/(float)1000);
	
	FILE * fpP = freopen("results.txt", "w", stdout);
	int ii, jj;
	for(ii=0;ii<SAMPLEPOINTS/DIVIDE;ii++){
		for(jj=0;jj<targetLengthJoint;jj++){
			int post = ii*targetLengthJoint+jj;
			printf("%lf  ",jointProbs[post]);
		}
			printf("\n");
	}
	fclose(fpP);
	free(jointProbs);
	printf("prob joint done\n");

	cudaFree(dev_waveDataProbs);
	cudaFree(dev_hammingDataProbs);
	cudaFree(dev_banwWidthWave);
	cudaFree(dev_banwWidthHamming);
	cudaFree(dev_hammingMat);
	cudaFree(dev_waveData);
	cudaFree(dev_targetWave);
	cudaFree(dev_targetHamming);
	cudaFree(dev_jointProbs);
	
	/*** NO ERRORS!! ***/
	
	/*****************************************************************************************
	 ******************************************************************************************/	
	
	
	/****************** Calling normaliseWaveDataKernel on waveData ********************/
	/*
	//Store wave data in device
	//double * dev_waveData;
	//cudaMalloc((void **)&dev_waveData, SAMPLEPOINTS * SAMPLES*sizeof(double));
	//cudaMemcpy(dev_waveData, waveData, SAMPLEPOINTS * SAMPLES*sizeof(double), cudaMemcpyHostToDevice);
	
	//Store normalised wave data
	unsigned int * waveDataNormalised = (unsigned int *)malloc(sizeof(unsigned int) * SAMPLEPOINTS * SAMPLES);
	unsigned int * dev_waveDataNormalised;
	cudaMalloc((void **)&dev_waveDataNormalised, SAMPLEPOINTS * SAMPLES*sizeof(unsigned int));
	
	//Store firstNumState
	unsigned int * dev_firstNumState;
	cudaMalloc((void **)&dev_firstNumState, sizeof(unsigned int)*SAMPLEPOINTS);
	unsigned int * firstNumState = (unsigned int *)malloc(sizeof(unsigned int) * SAMPLEPOINTS);
	
	dim3 numBlocksNorm(100, 1, 1);
	dim3 numThreadsNorm(1000, 1, 1);
	
	normaliseWaveDataKernel<<<numBlocksNorm, numThreadsNorm>>>(dev_waveData, dev_waveDataNormalised, dev_firstNumState);
	cudaDeviceSynchronize();
	
	cudaMemcpy(waveDataNormalised, dev_waveDataNormalised, SAMPLEPOINTS*SAMPLES*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(firstNumState, dev_firstNumState, SAMPLEPOINTS*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	
	cudaFree(dev_waveData);
	//cudaFree(dev_waveDataNormalised);
	
	printf("normalise wave done\n");*/
	
	/***********************Calling normalisingHammingKernel*********************************/
	/*
	//Hamming data
	double * dev_hamming;
	cudaMalloc((void **)&dev_hamming, KEYS*KEYBYTES*SAMPLES*sizeof(double));
	cudaMemcpy(dev_hamming, hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double), cudaMemcpyHostToDevice);
	
	//Result matrices
	unsigned int * hammingNormalised = (unsigned int *)malloc(sizeof(unsigned int) * KEYS*KEYBYTES*SAMPLES);
	unsigned int * dev_hammingNormalised;
	cudaMalloc((void **)&dev_hammingNormalised, KEYS*KEYBYTES*SAMPLES*sizeof(unsigned int));

	//Store secondNumState
	unsigned int * dev_secondNumState;
	cudaMalloc((void **)&dev_secondNumState, sizeof(unsigned int)*KEYS*KEYBYTES);
	unsigned int * secondNumState = (unsigned int *)malloc(sizeof(unsigned int) * KEYS*KEYBYTES);
	
	dim3 numBlocksHamming(256, 1, 1);
	dim3 numThreadsHamming(16, 1, 1);
	
	normaliseHammingKernel<<<numBlocksHamming, numThreadsHamming>>>(dev_hamming, dev_hammingNormalised, dev_secondNumState);
	cudaDeviceSynchronize();
	
	cudaMemcpy(hammingNormalised, dev_hammingNormalised, KEYS*KEYBYTES*SAMPLES*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(secondNumState, dev_secondNumState, KEYS*KEYBYTES*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	
	cudaFree(dev_hamming);
	printf("Hamming normalised, dev_hamming removed\n");*/
	
	/***********Sorting***************/
	/*
	double *dev_sortedMI;
	cudaMalloc((void **)&dev_sortedMI, KEYBYTES*KEYS*sizeof(double)); 
	cudaMemset(dev_sortedMI, 0, KEYBYTES*KEYS*sizeof(double));
	double *sortedMI = (double *)malloc(KEYBYTES*KEYS*sizeof(double));
	
	dim3 numBlockSort(1,1,1);
	dim3 threadsPerBlockSort(16,1,1);
	
	for (i=0;i<NUMOFKEYS;i++){
		simpleSort<<<numBlockSort, threadsPerBlockSort>>>(dev_MIvals, 100000, dev_sortedMI, i); checkCudaError();
	}
	cudaMemcpy(sortedMI, dev_sortedMI, KEYBYTES*KEYS*sizeof(double), cudaMemcpyDeviceToHost); checkCudaError();
	
	
	
	FILE * fpP = freopen("results_MISorted.txt", "w", stdout);
	int pp1,pp2;
	for(pp1=0;pp1<NUMOFKEYS;pp1++){
		for(pp2=0;pp2<KEYBYTES;pp2++){
			printf("%lf  ",sortedMI[pp1*KEYBYTES + pp2]);
		}
		printf("\n");
	}
	fclose(fpP);*/
	/****************************************************************************************************************************/
	//free(firstStateProbs);
	//free(secondStateProbs);
	//free(jointProbs);
	
	//firstStateProbs = NULL;
	//secondStateProbs = NULL;
	//jointProbs = NULL;
	
	free(plainTexts);
	free(hammingMat);
	free(waveData);
	//free(waveDataNormalised);
	//free(hammingNormalised);
	//free(MIvals);
	//free(sortedMI);
	
	/*cudaFree(dev_firstStateProbs);
	cudaFree(dev_secondStateProbs);
	
	cudaFree(dev_jointProbs);
	
	cudaFree(dev_waveDataNormalised);
	printf("dev_waveDataNormalised done\n");

	cudaFree(dev_hammingNormalised);
	printf("dev_hammingNormalised done\n");

	cudaFree(dev_firstNumState);
	printf("dev_firstNumState done\n");

	cudaFree(dev_secondNumState);
	printf("dev_secondNumState done\n");
	
	cudaFree(dev_MIvals);
	printf("dev_secondNumState done\n");
	printf("*********************END OF PROGRAM*********************\n");
	
	cudaFree(dev_sortedMI);*/
	
	return 0;
}
