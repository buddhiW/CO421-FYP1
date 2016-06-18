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
#include "data.cuh"
#include "arrayOperations.cuh"
#include "probOps.cuh"
#include "helpers.cuh"

#define PLAINFILE "plain.txt"
//#define WAVEFILE "wave.txt"
#define WAVEFILE "/home/e11444/CO421/test2/wave.txt"

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

//__global__ void makeHamming(unsigned int *plaintxt, double *hammingRes, int n){//n is the size of samples 
	
	//unsigned int key = 0x01;
	//int index = blockIdx.x*KEYBYTES + threadIdx.x;
	
	//if (index<n){
		//byte temp = hamming(plaintxt,index, key);
		//hammingRes[index] = temp;
	//}
//} 

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
//	Find maximum of each MI set
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

int arrayMax(unsigned int * data, int length){
	
	int max = data[0];
	for(int i=0; i<length; i++){
		if(data[i]>max)
			max = data[i];
	}
	return max;
}

int main(){

	//cudaSetDevice(1);
	


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
	
	cudaMemcpy(hammingMat, dev_hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double), cudaMemcpyDeviceToHost);
	
	cudaFree(dev_hammingMat);
	cudaFree(dev_plainTexts);
	printf("make hamming done\n");
	
	/****************** Calling normaliseWaveDataKernel on waveData ********************/

	//Store wave data in device
	double * dev_waveData;
	cudaMalloc((void **)&dev_waveData, SAMPLEPOINTS * SAMPLES*sizeof(double));
	cudaMemcpy(dev_waveData, waveData, SAMPLEPOINTS * SAMPLES*sizeof(double), cudaMemcpyHostToDevice);
	
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
	
	printf("normalise wave done\n");
	
	/***********************Calling normalisingHammingKernel*********************************/
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
	printf("Hamming normalised, dev_hamming removed\n");
	
	/***************************************Calling findProbWave******************************************/
	int maxOfFirst = arrayMax(firstNumState, SAMPLEPOINTS);
	
	double *dev_firstStateProbs;
	cudaMalloc((void **)&dev_firstStateProbs, SAMPLEPOINTS*maxOfFirst*sizeof(double));	checkCudaError();
	
	double *firstStateProbs = (double *)malloc(sizeof(double)*SAMPLEPOINTS*maxOfFirst);

	dim3 numBlocksProbsWave(16, 16, 1);
	dim3 threadsPerBlocksProbsWave(512, 1, 1);
	
	findProbsWave<<<numBlocksProbsWave, threadsPerBlocksProbsWave>>>(dev_waveDataNormalised, dev_firstStateProbs, dev_firstNumState, maxOfFirst);		checkCudaError(); 
	cudaDeviceSynchronize();
	/***************************************Calling findProbsHamming******************************************/
	int maxOfSecond = arrayMax(secondNumState, KEYS*KEYBYTES);
	
	double *dev_secondStateProbs;
	
	cudaMalloc((void **)&dev_secondStateProbs, KEYS*KEYBYTES*maxOfSecond*sizeof(double));	checkCudaError();
	double *secondStateProbs = (double *)malloc(sizeof(double)*KEYS*KEYBYTES*maxOfSecond);
	
	//we need 16 parallel operations
	dim3 numBlocksProbsHamming(256, 1, 1); 
	dim3 threadsPerBlocksProbsHamming(16,1,1);
	
	findProbsHamming<<<numBlocksProbsHamming, threadsPerBlocksProbsHamming>>>(dev_hammingNormalised, dev_secondStateProbs, dev_secondNumState, maxOfSecond);	checkCudaError();
	cudaDeviceSynchronize();
	
	printf("Hamming probabilities calculated\n");
	cudaMemcpy(secondStateProbs, dev_secondStateProbs, KEYS*KEYBYTES*maxOfSecond*sizeof(double), cudaMemcpyDeviceToHost);
	
	//FILE * fpP = freopen("results_prob_hamming.txt", "w", stdout);
	//printMatDouble(secondStateProbs,100,maxOfSecond*KEYBYTES);
	//fclose(fpP);
	/*************************************************************Calling joint Probs***************************************************/
	int maxOfJoint = maxOfFirst*maxOfSecond;
	printf("Max of first: %d\n",maxOfFirst);
	printf("Max of second: %d\n",maxOfSecond);
	printf("Max of joint: %d\n",maxOfJoint);

	double *dev_jointProbs;
	
	int threads = SAMPLEPOINTS/REPEAT;
	unsigned int index = 0;
	
	//cudaMalloc((void **)&dev_jointProbs, KEYBYTES*SAMPLEPOINTS*maxOfJoint*sizeof(double));
	cudaMalloc((void **)&dev_jointProbs, KEYBYTES*threads*maxOfJoint*sizeof(double)); 
	double *jointProbs = (double *)malloc(KEYBYTES*threads*maxOfJoint*sizeof(double));
	
	//double *dev_MIvals;
	//cudaMalloc((void **)&dev_MIvals, KEYBYTES*SAMPLEPOINTS*sizeof(double)); 
	//double *MIvals = (double *)malloc(KEYBYTES*SAMPLEPOINTS*sizeof(double));
	
	double *dev_MIvals;
	cudaMalloc((void **)&dev_MIvals, KEYS*KEYBYTES*SAMPLEPOINTS*sizeof(double)); 
	cudaMemset(dev_MIvals, 0, KEYS*KEYBYTES*SAMPLEPOINTS*sizeof(double));
	double *MIvals = (double *)malloc(KEYS*KEYBYTES*SAMPLEPOINTS*sizeof(double));
	
	//dim3 numBlocksJointProb(128,16,1);
	//dim3 threadsPerBlockJointProb(1024,1,1);
	dim3 numBlocksJointProb(32,16,1);
	dim3 threadsPerBlockJointProb(256,1,1); //OPTIMUM: 256 threads per block. Occupancy reduce for lesser value.
	
	/**findJointProbs occupancy is maximum**/
	
	cudaEvent_t start,stop;
	float elapsedtime;
	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	
	for(i=0; i<NUMOFKEYS; i++){
		int startPosition = 0;
		for(j=0; j<REPEAT; j++){
			//printf("This is the Start position: %d\n",startPosition);
			
			cudaMemset(dev_jointProbs, 0, KEYBYTES*threads*maxOfJoint*sizeof(double));
			
			
			findJointProbs<<<numBlocksJointProb, threadsPerBlockJointProb>>>(dev_MIvals, dev_waveDataNormalised, dev_hammingNormalised, dev_firstNumState, 
						dev_secondNumState, dev_firstStateProbs, dev_secondStateProbs, dev_jointProbs, 
						maxOfFirst, maxOfSecond, i, startPosition); 	
			cudaDeviceSynchronize();
			checkCudaError();
			//if(i==REPEAT-1)	
			//	cudaMemcpy(jointProbs, dev_jointProbs, (KEYBYTES*threads*maxOfJoint)*sizeof(double), cudaMemcpyDeviceToHost); checkCudaError();
			
			startPosition = startPosition + threads;
			//index = index + (KEYBYTES*threads*maxOfJoint);
			//printf("index: %ul\n", index);
		}
	}
	
	printf("MI cal done\n");
	cudaMemcpy(MIvals, dev_MIvals, KEYBYTES*SAMPLEPOINTS*sizeof(double), cudaMemcpyDeviceToHost); checkCudaError();
	
	/*int oo,aa;
	for(oo=0;oo<2;oo++){
		for(aa=0;aa<20;aa++){
			printf("MI[%d] = %lf\n",(oo*SAMPLEPOINTS + aa),MIvals[oo*SAMPLEPOINTS + aa]);
		}
		printf("==================================================\n");
	}*/
	
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedtime,start,stop);
	fprintf(stdout,"Time spent for operation : %.10f seconds\n",elapsedtime/(float)1000);
	
	/***********Sorting***************/
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
	fclose(fpP);
	/****************************************************************************************************************************/
	free(firstStateProbs);
	free(secondStateProbs);
	free(jointProbs);
	
	firstStateProbs = NULL;
	secondStateProbs = NULL;
	jointProbs = NULL;
	
	free(plainTexts);
	free(hammingMat);
	free(waveData);
	free(waveDataNormalised);
	free(hammingNormalised);
	free(MIvals);
	free(sortedMI);
	
	cudaFree(dev_firstStateProbs);
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
	
	cudaFree(dev_sortedMI);
	
	return 0;
}
