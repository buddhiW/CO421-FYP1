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
__global__ void simpleSort(double *origMat, int cols, double *maxOfEach){
	
	int xx = blockIdx.x*blockDim.x + threadIdx.x;//16
	int yy = blockIdx.y*blockDim.y + threadIdx.y;//256
	
	int i;
	double max = 0.0;
	for (i=0;i<cols;i++){
		int pos = yy*KEYBYTES*SAMPLEPOINTS + xx*SAMPLEPOINTS + i;
		if (max<origMat[pos])
			max = origMat[pos];
	}
	maxOfEach[yy*KEYBYTES + xx] = max;
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
	
	
	unsigned int*dev_plainTexts;
	cudaMalloc((void **)&dev_plainTexts, KEYBYTES*SAMPLES*sizeof(unsigned int));	checkCudaError();
	cudaMemcpy(dev_plainTexts, plainTexts, KEYBYTES*SAMPLES*sizeof(unsigned int), cudaMemcpyHostToDevice);	checkCudaError();
	
	double *dev_hammingMat;
	cudaMalloc((void **)&dev_hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double));	checkCudaError();
	
	dim3 numBlocks(16,1); //Blocks in the grid
	dim3 numThreads(16,16); // Threads per block
	
	makeAllHamming<<<numBlocks,numThreads>>>(dev_plainTexts, dev_hammingMat, width);
	cudaDeviceSynchronize();	checkCudaError();
	
	cudaMemcpy(hammingMat, dev_hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double), cudaMemcpyDeviceToHost);	checkCudaError();
	
	cudaFree(dev_hammingMat);	checkCudaError();
	cudaFree(dev_plainTexts);	checkCudaError();
	printf("make hamming done\n");
	
	/****************** Calling normaliseWaveDataKernel on waveData ********************/

	//Store wave data in device
	double * dev_waveData;
	cudaMalloc((void **)&dev_waveData, SAMPLEPOINTS * SAMPLES*sizeof(double));	checkCudaError();
	cudaMemcpy(dev_waveData, waveData, SAMPLEPOINTS * SAMPLES*sizeof(double), cudaMemcpyHostToDevice);	checkCudaError();
	
	//Store normalised wave data
	unsigned int * waveDataNormalised = (unsigned int *)malloc(sizeof(unsigned int) * SAMPLEPOINTS * SAMPLES);
	unsigned int * dev_waveDataNormalised;
	cudaMalloc((void **)&dev_waveDataNormalised, SAMPLEPOINTS * SAMPLES*sizeof(unsigned int));	checkCudaError();
	
	//Store firstNumState
	unsigned int * dev_firstNumState;
	cudaMalloc((void **)&dev_firstNumState, sizeof(unsigned int)*SAMPLEPOINTS);	checkCudaError();
	
	unsigned int * firstNumState = (unsigned int *)malloc(sizeof(unsigned int) * SAMPLEPOINTS);
	
	dim3 numBlocksNorm(100, 1, 1);
	dim3 numThreadsNorm(1000, 1, 1);
	
	normaliseWaveDataKernel<<<numBlocksNorm, numThreadsNorm>>>(dev_waveData, dev_waveDataNormalised, dev_firstNumState);
	cudaDeviceSynchronize();	checkCudaError();
	
	cudaMemcpy(waveDataNormalised, dev_waveDataNormalised, SAMPLEPOINTS*SAMPLES*sizeof(unsigned int), cudaMemcpyDeviceToHost);	checkCudaError();
	cudaMemcpy(firstNumState, dev_firstNumState, SAMPLEPOINTS*sizeof(unsigned int), cudaMemcpyDeviceToHost);	checkCudaError();
	
	cudaFree(dev_waveData);	checkCudaError();
	//cudaFree(dev_waveDataNormalised);
	
	printf("normalise wave done\n");
	
	/***********************Calling normalisingHammingKernel*********************************/
	//Hamming data
	double * dev_hamming;
	cudaMalloc((void **)&dev_hamming, KEYS*KEYBYTES*SAMPLES*sizeof(double));	checkCudaError();
	cudaMemcpy(dev_hamming, hammingMat, KEYS*KEYBYTES*SAMPLES*sizeof(double), cudaMemcpyHostToDevice);	checkCudaError();
	
	//Result matrices
	unsigned int * hammingNormalised = (unsigned int *)malloc(sizeof(unsigned int) * KEYS*KEYBYTES*SAMPLES);
	unsigned int * dev_hammingNormalised;
	cudaMalloc((void **)&dev_hammingNormalised, KEYS*KEYBYTES*SAMPLES*sizeof(unsigned int));	checkCudaError();

	//Store secondNumState
	unsigned int * dev_secondNumState;
	cudaMalloc((void **)&dev_secondNumState, sizeof(unsigned int)*KEYS*KEYBYTES);
	unsigned int * secondNumState = (unsigned int *)malloc(sizeof(unsigned int) * KEYS*KEYBYTES);
	
	dim3 numBlocksHamming(256, 1, 1);
	dim3 numThreadsHamming(16, 1, 1);
	
	normaliseHammingKernel<<<numBlocksHamming, numThreadsHamming>>>(dev_hamming, dev_hammingNormalised, dev_secondNumState);
	cudaDeviceSynchronize();	checkCudaError();
	
	cudaMemcpy(hammingNormalised, dev_hammingNormalised, KEYS*KEYBYTES*SAMPLES*sizeof(unsigned int), cudaMemcpyDeviceToHost);	checkCudaError();
	cudaMemcpy(secondNumState, dev_secondNumState, KEYS*KEYBYTES*sizeof(unsigned int), cudaMemcpyDeviceToHost);	checkCudaError();
	
	cudaFree(dev_hamming);	checkCudaError();
	printf("Hamming normalised, dev_hamming removed\n");
	
	/***************************************Calling findProbWave******************************************/
	int maxOfFirst = arrayMax(firstNumState, SAMPLEPOINTS);
	
	double *dev_firstStateProbs;
	cudaMalloc((void **)&dev_firstStateProbs, SAMPLEPOINTS*maxOfFirst*sizeof(double));	checkCudaError();
	
	double *firstStateProbs = (double *)malloc(sizeof(double)*SAMPLEPOINTS*maxOfFirst);

	dim3 numBlocksProbsWave(16, 16, 1);
	dim3 threadsPerBlocksProbsWave(512, 1, 1);
	
	findProbsWave<<<numBlocksProbsWave, threadsPerBlocksProbsWave>>>(dev_waveDataNormalised, dev_firstStateProbs, dev_firstNumState, maxOfFirst);		checkCudaError(); 
	cudaDeviceSynchronize(); checkCudaError();
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
	/*************************************************************Calling joint Probs***************************************************/
	/*calculating joint sizes*/
	unsigned long *jointSizes = (unsigned long*)malloc((KEYS/32)*KEYBYTES*SAMPLEPOINTS*sizeof(unsigned long)); //for 8 keys @ a momment
	int i1,i2,i3;
	
	double *dev_MIvals;
	cudaMalloc((void **)&dev_MIvals, KEYS*KEYBYTES*SAMPLEPOINTS*sizeof(double));   checkCudaError();
	cudaMemset(dev_MIvals, 0, KEYS*KEYBYTES*SAMPLEPOINTS*sizeof(double));	checkCudaError();
	
	double *MIvals = (double *)malloc(KEYS*KEYBYTES*SAMPLEPOINTS*sizeof(double));
	
	unsigned long *dev_jointSizes;
	cudaMalloc((void **)&dev_jointSizes, (KEYS/32)*KEYBYTES*SAMPLEPOINTS*sizeof(unsigned long)); 
	
	cudaEvent_t start,stop;
	float elapsedtime;
	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	
	/*confirm the correctness of this*/
	dim3 block3d(16,2,4);
	dim3 grid3d(1,4,25000);
		
	int turn;
	for (turn=0;turn<8*32;turn+=8){
		
		unsigned long totalSize = 0l;
	
		for(i1=0;i1<(KEYS/32);i1++){
			for(i2=0;i2<KEYBYTES;i2++){
				for(i3=0;i3<SAMPLEPOINTS;i3++){
					totalSize += (long)firstNumState[i3]*(long)secondNumState[(turn+i1)*KEYBYTES + i2];
					if (totalSize>7466209652) printf("THis sucksssssss\n");
					jointSizes[i1*KEYBYTES*SAMPLEPOINTS + i2*SAMPLEPOINTS + i3] = totalSize;
				}
			}
		}
	
		unsigned long ttl = totalSize*8;
		//printf("This is the total size: %ld bytes\n",ttl);
		//free(jointSizes);
	
		cudaMemcpy(dev_jointSizes, jointSizes, (KEYS/32)*KEYBYTES*SAMPLEPOINTS*sizeof(unsigned long), cudaMemcpyHostToDevice);  checkCudaError();
	
		double *dev_jointProbs;
		cudaMalloc((void **)&dev_jointProbs, totalSize*sizeof(double));   checkCudaError();
		
		//double *jointProbs = (double *)malloc(totalSize*sizeof(double));
	
		
	
		/**findJointProbs occupancy is maximum**/
	
		
	
		//cudaMemset(dev_jointProbs, 0, (KEYS/4)*KEYBYTES*SAMPLES*sizeof(double));   checkCudaError(); //??
		cudaMemset(dev_jointProbs, 0, totalSize*sizeof(double));   checkCudaError();
		
		printf("Now starting joint probs and MI calculation for turn = %d\n",turn);
		
		findJointProbs2<<<grid3d, block3d>>>(dev_MIvals, dev_waveDataNormalised, dev_hammingNormalised, 
											 dev_firstNumState, dev_secondNumState, dev_firstStateProbs, dev_secondStateProbs, dev_jointProbs, 
											 dev_jointSizes, maxOfFirst, maxOfSecond, turn); 		cudaDeviceSynchronize();		checkCudaError();
	
		cudaFree(dev_jointProbs);
	}	
	printf("MI cal done\n");
	
	//cudaMemcpy(MIvals, dev_MIvals, KEYBYTES*SAMPLEPOINTS*sizeof(double), cudaMemcpyDeviceToHost); checkCudaError();
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedtime,start,stop);
	fprintf(stdout,"Time spent for operation : %.10f seconds\n",elapsedtime/(float)1000);
	
	
	free(firstStateProbs);
	free(secondStateProbs);
	
	
	firstStateProbs = NULL;
	secondStateProbs = NULL;
	
	free(plainTexts);
	free(hammingMat);
	free(waveData);
	free(waveDataNormalised);
	free(hammingNormalised);
	
	cudaFree(dev_firstStateProbs);
	cudaFree(dev_secondStateProbs);
	
	
	
	cudaFree(dev_waveDataNormalised);
	printf("dev_waveDataNormalised done\n");

	cudaFree(dev_hammingNormalised);
	printf("dev_hammingNormalised done\n");

	cudaFree(dev_firstNumState);
	printf("dev_firstNumState done\n");

	cudaFree(dev_secondNumState);
	printf("dev_secondNumState done\n");
	
	/***********Sorting***************/
	printf("sorting\n");
	double *dev_sortedMI;
	cudaMalloc((void **)&dev_sortedMI, KEYBYTES*KEYS*sizeof(double)); 
	
	double *sortedMI = (double *)malloc(KEYBYTES*KEYS*sizeof(double));
	
	dim3 block3dS(16,16,1);
	dim3 grid3dS(1,16,1);
	
	simpleSort<<<grid3dS, block3dS>>>(dev_MIvals, SAMPLEPOINTS, dev_sortedMI); 	cudaDeviceSynchronize(); 	checkCudaError();
	
	cudaMemcpy(sortedMI, dev_sortedMI, KEYBYTES*KEYS*sizeof(double), cudaMemcpyDeviceToHost); checkCudaError();
	
	FILE * fpP = freopen("results_MISorted.txt", "w", stdout);
	
	int pp1,pp2;
	
	for(pp1=0;pp1<KEYS;pp1++){
		for(pp2=0;pp2<KEYBYTES;pp2++){
			printf("%lf  ",sortedMI[pp1*KEYBYTES + pp2]);
		}
		printf("\n");
	}
	fclose(fpP);
	/****************************************************************************************************************************/
	
	cudaFree(dev_MIvals);
	cudaFree(dev_sortedMI);
	
	free(MIvals);
	free(sortedMI);
	
	printf("dev_MI_vals done\n");
	printf("*********************END OF PROGRAM*********************\n");
	return 0;
}
