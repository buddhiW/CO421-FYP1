//-------------------------------------------------------------------------
//	MIA based power analysis attack on the CPU
//
//	Version information -- working
//
//	Author			Details
//	MP				Initial version
//	BW				Modularised some code into functions
//					Moved plaiText matrix to heap	
//					Implemented Threaded version 2 (Sample point calculation is threaded)	
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include "arrayOperations.h"
#include "data.h"
#include "MI.h"

//number of sample points in one encryption
#define SAMPLEPOINTS 100000
#define SAMPLES 80
//defined for 128 bit AES
#define KEYBYTES 16
//#define KEYS 256
#define KEYS 256

#define THREADS 32

//argument struct for pthread
typedef struct args{
	double **waveData;
	unsigned int **plainText;
	//double **maxMI;
	unsigned int startPosition;
	int threadNumber;
}threadArgs;

//argument struct to find maximum
typedef struct max{
	unsigned int keyByte;
}maxArgs;

double maxMIArray[KEYS][KEYBYTES] = {0}; //Stores final result
double maxMIAll[THREADS][KEYS][KEYBYTES] = {0}; //Stores intermediate results of all threads

//from Hasindu's code
//calculates hamming weight of a 8 bit number
byte hammingweight(byte H){
	// Count the number of set bits
	byte dist=0;
	while(H){
		dist++; 
		H &= H - 1;
	}
	return dist;
}

/********************************************************************** SELECTION FUNCTION ****************************************************************************/
//find hamming weight for the selection function
byte hamming(unsigned int sample[], unsigned int n,unsigned int key) { //n is byteno sample is the line in sample text
    byte inter = (byte)sbox[sample[n] ^ key];
    byte dist = hammingweight(inter);		
	return dist;
}

//---------------------------------------------------------------------------------------
//	calHammingMat
//	This function calculates hamming matrix (KEYBYTES * SAMPLES) which contains hamming weight
//  of each sample text byte for a given key guess
//---------------------------------------------------------------------------------------
void calHammingMat(byte keyGuess, double hammingMat[KEYBYTES][SAMPLES], unsigned int **plainText){
	//double hammingMat[KEYBYTES][SAMPLES];
	int i,j;
	for (i=0;i<SAMPLES;i++){
		for(j=0;j<KEYBYTES;j++){
			hammingMat[j][i] = (double)hamming(plainText[i], j, keyGuess);	
		}
	}
	
}

//------------------------------------------------------------------------------------------
//	calHammingArray
//------------------------------------------------------------------------------------------
void calHammingArray(byte keyGuess, double hammingArray[SAMPLES], unsigned int **plainText, unsigned int keyByte){
	//double hammingMat[KEYBYTES][SAMPLES];
	int i,j;
	for (i=0;i<SAMPLES;i++){
		hammingArray[i] = (double)hamming(plainText[i], keyByte, keyGuess);	
	}

}

//------------------------------------------------------------------------
// maxMutualInfo 
// Serial implementation of MIA based power analysis attack
// maxMIArray: 16 * 256 array
// 			   nth column gives maximum MI values for each key guess for 
//				nth key byte
//------------------------------------------------------------------------
void maxMutualInfo(double ** waveData, unsigned int **plainText, double maxMIArray[KEYS][KEYBYTES]){
	int l,k,m;
	double maxMI = 0.0;
	int col;	
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	
	for(m=0; m<KEYS; m++){				
		calHammingMat(m, hammingMat, plainText); //Calculate hamming weights for that key guess		
		for (l=0;l<KEYBYTES;l++){
			maxMI = 0.0;		
			//printf("This is the %d hamming array\n",l);
			for(k=0;k<SAMPLEPOINTS;k++){
				double mi = MI(waveData[k], hammingMat[l], SAMPLES);
				//printf("MI of %d = %lf\n",k,mi);
				
				if (mi>maxMI){
					maxMI = mi;
					col = k;
				}
			}
			maxMIArray[m][l] = maxMI;
		}
	}
}



//Hasindu's implementation
/*void *thread_process(void *voidargs){
	int t,i,j;
	struct args *myargs=(struct args *)voidargs;
	int count=0;
	for(t=myargs->position;t<(myargs->position)+KEYBYTES*KEYS/THREADS;t++){
		i=t/KEYBYTES;
		j=t-t/KEYBYTES*KEYBYTES;
		corelation[i][j]=maxCorelation(myargs->wavedata, myargs->sample, i, j);
		if(myargs->position==0){
			fprintf(stderr,"%d of %d completed\n",count,(KEYBYTES*KEYS/THREADS));
		}
		count++;
	}
	//fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
	pthread_exit(0);
}*/


//--------------------------------------------------------------------------
//	threaded_mia
//	Parallel version 2 => based on sample points
//--------------------------------------------------------------------------
void * threaded_mia(void * args){
	
	threadArgs * myargs = (threadArgs *) args;
	double ** waveData = myargs -> waveData;
	unsigned int **plainText = myargs -> plainText;
	//double ** maxMIArray = myargs -> maxMI;
	unsigned int startPosition = myargs -> startPosition;
	int threadNumber = myargs -> threadNumber;
	
	unsigned int m, k, l;
	int index = 0;
	double maxMI = 0.0;
	
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	unsigned int endPosition;
	double result[KEYS][KEYBYTES];
	
	endPosition = startPosition + (SAMPLEPOINTS/THREADS);

	//printf("thread no: %d\n", threadNumber);
	for(m=0; m<KEYS; m++){	
				
		calHammingMat(m, hammingMat, plainText); //Calculate hamming weights for that key guess		
		
		for (l=0;l<KEYBYTES;l++){
			maxMI = 0.0;		
			for (k=startPosition; k<endPosition; k++ ){
				
				double mi = MI(waveData[k], hammingMat[l], SAMPLES);	
				
				if (mi>maxMI){
					maxMI = mi;
				}	
			}
			
			maxMIAll[threadNumber][m][l] = maxMI;
		}
	}
}

//----------------------------------------------------------------------------------
//	threaded_max
//	Calculate final result matrix 
//----------------------------------------------------------------------------------
void * threaded_max(void * keyByteArg){
	//unsigned int * keyByteAddr = (unsigned int *)keyByteArg;
	
	maxArgs * myMaxArgs = (maxArgs *) keyByteArg;
	unsigned int keyByte = myMaxArgs -> keyByte;
	unsigned int i,j;
	double maxMI;
	
	for(i=0;i<KEYS;i++){
		maxMI = 0.0;
		for(j=0;j<THREADS;j++){
			if(maxMIAll[j][i][keyByte] > maxMI){
				maxMI = maxMIAll[j][i][keyByte];
			}				
		}				
		maxMIArray[i][keyByte] = maxMI;
	}	
}

int main(){
	
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	int i,j;
	int l=0,k,m;
	int ret;
	clock_t start, end;
	time_t wallStart, wallEnd;
	
	
	int startPosition = 0;
	threadArgs * myArgs = NULL;
	maxArgs * myMaxArgs = NULL;
	
	//Assigning heap to plaintext matrix
	unsigned int **plainText=malloc(sizeof(unsigned int*)*SAMPLES);
	checkMalloc(plainText);
	for (i=0; i<SAMPLES; i++){
		plainText[i]=malloc(sizeof(unsigned int)*KEYBYTES);
		checkMalloc(plainText[i]);
	}
	
	//Reading plaintext file
	FILE * plainTextFile = fopen("../test2/plain.txt","r");
	char temp1[32];
	for(i=0; i<SAMPLES ;i++){
		fscanf(plainTextFile,"%s",temp1);
		for(j=0; j<KEYBYTES; j++){
			sscanf(&temp1[2*j],"%02X",&plainText[i][j]);		
		}
	}
	fclose(plainTextFile);	
	
	//byte keyGuess = 0x00;

	//readng and storing wave.txt
	double **waveData=malloc(sizeof(double*) * SAMPLEPOINTS);
	
	for (i=0; i<SAMPLEPOINTS; i++)
		waveData[i]=malloc(sizeof(double) * SAMPLES);

	
	FILE * waveStuff = fopen("../test2/wave.txt","r");
	float dat;
	
	
	printf("Started reading wave.txt file\n");
	
	start = clock();
	
	for(i=0; i<SAMPLES ;i++){
		for(j=0; j<SAMPLEPOINTS; j++){
			fscanf(waveStuff,"%f",&dat);
			waveData[j][i]=(double)dat*1000.0;
		}
	}
	end = clock();
	
	double micros = (double)end - (double)start;
    double millis = micros / 1000;
    
	printf("Reading time= %lf\n",millis);
	
	fclose(waveStuff);
	
	printf("Finished reading wave.txt file\n");
	
//////////////////////////////////////////////Reading stuff are done///////////////////////////////////	

	
	FILE *fp; //To print results to a file
	fp = freopen("results_new.txt", "a", stdout);
	printf("\n********************************************\n");
	printf("Started calculating MI\n");
	printf("KEYS: %d, SAMPLEPOINTS: %d, THEADS: %d SAMPLES: %d\n", KEYS, SAMPLEPOINTS, THREADS, SAMPLES);
	

	do{
		
		start = clock(); //This is CPU time
		wallStart = time(NULL); //This is wall clock time 
		
			pthread_t threads[THREADS];
			
			for (k = 0; k<THREADS; k++){
	
				myArgs = (threadArgs *)malloc(sizeof(threadArgs));
				myArgs -> waveData = waveData;
				myArgs -> plainText = plainText;
				myArgs -> startPosition = startPosition;
				myArgs -> threadNumber = k;
				ret=pthread_create(&threads[k],NULL,threaded_mia,(void*)(myArgs));
				
				if(ret == -1)
				{
					printf("Error creating thread: thread(%d)\n", i);
					break;
				}
					startPosition = startPosition + (SAMPLEPOINTS/THREADS);
			}
			
			for(i=0;i<THREADS;i++){
				int ret=pthread_join(threads[i],NULL);
				if(ret!=0){
					printf("Error joining threads\n");
					break;
				}
			}
		
		double maxMI = 0.0;
		unsigned int keyByte;
		pthread_t threads_max[KEYBYTES];
		
		for(i=0; i<KEYBYTES; i++){
			myMaxArgs = (maxArgs *)malloc(sizeof(maxArgs));
			myMaxArgs -> keyByte = i;
			ret=pthread_create(&threads_max[i],NULL,threaded_max,(void*)(myMaxArgs));
			if(ret == -1){
				printf("Error creating thread\n");
				break;
			}
			
		}
		
		for(i=0;i<KEYBYTES;i++){
			int ret=pthread_join(threads_max[i],NULL);
			if(ret!=0){
				printf("Error joining threads\n");
				break;
			}
		}
		printf("Finished calculating MI\n");
		wallEnd = time(NULL);
		end = clock();
		
		
		micros = (double)end - (double)start;
		millis = micros / 1000;
		printf("Wall clock time: %.2f s\n", (double)(wallEnd - wallStart));
		printf("CPU time= %lf ms\n",millis);
		
		
		for(i=0;i<KEYS;i++){
			for(j=0;j<KEYBYTES;j++){
				//printf("0x%02x: MI= %lf ", i, maxMIArray[i][j]);
				printf("%lf  ", maxMIArray[i][j]);
			}
			printf("\n");
		}
		fclose(fp);
	}while(0);
	
///////////////////////////////////Calculations done////////////////////////////////////////////	
	
	//clearing thread arguments
	free(myArgs);
	free(myMaxArgs);
	
	//clearing waveData matrix
	for (i=0; i<SAMPLEPOINTS; i++)
		free(waveData[i]);
	
	free(waveData);
	
	//clearing plaiText matrix
	for (i=0; i<SAMPLES; i++)
		free(plainText[i]);
		
	free(plainText);
	
	
	return 0;
}
