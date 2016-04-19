//-------------------------------------------------------------------------
//	MIA based power analysis attack on the CPU
//
//	Version information
//
//	Author			Details
//	MP				Initial version
//	BW				Modularised some code into functions
//					Moved plaiText matrix to heap		
//-------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include "arrayOperations.h"
#include "data.h"
#include "MI.h"

//number of sample points in one encryption
#define SAMPLEPOINTS 100000
#define SAMPLES 200
//defined for 128 bit AES
#define KEYBYTES 16
//#define KEYS 256
#define KEYS 256

//#define THREADS 8

//argument struct for pthread
typedef struct args{
	double **waveData;
	unsigned int **plainText;
	unsigned int keyByte;
}threadArgs;

double maxMIArray[KEYS][KEYBYTES] = {0};

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
//This will have to be modified if your selection function/intermediate values are different

//find hamming weight for the selection function
byte hamming(unsigned int sample[], unsigned int n,unsigned int key) { //n is byteno sample is the line in sample text
    byte inter = (byte)sbox[sample[n] ^ key];
    byte dist = hammingweight(inter);		
	return dist;
}

void calHammingMat(byte keyGuess, double hammingMat[KEYBYTES][SAMPLES], unsigned int **plainText){
	//double hammingMat[KEYBYTES][SAMPLES];
	int i,j;
	for (i=0;i<SAMPLES;i++){
		for(j=0;j<KEYBYTES;j++){
			hammingMat[j][i] = (double)hamming(plainText[i], j, keyGuess);	
		}
	}
	
}

//------------------------------------------------------------------------
// maxMutualInfo
// maxMIArray: 16 * 256 array
// 			   nth column gives maximum MI values for each key guess for 
//				nth key byte
//------------------------------------------------------------------------
void maxMutualInfo(double ** waveData, unsigned int **plainText, double maxMIArray[KEYS][KEYBYTES]){
	int l,k,m;
	double maxMI = 0.0;
	int col;	
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	//int keyGuesses = 8;
	//double miArr[KEYS][KEYBYTES] = {0};
	
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

void maxMiaThreaded(double ** waveData, unsigned int **plainText, unsigned int keyByte){
	int l,k,m;
	double maxMI = 0.0;
	int col;	
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	//int keyGuesses = 8;
	//double miArr[KEYS][KEYBYTES] = {0};
	
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

void * threaded_mia(void * args){
	
	threadArgs * myargs = (threadArgs *) args;
	double ** waveData = myargs -> waveData;
	unsigned int **plainText = myargs -> plainText;
	unsigned int keyByte = myargs -> keyByte;
	unsigned int m, k;
	double maxMI = 0.0;
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	
	for(k=0; k<KEYS; k++){
		maxMI = 0.0;				
		calHammingMat(k, hammingMat, plainText); //Calculate hamming weights for that key guess		
			
			for(m=0;m<SAMPLEPOINTS;m++){
				double mi = MI(waveData[m], hammingMat[keyByte], SAMPLES);
				
				if (mi>maxMI){
					maxMI = mi;
					
				}
			}
			
			maxMIArray[k][keyByte] = maxMI;
	}
	
}


int main(){
	//unsigned int plainText[SAMPLES][KEYBYTES] //move this to heap!! BW;
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	
	//byte keyGuess[256];
	//FILE * plainTexts = fopen("../test2/plain.txt","r");
	
	int i,j;
	int ret;
	
	pthread_t threads[KEYBYTES];
	//char temp1[32];
	
	/*for(i=0; i<SAMPLES ;i++){
		fscanf(plainTexts,"%s",temp1);
		for(j=0; j<KEYBYTES; j++){
			sscanf(&temp1[2*j],"%02X",&plainText[i][j]);		
		}
	}*/
	
	//--------------------------------------------
	
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
	//------------------------------------------------
	
	byte keyGuess = 0x00;
	
	//900000/9
	//readng and storing wave.txt
	double **waveData=malloc(sizeof(double*) * SAMPLEPOINTS);
	
	for (i=0; i<SAMPLEPOINTS; i++)
		waveData[i]=malloc(sizeof(double) * SAMPLES);

	
	FILE * waveStuff = fopen("../test2/wave.txt","r");
	float dat;
	
	clock_t start, end;
	time_t wallStart, wallEnd;
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

	printf("Started calculating MI\n");
	
	printf("KEYS: %d, SAMPLEPOINTS: %d\n", KEYS, SAMPLEPOINTS);
	
	int l=0,k,m;
	//int keyGuesses = 8;
	
	//maxMutualInfo(waveData, plainText, miArr);
	threadArgs * myArgs = NULL;
	do{
		start = clock(); //??
		wallStart = time(NULL);
		for (i=0;i<KEYBYTES;i++){
			myArgs = (threadArgs *)malloc(sizeof(threadArgs));
			myArgs -> waveData = waveData;
			myArgs -> plainText = plainText;
			myArgs -> keyByte = i;
			ret=pthread_create(&threads[i],NULL,threaded_mia,(void*)(myArgs));
			//clearing data structure
			
	
			if(ret == -1)
			{
				printf("Error creating thread: thread(%d)\n", i);
				break;
			}			
		}
		
		for(i=0;i<KEYBYTES;i++){
			int ret=pthread_join(threads[i],NULL);
			if(ret!=0){
				printf("Error joining threads\n");
				break;
			}
		}
		free(myArgs);
		
		printf("Finished calculating MI\n");
		wallEnd = time(NULL);
		end = clock();
		
		for(i=0; i<KEYS; i++){
			for(j=0;j<KEYBYTES;j++){
				printf("0x%02X: MI= %lf ", i, maxMIArray[i][j]);
			}
			printf("\n");
		}
	
		micros = (double)end - (double)start;
		millis = micros / 1000;
		printf("Wall clock time: %.2f s\n", (double)(wallEnd - wallStart));
		printf("Testing time= %lf\n",millis);
		
	}while(0);
	
///////////////////////////////////Calculations done////////////////////////////////////////////	
	
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
