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
#include <string.h>
#include <math.h>

#include "arrayOperations.h"
#include "data.h"
#include "MI.h"
//#include "helpers.h"
#include "probOps.h"
//number of sample points in one encryption
#define SAMPLEPOINTS 100000
#define SAMPLES 200
//defined for 128 bit AES
#define KEYBYTES 16
//#define KEYS 256
#define KEYS 16
#define TEMP 194

double sum_array[2] = {0};
double sum_2_array[2] = {0};
double bandwidth[2] = {0};
double sourceMat[2][SAMPLES] = {0};

//argument struct for pthread
typedef struct args{
	double **wavedata;
	unsigned int **plainText;
	unsigned int keyguess;
	unsigned int keybyte;
}threadArgs;


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

void calc_bandwidth(){ //Default bandwidth
	
	int i, j;
	
	for(i=0; i<2; i++){
		for(j=0; j<SAMPLES; j++){
		
			sum_array[i] += sourceMat[i][j];
			sum_2_array[i] += sourceMat[i][j] * sourceMat[i][j];
		
		}
		
		double x  = sum_array[i]/SAMPLES;
		double x2 = sum_2_array[i]/SAMPLES;
		double sigma = sqrt(x2 - (x*x));
		double b = sigma*(pow((3.0*SAMPLES/4.0),(-1.0/5.0)));
		bandwidth[i] = b;
	}
	
}



int main(){
	//unsigned int plainText[SAMPLES][KEYBYTES] //move this to heap!! BW;
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	
	//byte keyGuess[256];
	//FILE * plainTexts = fopen("../test2/plain.txt","r");
	
	int i,j;
	//char temp1[32];
	
	/*for(i=0; i<SAMPLES ;i++){
		fscanf(plainTexts,"%s",temp1);
		for(j=0; j<KEYBYTES; j++){
			sscanf(&temp1[2*j],"%02X",&plainText[i][j]);		
		}
	}*/
	
	//--------------------------------------------
	
	unsigned int **plainText=malloc(sizeof(unsigned int*)*SAMPLES);
	//checkMalloc(plainText);
	for (i=0; i<SAMPLES; i++){
		plainText[i]=malloc(sizeof(unsigned int)*KEYBYTES);
		//checkMalloc(plainText[i]);
	}
	
	//Reading plaintext file
	FILE * plainTextFile = fopen("/home/e11444/CO421/test2/plain.txt","r");
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
	
	/*for (i=0;i<SAMPLES;i++){
		for(j=0;j<KEYBYTES;j++){
			hammingMat[j][i] = (double)hamming(plainText[i], j, keyGuess);	
		}
	}*/
	//calHammingMat(keyGuess, hammingMat, plainText);

	//900000/9
	//readng and storing wave.txt
	double **waveData=malloc(sizeof(double*) * SAMPLEPOINTS);
	
	for (i=0; i<SAMPLEPOINTS; i++)
		waveData[i]=malloc(sizeof(double) * SAMPLES);

	
	FILE * waveStuff = fopen("/home/e11444/CO421/test2/wave.txt","r");
	float dat;
	
	clock_t start, end;
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
	//generating key guesses
	/*for (i=0;i<256;i++)
		keyGuess[i] = (byte)i;
	*/
	
	//double miArr[KEYBYTES];
	
	printf("Started calculating MI\n");
	
	start = clock();
	
	int l=0,k,m;
	//int keyGuesses = 8;
	double miArr[KEYS][KEYBYTES] = {0};
	
	//double **sourceMat=malloc(sizeof(double*) * 2);
	
	//for (i=0; i<2; i++)
		//sourceMat[i]=malloc(sizeof(double) * SAMPLES);
	
	double * joint = malloc(sizeof(double) * SAMPLES * SAMPLES);
	double * marginalWave = malloc(sizeof(double) * SAMPLES);
	//************ Temp test **************
	/*FILE * test1 = fopen("set1.txt", "r");
	FILE * test2 = fopen("set2.txt", "r");
	for(i=0; i<SAMPLES; i++){
	
		fscanf(test1,"%f",&dat);
		sourceMat[0][i] = (double)dat;
		fscanf(test2,"%f",&dat);
		sourceMat[1][i] = (double)dat;
	}
	fclose(test1);
	fclose(test2);*/
	//*************************************
	
	for(m=0; m<1; m++){				
		calHammingMat(m, hammingMat, plainText); //Calculate hamming weights for that key guess		
		for (l=0;l<1;l++){
			double maxMI = 0.0;
			int col;
			
			//printf("This is the %d hamming array\n",l);
			for(k=0;k<1;k++){
				
				memcpy(sourceMat[0], waveData[k], sizeof(double) * SAMPLES);
				memcpy(sourceMat[1], hammingMat[l], sizeof(double) * SAMPLES);
				
				calc_bandwidth();
				
				//for(i=0; i<SAMPLES; i++)
					//printf("%lf\n", hammingMat[l][i]);
				
				//printf("b1: %lf b2: %lf\n", bandwidth[0], bandwidth[1]);
				gauss_kde_joint(sourceMat, waveData[k], bandwidth[0], 200, hammingMat[l], bandwidth[1], 200, joint);
				gauss_kde_marginal(sourceMat[0], bandwidth[0], 200, marginalWave);
				
				//double mi = MI(waveData[k], hammingMat[l], SAMPLES);
				
				
				
				//printf("MI of %d = %lf\n",k,mi);
				
				/*if (mi>maxMI){
					maxMI = mi;
					col = k;
				}*/
			}
		//	miArr[m][l] = maxMI;
		}
	}
	
	//maxMutualInfo(waveData, plainText, miArr);
	
	printf("Finished calculating MI\n");
	
	end = clock();
	
	//clearing waveData matrix
	for (i=0; i<SAMPLEPOINTS; i++)
		free(waveData[i]);
	
	free(waveData);
	
	//clearing plaiText matrix
	for (i=0; i<SAMPLES; i++)
		free(plainText[i]);
		
	free(plainText);
	
	/*for(i=0; i<KEYS; i++){
		for(j=0;j<KEYBYTES;j++){
			printf("0x%02X: MI= %lf ", i, miArr[i][j]);
		}
		printf("\n");
	}*/
	
	FILE * f = freopen("test.txt", "w", stdout);
	for(i=0; i<SAMPLES*SAMPLES; i++){
		//for(j=0; j<SAMPLES; j++){
			printf("%lf ", joint[i]);
		//}
		//printf("\n");
	}
	
	printf("\n");
	for(i=0; i<SAMPLES; i++){
			printf("%lf ", marginalWave[i]);
	}



	//printf("Resulting MI:\n");
	//printArrayDouble(miArr,KEYBYTES);
	
	micros = (double)end - (double)start;
    millis = micros / 1000;
    
	printf("Testing time= %lf\n",millis);
	fclose(f);
	//printArrayDouble(miArr,100);
	return 0;
}
