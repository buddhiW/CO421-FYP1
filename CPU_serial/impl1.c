#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "arrayOperations.h"
#include "data.h"
#include "MI.h"

//number of sample points in one encryption
#define SAMPLEPOINTS 1000000
#define SAMPLES 200
//defined for 128 bit AES
#define KEYBYTES 16
#define KEYS 256

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

int main(){
	unsigned int plainText[SAMPLES][KEYBYTES];
	double hammingMat[KEYBYTES][SAMPLES];//16x200
	
	//byte keyGuess[256];
	FILE * plainTexts = fopen("../test2/plain.txt","r");
	
	int i,j;
	char temp1[32];
	
	for(i=0; i<SAMPLES ;i++){
		fscanf(plainTexts,"%s",temp1);
		for(j=0; j<KEYBYTES; j++){
			sscanf(&temp1[2*j],"%02X",&plainText[i][j]);		
		}
	}
	
	fclose(plainTexts);
	
	byte keyGuess = 0x00;
	
	for (i=0;i<SAMPLES;i++){
		for(j=0;j<KEYBYTES;j++){
			hammingMat[j][i] = (double)hamming(plainText[i], j, keyGuess);	
		}
	}

	//900000/9
	//readng and storing wave.txt
	double **waveData=malloc(sizeof(double*) * SAMPLEPOINTS);
	
	for (i=0; i<SAMPLEPOINTS; i++)
		waveData[i]=malloc(sizeof(double) * SAMPLES);

	
	FILE * waveStuff = fopen("../test2/wave.txt","r");
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
	
	int l=0,k;
	double miArr[KEYBYTES] = {0};
	for (l=0;l<KEYBYTES;l++){
		double maxMI = 0.0;
		int col;
		
		//printf("This is the %d hamming array\n",l);
		for(k=0;k<SAMPLEPOINTS;k++){
			double mi = MI(waveData[k], hammingMat[l], SAMPLES);
			
			//printf("MI of %d = %lf\n",k,mi);
			
			if (mi>maxMI){
				maxMI = mi;
				col = k;
			}
		}
		miArr[l] = maxMI;
	}
	printf("Finished calculating MI\n");
	
	end = clock();
	
	//clearing waveData matrix
	for (i=0; i<SAMPLEPOINTS; i++)
		free(waveData[i]);
	
	free(waveData);
	
	for(i=0;i<KEYBYTES;i++)
		printf("MI= %lf and keybyte= %d\n",miArr[i], i);
	//printf("Resulting MI:\n");
	//printArrayDouble(miArr,KEYBYTES);
	
	micros = (double)end - (double)start;
    millis = micros / 1000;
    
	printf("Testing time= %lf\n",millis);
	//printArrayDouble(miArr,100);
	return 0;
}
