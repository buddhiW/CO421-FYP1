/* 

	Header file for CUDA helper functions
	Include this file in your code if you want to use the function checkCudaError() 
	checkCudaError() function checks if the last cuda function call or kernel launch caused an error
	and if yes it will print the error message and will abort the program

	For example to check any errors due to memcpy use it as
	cudaMemcpy(arguments...);
	checkCudaError();

	For example to check any errors due to a kernel launch use it as		
	kernel<<<parameters>>>(arguments...);
	checkCudaError();
	
	Hasindu Gamaarachchi <hasindu2008@live.com>
	22 Dec 2015

*/

#define SAMPLEPOINTS 100000

#define SAMPLES 200

//defined for 128 bit AES
#define KEYBYTES 16
#define KEYS 256

#define NUMOFKEYS 1 //Number of keys to be calculated
//#define TEST 6250
#define REPEAT 16
#define WIDTH 20 //Range of one bin
//#define BINS 80
#define checkCudaError() { gpuAssert(__FILE__, __LINE__); }

/* check whether the last CUDA function or CUDA kernel launch is erroneous and if yes an error message will be printed
and then the program will be aborted*/
void gpuAssert(const char *file, int line);

/* Check whether a previous memory allocation was successful. If RAM is full usually the returned value is a NULL pointer.
For example if you allocate memory by doing 
int *mem = malloc(sizeof(int)*SIZE)
check whether it was successful by calling
isMemoryFull(mem) afterwards */
void checkAllocRAM(void *ptr);

/* This checks whether a file has been opened corrected. If a file opening failed the returned value is a NULL pointer
FOr example if you open a file using
FILE *file=fopen("file.txt","r");
check by calling isFileValid(file); */
void isFileValid(FILE *fp);
