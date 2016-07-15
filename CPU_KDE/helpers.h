/*
 * Author: Malin Prematilake
 * Date: 23.03.2016 00:28:11
 * Version: 1.0
 * Description: This code contains all the definitions related with calculating probability
*/
#ifndef HELPER_H
#define HELPER_H

#define SAMPLES 200

FILE* openFile(char* file, char *mode);
void closeFile(FILE * file);
void checkMalloc(void *ptr,int position);

#endif
