/*
 * Author: Malin Prematilake
 * Date: 23.03.2016 00:28:11
 * Version: 1.0
 * Description: This code contains all the definitions related with calculating probability
*/
#ifndef PROBABILITYOPS_H
#define PROBABILITYOPS_H

typedef struct jpState
{
  double *jointProbabilityVector;
  int numJointStates;
  double *firstProbabilityVector;
  int numFirstStates;
  double *secondProbabilityVector;
  int numSecondStates;
} JointProbabilityState;

typedef struct pState
{
  double *probabilityVector;
  int numStates;
} ProbabilityState;

JointProbabilityState* jointProbability(double *array1, double *array2, int length);
ProbabilityState probability(double *array, int Length);
__global__ void findProbsWave(unsigned int *dev_normaliseArray1, double * dev_firstStateProbs, unsigned int *dev_firstNumState, int length);
__global__ void findProbsHamming(unsigned int *dev_normaliseArray2, double * dev_secondStateProbs, unsigned int *dev_secondNumState, int length);
__global__ void findJointProbs(double *MIvals, unsigned int *normaliseWave, unsigned int *normaliseHamming, unsigned int *firstNumStates, 
								unsigned int *secondNumStates, double *waveProbs, double *hammingProbs, double *jointProbs, 
								int maxOfFirst, int maxOfSecond, int key, int startPosition);
								
__global__ void KDE_findProbsHamming(double *dev_hammingData, double * dev_hammingDataProbs,
										double *bandWidthArray, int targetLength, double * maxMin);
__global__ void KDE_findProbsWave(double *dev_waveData, double * dev_waveDataProbs, 
									double *bandWidthArray, int targetLength, double * maxMin);
__global__ void KDE_findJointProbs(double *dev_hammingData, double *dev_waveData, unsigned int targetLengthW, 
									unsigned int targetLengthH, double *jointProbs, double * bandwidthW, double * bandwidthH, 
									double * maxMinH, double * maxMinW, int turn);
__global__ void KDE_findJointProbs2(double *dev_hammingData, double *dev_waveData, unsigned int targetLengthW, 
									unsigned int targetLengthH, double * bandwidthW, double * bandwidthH, 
									double *dev_hammingTarget, double *dev_waveTarget, double *jointProbs,
									int turn, int key, int keyByte);

#endif
