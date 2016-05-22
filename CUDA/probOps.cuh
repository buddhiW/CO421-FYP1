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

#endif
