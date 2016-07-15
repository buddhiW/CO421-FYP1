/*
 * Author: Malin Prematilake
 * Date: 23.03.2016 00:28:11
 * Version: 1.0
 * Description: This code contains all the definitions related with calculating probability
*/
#ifndef PROBABILITYOPS_H
#define PROBABILITYOPS_H

#define SAMPLES	200

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

#endif
