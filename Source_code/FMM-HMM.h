#ifndef GIBBSSAMPLING_H
#define GIBBSSAMPLING_H
#include <stdint.h>
#define bool int
#define true 1
#define false 0

#define MAX_OUTPUT_LINE 15
//#define NEW_CLUSTER 1000000

// ***** Distance functions *****

void f_add(double *x, double *y, long *n, long *m, double *z);
double emission(double *observation, double *lamda, long nstate, long nmarker, long ntime, long nregion, long ncluster, long k, long w, long h, long t);
void forward(double *f, double *b, double *lambda, double *observation, long nstate, long nmarker, long ntime, long nregion, long ncluster);
void backward(double *back, double *b, double *lambda, double *observation, long nstate, long nmarker, long ntime, long nregion, long ncluster);
double log_likelihood(long nregion, long ncluster, long ntime, long *cluster, double *pie, double *forwar);
double f_max(double *a, long a_length, long *index);
long f_which(double *a, long a_length, double b);
double *rowSums(double **probability, long nregion, long ncluster);
void matrixtovector(double **matrix, double *vector, long nregion, long ncluster);
void colMeans(double **probability, double *pie, long nregion, long ncluster);
void HMM_Recursion(double *b, double *observation, double *lambda, long *cluster, double *forwar, double *backwar, long nstate, long ncluster, 
                   long nmarker, long ntime, long nregion,long *H);
double Edistance(double *a, double *b, long length);
double mEdistance(double **a, double **b, long nrow, long ncol);
void copy_vector(double *a, double *b, long length);
void copy_matrix(double **a, double **b, long nrow, long ncol);
void cluster_matrix(long nregion, long ncluster, long ntime, long *cluster, double *pie, double *forwar, double **probability);

void f_HMM(double *observation, double *b, double *lambda, int32_t *cluster, double *pie, double *forwar, double *backwar, 
double *probability, double *diff_r, double *log_likeli_r, double *pie_r, int32_t *nstate, int32_t *ncluster, int32_t *nmarker, int32_t *ntime, int32_t *nregion, 
int32_t *maxiteration, int32_t *stop, int32_t *nstep, double *ndistance, int32_t *H);

// Test functions



#endif









