#ifndef PHILIB_H
#define PHILIB_H

/* -- Libraries -- */
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
/*square lattice dimensions*/

#define FALSE 0
#define TRUE 1

/* a macro used in function pythag() */

static float sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)


     /* Macro definitions for integer arguments only */

#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))



typedef struct {
 
  double edge;
  double d_edge;

} matrix_element;

typedef struct {
	int nVertices;   /* number of vertices */
	matrix_element *** edgeMatrix; /* 2D matrix array of edges */
} adj_matrix;

int *cluster;
int *cluster_coarse;

/*Monte Carlo functions*/
double *init_phi_lattice_coarse(adj_matrix* Graph,  double**R, int nold);
double *init_phi_lattice(adj_matrix* Graph);
void init_observables(adj_matrix* Graph, double *phi_vet, double *E_xbond, double *M_xspin, float lambda,int flag, double **R, int nold);
double totalEnergy(adj_matrix* Graph, double *phi_vet, double E_xbond, float lambda, int flag ,double **R, int nold);
double totalMagnetization(adj_matrix* Graph, double *phi_vet, double M_xspin, int flag ,double **R, int nold);
double calcAutocor(int sampleSize, double* M_vet_abs, double M_sum_abs, float Temp, FILE *out2);
double update_metropolis(double**R, double *phi_vet, int nold, int x, int flag);
void metropolis(adj_matrix* Graph, double *phi_vet, int x,  float Temp, float lambda, int flag ,double **R, int nold);
void sweep (adj_matrix* Graph, double *phi_vet, double *E_xbond, double *M_xspin, FILE *out, float Temp, float lambda, int flag ,double **R, int nold);
void print_vet(adj_matrix* Graph, double *phi_vet);

 //wolff cluster algorithm functions
int clusterCheck(double *phi_vet, int x, int xnext, float Temp) ;
void growClusterNeg(adj_matrix* Graph, double *phi_vet, int x,  float Temp);
void growClusterPos(adj_matrix* Graph, double *phi_vet, int x,  float Temp);
void flipCluster(adj_matrix* Graph,double *phi_vet);
void wolf(adj_matrix* Graph, double *phi_vet, int x, float Temp);

//random number generation functions
double genrand(); 
void sgenrand(unsigned long);/* initializing random generator with a NONZERO seed */

// network function declarations  
adj_matrix* CreateAdjMatrix(double, int);
adj_matrix* CreateAdjMatrixNew(double prob1,double prob2,double prob3, int nVertices);
void ViewAdjMatrix(adj_matrix* Graph);
void DestroyAdjMatrix(adj_matrix* Graph);

// find eigenvalues and eigenvectors
double ** allocate_eigvec_matrix(adj_matrix* Graph);
double * allocate_double_vector(int n);
double ** allocate_double_matrix(int row, int col);
void  eigen (double **M, double *eigenv, int n);
void tqli(double *d, double *e, int n, double **z);
void tred2(double **a, int n, double *d, double *e);
double pythag(double a, double b);

//coarse graining functions
double find_max(double *vett, int n);
double find_min(double *vett, int n);
int find_max_int(int *vett, int n);
void dswap(double *a, double *b);
void iswap(int *a, int *b);
void sort(double[],int [], int, int );

int* find_group(double *vett, int n, int num_interv);
int adjust_group(int * v_group, int n, int num_group);
int search(int *vet, int n, int num);
void merge_group(int n,int nargs, ...);
int * merge_group2(int n, int * group_vect, int * merged_vect);
double ** allocate_R_matrix(int num_group, int nVertices, int *group_vect);
double ** allocate_Coarse_matrix(int num_group, int nVertices, adj_matrix* Graph, double **R);
double ** allocate_eigvec_matrix2(double ** Mat,int n);
double * calc_phi_Coarse(double **R, double *phi_vet, int num_group, int nVertices);
int count(int * v_group, int num_group, int num);
void DestroyDoubleMatrix(double ** Mat,int nrow);
void DestroyDoubleVector(double * vet);
int * Coarse_Graining(adj_matrix* Graph, FILE * info ,int);
void grado(adj_matrix* Graph, FILE *info);
double  grado2(adj_matrix* Graph, int flag);
int *k_means(double **data, int n, int m, int k, double t, double **centroids);
// input & output functions
void write_adj(adj_matrix* Graph, int flag);
adj_matrix* read_adj(int nVertices, int flag);
void write_edge(adj_matrix* Graph, int flag);
adj_matrix* read_edge(int nVertices, int flag);
void write_pajek_adj(adj_matrix* Graph, int *group_vect);
void write_pajek_coarse(adj_matrix* Graph, int nold, double **R );

#endif




