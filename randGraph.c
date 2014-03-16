#include "phi_lib.h"


/* -- Internal function defintions -- */


/* Function : CreateAdjMatrix

 Return value   adj_matrix* - pointer to adjacency matrix 
*/
adj_matrix* CreateAdjMatrix(double prob, int nVertices)
{
	int i, j, k;
	   
	adj_matrix* Graph;
        double c;
	srand(time(NULL));
	/* allocate the base structure */
	Graph = (adj_matrix*)malloc(sizeof(adj_matrix));
	if (Graph == NULL)
		return NULL;
	Graph->edgeMatrix = NULL;
	Graph->nVertices = nVertices;
        
	/* read the number of vertices 
        printf("\ninserisci il numero di vertici: ");
	scanf("%d", &nVertices);
        printf("\n");
	Graph->nVertices = nVertices;*/


	/* allocate the base array of the 2D adjacency matrix */
	Graph->edgeMatrix = (matrix_element***)malloc(sizeof(matrix_element**) * nVertices);
	if (Graph->edgeMatrix == NULL) {
		DestroyAdjMatrix(Graph);
		return NULL;
	}


	/* this is to initialize the array  */
	for (i = 0; i < nVertices; i++)
		Graph->edgeMatrix[i] = NULL;

	
	/* allocate each row of the 2D adjacency matrix */
	for (i = 0; i < nVertices; ++i) {

		Graph->edgeMatrix[i] = (matrix_element**)malloc(sizeof(matrix_element*) * nVertices);
		
		if (Graph->edgeMatrix[i] == NULL) {
			DestroyAdjMatrix(Graph);
			return NULL;
		}
	}
        
       /*printf("\n----sono prima di read matrix-----\n");*/
	/* allocate  each data structure in the matrix */
	for (i = 0; i < nVertices; i++){
		for (j = 0; j < nVertices; j++){
	               Graph->edgeMatrix[i][j] = (matrix_element*)malloc(sizeof(matrix_element));
                       
               }
                
             }
             
       /* initialize each matrix diagonal data with 0  */ 
       
       for (i = 0; i < nVertices; i++){
            Graph->edgeMatrix[i][i]->edge=0;
       }
            
       
       /* initialize each matrix non-diagonal data with adjacency  */
           k=0;
           for (i = 0; i < nVertices; i++){
		
                for (j = 0; j < k; j++){
	               
                       c=rand()/(double)RAND_MAX;
                       

                       if (c<prob){
                       Graph->edgeMatrix[i][j]->edge=1;
                       Graph->edgeMatrix[j][i]->edge=1;
                       
                        }

                       else{
                        Graph->edgeMatrix[i][j]->edge=0;
                        Graph->edgeMatrix[j][i]->edge=0;
                       }       
                  //printf("%d ", Graph->edgeMatrix[i][j]->edge);
               }
            k++;
           // printf("\n");    
             }
 printf("\n");    
	
	
	return Graph;
        
}



adj_matrix* CreateAdjMatrixNew(double prob1,double prob2, double prob3, int nVertices)
{
	int i, j, k;
	  
	adj_matrix* Graph;
        double c;
	srand(time(NULL));
	/* allocate the base structure */
	Graph = (adj_matrix*)malloc(sizeof(adj_matrix));
	if (Graph == NULL)
		return NULL;
	Graph->edgeMatrix = NULL;
	Graph->nVertices = nVertices;

        if (nVertices%8!=0){
              printf("The number of vertices isn t multiple of 8 \n");
              exit(1);
          }


	/* allocate the base array of the 2D adjacency matrix */
	Graph->edgeMatrix = (matrix_element***)malloc(sizeof(matrix_element**) * nVertices);
	if (Graph->edgeMatrix == NULL) {
		DestroyAdjMatrix(Graph);
		return NULL;
	}


	/* this is to initialize the array  */
	for (i = 0; i < nVertices; i++)
		Graph->edgeMatrix[i] = NULL;

	
	/* allocate each row of the 2D adjacency matrix */
	for (i = 0; i < nVertices; ++i) {

		Graph->edgeMatrix[i] = (matrix_element**)malloc(sizeof(matrix_element*) * nVertices);
		
		if (Graph->edgeMatrix[i] == NULL) {
			DestroyAdjMatrix(Graph);
			return NULL;
		}
	}
        
       /*printf("\n----sono prima di read matrix-----\n");*/
	/* allocate  each data structure in the matrix */
	for (i = 0; i < nVertices; i++){
		for (j = 0; j < nVertices; j++){
	               Graph->edgeMatrix[i][j] = (matrix_element*)malloc(sizeof(matrix_element));
                       
               }
                
             }
             
       /* initialize each matrix diagonal data with 0  */ 
       
       for (i = 0; i < nVertices; i++){
            Graph->edgeMatrix[i][i]->edge=0;
       }
            
       
       /* initialize each matrix non-diagonal data with adjacency  */
           
           for (i = 1; i < nVertices; i++){
		
                for (j =i ; j < nVertices; j++){
	               
                       c=rand()/(double)RAND_MAX;
                       

                       if (c<prob1){
                       Graph->edgeMatrix[i][j]->edge=1;
                       Graph->edgeMatrix[j][i]->edge=1;
                       
                        }

                       else{
                        Graph->edgeMatrix[i][j]->edge=0;
                        Graph->edgeMatrix[j][i]->edge=0;
                       }       
                  //printf("%d ", Graph->edgeMatrix[i][j]->edge);
               }
           
           // printf("\n");    
             }
         int z=0;
         
         if (prob2!=0){
             for (z =0; z <nVertices; z=z+nVertices/4){
                   // printf("z vale %d \n",z);           
              for (i =z; i < z+nVertices/4; i++){
		       // printf("i vale %d",i);
                for (j =i ; j < z+nVertices/4; j++){
	               
                       
                       Graph->edgeMatrix[i][j]->edge=0;
                       Graph->edgeMatrix[j][i]->edge=0;
                       c=rand()/(double)RAND_MAX;
                       if (c<prob2){
                       Graph->edgeMatrix[i][j]->edge=1;
                       Graph->edgeMatrix[j][i]->edge=1;
                       
                        }

                       else{
                        Graph->edgeMatrix[i][j]->edge=0;
                        Graph->edgeMatrix[j][i]->edge=0;
                       }       
                  //printf("%d ", Graph->edgeMatrix[i][j]->edge);
               }
           
            printf("\n");    
             }
           }

         //elimino i loop
        for (i = 0; i < nVertices; i++){
            Graph->edgeMatrix[i][i]->edge=0;
         }

          }//end if

         if (prob3!=0){
             for (z =0; z <nVertices; z=z+nVertices/8){
                    //printf("z vale %d \n",z);           
              for (i =z; i < z+nVertices/8; i++){
		       //printf("i vale %d",i);
                for (j =i ; j < z+nVertices/8; j++){
	               
                       
                       Graph->edgeMatrix[i][j]->edge=0;
                       Graph->edgeMatrix[j][i]->edge=0;
                       c=rand()/(double)RAND_MAX;
                       if (c<prob3){
                       Graph->edgeMatrix[i][j]->edge=1;
                       Graph->edgeMatrix[j][i]->edge=1;
                       
                        }

                       else{
                        Graph->edgeMatrix[i][j]->edge=0;
                        Graph->edgeMatrix[j][i]->edge=0;
                       }       
                  //printf("%d ", Graph->edgeMatrix[i][j]->edge);
               }
           
            printf("\n");    
             }
           }        

          }//end if

 printf("\n");    
	
	 //elimino i loop
        for (i = 0; i < nVertices; i++){
            Graph->edgeMatrix[i][i]->edge=0;
         }

	return Graph;
        
}

/* Function : ViewAdjMatrix */

void ViewAdjMatrix(adj_matrix* Graph)
{
	int i, j;
            
            printf("\n");    
	
	for (i = 0; i < Graph->nVertices; i++) {

		for (j = 0; j < Graph->nVertices; j++) {
			printf("%8.1lf ", Graph->edgeMatrix[i][j]->edge);
		}

		printf("\n");
	}
}


/* Function : DestroyAdjMatrix
*/
void
DestroyAdjMatrix(adj_matrix* Graph)
{
	int i,j;

	if (Graph != NULL) {

		/* destroy the matrix array and elements */
		if (Graph->edgeMatrix != NULL) {

			for (i = 0; i < Graph->nVertices; i++){
                             for (j = 0; j < Graph->nVertices; j++){   
				free(Graph->edgeMatrix[i][j]);
                           
                             }
                             free(Graph->edgeMatrix[i]);
                           }
			
			free(Graph->edgeMatrix);
		}

		/* destroy the graph structure */
		free(Graph);
	}
}

double ** allocate_eigvec_matrix(adj_matrix* Graph)
{


         int i,j,n;
         double **a;
         n= Graph->nVertices;

	//for (i = 0; i < n; ++i) {

		a = (double**)malloc(sizeof(double*) * n);
				
	//}
        
       	/* allocate  each data structure in the matrix */
	for (i = 0; i < n; i++){
		
	               a[i] = (double*)malloc(sizeof(double)* n);
                                                   
             }

         for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
	               a[i][j] = Graph->edgeMatrix[i][j]->edge;
                                             
               }
                printf("\n");
             }

    return a;
}
 

double ** allocate_eigvec_matrix2(double ** Mat,int n)
{


         int i,j;
         double **a;
         
	///for (i = 0; i < n; ++i) {

		a = (double**)malloc(sizeof(double*) * n);
				
	//}
        
       	/* allocate  each data structure in the matrix */
	for (i = 0; i < n; i++){
		
	               a[i] = (double*)malloc(sizeof(double)* n);
                                                   
             }

         for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
	               a[i][j] = Mat[i][j];
                                             
               }
                printf("\n");
             }

    return a;
}

double ** allocate_double_matrix(int row, int col)
{


         int i,j,n;
         double **a;
         

	//for (i = 0; i < row; ++i) {

		a = (double**)malloc(sizeof(double*) * row);		
	//}
        
       
	/* allocate  each data structure in the matrix */
	for (i = 0; i < row; i++){
		
	               a[i] = (double*)malloc(sizeof(double)* col);
                                      
             }


    return a;
}
double * allocate_double_vector(int n){
        
        double *v;

        v = (double*)malloc(sizeof(double) * n);
        return v;


}


/* Computes eigenvalues (and eigen vectors if desired) for	*
*  symmetric matices. 						*/
void  eigen ( double **M, double *eigenv, int n	)
{
    int   i,j;
    double *e, **eigvec;
    
    e = (double*)malloc(sizeof(double) * n);

    tred2(M,n,eigenv,e);
   
/* allocate each row of the 2D adjacency matrix */
	for (i = 0; i < n; ++i) {

		eigvec = (double**)malloc(sizeof(double*) * n);
				
	}
               
	/* allocate  each data structure in the matrix */
	for (i = 0; i < n; i++){
		
	      eigvec[i] = (double*)malloc(sizeof(double)* n);
                                   
             }

         for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
	               if(i==j)
                            eigvec[i][j]=1;
                       else eigvec[i][j]=0;               
               }
                
             }

    tqli(eigenv,e,n,M);

      for (j=0; j<n; j++) {
           printf("eigenvalue %8.3lf ->", eigenv[j]);
           
          for (i=0; i<n; i++) {
                  printf("%8.3lf, ", M[i][j]);
         }
            printf("\n\n"); 
       }

      // DestroyDoubleMatrix(eigvec, n);
        //DestroyDoubleVector(e);   
}


/*
    ** The function
    **                 tqli()
    ** determine eigenvalues and eigenvectors of a real symmetric
    ** tri-diagonal matrix, or a real, symmetric matrix previously
    ** reduced by function tred2[] to tri-diagonal form. On input,
    ** d[] contains the diagonal element and e[] the sub-diagonal
    ** of the tri-diagonal matrix. On output d[] contains the
    ** eigenvalues and  e[] is destroyed. If eigenvectors are
    ** desired z[][] on input contains the identity matrix. If
    ** eigenvectors of a matrix reduced by tred2() are required,
    ** then z[][] on input is the matrix output from tred2().
    ** On output, the k'th column returns the normalized eigenvector
    ** corresponding to d[k]. 
    ** The function is modified from the version in Numerical recipe.
    */

void tqli(double *d, double *e, int n, double **z)
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;

   for(i = 1; i < n; i++) e[i-1] = e[i];
   e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         if(m != l) {
            if(iter++ == 30) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag(g,1.0);
            g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop */
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */
   
    
 /*
    ** The function
    **                tred2()
    ** perform a Housholder reduction of a real symmetric matrix
    ** a[][]. On output a[][] is replaced by the orthogonal matrix 
    ** effecting the transformation. d[] returns the diagonal elements
    ** of the tri-diagonal matrix, and e[] the off-diagonal elements, 
    ** with e[0] = 0.
    ** The function is modified from the version in Numerical recipe.
    */

void tred2(double **a, int n, double *d, double *e)
{
   register int    l,k,j,i;
   double          scale,hh,h,g,f;

   for(i = n - 1; i > 0; i--) {
      l = i-1;
      h = scale= 0.0;
      if(l > 0) {
         for(k = 0; k <= l; k++)
            scale += fabs(a[i][k]);
            if(scale == 0.0)               // skip transformation
               e[i] = a[i][l];
            else {
            for(k = 0; k <= l; k++) {
               a[i][k] /= scale;          // used scaled a's for transformation
               h       += a[i][k]*a[i][k];
            }
            f       = a[i][l];
            g       = (f >= 0.0 ? -sqrt(h) : sqrt(h));
            e[i]    = scale*g;
            h      -= f * g;
            a[i][l] = f - g;
            f       = 0.0;

            for(j = 0;j <= l;j++) {
               a[j][i] = a[i][j]/h;       // can be omitted if eigenvector not wanted
               g       = 0.0; 
               for(k = 0; k <= j; k++) {
                  g += a[j][k]*a[i][k];
               }
               for(k = j+1; k <= l; k++)
                  g += a[k][j]*a[i][k];
               e[j]=g/h;
               f += e[j]*a[i][j];
            }
            hh=f/(h+h);
            for(j = 0; j <= l;j++) {
               f = a[i][j];
               e[j]=g=e[j]-hh*f;
               for(k = 0; k <= j; k++)
                  a[j][k] -= (f*e[k]+g*a[i][k]);
            }
         }  // end k-loop
      }  // end if-loop for l > 1
      else {
         e[i]=a[i][l];
      }
      d[i]=h;
   }  // end i-loop
   d[0]  = 0.0;
   e[0]  = 0.0;

         /* Contents of this loop can be omitted if eigenvectors not
	 ** wanted except for statement d[i]=a[i][i];
         */

   for(i = 0; i < n; i++) {
      l = i-1;
      if(d[i]) {
         for(j = 0; j <= l; j++) {
            g= 0.0;
            for(k = 0; k <= l; k++) {
               g += a[i][k] * a[k][j];
            }
            for (k = 0; k <= l; k++) {
               a[k][j] -= g * a[k][i];
            }
         }
      }
      d[i]    = a[i][i];
      a[i][i] = 1.0;
      for(j = 0; j <= l; j++)  {
         a[j][i]=a[i][j] = 0.0;
      }
   }
} // End: function tred2(), (C) Copr. 1986-92 Numerical Recipes Software )



double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
// End: function pythag(), (C) Copr. 1986-92 Numerical Recipes Software )%.



double find_max(double *vett, int n){
 int i;
 double max=-DBL_MAX;
             for(i=0; i<n; i++){
                   if (vett[i]>max)
                       max=vett[i];
              }

return max;
}

double find_min(double *vett, int n){
 int i;
 double min=DBL_MAX;

             for(i=0; i<n; i++){
                   if (vett[i]<min)
                       min=vett[i];
              }

return min;
}

////////////sort functions/////////////
void dswap(double *a, double *b)
{
  double t=*a; *a=*b; *b=t;
}

void iswap(int *a, int *b)
{
  int t=*a; *a=*b; *b=t;
}

void sort (double a[],int posit [], int lo, int hi)
{
    //  lo is the lower index, hi is the upper index
    //  of the region of array a that is to be sorted
    int i=lo, j=hi, h;
    double x=a[(lo+hi)/2];

    //  partition
    do
    {
        while (a[i]<x) i++; 
        while (a[j]>x) j--;
        if (i<=j)
        {
            //h=a[i]; a[i]=a[j]; a[j]=h;
            dswap(&a[i], &a[j]);
            iswap ( &posit[i], &posit[j]);
             i++; j--;
        }
    } while (i<=j);

    //  recursion
    if (lo<j) sort(a, posit,lo, j);
    if (i<hi) sort(a,posit, i, hi);
}




 int * find_group(double *vett, int n, int num_interv){

     int *v_group;
     
     int  j,i, numg=0;
     double max, min, g, l;
     v_group= (int*)malloc(sizeof(int) * n);


      min=find_min(vett, n);
      max=find_max(vett, n);
      printf("\nMAX = %8.3lf, MIN = %8.3lf  \n", max, min);
      l= (max-min)/num_interv;  
      printf("num intervals = %d, l = %8.3lf  \n", num_interv, l);
      
           g=max;
      //for(g=max; g>min; g=g-l)
          for(j=0; j<num_interv; j++){
             printf("\ng= %8.3lf  \n", g);    

            for(i=0; i<n; i++)
                if (vett[i]<= g) 
                   v_group[i]=numg;  

            numg++;
            g=g-l;
            
      }
          for(i=0; i<n; i++){
                   printf("%d, ", v_group[i]);
      }
      printf("\n");    
      
      return v_group;
}



int adjust_group(int * v_group, int n, int num_group){
          
          int i, j, flag,z;
          int max=-INT_MAX;
         
    
           for(i=0; i<num_group; i++){
                
               flag=search(v_group,n,i);
               
               while(flag==FALSE ){
                 
                 for(j=0; j<n; j++)
                        if (v_group[j]>i)
                             v_group[j]=v_group[j]-1;
                 
               
                 flag=search(v_group,n,i);
                 
                 }// while
               for(z=0; z<n; z++){
                   if (v_group[z]>max)
                     max=v_group[z]; //find max value of groups
                    }
                  
               if (i==max) return max+1;
               else max=-INT_MAX;
            }        
          
       
         
        return max+1;


}

int search(int *vet, int n, int num){
     int i; 

      for(i=0; i<n; i++)
           if (vet[i]==num)
                 return TRUE;
      
        return FALSE;   

}


void merge_group(int n, int narg, ...){
      
      int *vett;    
      char a[50];
      char buffer[50];
      char flag_vet[n][50]; 
      int i, j;
    va_list args;
	  long long int ret_vett[n];
  
 
   for(i=0; i<n; i++){
        
          va_start(args,narg); 
      for(j=0; j<narg; j++){
           
       vett=va_arg(args, int*); 
      
       sprintf( a, "%d", vett[i] ); 
      
       if (j==0)
          strcpy(buffer, a);
       else
           strcat(buffer, a);
       
      }
       va_end(args);
       strcpy(flag_vet[i], buffer);
       puts(flag_vet[i]);
       ret_vett[i]=atoi(flag_vet[i]);
    }
   //for(i=0; i<n; i++){printf("%d ", ret_vett[i]);}
   
   int count =0, temp;
   for(i=0; i<n; i++){
        if (ret_vett[i]>count-1){   
                temp=ret_vett[i];
                ret_vett[i]=count;
           for(j=i+1; j<n; j++){
               if(ret_vett[j]==temp )
                     ret_vett[j]=count;
            }
         }
      count++;
   }
 printf("\n--------nuovo ret vet--------\n ");
    for(i=0; i<n; i++){printf("%d ", ret_vett[i]);}

}

int * merge_group2(int n, int * group_vect, int * merged_vect){
         
      char a[20];
      char buffer[20];
      int i, j;
      int *ret_vett;

     ret_vett=(int*)malloc(sizeof(int)*n); 
  
   for(i=0; i<n; i++){
       ret_vett[i]= 0;
       group_vect[i]++;
       merged_vect[i]++;      
       sprintf( a, "%d", group_vect[i] ); 
       sprintf( buffer, "%d", merged_vect[i] );
      
       
       strcat(buffer, a);
       
       ret_vett[i]=atoi(buffer);
   } 

   
   int count =0, temp;
   for(i=0; i<n; i++){
        if (ret_vett[i]>count-1){   
                temp=ret_vett[i];
                ret_vett[i]=count;
           for(j=i+1; j<n; j++){
               if(ret_vett[j]==temp )
                     ret_vett[j]=count;
            }
         }
      count++;
   }
 
   return ret_vett;

}
double ** allocate_R_matrix(int num_group, int nVertices, int *group_vect){

          int i,j,z,c=0;   
          double ** a;
          /* allocate each row of the 2D adjacency matrix */
        
	//for (i = 0; i < num_group; i++) {

		a = (double**)malloc(sizeof(double*) * num_group);		
	//}
        
       
	/* allocate  each data structure in the matrix */
	for (i = 0; i < num_group; i++){
		
	               a[i] = (double*)malloc(sizeof(double)* nVertices);
                                      
             }

         for (i = 0; i < num_group; i++){//fisso la riga i -> indice del gruppo

                   c=count(group_vect, nVertices, i);// conto quanti nodi ci sono nel gruppo i
                   if (c==0){
                         printf("\n\nR MATRIX has zero component in %d row\n\n", i);
                        exit(1);} 
                 //printf("\ni = %d, c = %d 1/sqrt(c) = %f \n",i, c, 1/sqrt(c));                                    
		for (j = 0; j < nVertices; j++){// riempio le colonne della riga i
	               
                        
                             if(group_vect[j]==i)                         
                               a[i][j] = 1/sqrt(c);
                             else a[i][j]=0;
                      
                       
               }
                
             }   
          printf("\n");


     printf("\n\n---------------R MATRIX ---------------\n\n");
       for (i = 0; i < num_group; i++){
                  
                 for (j = 0; j < nVertices; j++){

                         printf("%5.2lf ",a[i][j]);
           }
                printf("\n\n");
         }
	               
          printf("\n");

        return a;
                
}


double ** allocate_Coarse_matrix(int num_group, int nVertices, adj_matrix* Graph, double **R){
         
       int a, b, i,j;   
       double **mat;
          
       mat = (double**)malloc(sizeof(double*) * num_group);		

              
	/* allocate  each data structure in the matrix */
	for (i = 0; i < num_group; i++){
		
	               mat[i] = (double*)malloc(sizeof(double)* nVertices);
                                      
             }

      
       for(i=0; i<num_group; i++){
                   for(j=0; j<num_group; j++){
                       mat[i][j]=0;
                       
                   }
            
          }
              
      
     for(a=0; a<num_group; a++){
              
           for(b=0; b<num_group; b++){ 
                   
               for(i=0; i<nVertices; i++){
                         
                   for(j=0; j<nVertices; j++){
                            
                       mat[a][b]+=R[a][i]*(Graph->edgeMatrix[i][j]->edge)*R[b][j];
                     
                   }
                }
            }
        }
      
     printf("-\n-----------Coarse Grained Matrix--------------\n\n");
     for(i=0; i<num_group; i++){
                   for(j=0; j<num_group; j++){
                       printf("%5.2lf ",mat[i][j]);
                   } printf("\n\n");
          }   
       printf("\n");


      return mat;

}

void DestroyDoubleMatrix(double ** Mat, int nrow)
{
	int i;

	if (Mat != NULL) {

			for (i = 0; i < nrow; i++){
                             
                             free(Mat[i]);
                           }
			
			free(Mat);
		
		}
}

void DestroyDoubleVector(double * vet)
{


	if (vet != NULL) {
		free(vet);
		
		}
}


int count(int * v_group, int nVertices, int num){
        int count=0; 
        int i;
        
        for (i = 0; i <nVertices ; i++) {
              
		if (v_group[i]==num)                     
                                count++;	
	}
        
        return count;

}

double * calc_phi_Coarse(double **R, double *phi_vet, int num_group, int nVertices){

      int i, j;
      double *phi_coarse;


      phi_coarse=allocate_double_vector(num_group);
              
              for(i=0; i<num_group; i++){
                   for(j=0; j<nVertices; j++) {
                            phi_coarse[i]+=R[i][j]*phi_vet[j];
                            
              }printf("%8.2lf ",phi_coarse[i]);   
          }printf("\n");   

           return phi_coarse;

}

int find_max_int(int *vett, int n){
 int i;
 int max=-INT_MAX;
             for(i=0; i<n; i++){
                   if (vett[i]>max)
                       max=vett[i];
              }

return max;
}


int * Coarse_Graining(adj_matrix* Graph, FILE *info, int nEigenvalues){


    double ** eigenvectors, * eigenvalues; 
    int *group_vect, *posit, *merged_vect;
    int num_group, nEigen,i,j, z, pos_max_eigv,l,t;
    double max_eigv;
    double *temp_eigv, *sorted_eigv;
    int letter;
    double **data=NULL;
    eigenvectors = allocate_eigvec_matrix(Graph);
    eigenvalues = allocate_double_vector(Graph->nVertices);
    merged_vect = (int*)malloc(sizeof(int) * Graph->nVertices); 

     eigen(eigenvectors,eigenvalues, Graph->nVertices);
         
     nEigen=nEigenvalues;
     
         
     // sort eigenvalues into ascending order and remember their positions (vector posit)
     sorted_eigv=allocate_double_vector(Graph->nVertices);
     posit= (int*)malloc(sizeof(int) * Graph->nVertices);
     for (i=0;i<Graph->nVertices; i++){
           sorted_eigv[i]=fabs(eigenvalues[i]);
           posit[i]=i;
     }


     //debug
         
    sort(sorted_eigv, posit, 0, Graph->nVertices-1);
    fprintf(info,"eigenvalues sorted by descending abs:\n");
    printf("\n-------autovalori  ordinati: ----------\n");
    for(i=0; i<Graph->nVertices; i++){
                  printf("%2.3lf, ", eigenvalues[posit[i]]);
                  fprintf(info,"%2.3lf, ", eigenvalues[posit[i]]);
            }
     fprintf(info,"\n\n");
    // end debug
     printf("\n\n\nChoose grouping alghorithm: Default (each number) or KMeans (1)\n");
     scanf("%d",&letter);  
     fprintf(info,"\n[preserved eigenvalue-->eigenvector of original network in descending order]\n" ); 
     fprintf(info,"[grouping informations]\n\n\n" );
      if (letter==1){
             
                 
                 data =allocate_double_matrix(Graph->nVertices, 1);
                      
                 for(i=0; i<nEigen; i++)  {  
                          max_eigv = sorted_eigv[Graph->nVertices-i-1];
                          pos_max_eigv=posit[Graph->nVertices-i-1];
                          printf("\n\nmax eigv = %2.3lf in position %d\n", max_eigv, pos_max_eigv);
                          fprintf(info,"%1.4lf--> ", max_eigv ); 
                             
                          temp_eigv = allocate_double_vector(Graph->nVertices);
                          for(z=0; z<Graph->nVertices; z++){
                               if (fabs(eigenvalues[z])==max_eigv && z==pos_max_eigv)
                                    for(j=0; j<Graph->nVertices; j++) {
                                       temp_eigv[j]=eigenvectors[j][z];
                                       fprintf(info,"%1.4lf, ", temp_eigv[j] );
                             
                                     }                                                            
                             }
                                      
                          
                                  
                          for(l=0; l<Graph->nVertices; l++)                             
                             data[l][0]=temp_eigv[l];
                                       
                                      
                                           
                          printf("\n\n\nInsert number of groups for the %d` biggest eigenvalue:\n", i );
                          scanf("%d",&num_group);
                          printf("\n");
                          fprintf(info,"\n\n number of groups for the %d` biggest eigenvalue = %d \n\n", i, num_group );  
                          group_vect = (int*)malloc(sizeof(int) * Graph->nVertices);
                          group_vect=k_means(data, Graph->nVertices, 1, num_group, 0.0001, NULL);
                           //num_group=adjust_group(group_vect, Graph->nVertices, num_group);
      
                         if(i==0){
                              
                           merged_vect=group_vect;
                           if (nEigen==1) return merged_vect;
                           }
                         else{
                               
                             merged_vect=merge_group2(Graph->nVertices, group_vect, merged_vect);
                  
                             } 
                         free(temp_eigv);
                         free(group_vect);
                      }
                   }                                                            
      else {   
        
          for(i=0; i<nEigen; i++)  {  
             //take the biggest eigenvalue and insert the number of groups for the relative eigenvector 
               
              max_eigv = sorted_eigv[Graph->nVertices-i-1];
              pos_max_eigv=posit[Graph->nVertices-i-1];
              printf("\n\nmax eigv = %2.3lf in position %d\n", max_eigv, pos_max_eigv);    
              fprintf(info,"%1.4lf--> ", max_eigv );
              temp_eigv = allocate_double_vector(Graph->nVertices);
              
              for(z=0; z<Graph->nVertices; z++){
                   if (fabs(eigenvalues[z])==max_eigv && z==pos_max_eigv)
                         for(j=0; j<Graph->nVertices; j++) {
                              temp_eigv[j]=eigenvectors[j][z];
                              fprintf(info,"%1.4lf, ", temp_eigv[j] );
                             
                          }  
               }
              
              printf("\n\n\nInsert number of groups for the %d` biggest eigenvalue:\n", i );
              scanf("%d",&num_group);
              printf("\n");
              fprintf(info,"\n\n number of groups for the %d` biggest eigenvalue = %d \n\n", i, num_group );  
              group_vect = (int*)malloc(sizeof(int) * Graph->nVertices);
              group_vect=find_group(temp_eigv,  Graph->nVertices,  num_group); 
              //num_group=adjust_group(group_vect, Graph->nVertices, num_group);
      
             if(i==0){
                  merged_vect=group_vect;
                  if (nEigen==1) return merged_vect;
                   
                   }
             else{
                 
                  merged_vect=merge_group2(Graph->nVertices, group_vect, merged_vect);                   
                    } 
            free(temp_eigv);
            //free(group_vect);
                     
            }//end for                                                                 
       }//end else
       
                    
       
      for (l=0; l<Graph->nVertices; l++)
                 printf("%d ", merged_vect[l]);
     
      DestroyDoubleMatrix(eigenvectors, Graph->nVertices);
      if (data!=NULL ) DestroyDoubleMatrix(data, Graph->nVertices);
      
      return merged_vect;
  

}


void grado(adj_matrix* Graph, FILE *info){

  int i,j,c; 
  float mean_k=0, mean_k_square=0;
  float avg, avg2, Tc;
  float *k_adj;
  k_adj=(float*)malloc(sizeof(float) * Graph->nVertices);

  for(i=0; i<Graph->nVertices; i++){
     k_adj[i]=0;   
     for(j=0; j<Graph->nVertices; j++) {
           
                k_adj[i]+=Graph->edgeMatrix[i][j]->edge;
               // printf("%1.1f + ",Graph->edgeMatrix[i][j]->edge);
             }
    //printf("= %f \n", k_adj[i]);
     
  }

  for(i=0; i<Graph->nVertices; i++){
         mean_k_square=mean_k_square+k_adj[i]*k_adj[i];
         mean_k=mean_k+k_adj[i]; 
        // printf("%1.5f ", k_adj[i]);
         
  }


  avg= mean_k/Graph->nVertices;
  avg2= mean_k_square/Graph->nVertices;
  Tc=avg2/avg;
  printf("\n<k> = %f , <k2> = %lf, <k2>/<k> = %lf, dim rete= %d\n", avg,avg2, Tc, Graph->nVertices);
   fprintf(info,"\n\n<k> = %f , <k2> = %lf, <k2>/<k> = %lf, network dimension = %d\n\n", avg,avg2, Tc, Graph->nVertices);
  free(k_adj);
 // scanf("%d",&c);

}


double  grado2(adj_matrix* Graph, int flag){
      int i,j; 
      double grado=0;                        
 if (flag==0){
     for(i=0; i<Graph->nVertices; i++)     
          for(j=0; j<Graph->nVertices; j++) 
           
                grado+=Graph->edgeMatrix[i][j]->edge;
    
     return grado/2;}
 else{ 
        for(i=0; i<Graph->nVertices; i++)     
          for(j=i; j<Graph->nVertices; j++) 
           
                grado+=Graph->edgeMatrix[i][j]->edge;
          return grado;
     }
               
}


void write_adj(adj_matrix* Graph, int flag){

  FILE *fp;
  int i,j;
  if (flag==0){

     if ((fp = fopen("adj_matrix", "w")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
  }
  else{
     if ((fp = fopen("coarse_adj_matrix", "w")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
   }
     
  fprintf(fp,"%d\n",Graph->nVertices);
  for (i=0; i<Graph->nVertices; i++){
     for (j=0; j<Graph->nVertices; j++){
         
          fprintf(fp,"%lf ",Graph->edgeMatrix[i][j]->edge);
         
     }fprintf(fp,"\n");
   }
  fclose(fp);

}


adj_matrix* read_adj(int nVertices, int flag){
   FILE *fp;
   int i,j;
   adj_matrix* Graph;
   double val;
   if (flag==0){

     if ((fp = fopen("adj_matrix", "r")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
  }
  else{
     if ((fp = fopen("coarse_adj_matrix", "r")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
  }
  fscanf(fp,"%d",&nVertices);
  if (nVertices==0){
         printf("Error: Network with 0 vertices\n");
         exit(1);
  }      
  Graph = CreateAdjMatrix(0.0, nVertices);
  if (Graph == NULL) {
        printf("error in allocating adjacency matrix\n");
        exit(1);
  }
  for (i=0; i<nVertices; i++)
        for (j=0; j<nVertices; j++){
            fscanf(fp,"%lf",&val);
            Graph->edgeMatrix[i][j]->edge=val;
        }
  fclose(fp);

return Graph;
}

void write_edge(adj_matrix* Graph, int flag){

   FILE *fp;
   int i,j;
    if (flag==0){

     if ((fp = fopen("edge_matrix", "w")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }  
   fprintf(fp,"%d\n", Graph->nVertices);
   for (i=0; i<Graph->nVertices; i++){
     for (j=i; j<Graph->nVertices; j++){
       if (Graph->edgeMatrix[i][j]->edge!=0)
          fprintf(fp,"%d %d\n",i,j);
      }//fprintf(fp,"\n");
    }

  }
  else{
     if ((fp = fopen("coarse_edge_matrix", "w")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
     fprintf(fp,"%d\n", Graph->nVertices);
    for (i=0; i<Graph->nVertices; i++){
     for (j=i; j<Graph->nVertices; j++){
       if (Graph->edgeMatrix[i][j]->edge!=0)
          fprintf(fp,"%d %d %lf\n",i,j, Graph->edgeMatrix[i][j]->edge );
      }//fprintf(fp,"\n");
    }
   }
   
   fclose(fp);

}

adj_matrix* read_edge(int nVertices, int flag){
  FILE *fp;
  int i,j;
  adj_matrix* Graph;
  int val1, val2;
  double val3;
  int max, temp;
   if (flag==0){

     if ((fp = fopen("edge_matrix", "r")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
    fscanf(fp, "%d", &max);
    nVertices=max;
    Graph = CreateAdjMatrix(0.0, nVertices);
    if (Graph == NULL) {
       printf("error in allocating adjacency matrix\n");
       exit(1);
     } 

    while (fscanf(fp, "%d %d", &val1, &val2) != EOF) {

          Graph->edgeMatrix[val1][val2]->edge=1.0;
          Graph->edgeMatrix[val2][val1]->edge=1.0;
 
    }

  }
  else{
     if ((fp = fopen("coarse_edge_matrix", "r")) == NULL)
     {
        printf("Can t open file adj_matrix\n");
        exit(1);
     }
    fscanf(fp, "%d", &max);
    nVertices=max;
    Graph = CreateAdjMatrix(0.0, nVertices);
    if (Graph == NULL) {
       printf("error in allocating adjacency matrix\n");
       exit(1);
    }

    while (fscanf(fp, "%d %d %lf", &val1, &val2, &val3) != EOF) {

          Graph->edgeMatrix[val1][val2]->edge=val3;
          Graph->edgeMatrix[val2][val1]->edge=val3;
 
    }
  }
  


  fclose(fp);

  return Graph;
}


void write_pajek_adj(adj_matrix* Graph, int *group_vect){

  FILE *fp;
  int i,j;

  if ((fp = fopen("pajek_net.net", "w")) == NULL)
  {
  printf("Can t open file pajek\n");
  exit(1);
  }
  fprintf(fp,"*Vertices %d\n", Graph->nVertices);
  
  for (i=0; i<Graph->nVertices; i++)
       fprintf(fp,"%d \"%d\" %lf %lf 0.0 x_fact 5.000000 y_fact 5.000000\n", i+1, i+1,(200*sin(group_vect[i])+10*sin(i)+210)/420, (200*cos(group_vect[i])+10*cos(i)+210)/420); 
  
  fprintf(fp,"*Matrix\n");
  for (i=0; i<Graph->nVertices; i++){
         for (j=0; j<Graph->nVertices; j++){
             fprintf(fp,"%1.6lf ",Graph->edgeMatrix[i][j]->edge);
      }fprintf(fp,"\n");
  }
fclose(fp);

}

void write_pajek_coarse(adj_matrix* Graph, int nold, double **R ){

  FILE *fp;
  int i,j, count;

  if ((fp = fopen("pajek_net_coarse.net", "w")) == NULL)
  {
  printf("Can t open file edge_matrix\n");
  exit(1);
  }

 
  fprintf(fp,"*Vertices %d\n", Graph->nVertices);
  
    for (i=0; i<Graph->nVertices; i++){ 
       count = 0;
       count= count_group_elements(R, i, nold); 
       fprintf(fp,"%d \"%d\" 0.0 0.0 0.0 x_fact %1.6lf y_fact %1.6lf\n", i+1, i+1, 2*(float)count/sqrt(count), 2*(float)count/sqrt(count));
      }
  

  fprintf(fp,"*Matrix\n");
  for (i=0; i<Graph->nVertices; i++){
         for (j=0; j<Graph->nVertices; j++){
            fprintf(fp,"%1.6lf ",Graph->edgeMatrix[i][j]->edge);
         }fprintf(fp,"\n");
  }
  fclose(fp);

}

/*****
** kmeans.c
** - a simple k-means clustering routine
** - returns the cluster labels of the data points in an array
** - here's an example
**   extern int *k_means(double**, int, int, int, double, double**);
**   ...
**   int *c = k_means(data_points, num_points, dim, 20, 1e-4, 0);
**   for (i = 0; i < num_points; i++) {
**      printf("data point %d is in cluster %d\n", i, c[i]);
**   }
**   ...
**   free(c);
** Parameters
** - array of data points (double **data)
** - number of data points (int n)
** - dimension (int m)
** - desired number of clusters (int k)
** - error tolerance (double t)
**   - used as the stopping criterion, i.e. when the sum of
**     squared euclidean distance (standard error for k-means)
**     of an iteration is within the tolerable range from that
**     of the previous iteration, the clusters are considered
**     "stable", and the function returns
**   - a suggested value would be 0.0001
** - output address for the final centroids (double **centroids)
**   - user must make sure the memory is properly allocated, or
**     pass the null pointer if not interested in the centroids
** References
** - J. MacQueen, "Some methods for classification and analysis
**   of multivariate observations", Fifth Berkeley Symposium on
**   Math Statistics and Probability, 281-297, 1967.
** - I.S. Dhillon and D.S. Modha, "A data-clustering algorithm
**   on distributed memory multiprocessors",
**   Large-Scale Parallel Data Mining, 245-260, 1999.
** Notes
** - this function is provided as is with no warranty.
** - the author is not responsible for any damage caused
**   either directly or indirectly by using this function.
** - anybody is free to do whatever he/she wants with this
**   function as long as this header section is preserved.
** Created on 2005-04-12 by
** - Roger Zhang (rogerz@cs.dal.ca)
** Modifications
** -
** Last compiled under Linux with gcc-3
*/


int *k_means(double **data, int n, int m, int k, double t, double **centroids)
{
   /* output cluster label for each data point */
   int *labels = (int*)calloc(n, sizeof(int));

   int h, i, j; /* loop counters, of course :) */
   int *counts = (int*)calloc(k, sizeof(int)); /* size of each cluster */
   double old_error, error = DBL_MAX; /* sum of squared euclidean distance */
   double **c = centroids ? centroids : (double**)calloc(k, sizeof(double*));
   double **c1 = (double**)calloc(k, sizeof(double*)); /* temp centroids */

   assert(data && k > 0 && k <= n && m > 0 && t >= 0); /* for debugging */

   /****
   ** initialization */

   for (h = i = 0; i < k; h += n / k, i++) {
      c1[i] = (double*)calloc(m, sizeof(double));
      if (!centroids) {
         c[i] = (double*)calloc(m, sizeof(double));
      }
      /* pick k points as initial centroids */
      for (j = m; j-- > 0; c[i][j] = data[h][j]);
   }

   /****
   ** main loop */

   do {
      /* save error from last step */
      old_error = error, error = 0;

      /* clear old counts and temp centroids */
      for (i = 0; i < k; counts[i++] = 0) {
         for (j = 0; j < m; c1[i][j++] = 0);
      }

      for (h = 0; h < n; h++) {
         /* identify the closest cluster */
         double min_distance = DBL_MAX;
         for (i = 0; i < k; i++) {
            double distance = 0;
            for (j = m; j-- > 0; distance += pow(data[h][j] - c[i][j], 2));
            if (distance < min_distance) {
               labels[h] = i;
               min_distance = distance;
            }
         }
         /* update size and temp centroid of the destination cluster */
         for (j = m; j-- > 0; c1[labels[h]][j] += data[h][j]);
         counts[labels[h]]++;
         /* update standard error */
         error += min_distance;
      }

      for (i = 0; i < k; i++) { /* update all centroids */
         for (j = 0; j < m; j++) {
            c[i][j] = counts[i] ? c1[i][j] / counts[i] : c1[i][j];
         }
      }

   } while (fabs(error - old_error) > t);

   /****
   ** housekeeping */

   for (i = 0; i < k; i++) {
      if (!centroids) {
         free(c[i]);
      }
      free(c1[i]);
   }

   if (!centroids) {
      free(c);
   }
   free(c1);

   free(counts);

   return labels;
}

