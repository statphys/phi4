#include "phi_lib.h"

//initialize lattice
double *init_phi_lattice_coarse(adj_matrix* Graph,  double**R, int nold){

int i;
double c, phi_square=0, prod_scal=0, *temp, sum ;
double * phi_vet, *phi_coarse;
phi_vet = allocate_double_vector(nold);
phi_coarse = allocate_double_vector(Graph->nVertices);
sgenrand(time(NULL)); /*initializing seed for random number generator*/

/*initializing  phi_lattice to random configuration between 1.5 and -1.5 */

for (i=0; i<nold; i++){
            
         phi_vet[i]=(3*genrand()-1.5);
   //DEBUG
   // printf( "%f ", phi_vet[i]);

   }
  
  phi_coarse=calc_phi_Coarse(R, phi_vet, Graph->nVertices, nold);

  return phi_coarse;

}// end init



double *init_phi_lattice(adj_matrix* Graph){

int i;

double * phi_vet;
phi_vet = allocate_double_vector(Graph->nVertices);

sgenrand(time(NULL)); /*initializing seed for random number generator*/

/*initializing  phi_lattice to random configuration between 1.5 and -1.5 */

for (i=0; i<Graph->nVertices; i++){
            
         phi_vet[i]=3*genrand()-1.5;
   //DEBUG
    //printf( "%f ", phi_vet[i]);

   }
  

  return phi_vet;

}// end init




/*counting the initial total energy per site and the total magnetization per spin */
void init_observables(adj_matrix* Graph, double *phi_vet, double *E_xbond, double *M_xspin, float lambda, int flag, double **R, int nold){

int i,j, count;
double c, phi_square=0, prod_scal=0, *temp, sum=0 ;


i=0;
j=0;
for (i=0; i<Graph->nVertices; i++){
    
     *M_xspin+=phi_vet[i];
    } 
    

printf( "\nMxSPIN iniziale = %f\n ", *M_xspin);


temp = allocate_double_vector(Graph->nVertices);

//initalize temp vector
for (i=0; i<Graph->nVertices; i++){

       temp[i]=0;

}

// vett_trasposto*matrix = temp

for(j=0; j<Graph->nVertices; j++){
for (i=0; i<Graph->nVertices; i++) {

      temp[j] += Graph->edgeMatrix[i][j]->edge * phi_vet[i];

}
}
// temp*vett
for (i=0; i<Graph->nVertices; i++) {

      prod_scal += temp[i]*phi_vet[i];

}
    printf( "\nprod scal = %f, n = %d\n ", prod_scal, Graph->nVertices);

   
   for (i=0; i<Graph->nVertices; i++) {  
    
     phi_square=phi_vet[i]*phi_vet[i];
     
     
      if (flag==0)
         sum+= lambda*(phi_square-1)*(phi_square-1);
      else {
          count = 0;
          count= count_group_elements(R, i, nold);
           
          if(count!=0)
               sum+= lambda*(phi_square-count)*(phi_square-count)/count;
           
           }

    }

    *E_xbond = -prod_scal/2+sum;

   printf( "\nExBond iniziale = %f\n ", *E_xbond);
/*initial energy per bond and magnetization per spin*/

   free(temp);
//getchar();
}// end calc init observable/*end function*/





double totalEnergy(adj_matrix* Graph, double *phi_vet, double E_xbond, float lambda, int flag, double **R, int nold){
   
    int i,j, count;
   double phi_square=0, *temp, prod_scal, sum=0, term=0;
   
   temp = (double *)malloc(sizeof(double)*Graph->nVertices);

//initalize temp vector
for (i=0; i<Graph->nVertices; i++){

       temp[i]=0;
      // printf( "%f, ", phi_vet[i]);  
}

// vett_trasposto*matrix = temp

for(j=0; j<Graph->nVertices; j++){
for (i=0; i<Graph->nVertices; i++) {

      temp[j] += Graph->edgeMatrix[i][j]->edge * phi_vet[i];

}
}
// temp*vett
for (i=0; i<Graph->nVertices; i++) {

      prod_scal += temp[i]*phi_vet[i];

}
   E_xbond = E_xbond -prod_scal/2;
  // printf( "\nprod scal = %f\n\n (phi2 -1) \n", E_xbond); 
  for (i=0; i<Graph->nVertices; i++) {  
           phi_square=phi_vet[i]*phi_vet[i];
           term+=  lambda*(phi_square-1)*(phi_square-1);
          // printf( "%f, ",  lambda*(phi_square-1)*(phi_square-1));
     if (flag==0)
         E_xbond+= lambda*(phi_square-1)*(phi_square-1);
      else {
          count = 0;
          count= count_group_elements(R, i, nold);
          if (count!=0)
               E_xbond+= lambda*(phi_square-count)*(phi_square-count)/count;
           }
    }

     // printf( "\nterm = %f\n ", term);
     //printf( "\nE xbond = %f\n ", E_xbond); //getchar();
   free(temp);
return E_xbond;
}


double totalMagnetization(adj_matrix* Graph, double *phi_vet, double M_xspin,int flag, double **R, int nold){
     
 int i,j, count;
 
 for (i=0; i<Graph->nVertices; i++){
     if (flag==0)
          M_xspin+=phi_vet[i];
     else {
           count = 0;
           count= count_group_elements(R, i, nold);
           M_xspin+=phi_vet[i]/sqrt(count);

          }
           
    } 
return M_xspin;
}



double update_metropolis(double**R, double *phi_vet, int nold, int x, int flag){
        
        double phi_new=0;
        int i, count;

        
        if (flag==0){
                phi_new=phi_vet[x]+ (3*genrand()-1.5);
                return phi_new;
        }
        else{ 
          
          
            count = 0;
            //phi_new=phi_vet[x]+ (3*genrand()-1.5);
                for (i=0; i<nold; i++){
                     if (R[x][i]!=0.0)
                            count++;}
           if(count!=0.0)
               phi_new=phi_vet[x]+ (3*genrand()-1.5)/count;///////modificato
            //printf( "\nsono dentro update phi[%d] = %f, phi_new = %f, count = %d\n ", x, phi_vet[x], phi_new, count );
           //getchar();

             return phi_new;
            } 
 }


int count_group_elements(double ** R, int index, int nold){

            int i, count = 0;
            
                for (i=0; i<nold; i++)
                     if (R[index][i]!=0.0)
                              count++;
             return count;

}


void metropolis(adj_matrix* Graph, double *phi_vet, int x,  float Temp, float lambda, int flag, double **R, int nold){
    
    
    double phi_square_new, phi_square, phi_new, Enew , Eold;
    double delta, sum=0, c,d;
    int i, count; 
    

     phi_new=update_metropolis(R, phi_vet, nold, x, flag);
     phi_square_new=phi_new*phi_new;
    
    //DEBUG
    //printf( "\nsono dentro sweep->metropolis phi[%d] = %f, phi_new = %f\n ", x, phi_vet[x], phi_new );
    //getchar(); 
  for (i=0; i<Graph->nVertices; i++){
     
     if (Graph->edgeMatrix[x][i]->edge!=0){
           sum+=phi_vet[i]; 
           }
     }
     //DEBUG
     //printf( "\nsono dentro sweep->metropolis sum dei vicini del nodo %d e sum = %f\n ", x, sum);
     
     phi_square=phi_vet[x]*phi_vet[x];
     
     if (flag==0) {
            Eold=-sum*phi_vet[x]+lambda*(phi_square-1)*(phi_square-1);
            Enew=-sum*phi_new+lambda*(phi_square_new-1)*(phi_square_new-1);
      }
     else { 

          count = 0;
          count= count_group_elements(R, x, nold);
          Eold=-sum*phi_vet[x]+lambda*(phi_square-count)*(phi_square-count)/count;
          Enew=-sum*phi_new+lambda*(phi_square_new-count)*(phi_square_new-count)/count;
         }
 
      
     delta = Enew-Eold; 
 
     if (delta<=0) { /*the move is accepted*/
               
                
                phi_vet[x]=phi_new;
             //DEBUG
            // printf( "\nsono dentro sweep->metropolis flip accettato phi[%d] = %f, phi_new = %f\n ", x, phi_vet[x], phi_new );
              //printf( "\ndelta <0 cioe delta = %lf\n ", delta);  
      }


     else {   /*the move is accepted*/
        c=genrand();

        d=exp(-(1/Temp)*delta);
          if (c< exp(-(1/Temp)*delta)){

            
               phi_vet[x]=phi_new;
              //DEBUG
             //printf( "\nsono dentro sweep->metropolis flip accettato phi[%d] = %f, phi_new = %f\n ", x, phi_vet[x], phi_new );
              //printf( "\n%f < %f \n ",c, d ); 
              }                                         
         } 

        //getchar();
}//end function metropolis




void sweep (adj_matrix* Graph, double *phi_vet, double *E_xbond, double *M_xspin, FILE *out, float Temp, float lambda, int flag, double **R, int nold){

int i=0, x;
double grado;

sgenrand(time(NULL)); 

//metropolis sweep
for (i=0; i<(Graph->nVertices)*5; i++){ 
     
       
     x=(genrand())*Graph->nVertices; 
          
     metropolis(Graph, phi_vet, x, Temp, lambda, flag, R, nold);

     /*printf("(x,y) = (%d,%d) ", x, y);*/
     

 }//end for - sweep


    system ("clear");
    //print_lattice(phi_lattice);
    printf( "Temp = %f\n", Temp);    
    printf( "finiti i primi 5 sweep metropolis  \n");
    //for(i=0; i<Graph->nVertices; i++){
      //       printf( "%8.3lf ", phi_vet[i]);  
        //    }  printf( "\n");
    
    
   // getchar();

   x=(genrand())*Graph->nVertices; ; 
   

  wolf(Graph, phi_vet, x, Temp);
 // printf( "sono uscito da wolf \n"); getchar();


*E_xbond=totalEnergy(Graph, phi_vet, *E_xbond, lambda, flag, R, nold);
 *M_xspin=totalMagnetization(Graph, phi_vet, *M_xspin, flag, R, nold) ;

grado=grado2(Graph, flag);
*E_xbond=*E_xbond/(grado);
*M_xspin=*M_xspin/(Graph->nVertices);

 //fprintf(out, "%f %f\n", *E_xbond, *M_xspin);
   
    system ("clear");
     print_vet(Graph, phi_vet);
    printf( "Temp = %f", Temp);
    //getchar();    
}/*end function*/




void print_vet(adj_matrix* Graph, double *phi_vet){

int i, nmeno=0, npiu=0;

printf("\n");
for (i=0; i<Graph->nVertices; i++){
    
          

      if (phi_vet[i]>0.0 )   {npiu++; printf( "# "); }
      else if (phi_vet[i]<=0.0) {nmeno++; printf("- "); }
       
   
   
}

printf("\n");
printf("meno = %d piu = %d  punto à¼š\n", nmeno, npiu);



}/*end function */


// -----------------------------------------------------------------
// -----------------------------------------------------------------
// Wolff method 



// Grows positive cluster from a specified site 
void growClusterPos(adj_matrix* Graph, double *phi_vet, int x,  float Temp) {
    
    int i;
    
    for (i=0; i<Graph->nVertices; i++){
     
     if (Graph->edgeMatrix[x][i]->edge!=0.0)
          if (phi_vet[i] > 0 && clusterCheck(phi_vet, x, i, Temp)==TRUE)
               growClusterPos(Graph, phi_vet, i, Temp);
      
    }  
}

// Grows negative cluster from a specified site 
void growClusterNeg(adj_matrix* Graph, double *phi_vet, int x,  float Temp) {
    
     int i;
    
    for (i=0; i<Graph->nVertices; i++){
     
     if (Graph->edgeMatrix[x][i]->edge!=0.0)
          if (phi_vet[i] < 0 && clusterCheck(phi_vet, x, i, Temp)==TRUE)
               growClusterNeg(Graph, phi_vet, i, Temp);
      
    }  
}


void flipCluster(adj_matrix* Graph,double *phi_vet) {

   int i;

   for (i=0; i<Graph->nVertices; i++){
      
          if (cluster[i]==1)
                phi_vet[i]=-phi_vet[i];      
   } 
          
    
}

int clusterCheck(double *phi_vet, int x, int xnext, float Temp) {
    
    double pAdd;
   if (cluster[xnext]==1)
        return FALSE;       // Already in cluster
    
    pAdd = 1 -exp(-(2/Temp)*phi_vet[x]*phi_vet[xnext]);
    if (genrand() < pAdd) {
        cluster[xnext]=1;
        return TRUE;
    }
    return FALSE;
}


void wolf(adj_matrix* Graph, double *phi_vet, int x, float Temp){

    
    int i,j;
    
   // initialize cluster flag matrix     
   for (i=0; i<Graph->nVertices; i++){
         
       cluster[i]=0;  

     }
   // sire (x,y) belong to the cluster
    cluster[x]=1;
    //printf( "sono dentro wolf ho scelto il sito, (%d)\n", x);
    //printf( " psi vet(%d) = %8.3lf\n", x, phi_vet[x]);
    //getchar();
  
    if (phi_vet[x] > 0){
         growClusterPos(Graph, phi_vet, x, Temp);
        // printf( "ho applicato grow positiv\n");
        
         print_vet(Graph, phi_vet);
         //getchar();
         
        }
    else{
         growClusterNeg(Graph, phi_vet, x, Temp);
       //printf( "ho applicato grow negative\n");
         
        print_vet(Graph, phi_vet); 
        
       //getchar();  
  }
    
     
     flipCluster(Graph, phi_vet);
      //printf( "ho flippato il cluster\n");
         //getchar();
       // print_vet(Graph, phi_vet); 
        
   
    
}
// -----------------------------------------------------------------

