#include "phi_lib.h"

int main(void) {

/*MC varable declarations*/
int step;
int i,t, j, z, n,x,b;
char s;
FILE  *data, *out, *out_new, *fp, *info;
double E_xbond=0, M_xspin=0, mag;
float T;
char name[10];
double E_sum,E2_sum, M_sum, M_sum_abs, M_Mit_abs, M2_sum, E_Mit,M_Mit,E2_Mit,M2_Mit;
 double sheat,suscep;
 double sigmam2;
int equilibration_time=50;
float lambda=10;
int steps_per_site=1000;

/* Network variable declarations*/
adj_matrix* Graph;
double prob1, prob2, prob3, prob;
double *phi_vet;
int nVertices;
char letter='o';
int letter1,h;
double val, autocorTime=0;

/* Coarse Grained Network variable declarations*/

 int *group_vect=NULL, nEigenvalues;
 int num_group, nold;
 double **R, **A_Coarse;
 adj_matrix* Graph_Coarse;
 double *eigvalues_Coarse, *sorted_eigv;
 double **eigvectors_Coarse,max_eigv;
 int *posit, pos_max_eigv;
if ((info = fopen("info.txt", "w")) == NULL)
  {
  printf("Can t open file info\n");
  exit(1);
  }
fprintf(info,"\n--------------------Networks informations-----------\n\n");

while( letter!='F'  & letter!='f' & letter!='C' & letter!='c' ){
     
     printf("Read Adj Matrix from File (F), Create Adj Matrix (C), Exit (E)\n");
     scanf("%c", &letter);
     switch (letter){
     case 'C': case 'c': 
              printf("Create totally random network (all numbers), create grouped random network (1)\n");
              scanf("%d", &letter1);
              switch (letter1){
              case 1: 
                 printf("\nInsert number of vertices: (n must be 32, 64, 128, 256, ...)");
                 scanf("%d", &nVertices);
                 printf("\n");
                 if(nVertices%8!=0){
                       printf("n isn t divisible for 8\n ");
                       exit(1); 
                 }
                 printf("Insert low probability : ");
                 scanf("%lf", &prob1);
                 printf("Insert middle  probability : ");
                 scanf("%lf", &prob2);
                 printf("Insert high probability : ");
                 scanf("%lf", &prob3);
                
                 /* read the number of vertices */
                 
                 Graph = CreateAdjMatrixNew(prob1, prob2, prob3,nVertices);
                 if (Graph == NULL) {
                       printf("error in allocating adjacency matrix\n");
                       exit(1);
                 }   
                 write_adj(Graph, 0);                
                 ViewAdjMatrix(Graph);
                 
                 break;
             
              default:
                 printf("Insert probability : ");
                 scanf("%lf", &prob);
                 printf("\nInsert number of vertices: ");
                 scanf("%d", &nVertices);
                 printf("\n");

                 Graph = CreateAdjMatrix(prob,nVertices);
                 if (Graph == NULL) {
                       printf("error in allocating adjacency matrix\n");
                       exit(1);
                 } 
                 write_adj(Graph, 0);
                 write_edge(Graph, 0);
                 ViewAdjMatrix(Graph);
                 break;             
              }   
            break;
      case 'F': case 'f': 
              Graph=read_edge(nVertices, 0);    
              ViewAdjMatrix(Graph);
             
              break;


      case 'E':
      case 'e': printf("Exiting program... \n"); exit(1); break;
                 

    }//end switch
}//end while

grado(Graph, info);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                         COARSE GRAINING NETWORK                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
         
     printf("\nInsert number of eigenvalues you want to preserve: \n", nEigenvalues );
     scanf("%d",&nEigenvalues);
     printf("\n");
     if (nEigenvalues==0){
         printf("The net can not be Coarse Grained ");
      exit(1);}
          
      fprintf(info,"\n------The number of preserved eigenvalues is %d------\n\n", nEigenvalues );

      group_vect=Coarse_Graining(Graph, info, nEigenvalues);
      if (group_vect==NULL){
             printf("The net can not be Coarse Grained ");
             exit(1);}
      //printf("\n---Sono nel Main \n");
      for (i=0; i<Graph->nVertices; i++)
            printf("%d ", group_vect[i]);
      printf("\n");
      num_group=find_max_int(group_vect, Graph->nVertices);
      num_group++;
      num_group=adjust_group(group_vect, Graph->nVertices, num_group);
      printf("\nnum grup  %d \n", num_group);
      fprintf(info,"\n------The final number of groups is %d------\n\n", num_group );
      
      for (i=0; i<Graph->nVertices; i++){
            printf("%d ", group_vect[i]);
            fprintf(info,"%d, ", group_vect[i] );}
      printf("\n");
      write_pajek_adj(Graph,group_vect);
      R=allocate_R_matrix(num_group,Graph->nVertices, group_vect);
      
         
      A_Coarse=allocate_Coarse_matrix(num_group, Graph->nVertices, Graph, R);
      
      posit=(int*)malloc(sizeof(int) * num_group);
  
      
      eigvectors_Coarse = allocate_eigvec_matrix2(A_Coarse, num_group);
      eigvalues_Coarse = allocate_double_vector(num_group);
      eigen(eigvectors_Coarse,eigvalues_Coarse, num_group);
      
      sorted_eigv=allocate_double_vector(num_group);
      posit= (int*)malloc(sizeof(int) * num_group);
      for (i=0;i<num_group; i++){
           sorted_eigv[i]=fabs(eigvalues_Coarse[i]);
           posit[i]=i;
     }


     //debug
    fprintf(info,"\n\neigenvalues sorted by descending abs:\n");    
    sort(sorted_eigv, posit, 0, num_group-1);
    printf("\n-------autovalori  ordinati: ----------\n");
    for(i=0; i<num_group; i++){
         fprintf(info,"%2.3lf, ", eigvalues_Coarse[posit[i]]);   
         printf("%2.3lf, ", eigvalues_Coarse[posit[i]]);
    }
              
     fprintf(info,"\n\n");
    fprintf(info,"\n[preserved eigenvalue-->eigenvector of original network in descending order]\n\n" );  
    for(i=0; i<nEigenvalues; i++)  {  
            max_eigv = sorted_eigv[num_group-i-1];
            pos_max_eigv=posit[num_group-i-1];            
            fprintf(info,"%1.4lf--> ", max_eigv ); 
            
           for(z=0; z<Graph->nVertices; z++){
                 if (fabs(eigvalues_Coarse[z])==max_eigv && z==pos_max_eigv)
                      for(j=0; j<num_group; j++) {
                                      
                           fprintf(info,"%1.4lf, ", eigvectors_Coarse[j][z] );
                             
                      }                                                            
                  }
              fprintf(info,"\n\n");                        
            }  
    // copy Coarse matrix into a Graph
      
    Graph_Coarse = CreateAdjMatrix(0, num_group);
    for(i=0; i<num_group; i++){
        for(j=0; j<num_group; j++) {
              Graph_Coarse->edgeMatrix[i][j]->edge=A_Coarse[i][j];
                            
        }printf("\n\n");   
    }printf("\n");   
      
    grado(Graph_Coarse, info);     
    printf("\n----------THE END OF COARSE---------- \n");
    write_adj(Graph_Coarse, 1);
    write_edge(Graph_Coarse, 1);  
    
    nold=Graph->nVertices;
    write_pajek_coarse(Graph_Coarse,nold, R ); 
    printf("\n");

    fclose(info);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                         MONTE CARLO SIMULATION                                          //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
// allocate memory for vector

cluster = (int *)malloc(sizeof(int)*Graph->nVertices);
double M_vet_abs[steps_per_site];
printf("\nInsert number of step: ");
scanf("%d",&step);
printf("\n");


printf("insert a char to start the simulation\n");
//getchar();


if ((data = fopen("measurement", "a")) == NULL)
   {
      printf("Can t open file measurement");
      exit(1);
   }

 for (x = 1; x <= 100; x += 1)
 {
   
  T = x/10.0;
  E_sum = 0;
  M_sum = 0;
  E2_sum = 0;
  M2_sum = 0;
  M_sum_abs=0;
   E_xbond=0;
    M_xspin=0;
  phi_vet=init_phi_lattice(Graph);
  init_observables(Graph,phi_vet,&E_xbond, &M_xspin, lambda,0, R, 0);
  //print_lattice(lattice);
  printf( "Temp = %f\n", T); 
  //getchar();
   

   //sprintf(name, "data_%.1f", T);
   
   //if ((out = fopen(name, "a")) == NULL)
   //{
   //   printf("Can t open file %sn", name);
   //   exit(1);
   //}



  for (i=0;i<equilibration_time;i++)
  {
   sweep(Graph, phi_vet, &E_xbond, &M_xspin, out, T, lambda,0, R, 0);
  }
  
 //printf("Ho equilibrato, Temp = %f, M_xspin all equilibrio = %lf \n", T, M_xspin);
  //getchar();
  
  n=1;
  for (i=0;i<steps_per_site;i++)
  {
    
     sweep(Graph, phi_vet, &E_xbond, &M_xspin, out, T, lambda,0, R, 0);
    
    n++;
    M_sum += M_xspin;
    M_sum_abs += fabs(M_xspin);
    M_vet_abs[i]=fabs(M_xspin);
    E_sum += E_xbond;
    E2_sum += E_xbond*E_xbond;
    M2_sum +=  M_xspin*M_xspin;
   
 
  }
   
   M_Mit = M_sum / (steps_per_site);
   M_Mit_abs = M_sum_abs / (steps_per_site);
   M2_Mit = M2_sum / (steps_per_site);
   E_Mit = E_sum / (steps_per_site);
   E2_Mit = E2_sum / (steps_per_site);
   
   //sprintf(name, "data_corr_%.1f", T);   
   //if ((out_new = fopen(name, "a")) == NULL)
   //{
     // printf("Can t open file %sn", name);
      //exit(1);
   //}

   //autocorTime=calcAutocor(steps_per_site, M_vet_abs, M_Mit_abs, T, out_new);
   
   sigmam2 = (1/(steps_per_site-1)*(M2_Mit-M_Mit_abs*M_Mit_abs));
    
   suscep = (1/T)*(Graph->nVertices)*(M2_Mit-M_Mit_abs*M_Mit_abs);
   sheat = (1/(T*T))*(Graph->nVertices)*(E2_Mit-E_Mit*E_Mit);
  fprintf(data,"%f %.6le %.6le %.6le %.6le %.6le %.6le %.6le %.6le %.6le %6le\n", T,M_Mit,M_Mit_abs,M2_Mit,E_Mit, E2_Mit, sigmam2,E_xbond,suscep,sheat,    autocorTime);
 
 }
//fclose (out_new);
//fclose(out);
fclose(data);
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                         MONTE CARLO COARSE SIMULATION                                   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

printf("\n\n\n----------------II simulation--------------\n\n\n");


// allocate memory for vector
double *psi_vet;

free(cluster);
cluster = (int *)malloc(sizeof(int)*Graph_Coarse->nVertices);

//printf("\nGet char to start coarse simulation: ");
//scanf("%d",&step);


sprintf(name, "measurement_coarse_N%d_ng%d", nold, Graph_Coarse->nVertices );

if ((data = fopen(name, "a")) == NULL)
   {
      printf("Can t open file measurement");
      exit(1);
   }

 for (x = 1; x <= 120; x += 1)
 {
   
  T = x/10.0;
  E_sum = 0;
  M_sum = 0;
  E2_sum = 0;
  M2_sum = 0;
  M_sum_abs=0;
   E_xbond=0;
    M_xspin=0;
 
  psi_vet=init_phi_lattice_coarse(Graph_Coarse,R, nold);
  //for (i=0;i<equilibration_time;i++){
  //  printf("%5.3f /n", psi[i])
   //  }

  printf("Sono prima di init observables");
   
  init_observables(Graph_Coarse,psi_vet,&E_xbond, &M_xspin, lambda,1, R, nold);
  //print_lattice(lattice);
  printf( "Temp = %f\n", T); 
  //getchar();
   

  // sprintf(name, "coarse_data_%.1f", T);
   
  // if ((out = fopen(name, "a")) == NULL)
  // {
  //    printf("Can t open file %sn", name);
  //    exit(1);
  // }



  for (i=0;i<equilibration_time;i++)
  {
   sweep(Graph_Coarse, psi_vet, &E_xbond, &M_xspin, out, T, lambda,1, R, nold);
  }
  
 
  for (i=0;i<steps_per_site;i++)
  {
   
     sweep(Graph_Coarse, psi_vet, &E_xbond, &M_xspin, out, T, lambda,1, R, nold);
    
   
    M_sum += M_xspin;
    M_sum_abs += fabs(M_xspin);
    E_sum += E_xbond;
    E2_sum += E_xbond*E_xbond;
    M2_sum +=  M_xspin*M_xspin;
   
 
  }
   
   M_Mit = M_sum / (steps_per_site);
   M_Mit_abs = M_sum_abs / (steps_per_site);
   M2_Mit = M2_sum / (steps_per_site);
   E_Mit = E_sum / (steps_per_site);
   E2_Mit = E2_sum / (steps_per_site);
 
  
   sigmam2 = (1/(steps_per_site-1)*(M2_Mit-M_Mit_abs*M_Mit_abs));
     
   suscep = (1/T)*(Graph_Coarse->nVertices)*(M2_Mit-M_Mit_abs*M_Mit_abs);
   sheat = (1/(T*T))*(Graph_Coarse->nVertices)*(E2_Mit-E_Mit*E_Mit);
  fprintf(data,"%f %.6le %.6le %.6le %.6le %.6le %.6le %.6le %.6le %.6le\n", T,M_Mit,M_Mit_abs,M2_Mit,E_Mit, E2_Mit, sigmam2,E_xbond,suscep,sheat);
 }
 

//fclose(out);
fclose(data);

DestroyDoubleMatrix(R, num_group);
DestroyDoubleMatrix(A_Coarse, num_group);


DestroyDoubleVector((double *)group_vect);

DestroyAdjMatrix(Graph); 

DestroyAdjMatrix(Graph_Coarse); 
return 0;

}/*end main*/
