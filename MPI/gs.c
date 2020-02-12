#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */

/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */

/********************************/


/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;
  
  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);
    
    for(j = 0; j < num; j++)
       if( j != i)
     sum += fabs(a[i][j]);
       
    if( aii < sum)
    {
      printf("The matrix will not converge.\n");
      exit(1);
    }
    
    if(aii > sum)
      bigger++;
    
  }
  
  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}

/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;  
 
  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
  {
    printf("Cannot allocate a!\n");
    exit(1);
  }

 for(i = 0; i < num; i++) 
  {
    a[i] = (float *)malloc(num * sizeof(float)); 
    if( !a[i])
    {
        printf("Cannot allocate a[%d]!\n",i);
        exit(1);
    }
  }
 
 x = (float *) malloc(num * sizeof(float));
 if( !x)
  {
    printf("Cannot allocate x!\n");
    exit(1);
  }


 b = (float *) malloc(num * sizeof(float));
 if( !b)
  {
    printf("Cannot allocate b!\n");
    exit(1);
  }

 /* Now .. Filling the blanks */ 

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
    fscanf(fp,"%f ", &x[i]);
 
 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i][j]);
   
   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
 }
 
 fclose(fp); 

}


/************************************************************/
// this functions computes the new Xs
float* compute_new_xs(int start_range, int end_range, int work_per_pro){
    float* new_Xs = (float *) malloc(work_per_pro * sizeof(float));
    int counter = 0;
    for (int row = start_range; row < end_range; row++)
    {
        for (int col = 0; col < num; col++){
            if(col==0){
                new_Xs[counter] = b[row];
            }
            if (col != row){
                new_Xs[counter] = new_Xs[counter] - (a[row][col] * x[col]);
            }
        }
        new_Xs[counter] = new_Xs[counter]/a[row][row];
        counter++;
    }
    return new_Xs;

}

int main(int argc, char *argv[])
{

     int nit = 0; /* number of iterations */
     FILE * fp;
     char output[100] ="";
      
     if( argc != 2)
     {
       printf("Usage: ./gsref filename\n");
       exit(1);
     }
      
     /* Read the input file and fill the global data structure above */ 
     get_input(argv[1]);
     
     /* Check for convergence condition */
     /* This function will exit the program if the coffeicient will never converge to 
      * the needed absolute error. 
      * This is not expected to happen for this programming assignment.
      */

    /* Check for convergence condition */
    check_matrix();

    // my code starts here
    int comm_size;     
    int my_rank;        
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    
    // Since we don't want to use more processes than we have to, let's
    // determine how many we should use based on how many the user enters
    int num_of_pro = 0;
    if((num/comm_size)!=0){
        num_of_pro = comm_size;
    }
    else{
        num_of_pro = num;
    }

    
    // number of equations to use per process 
   int elements_to_use = 0;
   if(num%num_of_pro==0){
        elements_to_use = num/num_of_pro;
   }
   else{
        elements_to_use = (num/num_of_pro)+1;
   }
   
   
   // create arrays to hold the start and end positions that each process
   // will use in order to compute the different equations 
   int* start_per_pro = (int *) malloc(num_of_pro * sizeof(int));
   int* end_per_pro = (int *) malloc(num_of_pro * sizeof(int));

   // initialize arrays
   for(int x=0;x<num_of_pro;x++){
        start_per_pro[x] = 0;
        end_per_pro[x] = 0;
   }

   // fill arrays with correct start and end pos
   int pos=0;
   for(int x=0;x<num_of_pro;x++){
        if((pos+elements_to_use)<num){
            start_per_pro[x] = pos;
            end_per_pro[x] = pos+elements_to_use;
            pos+=elements_to_use;
        }
        else{
            start_per_pro[x] = pos;
            end_per_pro[x] = num;
            break;
        }
   }


   // finally give each process its start and end range 
   int start_range = 0;
   int end_range = 0;
   if(my_rank<num_of_pro){
    start_range = start_per_pro[my_rank];
    end_range = end_per_pro[my_rank];
   }
   else{
    start_range = my_rank;
   }

   // determine the amount of work each process will do
   int work_per_pro;
   if(my_rank<num_of_pro){
    work_per_pro = end_range - start_range;
    if(work_per_pro==0){
        work_per_pro = 1;
    }
   }
   else{
    work_per_pro = 1;
   }

   // establish arrays to be used in computation
    float* new_Xs = (float *) malloc(work_per_pro * sizeof(float));
    float* final_new;
    if(num < comm_size){
        final_new = (float *) malloc(comm_size * sizeof(float));
    }
    else{
        final_new = (float *) malloc(num * sizeof(float));
    }

    // compute the new Xs
    bool continue_loop = true;
    while (continue_loop){   
        nit++;
        if (my_rank < num_of_pro){
            // compute Xs
            new_Xs = compute_new_xs(start_range,end_range,work_per_pro);
        } 

        MPI_Allgather(new_Xs, work_per_pro, MPI_FLOAT, final_new, elements_to_use, MPI_FLOAT, MPI_COMM_WORLD);

        // determine whether the loop can end or not based on whether the error
        // is low enough
        int correct = 0;
        for (int i = 0; i < num; ++i){
            if ((fabs((final_new[i] - x[i])/(final_new[i]))) <= err){
                correct++;
            }
            x[i] = final_new[i];
        }

        if (correct == num){
            continue_loop=false;
        }
    } 

   
    free(final_new);
    free(new_Xs);
    free(start_per_pro);
    free(end_per_pro);

    // end mpi part
    MPI_Finalize();

    /* Writing results to file */
    sprintf(output,"%d.sol",num);
    fp = fopen(output,"w");
    if(!fp)
    {
        printf("Cannot create the file %s\n", output);
        exit(1);
    }
        
    for(int i = 0; i < num; i++){
       fprintf(fp,"%f\n",x[i]);
   }
     
    printf("total number of iterations: %d\n", nit);
     
    fclose(fp);
    
    exit(0);

} 