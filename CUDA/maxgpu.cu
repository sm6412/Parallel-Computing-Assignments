#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

__global__ void getmaxcu(unsigned int *, unsigned int *);
unsigned int getmax(unsigned int *, unsigned int);


__global__ 
void getmaxcu(unsigned int * gpu_numbers, unsigned int * gpu_max){
  extern __shared__ unsigned int block_data[];
  unsigned int tid = threadIdx.x;
  block_data[tid] = gpu_numbers[(blockIdx.x*blockDim.x) + threadIdx.x];
  __syncthreads();

  // reduce 
  for(unsigned int offset=1; offset < blockDim.x; offset *= 2) {
    int compare_val = 2*offset; 
    if (tid % (compare_val) == 0) {
      if(block_data[tid] < block_data[tid + offset]){
        block_data[tid] = block_data[tid + offset];
      }
    
    }
    __syncthreads();
  }
  // write block max to list of maxes 
  if (tid == 0){
    gpu_max[blockIdx.x] = block_data[0];
  } 
  __syncthreads();
}


int main(int argc, char *argv[])
{
  unsigned int size = 0;  // The size of the array
  unsigned int i;  // loop index
  unsigned int * numbers; //pointer to the array
	unsigned int * gpu_numbers; 
  unsigned int * host_max;
	unsigned int * gpu_max;
  unsigned int * gpu_max2;


    
  if(argc !=2)
  {
      printf("usage: maxseq num\n");
      printf("num = size of the array\n");
      exit(1);
  }
   
  size = atol(argv[1]);

  numbers = (unsigned int *)malloc(size * sizeof(unsigned int));
  if( !numbers )
  {
      printf("Unable to allocate mem for an array of size %u\n", size);
      exit(1);
  }    

  srand(time(NULL)); // setting a seed for the random number generator
  // Fill-up the array with random numbers from 0 to size-1 
  for( i = 0; i < size; i++)
      numbers[i] = rand()  % size;

  // print sequential answer 
  //printf("Correct answer is: %u\n", getmax(numbers, size));



	// allocate gpu memory for array of randomly generated numbers 
  int allocation_size = size * sizeof(unsigned int);
  cudaMalloc((void**)&gpu_numbers,allocation_size);
  cudaMemcpy(gpu_numbers,numbers,allocation_size,cudaMemcpyHostToDevice);
   
   // get block num, and threads per block
  int threads_per_block = 1024;
  int num_of_blocks = (int)ceil(size/(double)threads_per_block);

  // allocate gpu memory for max number 
  cudaMalloc((void**)&gpu_max,num_of_blocks*sizeof(unsigned int));
  
  // find maxes of each block
  getmaxcu<<<num_of_blocks,threads_per_block,threads_per_block*sizeof(unsigned int)>>>(gpu_numbers,gpu_max);
  
  // while there is still more than 1 block, continue
  // to reduce maxes
  while(num_of_blocks>1){
    // get new number of blocks
    num_of_blocks = (int)ceil(num_of_blocks/(double)threads_per_block);
    cudaMalloc((void**)&gpu_max2,num_of_blocks*sizeof(unsigned int));

    // rerun kernel
    getmaxcu<<<num_of_blocks,threads_per_block,threads_per_block*sizeof(unsigned int)>>>(gpu_max,gpu_max2); 

    // move over data
    cudaMalloc((void**)&gpu_max,num_of_blocks*sizeof(unsigned int));
    cudaMemcpy(gpu_max,gpu_max2,num_of_blocks*sizeof(unsigned int),cudaMemcpyDeviceToDevice);
  }

  // allocate memory for max on host
  host_max = (unsigned int *)malloc(num_of_blocks * sizeof(unsigned int));
  // copy max from device to host
  cudaMemcpy(host_max,gpu_max,num_of_blocks * sizeof(unsigned int),cudaMemcpyDeviceToHost);

  // display max
  printf(" The maximum number in the array is: %u\n", host_max[0]);

  // free memory
  cudaFree(gpu_numbers);
  cudaFree(gpu_max);
  free(numbers);
  free(host_max);
  exit(0);
}

unsigned int getmax(unsigned int num[], unsigned int size)
{
  unsigned int i;
  unsigned int max = num[0];

  for(i = 1; i < size; i++)
  if(num[i] > max)
     max = num[i];

  return( max );

}




