#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

// determine how long it takes to generate prime numbers
double tstart = 0.0;
double tend = 0.0;
double ttaken;
FILE * fp;
char output[100] ="";

int main(int argc, char* argv[]) {
	// read in command line arguments
	// N is a positive number bigger than 2 and less than 100,000
	int N = atoi(argv[1]);
	// positive number of threads, does not exceed 100
	int t = atoi(argv[2]);

	// start code to generate prime numbers
	tstart = omp_get_wtime();

	// array to hold all numbers up to N
	int *N_array = (int *) malloc((N-1) * sizeof(int));

	// create array to hold whether the number in N_array is prime or not
	bool *not_prime = (bool *) malloc((N-1) * sizeof(bool));
	
	// populate array that contains N numbers and
	// whether they are prime or not
	int count = 2;
	for (int i = 0; i < N-1; ++i)
	{
		N_array[i] = count;
		if(i == 0){
			not_prime[i] = false;
		}
		else{
			not_prime[i] = true;
		}
		count++;
		
	}

	// gather number that will mark when prime number generation must stop 
	int end_condition = (N+1)/2;
	
	# pragma omp parallel num_threads(t)
	for (int i = 0; i < N-1; ++i)
	{
		int current_num = N_array[i];
		// if the number is prime and larger than the stop number, stop
		// prime number generation 
		if((not_prime[i] == false) && (current_num > end_condition)){
			break;
		}
		// if number is prime, divide every number by it 
		else if(not_prime[i] == false){
			// assign a number to each thread 
			# pragma omp for
			for (int j = i+1; j < N-1; j++){
                int compare_num = N_array[j];
                // if the prime number is 2, initially divide every 
                // number by it to determine current non-prime numbers
                if(current_num == 2){
                	if(compare_num%current_num != 0){
                		not_prime[j] = false;
	                }
				}
				// if the number is not 2, divide all current prime
				// numbers by it to determine which are actually not
				// prime
                else{
                	if(not_prime[j] == false){
                		if(compare_num%current_num == 0){
                			// make the number not prime
                			not_prime[j] = true;
		                }
					}// end if 
				}// end else 

			}
		}
	}



	

	// reveal how much time it took to generate prime numbers
	ttaken = omp_get_wtime() - tstart;
	printf("Time taken for the main part: %f\n",ttaken );


	
	/** Writing results to file **/
  	sprintf(output,"%d.txt",N);
  	fp = fopen(output,"w");
  	if(!fp)
  	{
    	printf("Cannot create the file %s\n", output);
    	exit(1);
  	}

  	int position = 1;
  	for (int i = 0; i < N-1; ++i)
	{
		if(not_prime[i] == false){
			int prime_num = N_array[i];
			int diff = 0;
			if(prime_num != 2){
				int back_check_pos = i-1;
				bool is_prime = not_prime[back_check_pos];
				while(is_prime != false){
					back_check_pos--;
					is_prime = not_prime[back_check_pos];
				}
				diff = prime_num - N_array[back_check_pos];
			}

			fprintf(fp,"%d, %d, %d\n",position,N_array[i],diff);
			position++;
		}
	}

	// free memory
	free(N_array);
	free(not_prime);

	exit(0);
}