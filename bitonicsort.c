#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

struct timeval startwtime, endwtime;
double seq_time;

int thread_count;

int searchNearestN(int n);

void print(int * a, int n);
void printToFile(int * a, int n, char * file);
void test(int * a, int n);
void rng(int * a, int n, int m);
void fillDummy(int * a, int from, int to, int dummy);

void exchange(int * a, int i, int j);
void impBitonicSort(int * a, int n);
void impBitonicSortParallel(int * a, int n);

int main(int argc, char **argv) {
  if (argc != 3) {
  // if (argc != 2) {
    printf("Usage: %s n\n  where n is problem size\n", argv[0]);
    exit(1);
  }

  thread_count = strtol(argv[2], NULL, 10);
  // thread_count = 4;
  int N = atoi(argv[1]);
  int N2 = searchNearestN(N);
  int * a = (int *) malloc(N2 * sizeof(int));
  double time_parallel = 0.0;
  double time_serial = 0.0;

  rng(a, N, N2);
  printToFile(a, N, "unsorted");

  // Parallel
  printf("Sort Parallel\n");
  for (int i = 0; i < 3; i++) {
	  rng(a, N, N2);
	  // print(a, N2);

	  gettimeofday (&startwtime, NULL);
	  impBitonicSortParallel(a, N2);
	  gettimeofday (&endwtime, NULL);

	  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
			      + endwtime.tv_sec - startwtime.tv_sec);
	  time_parallel += seq_time;

	  printf("%d - Imperative Parallel wall clock time = %f\n", i, seq_time);
	  test(a, N);
  }
  printf("Avg time %f\n\n", time_parallel/3);
  printToFile(a, N, "parallelsorted");
  // Parallel
  printf("Sort Serial\n");
  for (int i = 0; i < 3; i++) {
	  rng(a, N, N2);
	  // print(a, N2);

	  gettimeofday (&startwtime, NULL);
	  impBitonicSort(a, N2);
	  gettimeofday (&endwtime, NULL);

	  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
			      + endwtime.tv_sec - startwtime.tv_sec);
	  time_serial += seq_time;

	  printf("%d - Imperative wall clock time = %f\n", i, seq_time);
	  test(a, N);
  }
  printf("Avg time %f\n", time_serial/3);

  printToFile(a, N, "serialsorted");
}

// nilai pangkat 2 terkecil yang lebih besar dari n
int searchNearestN(int n){
  int x = 1;
  while (x < n)
    x *= 2;
  return x;
}

/** procedure rng() : random number generator**/
void rng(int* a, int n, int m) {
    int seed = 13515112; // Ganti dengan NIM anda sebagai seed.
    int max = -9999;
    srand(seed);
    for(long i = 0; i < n; i++) {
        a[i] = (int)rand();
        if (max < a[i]) {
        	max = a[i];
        }
    }
    fillDummy(a, n, m, max+1);
}

// fill array with dummy integer until 2^n
void fillDummy(int * a, int from, int to, int dummy){
	for (long i = from; i < to; ++i) {
		a[i] = dummy;
	}
}

void printToFile(int * a, int n, char * file) {
  FILE * fp;
  char * namefile = (char *) malloc(50*sizeof(char));
  sprintf(namefile, "data/%s.txt", file);
  fp = fopen(namefile, "w+");
  for (int i = 0; i < n; ++i){
   fprintf(fp, "%d\n", a[i]); 
  }  
  fclose(fp);
}

/** procedure  print() : print array elements **/
void print(int * a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}

/** procedure test() : verify sort results **/
void test(int * a, int n) {
  int pass = 1;
  int i;
  for (i = 1; i < n; i++) {
    pass &= (a[i-1] <= a[i]);
  }
  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}

// swap
void exchange(int * a,int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}

// imperative parallel version of bitonic sort 
void impBitonicSort(int * a, int n) {
  int i,j,k;
  // dari 2 sampe 
  for (k=2; k<=n; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=0; i<n; i++) {
        int ij=i^j;
        if ((ij)>i) {
          if ((i&k)==0 && a[i] > a[ij]) 
              exchange(a,i,ij);
          if ((i&k)!=0 && a[i] < a[ij])
              exchange(a,i,ij);
        }
      }
    }
  }
}

// imperative parallel version of bitonic sort 
void impBitonicSortParallel(int * a, int n) {
  int i,j,k;
  // #pragma omp parallel num_threads(thread_count)
  for (k=2; k<=n; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
	  #pragma omp parallel for num_threads(thread_count)
    // #pragma omp for
  	  for (i=0; i<n; i++) {
  	    int ij=i^j;
  	    if ((ij)>i) {
  	      if ((i&k)==0 && a[i] > a[ij]) 
  	          exchange(a,i,ij);
  	      if ((i&k)!=0 && a[i] < a[ij])
  	          exchange(a,i,ij);
        }
      }
    }
  }
}