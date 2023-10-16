#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main(int argc, char *argv[]) {
   int n_size = 0, n_threads = 0;

   for (int i = 1; i < argc; i++) {		
      if (strncmp(argv[i], "-l", 2) == 0) 
         n_size = atoi(argv[i] + 2);
		
      else if (strncmp(argv[i], "-t", 2) == 0) 
         n_threads = atoi(argv[i] + 2);      
   }
	
	assert(n_size != 0 && n_threads != 0);
	
	
   return 0;
}