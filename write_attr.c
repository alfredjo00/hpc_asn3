# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

int
write_header(FILE *file, int n_size, int max_color_val){

    char str[20];
    fprintf(file, "%s\n%i %i\n%i\n", "P3", n_size, n_size, max_color_val);
    

    return 0;
}

int
write_row(FILE *file, int* attractor, int n_size, int n_degree){

    int str_len = 6;

    char *row_str = (char*) malloc(n_size*str_len*sizeof(char));

    for (int ix = 0, jx = 0; ix < n_size*str_len; ix+=str_len, ++jx){
        row_str[ix]         = attractor[jx] + 48;
        row_str[ix+1]       = 32;
        row_str[ix+2]       = attractor[jx] + 48;
        row_str[ix+3]       = 32;
        row_str[ix+4]       = attractor[jx] + 48;
        row_str[ix+5]       = 10;
    }
    row_str[n_size*str_len] = 0; // null termination


    fwrite( row_str, sizeof(char), n_size*str_len, file);
    free(row_str);
    
    return 0;
}


int
main(int argc, char* argv[])
{
    srand(time(NULL));

    // initialize working variables
    int n_size = 1000;       // length of row, n.o. pixels in row
    // int str_len = 6;    // length of color string
    int *attractor = (int*) malloc(n_size*sizeof(int)); // [n_size];   // array of attractor indicies
    int n_degree = 10;         // degree of polynomial

    // create filename
    char degree[] = "2";
    char f_name[] = "newton_attractors_x";
    strcat(f_name, degree);
    strcat(f_name, ".ppm");
    //

    // open file
    FILE *file_attr = fopen( f_name, "w" );
    if (file_attr == NULL){
        printf("Error opening file\n");
        return -1;
    }
    //

    // create temporary attractors
    for ( int ix = 0; ix < n_size; ++ix)
        attractor[ix] = (int) rand() % 10;

    int r;

    r = write_header(file_attr, n_size, 9);

    r = write_row(file_attr, attractor, n_size, n_degree);

    fclose(file_attr);

    return 0;

}