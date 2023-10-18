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
write_attr(FILE *file, int* attractor, int n_size, int n_degree){

    int color_str_len = 12;

    // color scheme
    char c_0[] = "100 100 100\n";
    char c_1[] = "100 100 255\n";
    char c_2[] = "100 255 100\n";
    char c_3[] = "100 255 255\n";
    char c_4[] = "255 100 100\n";
    char c_5[] = "255 100 255\n";
    char c_6[] = "255 255 100\n";
    char c_7[] = "255 255 255\n";
    char c_8[] = "100 150 250\n";
    char c_9[] = "250 150 100\n";

    char colors[140]; // = (char*) malloc(120*sizeof(char));
    sprintf(colors, "%s%s%s%s%s%s%s%s%s%s", c_0, c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9);
    // printf("%s", colors);

    char *row_str = (char*) malloc(n_size*color_str_len*sizeof(char));

    for ( size_t ix = 0, jx = 0; jx < n_size; ix += color_str_len, ++jx ){
        memcpy( row_str + ix, colors + attractor[jx]*color_str_len, color_str_len);
    }
    
    fwrite( row_str, sizeof(char), n_size*color_str_len, file);
    free(row_str);
    
    return 0;
}


int
main(int argc, char* argv[])
{
    srand(time(NULL));

    // initialize working variables
    int n_size = 1000;                                      // length of row, n.o. pixels in row
    int *attractor = (int*) malloc(n_size*sizeof(int));
    int n_degree = 10;                                      // degree of polynomial
    int max_color_val = 255;

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

    r = write_header(file_attr, n_size, max_color_val);

    r = write_attr(file_attr, attractor, n_size, n_degree);

    fclose(file_attr);
    free(attractor);
    return 0;

}