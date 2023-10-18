# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <assert.h>

int
write_header(FILE *file, int n_size, int max_color_val){

    fprintf(file, "%s\n%i %i\n%i\n", "P3", n_size, n_size, max_color_val);
    
    return 0;
}

int
write_conv(FILE *file, int* convergence, int n_size){

    int color_str_len = 20;
    int row_str_len_sum = 0;
    int offset = 0;
    char *row_str = (char*) malloc(n_size*color_str_len*sizeof(char));
    char colors[10000];
    colors[0] = '\0';

    for (int ix = 0; ix <= 254; ix += 2) {
        // Create a line with the current values
        char line[14];  // Make sure the buffer is large enough
        snprintf(line, sizeof(line), "%i %i %i\n", ix, ix, ix);

        // Append the line to the result string
        strcat(colors, line);
    }

    // for ( int ix = 1200; ix < 1300; ++ix ){
    //     if ( colors[ix] == NULL ){
    //         printf("colors[%i] == NULL\n", ix);
    //     }
    // }
    // printf("%i\n", colors[1370]);

    for ( int ix = 0, jx = 0; jx < n_size; ix += color_str_len, ++jx ){
        if ( convergence[jx] <= 5 ){
            color_str_len = 6;
            row_str_len_sum += color_str_len;
            offset = 0;
        }
        else if ( convergence[jx] <= 50 ){
            color_str_len = 9;
            row_str_len_sum += color_str_len;
            offset = -3*5;
        }
        else {
            color_str_len = 12;
            row_str_len_sum += color_str_len;
            offset = -6*5 - 3*45;
        }

        memcpy( row_str + ix, colors + offset + color_str_len*(convergence[jx]-1), color_str_len);

        // if ( offset + color_str_len*(convergence[jx]-1) < 0 ) {
        //     printf("ix = %i\n", ix);
        //     printf("jx = %i\n", jx);
        //     printf("colors[%i] = %c\n", offset + color_str_len*(convergence[jx]-1), colors[offset + color_str_len*(convergence[jx]-1)]);
        //     printf("color_str_len = %i\n", color_str_len);
        //     printf("offset = %i\n", offset);
        //     printf("convergence[%i] = %i\n", jx, convergence[jx]);
        // }
        
    }

    fwrite( row_str, sizeof(char), row_str_len_sum, file);
    free(row_str);
    
    return 0;
}


int
main(int argc, char* argv[])
{
    srand(time(NULL));

    // initialize working variables
    int n_size = 1000;                                      // length of row, n.o. pixels in row
    int *convergence = (int*) malloc(n_size*sizeof(int));   
    int max_color_val = 255;

    // create filename
    char degree[] = "2";
    char f_name[] = "newton_convergence_x";
    strcat(f_name, degree);
    strcat(f_name, ".ppm");
    //

    // open file
    FILE *file_conv = fopen( f_name, "w" );
    if (file_conv == NULL){
        printf("Error opening file\n");
        return -1;
    }
    //

    // create temporary attractors
    for ( int ix = 0; ix < n_size; ++ix){
        convergence[ix] = (int) rand() % 128 + 1;
    }
    // for ( int ix = 0; ix < 1000; ++ix )
    //     printf("%i\n", convergence[ix]);

    int r;
    r = write_header(file_conv, n_size, max_color_val);
    r = write_conv(file_conv, convergence, n_size);

    fclose(file_conv);
    free(convergence);

    return 0;

}