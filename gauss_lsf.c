#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define NMAX 8200
#define NDIM 3


// ---Prototype functions defined here---
void multiply_matrix(float mat_a[][NMAX], float mat_b[][NMAX], float mat_c[][NMAX], int rows1, int cols1, int rows2, int cols2);
void inverse(int n, float a[][NMAX], float inv_a[][NMAX]);
void gauss_fit(int n, float* x, float* y, float* a);
void print_mat(float array[][NMAX], int rows, int cols);


int main()
{
    FILE *data;
    data = fopen("data.txt","r");
    if (data==NULL){
        printf("Error opening data file!\n");
        exit (0);
    }

    int n = 0;
    char line[NMAX];
    while (fgets(line, NMAX, data))
    {
        n++;
    }
    fseek(data, 0, SEEK_SET);

    // Because the code reads the number of data points first, it can safely allocate the necessary memory before it scans in all the data in.

    float* x = (float*)calloc(n, sizeof(float));
    float* y = (float*)calloc(n, sizeof(float));
    float* d_y = (float*)calloc(n, sizeof(float));

    printf("Number of data points in file = %i \n", n);

    int n_inputs = 0;
    // Data scanned into 1D arrays
    for (int i = 0; i < n; i++){
        n_inputs = fscanf(data, "%f %f %f", &x[i], &y[i], &d_y[i]);
        if (n_inputs != 3){
            printf("Error in file---Not enough entries in file\n\n");
            exit(0);
        }
    }
//------------------------------------------End of File Reading/ set up-------------------------------------------------------------
    float A[NDIM];
    gauss_fit(n, x, y, A);

    for (int i = 0; i < NDIM; i++){
        printf("%.2f\n", A[i]);
    }

    printf("\n--End--");
    fclose(data);
    free(x);
    free(y);
    free(d_y);

    return 0;
}

// double * function(){
//     static double arr[];
//     return arr;
// }

void gauss_fit(int n, float* x, float* y, float a[NDIM])
{
    float X[NMAX][NMAX];
    float inv_X[NMAX][NMAX];
    float Y[NMAX][NMAX];
    float A[NMAX][NMAX];
    float x_sum = 0;
    float x_sqr = 0;
    float x_cb = 0;
    float x_qd = 0;

    for (int i=0; i<n; i++){
        Y[i][0] = y[i];
        x_sum = x_sum + x[i];
        x_sqr = x_sqr + x[i]*x[i];
        x_cb = x_cb + x[i]*x[i]*x[i];
        x_qd = x_qd + x[i]*x[i]*x[i]*x[i];
        printf("%.2f", Y[i][0]);
    }

    X[0][0] = n; X[1][0] = x_sum; X[2][0] = x_sqr;
    X[0][1] = x_sum; X[1][1] = x_sqr; X[2][1] = x_cb;
    X[0][2] = x_sqr; X[1][2] = x_cb; X[2][2] = x_qd;

    print_mat(X, NDIM, NDIM);

    inverse(NDIM, X, inv_X);
    multiply_matrix(inv_X, Y, A, NDIM, NDIM, NDIM, 1);

    for (int i=0; i<NDIM; i++){
        a[i] = A[i][0];
    }
}

void multiply_matrix(float mat_a[][NMAX], float mat_b[][NMAX], float mat_product[][NMAX], int rows1, int cols1, int rows2, int cols2) // argument has to take 3 matrices, multiplied and product, and sizes
{                                                                                                                                     // the product matrix is where the resultant data is stored and then used
    int i,j,k;

    for (i = 0; i < rows1; ++i) {
        for (j = 0; j < cols2; ++j) {
                mat_product[i][j] = 0;  // initialises to zero.
            for (k = 0; k < cols1; ++k) {
                mat_product[i][j] += mat_a[i][k] * mat_b[k][j];
            }
        }
    }
}

void inverse(int n, float a[][NMAX], float inv_a[][NMAX]) // Gauss-Jordan elimination. Creates an identity matrix augment [ A | I ]. Then [A-1 A | A-1 I ] by definition returns [ I | A-1 ]
{                                                         // by row reduction, i.e. by making the A matrix into an identity matrix, the operations in order to do this are necessarily equivalent
                                                          // to the inverse of A.
    int i,j,k;
    float ratio;

    if(n < 2){
        printf("Error! Not enough dimensions to inverse\n");
        exit (0);
    }

    // Augment Identity Matrix of Order n - Creating Augmented matrix by looping through k<2*n.
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i==j){
            a[i][j+n] = 1;
            }
            else {
                a[i][j+n] = 0;
            }
        }
    }
    // Gauss-Jordan Elimination
    for (i = 0; i < n; i++){
        if (a[i][i] == 0.0){
            printf("\nZero Error in Matrix - Cannot inverse\n");
            exit (0);
        }
        for (j = 0; j < n; j++){
            if (i!=j){
                ratio = a[j][i] / a[i][i];
                for (k = 0; k < 2*n; k++){
                    a[j][k] = a[j][k] - ratio * a[i][k];
                }
            }
        }
    }
    // Row Operation to Make the A matrix Diagonal and equivalent to identity matrix.
    for (i = 0; i < n; i++){
        for (j = n; j < 2*n; j++){
            a[i][j] = a[i][j] / a[i][i];
        }
    }
    //
    for (i = 0; i < n; i++){
        for (j = n; j < 2*n; j++){
            inv_a[i][j-n] = a[i][j];
        }
    }
}

void print_mat(float array[][NMAX], int rows, int cols)
{
    int i;
    int j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            printf("%.2f ", array[i][j]);
        }
        printf("\n");
    }
}