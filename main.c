#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define NMAX 8200


// ---Prototype functions defined here---
void setup();
void menu();
void choose();
void switch_menu(int n, int n_params, float* x, float* y, float* d_y);
void parameters (int n, int n_params, float* x, float* y, float* d_y);
void chi_squared(int n, int n_params, float* x, float* y, float* d_y);
void multiply_matrix(float mat_a[][NMAX], float mat_b[][NMAX], float mat_c[][NMAX], int rows1, int cols1, int rows2, int cols2);
void transpose(int mrows, int ncols, float matrix[][NMAX], float matrix_T[][NMAX]);
void inverse(int n, float a[][NMAX], float inv_a[][NMAX]);
void print_mat(float matrix[][NMAX], int rows, int cols);
void print_data(int n, float* x, float* y, float* d_y);
void linear_menu(int n, int n_params, float* x, float F[NMAX][NMAX]);


int main()
{
    setup();

    FILE *data;
    data = fopen("data.txt","w");
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

    int n_params;
    printf("Please enter the order of the linear equation you would like to fit for the data file:\t");
    scanf("%i", &n_params);

//------------------------------------------All of the functions are referenced in this menu-------------------------------------
    switch_menu(n, n_params, x, y, d_y);

    printf("\n");
    fclose(data);
    free(x);
    free(y);
    free(d_y);

    return(0);
}

void switch_menu(int n, int n_params, float* x, float* y, float* d_y)
{
    int ex = 1;
    char val[NMAX];
    int option;
    float a_arr[NMAX][NMAX];

    while (ex > 0){

        menu(); // Just text


        // sanity check to make sure input is integer.
        scanf("%s", val);
        option = atoi(val);

        // switch command to create a menu.
        switch(option){
            case 1:
                print_data(n, x, y, d_y);

                choose(); // Just text
                scanf("%i", &ex); // If number contradicts while argument, breaks menu. any input less than 0 will work to break.
                break;

            case 2:
                // function to calculate and print the parameters.
                parameters(n, n_params, x, y, d_y);

                choose();
                scanf("%i", &ex);
                break;

            case 3:
                // function to calculate and print the chi squared value.
                chi_squared(n, n_params, x, y, d_y);

                choose();
                scanf("%i", &ex);
                break;

            case -1:
                printf("\nAre you sure?\n");
                choose();
                scanf("%d", &ex);
                break;

            //default if input does not match any option. if user enters a string, atoi() of a string will return 0 which is not an option leading to here.
            default:
                printf("\nNo Option Chosen. %s is not an option\n", val);
                printf("Enter 1 to try again\n");
                printf("Enter -1 to close\n");
                scanf("%d", &ex);
        }

    }
}

void setup()
{
    printf("Hello, welcome to the fitting programm.\n");
    printf("Please make sure you're file is in the format x y d_y in column format, "
           "with each variable separated with a single white space, and that it is named DataExperiment.\n");
}

void menu()
{
    printf("\n---Menu---\n"
    "- 1. Read data.\n"
    "- 2. Calculate parameters for a linear function and associated standard deviation.\n"
    "- 3. Calculate chi squared and associated standard deviation.\n"
    "- -1. Exit.\n"
    "- Please enter option - ");
}

void choose(){
    printf("\nEnter 1 to return to menu\n");
    printf("Enter -1 to close\n");
}

void print_data(int n, float* x, float* y, float* d_y){ //Prints the data if the user chooses to see it again in order to fit the equations.

    printf("x\ty\td_y\n");
    printf("----------------\n");
    for (int i = 0; i < n; i++){
        printf("%.3f\t%.3f\t%.3f\n", x[i], y[i], d_y[i]);
    }
    printf("----------------\n");
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

void transpose(int mrows, int ncols, float matrix[][NMAX], float matrix_T[][NMAX]) // takes the size of matrix and stores data in matrix_T by looping through and assigning to opposite i and j.
{
    for (int i = 0; i < mrows; i++){
        for (int j = 0; j < ncols ; j++){
            matrix_T[j][i] = matrix[i][j];
        }
    }
}


void inverse(int n, float a[][NMAX], float inv_a[][NMAX]) // Gauss-Jordan elimination. Creates an identity matrix augment [ A | I ]. Then [A-1 A | A-1 I ] by definition returns [ I | A-1 ]
{                                                         // by row reduction, i.e. by making the A matrix into an identity matrix, the operations in order to do this are necessarily equivalent
                                                          // to the inverse of A.
    int i,j,k;
    float ratio;

    if(n < 2){
        printf("Error! Not enough dimensions\n");
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
            printf("%.1f, %i", a[i][i], i);
            printf("\nZero Error in Matrix - Cannot Inverse\n");
            exit(0);
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
void print_mat(float matrix[][NMAX], int rows, int cols)
{
    printf("\nMat\n");
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            printf("%.3f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void parameters (int n, int n_params, float* x, float* y, float* d_y) // This combines all the functions in the order needed to produce the parameters.
{
    float a_arr[NMAX][NMAX];                                                  // Performs the matrix operations: (Ft * aSigma * F)^-1 * Ft * aSigma * y = a
    float F_TxSigma[NMAX][NMAX];
    float inv_Cov[NMAX][NMAX];
    float Cov[NMAX][NMAX];
    float CovxF_T[NMAX][NMAX];
    float proto_arr[NMAX][NMAX];
    float sig[NMAX][NMAX];
    float err[NMAX];

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if(i==j){
            sig[i][j] = d_y[i];
            }
            else{
                sig[i][j] = 0;
            }
        }
    }

    float asig[NMAX][NMAX];
    inverse(n, sig, asig);

    float F[NMAX][NMAX];
    float F_T[NMAX][NMAX];

    // Can fit for linear equations of up to order x^n_parameters-1 (Where no of parameters is limited to 5, can be set to more, just a few lines of code)
    linear_menu(n, n_params, x, F);

    float y_arr[NMAX][NMAX];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < 1; j++){
            y_arr[i][j] = y[i];
        }
    }
    transpose(n, n_params, F, F_T);

    // (Ft * aSigma * F)

    multiply_matrix(F_T, asig, F_TxSigma, n_params, n, n, n);
    multiply_matrix(F_TxSigma, F, inv_Cov, n_params, n, n, n_params);
    inverse(n_params, inv_Cov, Cov);
    printf("\n--- Covariance Matrix ---\n");
    print_mat(Cov, n_params, n_params);
    // At this point you have the covariance matrix

    //(Ft * aSigma * F)^-1 * Ft
    multiply_matrix(Cov, F_T, CovxF_T, n_params, n_params, n_params, n);
    multiply_matrix(CovxF_T, asig, proto_arr, n_params, n, n, n);

    //(Ft * aSigma * F)^-1 * Ft * aSigma * y
    multiply_matrix(proto_arr, y_arr, a_arr, n_params, n, n, 1);

    printf("\n---Parameter Fitted Linear Equation---\n\n");
    for (int i = 0; i < n_params; i++){
        for (int j = 0; j < 1; j++){
            printf("a[%i] = %.5f\t", i, a_arr[i][j]);
        }
    }

    // Calculates the error of the fitted parameters using the square-rooted diagonal elements of the covariance matrix.

    for (int i = 0; i < n_params; i++){
        for (int j = 0; j < n_params; j++){
            if(i==j){
                err[i] = Cov[i][j];
            }
        }
    }
    printf("\n");
    for (int i = 0; i < n_params; i++){
        printf("\nError in a[%i] = %.3f\n", i, sqrt(err[i]));
    }

    // ----------------------------Printing to file----------------------------

    FILE *clear;
    clear = fopen("DataFittingProject.txt", "w");
    fclose(clear);

    FILE *results;
    results = fopen("DataFittingProject.txt","a");
    if (results==NULL){
        printf("Error opening results file!\n");
        exit (0);
    }
    fprintf(results, "\n---Parameter Fitted Linear Equation---\n\n");
    for (int i = 0; i < n_params; i++){
        for (int j = 0; j < 1; j++){
            fprintf(results, "a[%i] = %.5f\t", i, a_arr[i][j]);
        }
    }
    for (int i = 0; i < n_params; i++){
        fprintf(results, "\nError in a[%i] = %.5f ", i, sqrt(err[i]));
    }
    fprintf(results, "\n");
}

//--------------------------------------------------------------------CHI SQUARED FUNCTION-------------------------------------

void chi_squared(int n, int n_params, float* x, float* y, float* d_y)
{
    float F[NMAX][NMAX];
    float Fxa[NMAX][NMAX];
    float input_arr[NMAX][NMAX];
    float H[NMAX];

    // Can fit for linear equations of up to order x^n_parameters-1 (Where no of parameters is limited to 5, can be set to more, just a few lines of code) and sin/cos and exp graphs.

    linear_menu(n , n_params, x, F);

    //--------------------------------------------END OF INITIALISING THE MATRICES------------------------------------------------

    printf("\nPlease input guess parameters:\n");
    for (int i = 0; i < n_params; i++){
        printf("Enter a[%i]: ", i);
        scanf("%f", &input_arr[i][0]);
    }
    multiply_matrix(F, input_arr, Fxa, n, n_params, n_params, 1);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < 1; j++){
            H[i] = Fxa[i][0];
        }
    }

    float X2 = 0.0;
    for (int j = 0; j < n; j++){
        X2 = X2 + ( ((y[j] - H[j]) * (y[j] - H[j])) / (d_y[j] * d_y[j]) );
    }

    // Uncertainty in chi squared.
    float d_chi, d_chisqr, chi;
    d_chi = 2 * ( n - n_params );
    d_chisqr = sqrt(d_chi);
    chi = n - n_params;

    printf("\nChi Squared = %.5f\n"
           "\nExpected Chi Squared = %.5f +/- %.5f\n", X2, chi, d_chisqr);

    // -----------Printing to file--------------
    FILE *results;
    results = fopen("DataFittingProject.txt","a");
    if (results==NULL){
        printf("Error opening results file!\n");
        exit (0);
    }
    fprintf(results, "\nChi Squared = %.5f\n"
                     "\nExpected Chi Squared = %.5f +/- %.5f\n", X2, chi, d_chisqr);

    fclose(results);
}

void linear_menu(int n, int n_params, float* x, float F[NMAX][NMAX])
{
    char val[NMAX];
    int option;

    printf("\n\n---Choose your linear equation---\n\n"
    "- 1. Power series of x. a[0] + a[1]x + a[2]x^2 + a[n]x^n\n"
    "- 2. Linear equation of e^(x^n). a[0] + a[1]*e^x + a[n]*e^x^n\n"
    "- 3. Cosine function. a[0] + a[1]*cos(x)\n"
    "- 4. Sine function. a[0] + a[1]*sin(x)\n"
    "- 0. Exit.\n"
    "- Please enter option - ");

    scanf("%s", val);
    option = atoi(val);

    switch(option){
        case 1:
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n_params; j++){
                    F[i][0] = 1;
                    F[i][1] = x[i];
                    F[i][2] = x[i]*x[i];
                    F[i][3] = x[i]*x[i]*x[i];
                    F[i][4] = x[i]*x[i]*x[i]*x[i];
                }
            }
            break;

        case 2:
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n_params; j++){
                    F[i][0] = 1;
                    F[i][1] = exp(x[i]);
                    F[i][2] = exp(x[i]*x[i]);
                    F[i][3] = exp(x[i]*x[i]*x[i]);
                    F[i][4] = exp(x[i]*x[i]*x[i]*x[i]);
                }
            }
            break;

        case 3:
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n_params; j++){
                    F[i][0] = 1;
                    F[i][1] = cosf(x[i]);
                }
            }
            break;

        case 4:
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n_params; j++){
                    F[i][0] = 1;
                    F[i][1] = sinf(x[i]);
                }
            }
            break;

        default:
            printf("No option chosen. Start again.\n\n");
            break;
            }
}
