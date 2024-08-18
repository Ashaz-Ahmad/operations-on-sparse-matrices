#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix) 
{
    FILE *file = fopen(filename, "r");
    if (file == NULL) 
    {
        fprintf(stderr, "Failed to open file %s\n", filename);
        return;
    }

    char line[1500];

    while (fgets(line, sizeof(line), file)) 
    {
        if (line[0] != '%') break; //skip lines that are comments (any line starting with %)
    }

    //Read the matrix dimensions and number of non-zero elements
    sscanf(line, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->row_ptr = (int *)calloc(matrix->num_rows + 1, sizeof(int));

    if (matrix->csr_data == NULL || matrix->col_ind == NULL || matrix->row_ptr == NULL) //check if memory allocation was successfull
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        fclose(file);
        return;
    }

    for (int i = 0; i <= matrix->num_rows; i++) 
    {
        matrix->row_ptr[i] = 0; //initialize row_ptr
    }

    int row, col;
    double matrix_value; //initialize row, col, and value with random values for now
    //read the data entries and populate row_ptr
    while (fscanf(file, "%d %d %lf", &row, &col, &matrix_value) == 3) 
    {
        row-=1;  //Convert to 0 based index
        matrix->row_ptr[row + 1]+=1;
    }

    //Accumulate the row pointers
    for (int i = 1; i <= matrix->num_rows; i++) 
    {
        matrix->row_ptr[i] += matrix->row_ptr[i - 1];
    }

    fseek(file, 0, SEEK_SET); //fseek() basically allows us to jump to another place in the file that we want to read
    //in our case, we jump to reading from the beginning of the same file again

    while (fgets(line, sizeof(line), file)) 
    {
        if (line[0] != '%') break; //skip comments again
    }

    int *temp_row_ptr = (int *)malloc((matrix->num_rows + 1) * sizeof(int)); //create temporary row pointerss
    if (temp_row_ptr == NULL)
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        fclose(file);
        return;
    }
    
    memcpy(temp_row_ptr, matrix->row_ptr, (matrix->num_rows + 1) * sizeof(int)); //copy memory from row_ptr to temp_row_ptr

    while (fscanf(file, "%d %d %lf", &row, &col, &matrix_value) == 3) //read the data and assign values for csr_data and col_ind
    {
        row-=1;  
        col-=1;  //Convert both rows and cols to 0 based index
        int index = temp_row_ptr[row]++;
        matrix->csr_data[index] = matrix_value;
        matrix->col_ind[index] = col;
    }

    free(temp_row_ptr);
    fclose(file);
}

CSRMatrix addCSRMatrices(const CSRMatrix *A, const CSRMatrix *B) // This function adds two CSR matrices
{
    if (A->num_rows != B->num_rows || A->num_cols != B->num_cols)
    {
        fprintf(stderr, "Error: Matrices dimensions do not match for addition.\n");
        exit(EXIT_FAILURE);
    }

    CSRMatrix C;              // initialize C
    C.num_rows = A->num_rows; // Sum of matrices will have same dimensions as A or B, doesn't matter which one you use for rows and cols
    C.num_cols = A->num_cols;

    int max_non_zeros = A->num_non_zeros + B->num_non_zeros; // max number of non zero entries in matrix C will be the sum of non zeroes in A + non zeroes in B

    C.row_ptr = (int *)calloc(C.num_rows + 1, sizeof(int));
    C.csr_data = (double *)malloc(max_non_zeros * sizeof(double));
    C.col_ind = (int *)malloc(max_non_zeros * sizeof(int));
    int *col_tracker = (int *)malloc(C.num_cols * sizeof(int));

    if (col_tracker == NULL || C.row_ptr == NULL || C.csr_data == NULL || C.col_ind == NULL)
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(C.csr_data);
        free(C.col_ind);
        free(C.row_ptr); 
        free(col_tracker);
        exit(EXIT_FAILURE);
    }

    memset(col_tracker, -1, C.num_cols * sizeof(int)); // memset is used to fill blocks of memory. In this case, we put -1 in the blocks
    // of memory that are being pointed at by the column_marker pointer. We fill in C.num_cols * sizeof(int) bytes of memory

    // we will use column_marker to more efficiently perform addition below

    int counter = 0;
    for (int i = 0; i < C.num_rows; i++) // iterate over all rows in matrix C (A and B also have the same amount of rows)
    {
        C.row_ptr[i] = counter;

        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) // iterate over all non-zero entries in the i-th row of A
        {
            int col = A->col_ind[j]; // this chunk of code in the for loop copies all non-zero elements from A into C
            C.csr_data[counter] = A->csr_data[j];
            C.col_ind[counter] = col;
            col_tracker[col] = counter; // before we set all entries in column_marker to -1.
            // Now we change them from -1 to the value of the counter for any column in which there is a value of A copied to C
            counter+=1;
        }

        for (int j = B->row_ptr[i]; j < B->row_ptr[i + 1]; j++) // iterate over all non-zero entries in the i-th row of B
        {
            int col = B->col_ind[j];
            if (col_tracker[col] != -1) // any column where column_marker[col] != -1 is a column where a value from A was copied into C
            {
                C.csr_data[col_tracker[col]] += B->csr_data[j]; // so, we perform addition on these entries
            }
            else
            {
                C.csr_data[counter] = B->csr_data[j]; // this chunk of code in the else block creates new entries in C and copies the values from B if needed
                C.col_ind[counter] = col;
                col_tracker[col] = counter;
                counter+=1;
            }
        }

        // the 2 for loops below are used to reset the entries in column_marker to -1 since they will vary row by row
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            int col = A->col_ind[j];
            col_tracker[col] = -1;
        }

        for (int j = B->row_ptr[i]; j < B->row_ptr[i + 1]; j++)
        {
            int col = B->col_ind[j];
            col_tracker[col] = -1;
        }
    }

    C.row_ptr[C.num_rows] = counter;
 
    // Allocate new arrays that will be used to filter out any zeroes that got added to the values during computation
    double *correct_csr_data = (double *)malloc(counter * sizeof(double));
    int *correct_col_ind = (int *)malloc(counter * sizeof(int));
    int *correct_row_ptr = (int *)calloc(C.num_rows + 1, sizeof(int));
    if (correct_csr_data == NULL || correct_col_ind == NULL || correct_row_ptr == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(C.row_ptr);
        free(C.csr_data);
        free(C.col_ind);
        free(col_tracker);
        exit(EXIT_FAILURE);
    }
 
    int correct_num_values = 0;
    for (int i = 0; i < C.num_rows; i++) 
    {
        correct_row_ptr[i] = correct_num_values;
 
        for (int j = C.row_ptr[i]; j < C.row_ptr[i + 1]; j++) 
        {
            if (C.csr_data[j] != 0) 
            {
                correct_csr_data[correct_num_values] = C.csr_data[j];
                correct_col_ind[correct_num_values] = C.col_ind[j];
                correct_num_values+=1;
            }
        }
    }
 
    correct_row_ptr[C.num_rows] = correct_num_values;
 
    free(C.csr_data);
    free(C.col_ind);
    free(C.row_ptr);
    free(col_tracker);
 
    C.csr_data = correct_csr_data;
    C.col_ind = correct_col_ind;
    C.row_ptr = correct_row_ptr;
    C.num_non_zeros = correct_num_values;
 
    return C;
}

CSRMatrix subtractCSRMatrices(const CSRMatrix *A, const CSRMatrix *B) // This function subtracts two CSR matrices
{
    // The logic and implementation of this function is identical to the addCSRMatrices function, except there is subtraction in this function
    // I will only highlight the changes from addCSRMatrices.

    if (A->num_rows != B->num_rows || A->num_cols != B->num_cols)
    {
        fprintf(stderr, "Error: Matrices dimensions do not match for subtraction.\n");
        exit(EXIT_FAILURE);
    }

    CSRMatrix C;
    C.num_rows = A->num_rows;
    C.num_cols = A->num_cols;

    int max_non_zeros = A->num_non_zeros + B->num_non_zeros;

    C.row_ptr = (int *)calloc(C.num_rows + 1, sizeof(int));
    C.csr_data = (double *)malloc(max_non_zeros * sizeof(double));
    C.col_ind = (int *)malloc(max_non_zeros * sizeof(int));
    int *col_tracker = (int *)malloc(C.num_cols * sizeof(int));

    if (col_tracker == NULL || C.row_ptr == NULL || C.csr_data == NULL || C.col_ind == NULL)
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(C.csr_data);
        free(C.col_ind);
        free(C.row_ptr); 
        free(col_tracker);
        exit(EXIT_FAILURE);
    }
    memset(col_tracker, -1, C.num_cols * sizeof(int));

    int counter = 0;
    for (int i = 0; i < C.num_rows; i++)
    {
        C.row_ptr[i] = counter;

        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            int col = A->col_ind[j];
            C.csr_data[counter] = A->csr_data[j];
            C.col_ind[counter] = col;
            col_tracker[col] = counter;
            counter+=1;
        }

        for (int j = B->row_ptr[i]; j < B->row_ptr[i + 1]; j++)
        {
            int col = B->col_ind[j];
            if (col_tracker[col] != -1)
            {
                C.csr_data[col_tracker[col]] -= B->csr_data[j]; // Same as before but subtract entries in matrix C from those in B (instead of add)
            }
            else
            {
                C.csr_data[counter] = -B->csr_data[j]; // here we copy over the negative values of entries in B to C instead of positive
                // (since if there is no entry in the C data, then is must be 0 and we would have 0 - <Entry in B> = -<Entry in B>)
                C.col_ind[counter] = col;
                col_tracker[col] = counter;
                counter+=1;
            }
        }

        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            int col = A->col_ind[j];
            col_tracker[col] = -1;
        }
        for (int j = B->row_ptr[i]; j < B->row_ptr[i + 1]; j++)
        {
            int col = B->col_ind[j];
            col_tracker[col] = -1;
        }
    }

    C.row_ptr[C.num_rows] = counter;
 
    double *correct_csr_data = (double *)malloc(counter * sizeof(double));
    int *correct_col_ind = (int *)malloc(counter * sizeof(int));
    int *correct_row_ptr = (int *)calloc(C.num_rows + 1, sizeof(int));
    if (correct_csr_data == NULL || correct_col_ind == NULL || correct_row_ptr == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(C.row_ptr);
        free(C.csr_data);
        free(C.col_ind);
        free(col_tracker);
        exit(EXIT_FAILURE);
    }
 
    int correct_num_values = 0;
    for (int i = 0; i < C.num_rows; i++) 
    {
        correct_row_ptr[i] = correct_num_values;
 
        for (int j = C.row_ptr[i]; j < C.row_ptr[i + 1]; j++) 
        {
            if (C.csr_data[j] != 0) 
            {
                correct_csr_data[correct_num_values] = C.csr_data[j];
                correct_col_ind[correct_num_values] = C.col_ind[j];
                correct_num_values++;
            }
        }
    }
 
    correct_row_ptr[C.num_rows] = correct_num_values;
 
    free(C.csr_data);
    free(C.col_ind);
    free(C.row_ptr);
    free(col_tracker);
 
    C.csr_data = correct_csr_data;
    C.col_ind = correct_col_ind;
    C.row_ptr = correct_row_ptr;
    C.num_non_zeros = correct_num_values;
 
    return C;
}

CSRMatrix multiplyCSRMatrices(const CSRMatrix *A, const CSRMatrix *B)
{
    if (A->num_cols != B->num_rows)
    {
        fprintf(stderr, "Error: Matrices dimensions do not match for multiplication.\n");
        exit(EXIT_FAILURE);
    }
 
    CSRMatrix C; //initialize C
    C.num_rows = A->num_rows; 
    C.num_cols = B->num_cols; //When you multiply matrices A and B, matrix C (the product) will have the number of rows in A and the number of columns in B
    
    int max_num_values = 160000000; // A->num_non_zeros * B->num_non_zeros; this is the maximum amount of non-zero entries the product matrix C can have
                                   //The line above was commented and replaced because it was causing the program to fail memory allocation because the arrays would be way too big
                                   //Using 160million is big enough and does not cause memory allocation to fail
    
    C.row_ptr = (int *)calloc(C.num_rows + 1, sizeof(int)); 
    C.csr_data = (double *)malloc(max_num_values * sizeof(double));
    C.col_ind = (int *)malloc(max_num_values * sizeof(int));
    int *column_marker = (int *)malloc(C.num_cols * sizeof(int));

    if (column_marker == NULL || C.row_ptr == NULL || C.csr_data == NULL || C.col_ind == NULL)
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(C.row_ptr);
        free(C.csr_data);
        free(C.col_ind);
        free(column_marker);
        exit(EXIT_FAILURE);
    }

    memset(column_marker, -1, C.num_cols * sizeof(int)); //Same concept as before in add and subtract functions using memset
 
    // again, column_marker will be used to help us perform multiplication efficiently 

    int counter = 0;
    for (int i = 0; i < A->num_rows; i++) // Iterate over the rows of matrix A
    {
        C.row_ptr[i] = counter;

        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) // iterate over all non-zero elements in the rows of A
        {
            int a_col = A->col_ind[j];
            double a_value = A->csr_data[j]; // store the values and columns where those values occur

            for (int k = B->row_ptr[a_col]; k < B->row_ptr[a_col + 1]; k++) // iterate over the non-zero elements in B which are in a column corresponding to the values in the rows of A
            {
                int b_col = B->col_ind[k];
                double b_value = B->csr_data[k];
                if (column_marker[b_col] < C.row_ptr[i]) // check if the current column of B is already present in C
                {
                    column_marker[b_col] = counter; // This is executed if the column of B is not present in C
                    C.col_ind[counter] = b_col;
                    C.csr_data[counter] = a_value * b_value;
                    counter+=1; // it will add the entry to C and the product of the values
                }
                else
                {
                    C.csr_data[column_marker[b_col]] += a_value * b_value; // this is executed if the column of B is already present in C
                }
            }
        }
        for (int m = A->row_ptr[i]; m < A->row_ptr[i + 1]; m++) 
        {
            int a_col = A->col_ind[m];
 
            for (int n = B->row_ptr[a_col]; n < B->row_ptr[a_col + 1]; n++) 
            {
                int b_col = B->col_ind[n];
                column_marker[b_col] = -1;
            }
        }
    }

    C.row_ptr[C.num_rows] = counter;

    //This chunck of code has the same concept in using a filter like in the add and subtract functions

    double *correct_csr_data = (double *)malloc(counter * sizeof(double));
    int *correct_col_ind = (int *)malloc(counter * sizeof(int));
    int *correct_row_ptr = (int *)calloc(C.num_rows + 1, sizeof(int));

    if (correct_csr_data == NULL || correct_col_ind == NULL || correct_row_ptr == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(C.row_ptr);
        free(C.csr_data);
        free(C.col_ind);
        free(column_marker);
        exit(EXIT_FAILURE);
    }
    
    int correct_num_values = 0; 
    for (int i = 0; i < C.num_rows; i++)
    {
        correct_row_ptr[i] = correct_num_values; 
        for (int j = C.row_ptr[i]; j < C.row_ptr[i + 1]; j++)  
        {
            if (C.csr_data[j] != 0) 
            {
                correct_csr_data[correct_num_values] = C.csr_data[j]; 
                correct_col_ind[correct_num_values] = C.col_ind[j]; 
                correct_num_values++; 
            }
        }
    }
 
    correct_row_ptr[C.num_rows] = correct_num_values; 

    free(C.csr_data);
    free(C.col_ind);
    free(C.row_ptr);
    free(column_marker);
 
    C.csr_data = correct_csr_data; 
    C.col_ind = correct_col_ind; 
    C.row_ptr = correct_row_ptr; 
    C.num_non_zeros = correct_num_values; 
 
    return C; 
}

CSRMatrix transposeCSRMatrix(const CSRMatrix *A) //This function takes the transpose of a matrix
{
    CSRMatrix At;
    At.num_rows = A->num_cols;
    At.num_cols = A->num_rows; //When you take the transpose of a matrix, the rows and columns are switched
    At.num_non_zeros = A->num_non_zeros; //Amount of non-zero entries stay the same after you take the transpose

    At.row_ptr = (int *)calloc(At.num_rows + 1, sizeof(int));
    At.csr_data = (double *)malloc(At.num_non_zeros * sizeof(double));
    At.col_ind = (int *)malloc(At.num_non_zeros * sizeof(int));

    if (At.row_ptr == NULL || At.csr_data == NULL || At.col_ind == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(At.row_ptr);
        free(At.col_ind);
        free(At.csr_data);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < A->num_non_zeros; i++) 
    {
        At.row_ptr[A->col_ind[i] + 1]++; //count the number of entries in each column of the original matrix A
    }

    for (int i = 1; i <= At.num_rows; i++) 
    {
        At.row_ptr[i] += At.row_ptr[i - 1]; //Cumulative sum to convert counts to row_ptr
    }

    int *current_pos = (int *)malloc(At.num_rows * sizeof(int)); //Create a temporary array to track the current position in each row of At
    if (current_pos == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(At.row_ptr);
        free(At.csr_data);
        free(At.col_ind);
        exit(EXIT_FAILURE);
    }
    memcpy(current_pos, At.row_ptr, At.num_rows * sizeof(int)); //memcpy() is used to copy a block of memory from one location to another
    //in this case, we are copying memory from At.row_ptr to current_pos and the size of memory copied is At.num_rows * sizeof(int) bytes

    // Fill the csr_data and col_ind arrays for At
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) 
        {
            int col = A->col_ind[j];
            int destination = current_pos[col];

            At.csr_data[destination] = A->csr_data[j];
            At.col_ind[destination] = i;

            current_pos[col]++;
        }
    }

    free(current_pos);
    return At;
}

void printCSRMatrix(const CSRMatrix *matrix)
{
    printf("Number of non-zeros: %d\n", matrix->num_non_zeros); 

    printf("Row Pointer: ");
    for (int i = 0; i <= matrix->num_rows; i++) 
    {
        printf("%d ", matrix->row_ptr[i]); 
    }
    printf("\n");

    printf("Column Index: ");
    for (int i = 0; i < matrix->num_non_zeros; i++) 
    {
        printf("%d ", matrix->col_ind[i]); 
    }
    printf("\n");

    printf("Values: "); 
    for (int i = 0; i < matrix->num_non_zeros; i++) 
    {
        printf("%0.4lf ", matrix->csr_data[i]); 
    }
    printf("\n");
}

void freeCSRMatrix(CSRMatrix *matrix)
{
    free(matrix->csr_data);
    free(matrix->col_ind);
    free(matrix->row_ptr);
    
    matrix->num_non_zeros = 0;
    matrix->num_rows = 0;
    matrix->num_cols = 0;
}