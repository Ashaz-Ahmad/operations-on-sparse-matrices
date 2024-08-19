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
        if (line[0] != '%') break; 
    }

   
    sscanf(line, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->row_ptr = (int *)calloc(matrix->num_rows + 1, sizeof(int));

    if (matrix->csr_data == NULL || matrix->col_ind == NULL || matrix->row_ptr == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        fclose(file);
        return;
    }

    for (int i = 0; i <= matrix->num_rows; i++) 
    {
        matrix->row_ptr[i] = 0; 
    }

    int row, col;
    double matrix_value; 
    while (fscanf(file, "%d %d %lf", &row, &col, &matrix_value) == 3) 
    {
        row-=1;  
        matrix->row_ptr[row + 1]+=1;
    }

    
    for (int i = 1; i <= matrix->num_rows; i++) 
    {
        matrix->row_ptr[i] += matrix->row_ptr[i - 1];
    }

    fseek(file, 0, SEEK_SET); 

    while (fgets(line, sizeof(line), file)) 
    {
        if (line[0] != '%') break; 
    }

    int *temp_row_ptr = (int *)malloc((matrix->num_rows + 1) * sizeof(int)); 
    if (temp_row_ptr == NULL)
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        fclose(file);
        return;
    }
    
    memcpy(temp_row_ptr, matrix->row_ptr, (matrix->num_rows + 1) * sizeof(int)); 

    while (fscanf(file, "%d %d %lf", &row, &col, &matrix_value) == 3)
    {
        row-=1;  
        col-=1; 
        int index = temp_row_ptr[row]++;
        matrix->csr_data[index] = matrix_value;
        matrix->col_ind[index] = col;
    }

    free(temp_row_ptr);
    fclose(file);
}

CSRMatrix addCSRMatrices(const CSRMatrix *A, const CSRMatrix *B) 
{
    if (A->num_rows != B->num_rows || A->num_cols != B->num_cols)
    {
        fprintf(stderr, "Error: Matrices dimensions do not match for addition.\n");
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
                C.csr_data[col_tracker[col]] += B->csr_data[j]; 
            }
            else
            {
                C.csr_data[counter] = B->csr_data[j];
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

CSRMatrix subtractCSRMatrices(const CSRMatrix *A, const CSRMatrix *B) 
{

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
                C.csr_data[col_tracker[col]] -= B->csr_data[j]; 
            }
            else
            {
                C.csr_data[counter] = -B->csr_data[j]; 
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
 
    CSRMatrix C; 
    C.num_rows = A->num_rows; 
    C.num_cols = B->num_cols;     
    int max_num_values = 160000000; 
    
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

    memset(column_marker, -1, C.num_cols * sizeof(int)); 

    int counter = 0;
    for (int i = 0; i < A->num_rows; i++) 
    {
        C.row_ptr[i] = counter;

        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            int a_col = A->col_ind[j];
            double a_value = A->csr_data[j]; 

            for (int k = B->row_ptr[a_col]; k < B->row_ptr[a_col + 1]; k++) 
            {
                int b_col = B->col_ind[k];
                double b_value = B->csr_data[k];
                if (column_marker[b_col] < C.row_ptr[i])
                {
                    column_marker[b_col] = counter; 
                    C.col_ind[counter] = b_col;
                    C.csr_data[counter] = a_value * b_value;
                    counter+=1;
                }
                else
                {
                    C.csr_data[column_marker[b_col]] += a_value * b_value; 
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

CSRMatrix transposeCSRMatrix(const CSRMatrix *A) 
{
    CSRMatrix At;
    At.num_rows = A->num_cols;
    At.num_cols = A->num_rows; 
    At.num_non_zeros = A->num_non_zeros; 

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
        At.row_ptr[A->col_ind[i] + 1]++; 
    }

    for (int i = 1; i <= At.num_rows; i++) 
    {
        At.row_ptr[i] += At.row_ptr[i - 1]; 
    }

    int *current_pos = (int *)malloc(At.num_rows * sizeof(int)); 
    if (current_pos == NULL) 
    {
        fprintf(stderr, "Failed to allocate memory.\n");
        free(At.row_ptr);
        free(At.csr_data);
        free(At.col_ind);
        exit(EXIT_FAILURE);
    }
    memcpy(current_pos, At.row_ptr, At.num_rows * sizeof(int)); 
    
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