#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct 
{
    double *csr_data;   
    int *col_ind;       
    int *row_ptr;       
    int num_non_zeros;  
    int num_rows;       
    int num_cols;       
} CSRMatrix;


void ReadMMtoCSR(const char *filename, CSRMatrix *matrix);

CSRMatrix addCSRMatrices(const CSRMatrix *A, const CSRMatrix *B);
CSRMatrix subtractCSRMatrices(const CSRMatrix *A, const CSRMatrix *B);
CSRMatrix multiplyCSRMatrices(const CSRMatrix *A, const CSRMatrix *B);
CSRMatrix transposeCSRMatrix(const CSRMatrix *A);
void printCSRMatrix(const CSRMatrix *matrix);
void freeCSRMatrix(CSRMatrix *matrix);


#endif
