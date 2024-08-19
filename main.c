#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "functions.h"

int main(int argc, char *argv[]) 
{
	if (argc < 2 || argc > 5)
	{
		fprintf(stderr, "Incorrect format. Please try again with one of the following formats: ./main <file1.mtx> or ./main <file1.mtx> <file2.mtx> <operation> <print> or ./main <file1.mtx> tranpose <print>\n");
		return EXIT_FAILURE;
	}

	const char *filename = argv[1]; 
	CSRMatrix A;
	ReadMMtoCSR(filename, &A);

	clock_t start_time, end_time;
	double cpu_time_used;
	start_time = clock();

	if (argc == 2)
	{
		printCSRMatrix(&A);
		return EXIT_SUCCESS;
	}

	if (argc == 3)
	{
		fprintf(stderr, "Incorrect format. Please try again with one of the following formats: ./main <file1.mtx> or ./main <file1.mtx> <file2.mtx> <operation> <print> or ./main <file1.mtx> tranpose <print>\n");
		return EXIT_FAILURE;
	}
	

	if (argc == 4 && (strcmp(argv[2], "transpose") == 0)) 
	{
        CSRMatrix AT = transposeCSRMatrix(&A);
        if (atoi(argv[3]) == 1)
		{
            printf("Matrix A:\n"); 
            printCSRMatrix(&A); 
			printf("\n");
            printf("Matrix At:\n"); 
            printCSRMatrix(&AT);
			printf("\n");
        }
        freeCSRMatrix(&A); 
        freeCSRMatrix(&AT); 
        return EXIT_SUCCESS;
    }

	const char *filename_2 = argv[2]; 
	CSRMatrix B; 
	ReadMMtoCSR(filename_2, &B); 
	CSRMatrix C; 

	const char *operation = argv[3];

	if (argc == 5) 
	{
		if (strcmp(operation, "addition") == 0) 
		{
			C = addCSRMatrices(&A, &B); 
		} 
		else if (strcmp(operation, "subtraction") == 0) 
		{
			C = subtractCSRMatrices(&A, &B);
		} 
		else if (strcmp(operation, "multiplication") == 0)
		{
			C = multiplyCSRMatrices(&A, &B); 
		} 
		else 
		{
			fprintf(stderr, "Unsupported operation. Please use one of the following: addition, subtraction, multiplication, transpose.\n");
			freeCSRMatrix(&A);
			freeCSRMatrix(&B);
			return EXIT_FAILURE;
		}
	}

    if (atoi(argv[4]) == 1) 
	{
        printf("Matrix A:\n");
        printCSRMatrix(&A);
		printf("\n");
        printf("Matrix B:\n");
        printCSRMatrix(&B);
		printf("\n");
        printf("Matrix C (result):\n");
        printCSRMatrix(&C);
		printf("\n");
    }

    freeCSRMatrix(&A);
    freeCSRMatrix(&B);
    freeCSRMatrix(&C);

	end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time: %f seconds\n", cpu_time_used);

    return EXIT_SUCCESS;
}

