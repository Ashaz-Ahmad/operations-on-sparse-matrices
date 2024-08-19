# Sparse Matrix Arithmetic Operations

This project implements various arithmetic operations on sparse matrices using the **Compressed Sparse Row (CSR)** format. It efficiently handles sparse matrices stored in **Matrix Market (.mtx)** format, performing operations such as **addition**, **subtraction**, **multiplication**, and **transpose**.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Usage](#usage)


## Introduction

This project is an implementation of arithmetic operations on sparse matrices, which are matrices consisting mostly of zeros. Using the CSR format allows for efficient memory use and faster computations, which are crucial for large-scale problems in scientific computing and machine learning.

## Features

- **Matrix operations:** Addition, Subtraction, Multiplication, Transpose
- **Matrix Market File Handling:** Reads sparse matrices from `.mtx` files
- **Memory-efficient computations:** Implements the CSR format
- **Command-line interface:** Supports operations through command-line arguments

## Getting Started

### Prerequisites

- A C compiler (I have used GCC in my Makefile)
- Make utility
- Git (optional, for version control)
- A Unix-based OS (Linux/Mac) or Windows with a POSIX environment (e.g., WSL or MinGW)
- Matrices in the .mtx file format (I have included just one in this repo but you can download many more from the following website: https://sparse.tamu.edu/)

### Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/Ashaz-Ahmad/operations-on-sparse-matrices.git
    cd operations-on-sparse-matrices
    ```

2. Compile the project using `make`:

    ```bash
    make
    ```

    This will compile the program and generate the executable.

### Usage

The program accepts the following command-line arguments:

```bash
./main <file1.mtx>
./main <file1.mtx> <file2.mtx> <operation> <print>
./main <file1.mtx> tranpose <print>
```

The first command listed above is used to print the contents of any matrix in the .mtx format to the terminal. The second command listed above is used to perform an operation (addition, subtraction, or multiplication) on any two compatible matrices in the .mtx format. The third and final command listed above is used to take the transpose of a single matrix in the .mtx format. The `<print>` option can be either 0 or 1. 1 indicates that the user would like to print out all matrices to the terminal. 0 indicates that nothing will be printed out to the terminal except CPU time.