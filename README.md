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

This project is an implementation of arithmetic operations on sparse matrices, which are matrices mostly consisting of zero elements. Using the CSR format allows for efficient memory use and faster computations, which are crucial for large-scale problems in scientific computing and machine learning.

## Features

- **Matrix operations:** Addition, Subtraction, Multiplication, Transpose
- **Matrix Market File Handling:** Reads sparse matrices from `.mtx` files
- **Memory-efficient computations:** Implements the CSR format
- **Command-line interface:** Supports operations through command-line arguments

## Getting Started

### Prerequisites

- A C compiler (e.g., GCC)
- Make utility
- Git (optional, for version control)
- A Unix-based OS (Linux/Mac) or Windows with a POSIX environment (e.g., WSL or MinGW)

### Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/yourusername/sparse-matrix-arithmetic.git
    cd sparse-matrix-arithmetic
    ```

2. Compile the project using `make`:

    ```bash
    make
    ```

    This will compile the program and generate the executable.

### Usage

The program accepts the following command-line arguments:

```bash
./main <file1.mtx> <file2.mtx> <operation> <print>
./main <file1.mtx>
