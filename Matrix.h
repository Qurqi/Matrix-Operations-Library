// Matrix.h
//
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// 
// Defines functions for all matrix operations
// neccessary to create a machine learning model
// 
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// 
// Functions include:
// 
// || Display Matrix (DM)
// || Matrix Addition (MA)
// || Matrix Subtraction (MS)
// || Matrix Multiplication (MM)
// || Matrix Scalar Multiplication (SM)
// || Matrix Inversion (MI)
// || Matrix Transposition (MT)
// || Matrix Eigenvalue (MVal)
// || Matrix Eigenvector (MVec)
// || More To come...
// 
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// 
// How do the functions work ?
// 
// - Matrices are formatted as float type vectors
// 
// - Functions take Matrix or Matrices and float type Scalar ,if 
//	 neccessary, as input.
// 
// - Functions output Matrix, Column vector or Scalar of 
//	 type float.
// 
//	-Assume all matrices are square
// 
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- 
// 
// 

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
using namespace std;

// Display Matrix

void DM(vector<float> m1);

// Matrix Addition 

vector<float> MA(vector<float> matrix1, vector<float> matrix2);

// Matrix Subtraction (MS)

vector<float> MS(vector<float> matrix1, vector<float> matrix2);

// Matrix Multiplication (MM)

vector<float> MM(vector<float> matrix1, vector<float> matrix2);

// Matrix Scalar Multiplication (SM)

vector<float> SM(vector<float> matrix1, float scale);

// Matrix Inversion (MI)

vector<float> MI(vector<float> matrix1);

// Matrix Transposition (MT)

vector<float> MT(vector<float> matrix1);

// Matrix Eigenvalue (MVal)

vector<float> MVal(vector<float> matrix1);

// Matrix Eigenvector (MVec)

vector<float> MVec(vector<float> matrix1);

#endif