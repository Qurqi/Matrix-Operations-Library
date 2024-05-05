#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
class Matrix
{
public:

	//Constructors

	Matrix(float dim_row, float dim_col); 

	Matrix(float dim_row, float dim_col, vector<float> new_num);

	// Display Matrix

	void DM();

	Matrix RREF();

	Matrix Row_Dec(float row_num);

	Matrix Col_Dec(float col_num);

	Matrix Col_Adj(Matrix m1);

	// Matrix Addition 

	Matrix MA(Matrix m1);

	// Matrix Subtraction (MS)

	Matrix MS(Matrix m1, Matrix m2);

	// Matrix Multiplication (MM)

	Matrix MM(Matrix m1, Matrix m2);

	// Matrix Scalar Multiplication (SM)

	void SM(float scale);

	// Matrix Inversion (MI)

	Matrix MI(Matrix m1);

	// Matrix Transposition (MT)

	Matrix MT();

	// Matrix Eigenvalue (MVal)

	//Matrix MVal(Matrix m1);

	// Matrix Eigenvector (MVec)

	//Matrix MVec(Matrix m1);

private:

	float col;
	float row;

	vector<float> num;

};


Matrix::Matrix(float dim_row, float dim_col) {

	row = dim_row;
	col = dim_col;
}

Matrix::Matrix(float dim_row, float dim_col, vector<float> new_num) {

	row = dim_row;
	col = dim_col;

	vector<float>::iterator f1;

	for (f1 = new_num.begin(); f1 != new_num.end(); ++f1) {

		num.push_back((*f1));
	}
}



// Display Matrix
// 
// Intakes square matrix of dim N x N
// prints out N elements of a row per line
// surrounds matrix with square brackets
// 

void Matrix::DM() {

	//Get matrix dimension

	for (float i = 0; i < (row * col); (i + col)) {

		cout << "[ ";

			for (float j = 0; j < col; j++){
				
				cout << num.at(j+i) << " ";
			}

			cout << "]\n";

			i = i + col; //increment to the second row
				
	}

	cout << "\n";
}


//Matrix Row reduction
// 
// 
// 
// 
/*Matrix Matrix::RREF() {

	//scalar to reduce diagonal values to 1
	//scale always equals the inverse of the number it is multiplying
	float scale = 0;

	vector<float>::iterator i1;

	i1 = num.begin();

	//Set scale equal to the inverse of the pivot position of the respective column
	scale = (1 / (*i1));

	//multiply whole row by the scalar value

	//subtract pivot position of column-now equal to 1- multiplied by value of matrix position you are trying to eliminate

}*/

//Matrix Row_Dec
//
// Extracts specified row from object matrix
// 
// push_back row elements from 0+(row# - 1)(row) - row+(row# - 1)(row) into new matrix
// 
//

Matrix Matrix::Row_Dec(float row_num) {

	Matrix mR(1, col);

	float index = (row_num - 1) * col;

	for (float i = (0+index); i < (col+index); i++) {

		mR.num.push_back(num.at(i));
	}

	return mR;
}

//Matrix Col_Dec
//
// Extracts specified col from object matrix
// 
//
Matrix Matrix::Col_Dec(float col_num) {

	Matrix mR(row, 1);

	float index = (col_num - 1) * col;

	for (float i = (0 + index); i < (row + index); i++) {

		mR.num.push_back(num.at(i*col));
	}

	return mR;
}

//Matrix Col_Adj
//
// Adjoins argument matrix to the right side of object matrix
// 
// Break matrices into row vectors, push_back m1 row onto object matrix row.
// Repeat this for each row.
// Update dimensions
//
Matrix Matrix::Col_Adj(Matrix m1) {

	Matrix mR(row, col, num), mD(row, col), mC(m1.row, m1.col), mU(row, (col+m1.col)); // mR created to use member function on object matrix
	//mD and mC are dummy matrices to store row info without overwriting orignal matrix info

	vector<float>::iterator f1, f2;

	for (float i = 1; i <= m1.row ; i++) { //Extract rows from mD and mC
		
		mD = mR.Row_Dec(i);
		mC = m1.Row_Dec(i);

		for (f1 = mD.num.begin(); f1 != mD.num.end(); ++f1) {

			mU.num.push_back((*f1));
		}
		for (f2 = mC.num.begin(); f2 != mC.num.end(); ++f2) {

			mU.num.push_back((*f2));
		}
	}

	return mU;
}

//Matrix Row_Adj
//
// Adjoins argument matrix to bottom of object matrix
// push_back elements of argument onto object. return object.
//


// Matrix Addition
//
//	- Take in Vectors of defined size
//	- Define four iterators
//	- Point to beginning and end of both vectors
//	- Iterate through vectors and add each index value
//	- Store result of addition in new vector
// 
//

Matrix Matrix::MA(Matrix m1) {
	
	if (!((col == m1.col) && (row == m1.row))) {
		cout << "Matrices must be of equal dimension to be summed! " << "\n";
		exit(1);
	}
	Matrix mR(row,m1.row);
	
	vector<float>::iterator f1;
	vector<float>::iterator b1;
	
	vector<float>::iterator f2;
	vector<float>::iterator b2;

	for (f1 = m1.num.begin(), f2 = num.begin(); f1 != m1.num.end(), f2 != num.end(); ++f1, ++f2) {

		mR.num.push_back((*f1) + (*f2));
	}

	return mR;
}

// Matrix Subtraction (MS)

Matrix Matrix::MS(Matrix m1, Matrix m2) {

	if (!((m1.col == m2.col) && (m1.row == m2.row))) {
		cout << "Matrices must be of equal dimension to be summed! " << "\n";
		exit;
	}
	Matrix mR(m1.row, m1.row);

	vector<float>::iterator f1;
	vector<float>::iterator b1;

	vector<float>::iterator f2;
	vector<float>::iterator b2;

	for (f1 = m1.num.begin(), f2 = m2.num.begin(); f1 != m1.num.end(), f2 != m2.num.end(); ++f1, ++f2) {

		mR.num.push_back((*f1) - (*f2));
	}

	return mR;
}

// Matrix Multiplication (MM)
// Assume NxN matrix

Matrix Matrix::MM(Matrix m1, Matrix m2) {

	if (!(m1.col == m2.row)) {
		cout << "The column dimension of first matrix must equal the row dimension of the second matrix to multiply them! " << "\n";
		exit;
	}

	Matrix mR(m1.row, m2.col);

	vector<float>::iterator f1;
	vector<float>::iterator f2;

	f1 = m1.num.begin();
	f2 = m2.num.begin();

	
	float k = 0;

	for (float g = 0; g < m2.row; g++) {

		for (float j = 0; j < m1.col; j++) {

			k = 0;
				//Increments col dimension of second matrix
				for (float i = 0; i < m2.col; i++) {

					k += (*(f1 + j + (m1.col * i))) * (*(f2 + (m2.col * g) + i));
				}

				mR.num.push_back(k);

			}
		}

	return mR;
}

// Matrix Scalar Multiplication (SM)

void Matrix::SM(float scale) {

	vector<float>::iterator f1;

	for (f1 = num.begin(); f1 != num.end(); ++f1) {

		(*f1) *= scale ;
	}

}

// Matrix Inversion(MI)

Matrix Matrix::MI(Matrix m1) {

	if (!(row == col)) {

		cout << "Non-square matrices cannot be inverted!" << "\n";
		exit(1);
	}
	
	Matrix mR(m1.row, m1.row);

	vector<float>::iterator f1;

	for (f1 = m1.num.begin(); f1 != m1.num.end(); ++f1) {

	}

	return mR;
}

// Matrix Transposition (MT)

Matrix Matrix::MT() {

	//Change dimension
	Matrix mR(col, row);

	//Take each element column by column and push_back into new matrix

	for (float i = 0; i < col; i++) { // Increment horizontally through the matrix elements

		for (float j = 0; j < row; j++) { // Increment vertically through the matrix elements

			mR.num.push_back(num.at((row*j) + i));
			
		}
	}
	return mR;

}

// Matrix Eigenvalue (MVal)
/*
vector<float> Matrix::MVal(vector<float> matrix1) {


}

// Matrix Eigenvector (MVec)

vector<float> Matrix::MVec(vector<float> matrix1) {


}
*/
int main() {

	vector<float> mat = {1,2,3,4};

	vector<float> mat1 = { 1,0,0,0 };
	
	float row = 4;

	float col = 1;

	float s = 1;

	Matrix m1(row,col,mat), m2(row, col, mat1);
	
	m1.DM();

	m1 = m1.Col_Adj(m2);

	m1.DM();

	return 0;
}

