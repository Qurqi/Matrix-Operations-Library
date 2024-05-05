#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
class Matrix
{
public:

	Matrix(float dim_row, float dim_col); 

	Matrix(float dim_row, float dim_col, vector<float> new_num);

	void DM();

	Matrix RREF();

	Matrix Row_Dec(float row_num);

	Matrix Col_Dec(float col_num);

	Matrix Row_Adj(Matrix m1);

	Matrix Col_Adj(Matrix m1);

	Matrix MA(Matrix m1);

	Matrix MS(Matrix m1);

	Matrix MM(Matrix m1);

	void SM(float scale);

	Matrix MI(Matrix m1);

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


// ****New Matrix Constructor
//
// Make new empty matrix of specified dimensions
// for combining elements of other matrices 
//  
Matrix::Matrix(float dim_row, float dim_col) {

	row = dim_row;
	col = dim_col;
}

// ****Duplicate Matrix Constructor
//
// Make new matrix of specified dimensions filled with new_num values
//
Matrix::Matrix(float dim_row, float dim_col, vector<float> new_num) {

	row = dim_row;
	col = dim_col;

	vector<float>::iterator f1;

	for (f1 = new_num.begin(); f1 != new_num.end(); ++f1) {

		num.push_back((*f1));
	}
}



// ****Display Matrix
// 
// Displays Object matrix surrounded by square brackets
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


// ****Matrix Row reduction
// 
// Reduces NxN matrix into RREF form if possible
// 
// Returns Matrix in RREF form
// 
Matrix Matrix::RREF() {

	//scalar to reduce diagonal values to 1
	//scale always equals the inverse of the number it is multiplying
	float scale = 0;

	Matrix m1(row, col, num), mR(row, col), mD(row, col);

	vector<float>::iterator i1;

	mR = m1.Row_Dec(1);

	scale = (1 / (mR.num.at(0)));

	mR.SM(scale);


	i1 = num.begin();

	//Set scale equal to the inverse of the pivot position of the respective column
	scale = (1 / (*i1));

	//multiply whole row by the scalar value

	//subtract pivot position of column-now equal to 1- multiplied by value of matrix position you are trying to eliminate

}

// ****Matrix Row_Dec
//
// Extracts row_num'nth row from object matrix
// 
// Returns a new matrix of dimension 1xN
// 
Matrix Matrix::Row_Dec(float row_num) {

	Matrix mR(1, col);

	float index = (row_num - 1) * col;

	for (float i = (0+index); i < (col+index); i++) {

		mR.num.push_back(num.at(i));
	}

	return mR;
}

// ****Matrix Col_Dec
//
// Extracts col_num'nth column from object matrix. 
// 
// Returns a new matrix of dimension Mx1
// 
Matrix Matrix::Col_Dec(float col_num) { 
	
	if (col_num > col) {
		cout << "This Matrix only has " << col << " column(s)!" << "\n";
		exit(1);
	}

	Matrix mR(row, 1);

	for (float i = 0 ; (i+1) <= row; i++) {

		mR.num.push_back(num.at((col_num-1)+(col*i)));
	}

	return mR;
}

// ****Matrix Col_Adj
//
// Adjoins argument matrix m1 to the right side of object matrix
// 
// Returns matrix of dimension Mx(2N)
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

// ****Matrix Row_Adj
//
// Adjoins argument matrix m1 to bottom of object matrix
// 
// Returns matrix of dimension (2M)xN
//
Matrix Matrix::Row_Adj(Matrix m1) {

	Matrix mR(row, col, num), mD(row, col), mC(m1.row, m1.col), mU(row+m1.row, (col)); // mR created to use member function on object matrix
	//mD and mC are dummy matrices to store row info without overwriting orignal matrix info

	vector<float>::iterator f1, f2;

	for (float i = 1; i <= m1.row; i++) { //Extract rows from mD and mC

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

// ****Matrix Addition
//
// Adds NxN matrices.
// 
// Returns sum of matrices as NxN matrix 
// 
Matrix Matrix::MA(Matrix m1) {
	
	if (!((col == m1.col) && (row == m1.row))) {
		cout << "Matrices must be of equal dimension to be summed! " << "\n";
		exit(1);
	}
	Matrix mR(row,m1.row);
	
	vector<float>::iterator f1, f2;

	for (f1 = m1.num.begin(), f2 = num.begin(); f1 != m1.num.end(), f2 != num.end(); ++f1, ++f2) {

		mR.num.push_back((*f1) + (*f2));
	}

	return mR;
}

// ****Matrix Subtraction (MS)
//
// Subtracts NxN matrices
// 
// Subtracts argument matrix m1 from object matrix
// 
// Returns difference of matrices as NxN matrix
// 
Matrix Matrix::MS(Matrix m1) {

	if (!((col == m1.col) && (row == m1.row))) {
		cout << "Matrices must be of equal dimension to be summed! " << "\n";
		exit(1);
	}
	Matrix mR(row, m1.row);

	vector<float>::iterator f1, f2;

	for (f1 = m1.num.begin(), f2 = num.begin(); f1 != m1.num.end(), f2 != num.end(); ++f1, ++f2) {

		mR.num.push_back((*f1) - (*f2));
	}

	return mR;
}

// ****Matrix Multiplication (MM)
//
// Multiplies AxN and NxB Matrices
// 
// Returns AxB matrix product
//   
Matrix Matrix::MM(Matrix m1) {

	if (!(col == m1.row)) {
		cout << "The column dimension of first matrix must equal the row dimension of the second matrix to multiply them! " << "\n";
		exit;
	}

	Matrix mR(row, m1.col), mD(row, col, num);

	vector<float>::iterator f1, f2;

	f1 = mD.num.begin();
	f2 = m1.num.begin();

	
	float k = 0;

	for (float g = 0; g < m1.row; g++) {

		for (float j = 0; j < mD.col; j++) {

			k = 0;
				//Increments col dimension of second matrix
				for (float i = 0; i < m1.col; i++) {

					k += (*(f1 + j + (mD.col * i))) * (*(f2 + (m1.col * g) + i));
				}

				mR.num.push_back(k);

			}
		}

	return mR;
}

// ****Matrix Scalar Multiplication (SM)
//
// Scales MxN matrix by scalar scale
// 
// Returns scaled MxN matrix
// 
void Matrix::SM(float scale) {

	vector<float>::iterator f1;

	for (f1 = num.begin(); f1 != num.end(); ++f1) {

		(*f1) *= scale ;
	}

}

// ****Matrix Inversion (MI)
//
// Inverts NxN matrices if possible
// 
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

// ****Matrix Transpose (MT)
//
// Transposes MxN matrix
// 
// Returns NxM matrix 
//
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

// ****Matrix Eigenvalue (MVal)
/*
vector<float> Matrix::MVal(vector<float> matrix1) {


}

// Matrix Eigenvector (MVec)

vector<float> Matrix::MVec(vector<float> matrix1) {


}
*/

int main() {

	vector<float> mat = {1,2,3,4,5,6,7,8,9,10,11,12};

	vector<float> mat1 = { 1,0,0,0,0,1,0,0,0,0,1,0 };
	
	float row = 1;

	float col = 12;

	float s = 1;



	return 0;
}

