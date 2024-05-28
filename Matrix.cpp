#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
class Matrix
{
public:

	Matrix(float dim_row, float dim_col); 

	Matrix(float dim_col);

	Matrix(float dim_row, float dim_col, vector<float> new_num);

	void DM();

	Matrix RREF();

	Matrix Row_Dec(float row_num);

	Matrix Col_Dec(float col_num);

	Matrix Row_Adj(Matrix m1);

	void Row_Swap(float row_num1, float row_num2);

	Matrix Col_Adj(Matrix m1);

	Matrix Row_ExtRCT(float row_num);

	Matrix Row_ExtRCB(float row_num);

	float Ext_PV(float row_num);

	float Ext_RV(float row_num, float col_num);
	
	Matrix Row_Rep(float row_num, Matrix mR);

	Matrix Row_Sub(float row_num, Matrix mS);

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

// ****New Empty Matrix Constructor Row
//
// Make new empty matrix of specified dimensions
// for combining elements of other matrices 
//  
Matrix::Matrix(float dim_col) {

	row = 0;
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

Matrix Matrix::Row_ExtRCT(float row_num) {
	
	Matrix mT(1,col), mC(row_num - 1,col), mO(row,col,num);

	for (float i = 1; i < row_num; i++) {

		mT = mO.Row_Dec(i);

		for (float k = 0; k < mT.col; k++) {

			mC.num.push_back(mT.num.at(k));
		}
	}

	return mC;
}

Matrix Matrix::Row_ExtRCB(float row_num) {

	Matrix mT(1, col), mC(row-row_num, col), mO(row, col, num);

	for (float i = (row_num+1); i < (row+1); i++) {

		mT = mO.Row_Dec(i);

		for (float k = 0; k < mT.col; k++) {

			mC.num.push_back(mT.num.at(k));
		}
	}

	return mC;

}

void Matrix::Row_Swap(float row_num1, float row_num2) {

	vector<float>::iterator f1, f2, fT;
	
	vector<float> temp = { 0 };

	f1 = num.begin();
	f2 = num.begin();
	fT = temp.begin();

	f1 = (f1 + (col * (row_num1 - 1)) );

	f2 = (f2 + (col * (row_num2 - 1)) );

	for (float i = 0; i < col; i++) {

		*(fT) = *(f1);
		*(f1) = *(f2);
		*(f2) = *(fT);

		f1++;
		f2++;

	}


}

float Matrix::Ext_PV(float row_num) {

	if (row != col) {
		printf("Matrix not NxN! Only working for NxN matrices");
		exit(1);
	}

	vector<float>::iterator f1;
	float PV;

	f1 = num.begin();

	PV = *(f1 + (row_num - 1) * col + (row_num - 1));

	return PV;

}

float Matrix::Ext_RV(float row_num, float col_num) {

	if (row != col) {
		printf("Matrix not NxN! Only working for NxN matrices");
		exit(1);
	}

	vector<float>::iterator f1;
	float RV;

	f1 = num.begin();

	RV = *(f1 + (row_num - 1) * col + col_num);

	return RV;

}

Matrix Matrix::Row_Rep(float row_num, Matrix mR) {

	if (mR.col != col) {
		printf("Matrix column dimensions dont match!");
		exit(1);
	}

	Matrix mT(col), mT1(col);

	if (row_num == 1) {

		mT = Row_ExtRCB(1);

		mR = mR.Row_Adj(mT);

		return mR;
	}
	else {
		mT = Row_ExtRCT(row_num);
		mT1 = Row_ExtRCB(row_num);
		mR = mT.Row_Adj(mR);
		mR = mR.Row_Adj(mT1);
		return mR;
	}

}

// ****Matrix Row Subtraction operation
// 
// Extracts indicated row(arg1) from obj matrix and subtracts 
// given row (arg2) from extracted row.
// 
// Returns matrix of original dimension with indicated row replaced with 
// result of row subtraction
// 
Matrix Matrix::Row_Sub(float row_num, Matrix mS) {

	if (mS.col != col) {

		printf("The column dimension of the subtraction matrix and the object matrix do not match.\n\n Please try again.\n");
		exit(1);
	}

	Matrix mT(0,0), mC(row_num-1, mS.col), mC1(row-row_num, mS.col);

	mT = Row_Dec(row_num);

	mT = mT.MS(mS);

	if (row_num > 1) {

		mC = Row_ExtRCT(row_num);

		mT = mC.Row_Adj(mT);		
	}
	
	mC1 = Row_ExtRCB(row_num);
	mT = mT.Row_Adj(mC1);

	return mT;


}

// ****Matrix Row reduction
// 
// Reduces NxN matrix into RREF form if possible
// 
// Returns Matrix in RREF form
// 
Matrix Matrix::RREF() {

	float scal;

	Matrix mO(row, col, num), mT(col), mT1(col);

	for (float i = 1; i < row +1 ; i++) {

		scal = mO.Ext_PV(i);

		if (scal != 0) {
			scal = (1 / scal);
			mT = mO.Row_Dec(i);
			mT.SM(scal);
			mO = mO.Row_Rep(i, mT);
		}
		else if (scal == 0) {
			
			mO.Row_Swap(i, (i + 1));

			scal = mO.Ext_PV(i);
			scal = (1 / scal);
			mT = mO.Row_Dec(i);
			mT.SM(scal);
			mO = mO.Row_Rep(i, mT);
		}

		for (float j = (i + 1); j < row + 1; j++) {

			scal = mO.Ext_RV(j, (i-1));
			mT1 = mT;
			mT1.SM(scal);
			mO = mO.Row_Sub(j, mT1);

		}
	}

	return mO;

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

	Matrix mR(row, col, num), mU(row+m1.row, col); // mR created to use member function on object matrix
	//mD and mC are dummy matrices to store row info without overwriting orignal matrix info

	vector<float>::iterator f1, f2;

		for (f1 = mR.num.begin(); f1 != mR.num.end(); ++f1) {

			mU.num.push_back((*f1));
		}
		for (f2 = m1.num.begin(); f2 != m1.num.end(); ++f2) {

			mU.num.push_back((*f2));
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
		cout << "Matrices must be of equal dimension to be subtracted! " << "\n";
		exit(1);
	}
	Matrix mR(row, col);

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

	vector<float> mat = {1,2,3,1,2,3,7,8,9};

	vector<float> mat1 = { 4,5,6 };
	
	float row = 3;

	float col = 3;

	float s = 2;

	Matrix m1(row, col, mat), m2(1, col, mat1);

	m1.RREF();

	m1.DM();



	return 0;
}

