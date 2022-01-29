#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double p(double x) {
	return -1.0 / (powl(x, 2) * (x + 1));
	//return 2.0 / x;
}

double q(double x) {
	return -2.0 / (powl(x, 2) * (x + 1));
	//return 2.0 / powl(x, 2);
}

double f(double x) {
	return 1.0 / ( powl(x, 4) * (x + 1) );
}

double *GetEvenGrid(int numOfDots, double a, double b, double *delta) {
	*delta = (b - a) / (double)(numOfDots);
	double *grid = new double[numOfDots + 1];
	grid[0] = a;
	for (int i = 1; i < numOfDots + 1; i++) {
		grid[i] = grid[i - 1] + (*delta);
	}
	return grid;
}

double EuclideanNorm(double *vector, int lenght) {
	double sum = 0;
	for (int i = 0; i < lenght; i++) {
		sum += powl(vector[i], 2);
	}
	sum = sqrt(sum);
	return sum;
}

void FiniteDifferenceMethod(double xa, double ya, double xb, double yb, int n, double **A, double *b, double *delta /* out */, double **grid /* out */) {
	//define the grid
	*grid = GetEvenGrid(n, xa, xb, delta);

	//fill b
	b[0] = ya;

	for (int i = 1; i <= n - 1; i++) {
		b[i] = (2.0 * powl(*delta, 2)) / (2.0 + *delta * p( (*grid)[i] )) * f( (*grid)[i] );
	}

	b[n] = yb;

	//fill A matrix
	A[0][0] = 1;

	for (int i = 1; i < n; i++) {
		double coef1 = (2 - *delta * p( (*grid)[i] )) / (2 + *delta * p((*grid)[i])); 
		double coef2 = (2 * q( (*grid)[i]) * powl(*delta, 2) - 4) / (2 + *delta * p((*grid)[i]));

		//in each line we fill in only 3 values
		A[i][(i - 1) + 0] = coef1;
		A[i][(i - 1) + 1] = coef2;
		A[i][(i - 1) + 2] = 1;
	}

	
	A[n][n] = 1;
	return;
}

void initMatrixByZero(double **A, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = 0;
		}
	}
}

int main(void) {
	/*
		y'' + p(x) * y' + q(x) * y = f(x)

		=> Ay = b -->sweep method
	*/

	const int n = 19;
	double **A = new double* [n + 1];
	
	for (int i = 0; i < n + 1; i++) {
		A[i] = new double[n + 1];
	}
	initMatrixByZero(A, n + 1);

	double *b = new double[n + 1];

	double xa = 0.2;
	double ya = 6;
	double xb = 1;
	double yb = 2;
	double delta;
	double *grid = NULL;

	FiniteDifferenceMethod(xa, ya, xb, yb, n, A, b, &delta, &grid);

	//print
	{
		ofstream file;
		file.open("matr.txt");


		file << n + 1 << " " << n + 1 << endl;

		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				file << A[i][j] << " ";
			}
			file << endl;
		}

		file << endl;

		for (int i = 0; i < n + 1; i++) {
			file << b[i] << endl;
		}

		file.close();
	}

	{
		ofstream file;
		file.open("x.txt");

		file << n + 1 << endl;

		for (int i = 0; i < n + 1; i++) {
			file << grid[i] << endl;
		}

		file.close();
	}

	//clean
	for (int i = 0; i < n + 1; i++) {
		delete[] A[i];
	}
	delete[] A;
	delete[] b;
	delete[] grid;

	return 0;
}
