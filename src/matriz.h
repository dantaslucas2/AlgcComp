#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

using namespace std;

struct Matrix {
    private:
        int rows, cols;
        std::vector<std::vector<double>> mat;

    public:
        Matrix() : rows(0), cols(0) {}

        Matrix(int r, int c);

        void set_element(int row, int col, double val);

        double get_element(int row, int col) const;

        int get_rows() const;

        int get_cols() const;

        void print() const;

        Matrix operator+(const Matrix& m) const;

        Matrix operator-(const Matrix& m) const;

        Matrix operator*(const Matrix& m) const;

        Matrix& operator+=(const Matrix& m);

        Matrix& operator-=(const Matrix& m);

        Matrix& operator*=(const Matrix& m);

        //Task 1
        void LU_decomposition(Matrix& L, Matrix& U) const; 

        std::vector<double> solve_LU(const vector<double>& b) const;

        Matrix Cholesky_decomposition() const;

        std::vector<double> solve_Cholesky(const std::vector<double>& b, const Matrix& L) const;

        vector<double> jacobi_iterative_method(const vector<double>& b, double tolerance = 0.0000001, int max_iterations = 1000) const;

        vector<double> gauss_seidel_iterative_method(const vector<double>& b, double tolerance = 0.0000001, int max_iterations = 1000) const;

        //Task 2
        //Autovetores e Autovalores
        void power_method(double& lambda, std::vector<double>& v, double tol, int max_iterations = 1000000) const;

        void jacobi_method(std::vector<double>& eigenvalues, Matrix& eigenvectors, double tol, int max_iterations = 10000);
};

#endif // MATRIX_H
