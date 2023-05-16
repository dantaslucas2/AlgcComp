#include <iostream>
#include <vector>
#include "matriz.h"
#include <cmath>
#include <iomanip>

using namespace std;

Matrix::Matrix(int r, int c) {
    rows = r;
    cols = c;
    mat.resize(rows, vector<double>(cols, 0.0));
}

void Matrix::set_element(int row, int col, double val) {
    mat[row][col] = val;
}

double Matrix::get_element(int row, int col) const {
    return mat[row][col];
}

int Matrix::get_rows() const {
    return rows;
}

int Matrix::get_cols() const {
    return cols;
}

void Matrix::print() const {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << setw(9) << std::fixed << std::setprecision(5) << mat[i][j] << " ";
        }
        std::cout << endl;
    }
}

Matrix Matrix::operator+(const Matrix& m) const {
    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.mat[i][j] = mat[i][j] + m.mat[i][j];
        }
    }

    return result;
}

Matrix Matrix::operator-(const Matrix& m) const {
    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result.mat[i][j] = mat[i][j] - m.mat[i][j];
        }
    }

    return result;
}

Matrix Matrix::operator*(const Matrix& m) const {
    Matrix result(rows, m.cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < m.cols; j++) {
            for (int k = 0; k < cols; k++) {
                result.mat[i][j] += mat[i][k] * m.mat[k][j];
            }
        }
    }

    return result;
}

Matrix& Matrix::operator+=(const Matrix& m) {
    *this = *this + m;
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& m) {
    *this = *this - m;
    return *this;
}

void Matrix::LU_decomposition(Matrix& L, Matrix& U) const {
    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square for LU decomposition");
    }

    int n = rows;

    // Inicializa matriz L com elementos 0 na diagonal e elementos 1 no resto
    L = Matrix(n, n);
    for (int i = 0; i < n; i++) {
        L.set_element(i, i, 1.0);
    }

    // Inicializa matriz U com os elementos da matriz original
    U = *this;

    // Executa o algoritmo de decomposição LU
    for (int j = 0; j < n; j++) {
        for (int i = j+1; i < n; i++) {
            double factor = U.get_element(i, j) / U.get_element(j, j);
            L.set_element(i, j, factor);
            for (int k = j; k < n; k++) {
                U.set_element(i, k, U.get_element(i, k) - factor * U.get_element(j, k));
            }
        }
    }
}

std::vector<double> Matrix::solve_LU(const vector<double>& b) const {
    if (rows != cols || rows != b.size()) {
        throw std::invalid_argument("Invalid dimensions for solving linear system");
    }

    int n = rows;

    // Executa a decomposição LU da matriz
    Matrix L, U;
    LU_decomposition(L, U);

    // Resolve o sistema linear L*y = b usando substituição para frente
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L.get_element(i, j) * y[j];
        }
        y[i] = b[i] - sum;
    }

    // Resolve o sistema linear U*x = y usando substituição para trás
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i+1; j < n; j++) {
            sum += U.get_element(i, j) * x[j];
        }
        x[i] = (y[i] - sum) / U.get_element(i, i);
    }

    return x;
}

Matrix Matrix::Cholesky_decomposition() const {
    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square for Cholesky decomposition");
    }

    int n = rows;

    // Inicializa a matriz L com elementos 0
    Matrix Q(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i+1; j++) {
            double s = 0.0;
            for (int k = 0; k < j; k++) {
                s += Q.mat[i][k] * Q.mat[j][k];
            }
            if (i == j) {
                Q.mat[i][j] = std::sqrt(mat[i][i] - s);
            } else {
                Q.mat[i][j] = (mat[i][j] - s) / Q.mat[j][j];
            }
        }
        if (Q.mat[i][i] <= 0) {
            throw std::invalid_argument("Matrix is not positive definite");
        }
    }

    return Q;
}
std::vector<double> Matrix::solve_Cholesky(const std::vector<double>& b, const Matrix& L) const {
    if (rows != cols || rows != L.get_rows() || cols != b.size()) {
        throw std::invalid_argument("Invalid dimensions for solving linear system with Cholesky decomposition");
    }

    int n = rows;

    // Resolve o sistema linear L*y = b usando substituição para frente
    std::vector<double> y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L.get_element(i, j) * y[j];
        }
        y[i] = (b[i] - sum) / L.get_element(i, i);
    }

    // Resolve o sistema linear L^T*x = y usando substituição para trás
    std::vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i+1; j < n; j++) {
            sum += L.get_element(j, i) * x[j];
        }
        x[i] = (y[i] - sum) / L.get_element(i, i);
    }

    return x;
}

vector<double> Matrix::jacobi_iterative_method(const vector<double>& b, double tolerance, int max_iterations) const {
    vector<double> x(rows, 0.0);
    vector<double> x_new(rows, 0.0);

    int iterations = 0;
    double error = tolerance + 1.0;
    while (error > tolerance && iterations < max_iterations) {
        for (int i = 0; i < rows; i++) {
            double sum = 0.0;
            for (int j = 0; j < cols; j++) {
                if (j != i) {
                    sum += mat[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / mat[i][i];
        }

        error = 0.0;
        for (int i = 0; i < rows; i++) {
            error += abs(x_new[i] - x[i]);
            x[i] = x_new[i];
        }

        iterations++;
    }
    if (iterations >= max_iterations) {
        cerr << endl << "Warnning : Maximum number of iterations reached in Jacobi." << endl;
    }
    return x;
}
vector<double> Matrix::gauss_seidel_iterative_method(const vector<double>& b, double tolerance, int max_iterations) const {
    vector<double> x(rows, 0.0);

    int iterations = 0;
    double error = tolerance + 1.0;
    while (error > tolerance && iterations < max_iterations) {
        for (int i = 0; i < rows; i++) {
            double sum = 0.0;
            for (int j = 0; j < cols; j++) {
                if (j != i) {
                    sum += mat[i][j] * x[j];
                }
            }
            double new_x = (b[i] - sum) / mat[i][i];
            error = max(error, abs(new_x - x[i]));
            x[i] = new_x;
        }

        iterations++;
    }

    if (iterations >= max_iterations) {
        cerr << endl << "Warnning : Maximum number of iterations reached in Gauss-Seidel." << endl;
    }

    return x;
}

//Task 2
//Autovetores e Autovalores
void Matrix::power_method(double& lambda, std::vector<double>& v, double tol, int max_iterations) const {
    int n = rows;
    v = std::vector<double>(n, 1.0); // inicializa com um vetor de uns
    int k = 0;
    double norm_v = 0.0;
    double norm_v_old = 0.0;
    while (k < max_iterations) {
        // multiplica a matriz pela direita pelo vetor v
        std::vector<double> u(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                u[i] += mat[i][j] * v[j];
            }
        }
        // calcula a norma do vetor u
        norm_v = 0.0;
        for (int i = 0; i < n; i++) {
            norm_v += std::pow(u[i], 2);
        }
        norm_v = std::sqrt(norm_v);
        // normaliza o vetor u e calcula o autovalor correspondente
        for (int i = 0; i < n; i++) {
            v[i] = u[i] / norm_v;
        }
        lambda = 0.0;
        for (int i = 0; i < n; i++) {
            lambda += mat[i][i] * v[i];
        }
        // verifica o critério de parada
        if (std::abs(norm_v - norm_v_old) < tol) {
            break;
        }
        norm_v_old = norm_v;
        k++;
    }    
}
void Matrix::jacobi_method(std::vector<double>& eigenvalues, Matrix& eigenvectors, double tol, int max_iterations) {
    int n = rows;
    eigenvectors = Matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            eigenvectors.set_element(i, j, i == j ? 1.0 : 0.0);
        }
    }
    eigenvalues = std::vector<double>(n, 0.0);
    int k = 0;
    double off = 0.0;
    double off_old = 0.0;
    while (k < max_iterations) {
        // procura o maior elemento fora da diagonal
        int p = -1, q = -1;
        double max_off = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double off_ij = std::abs(mat[i][j]);
                if (off_ij > max_off) {
                    max_off = off_ij;
                    p = i;
                    q = j;
                }
            }
        }
        if (max_off < tol) {
            break; // termina se todos os elementos fora da diagonal forem menores que a tolerância
        }
        // calcula o ângulo de rotação
        double theta = 0.5 * std::atan2(2 * mat[p][q], mat[q][q] - mat[p][p]);
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        // atualiza a matriz e os vetores próprios
        for (int i = 0; i < n; i++) {
            double temp1 = mat[p][i] * cos_theta + mat[q][i] * sin_theta;
            double temp2 = -mat[p][i] * sin_theta + mat[q][i] * cos_theta;
            mat[p][i] = temp1;
            mat[q][i] = temp2;
            temp1 = mat[i][p] * cos_theta + mat[i][q] * sin_theta;
            temp2 = -mat[i][p] * sin_theta + mat[i][q] * cos_theta;
            mat[i][p] = temp1;
            mat[i][q] = temp2;
            temp1 = eigenvectors.get_element(i, p) * cos_theta + eigenvectors.get_element(i, q) * sin_theta;
            temp2 = -eigenvectors.get_element(i, p) * sin_theta + eigenvectors.get_element(i, q) * cos_theta;
            eigenvectors.set_element(i, p, temp1);
            eigenvectors.set_element(i, q, temp2);
        }
        // calcula a norma dos elementos fora da diagonal
        off = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                off += std::pow(mat[i][j], 2);
            }
        }
        off = std::sqrt(off);
        // verifica o critério de parada
        if (std::abs(off - off_old) < tol) {
            break;
        }
        off_old = off;
        k++;
    }
}