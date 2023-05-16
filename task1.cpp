#include "src/matriz.cpp"
#include <iostream>
#include<string>
#include <sstream>
#include <fstream>
#include "Data/load_data.cpp"

using namespace std;

void calculate_for_all_methods(Matrix A, vector<double> b){

    cout << "Matriz A: " << endl;
    A.print();
    Matrix L, U;
    A.LU_decomposition(L, U);
    Matrix Q = A.Cholesky_decomposition();

    cout << endl << "Matriz L: " << endl;
    L.print();
    cout << endl << "Matriz U: " << endl;
    U.print();
    cout << endl << "Matriz Q cholesky: " << endl;
    Q.print();

    cout << setw(7) << endl << "Vetor b: " << endl;
    for (int i = 0; i < b.size(); i++) {
        cout << b[i] << endl;
    }
    
    std::vector<double> x1_lu = U.solve_LU(L.solve_LU(b));
    std::cout << setw(7) << endl << "vetor x1 LU:" << endl;
    for (int i = 0; i < x1_lu.size(); i++) {
        cout << x1_lu[i] << endl;
    }

    std::vector<double> x1_ch = A.solve_Cholesky(b, Q);
    std::cout << setw(7) << endl << "vetor x1 cholesky:" << endl;
    for (int i = 0; i < x1_ch.size(); i++) {
        cout << x1_ch[i] << endl;
    }

    vector<double> jb1 = A.jacobi_iterative_method(b);
    std::cout  << endl << "vetor x1 jacobi:" << endl;
    for (int i = 0; i < jb1.size(); i++) {
        cout << jb1[i] << endl;
    }

    vector<double> gs1 = A.gauss_seidel_iterative_method(b);
    std::cout  << endl << "vetor x1 gauss seidel:" << endl;
    for (int i = 0; i < gs1.size(); i++) {
        cout << gs1[i] << endl;
    }
}
int main(int argc, char const *argv[])
{
    Matrix A = load_matriz_A();
    // vector<double> b1 = load_vector_1();
    // calculate_for_all_methods(A,b1);

    vector<double> b2 = load_vector_1();
    calculate_for_all_methods(A,b2);

    // vector<double> b3 = load_vector_1();
    // calculate_for_all_methods(A,b3);

    return 0;
}
