#include "src/matriz.cpp"
#include <iostream>
#include<string>
#include <sstream>
#include <fstream>
#include "Data/load_data.cpp"

using namespace std;

void calculate_for_all_methods(Matrix A, Matrix L, Matrix U, Matrix Q,vector<double> b, string name){

    cout << setw(7) << endl << "Vetor b " << name << " : " << endl;
    for (int i = 0; i < b.size(); i++) {
        cout << b[i] << endl;
    }
    
    std::vector<double> x1_lu = U.solve_LU(L.solve_LU(b));
    std::cout << setw(7) << endl << "vetor x " << name << " LU:" << endl;
    for (int i = 0; i < x1_lu.size(); i++) {
        cout << x1_lu[i] << endl;
    }

    std::vector<double> x1_ch = A.solve_Cholesky(b, Q);
    std::cout << setw(7) << endl << "vetor x " << name << " cholesky:" << endl;
    for (int i = 0; i < x1_ch.size(); i++) {
        cout << x1_ch[i] << endl;
    }

    vector<double> jb1 = A.jacobi_iterative_method(b);
    std::cout  << endl << "vetor x " << name << " jacobi:" << endl;
    for (int i = 0; i < jb1.size(); i++) {
        cout << jb1[i] << endl;
    }

    vector<double> gs1 = A.gauss_seidel_iterative_method(b);
    std::cout  << endl << "vetor x " << name << " gauss seidel:" << endl;
    for (int i = 0; i < gs1.size(); i++) {
        cout << gs1[i] << endl;
    }
}
int main(int argc, char const *argv[])
{
    Matrix A = load_matriz_A();
    A.print();
    cout << "Matriz A: " << endl;
    Matrix L, U;
    A.LU_decomposition(L, U);
    Matrix Q = A.Cholesky_decomposition();

    cout << endl << "Matriz L: " << endl;
    L.print();
    cout << endl << "Matriz U: " << endl;
    U.print();
    cout << endl << "Matriz Q cholesky: " << endl;
    Q.print();

    vector<double> b1 = load_vector_1();
    calculate_for_all_methods(A, L, U , Q, b1, "Vector_B_01");

    vector<double> b2 = load_vector_1();
    calculate_for_all_methods(A, L, U , Q, b2, "Vector_B_02");

    vector<double> b3 = load_vector_1();
    calculate_for_all_methods(A, L, U , Q, b3, "Vector_B_03");

    return 0;
}
