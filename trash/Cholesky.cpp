#include "src/matriz.cpp"
#include <iostream>

using namespace std;


int main(int argc, char const *argv[])
{
    // Matrix A(2, 2);
    // A.set_element(0, 0, 1.0);
    // A.set_element(0, 1, 2.0);
    // A.set_element(0, 2, 3.0);
    // A.set_element(1, 0, 2.0);
    // A.set_element(1, 1, 5.0);
    // A.set_element(1, 2, 2.0);
    // A.set_element(2, 0, 6.0);
    // A.set_element(2, 1, 2.0);
    // A.set_element(2, 2, 1.0);

    // Matrix B(2, 2);
    // B.set_element(0, 0, 5.0);
    // B.set_element(0, 1, 6.0);
    // B.set_element(1, 0, 7.0);
    // B.set_element(1, 1, 8.0);

    // Matrix C = A + B;
    // C.print();

    // std::vector<double> b = {6, -4, 27};
    // std::vector<int> p;
    // A.LU_decomposition(p);
    // std::vector<double> x = A.solve_LU(b, p);
    // //std::cout << x << endl;


// Cria uma matriz 3x3 e atribui alguns valores a ela
    // Cria uma matriz aleatória simétrica positiva definida
    Matrix A(3, 3);
    A.set_element(0, 0, 5.0);
    A.set_element(0, 1, 2.0);
    A.set_element(0, 2, 3.0);
    A.set_element(1, 0, 2.0);
    A.set_element(1, 1, 6.0);
    A.set_element(1, 2, 4.0);
    A.set_element(2, 0, 3.0);
    A.set_element(2, 1, 4.0);
    A.set_element(2, 2, 7.0);

    // Realiza a decomposição Cholesky
    Matrix Q = A.Cholesky_decomposition();
    std::vector<double> b {1,2,3};
    std::vector<double> x = A.solve_Cholesky(b, Q);
    // Imprime a matriz original e a matriz L
    std::cout << "Matriz A:" << std::endl;
    A.print();
    std::cout << "Matriz Q:" << std::endl;
    Q.print();
    std::cout << "Vetor X:" << std::endl;
    for (int i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }


    return 0;
}
