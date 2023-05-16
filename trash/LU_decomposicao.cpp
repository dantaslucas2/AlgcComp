#include "src/matriz.cpp"
#include <iostream>

using namespace std;


int main(int argc, char const *argv[])
{

// Cria uma matriz 3x3 e atribui alguns valores a ela
    Matrix A(3, 3);
    A.set_element(0, 0, 2.0);
    A.set_element(0, 1, 1.0);
    A.set_element(0, 2, 1.0);
    A.set_element(1, 0, 4.0);
    A.set_element(1, 1, 3.0);
    A.set_element(1, 2, 3.0);
    A.set_element(2, 0, 8.0);
    A.set_element(2, 1, 7.0);
    A.set_element(2, 2, 9.0);

    // Cria um vetor b com os termos independentes da equação linear
    std::vector<double> b = {1.0, 2.0, 3.0};

    // Resolve o sistema linear Ax = b usando decomposição LU
    try {
        Matrix L, U;
        A.LU_decomposition(L, U);

        cout<< "Matriz A:"<<endl;
        A.print();
        cout<< "Matriz L:"<<endl;
        L.print();
        cout<< "Matriz U:"<<endl;
        U.print();
        std::vector<double> x = U.solve_LU(L.solve_LU(b));

        // Imprime a solução do sistema
        for (int i = 0; i < x.size(); i++) {
            std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }
    }
    catch (const std::invalid_argument& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}
