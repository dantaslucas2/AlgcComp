#include "src/matriz.cpp"
#include <iostream>
#include<string>
#include <sstream>
#include <fstream>
#include "Data/load_data.cpp"

using namespace std;

void calculate_for_all_methods(Matrix A){

    cout << "Matriz A: " << endl;
    A.print();
    cout << "power metodo " << endl;

    double lambda; // armazenará o autovalor dominante
    vector<double> v; // armazenará o autovetor correspondente
    double tol = 0.00001; // tolerância para o critério de parada
    A.power_method(lambda, v, tol);

    cout << "Autovalor dominante: " << lambda << endl;
    cout << "Autovetor correspondente:" << endl;
    for (int i = 0; i < 10; i++) {
        cout << v[i] << endl;
    }

    cout << "jacobi metodo " << endl;
    std::vector<double> eigenvalues;
    Matrix eigenvectors;
    A.jacobi_method(eigenvalues, eigenvectors, tol);
    std::cout << "Autovalores:" << std::endl;
    for (int i = 0; i < eigenvalues.size(); i++) {
        std::cout << eigenvalues[i] << std::endl;
    }
    std::cout << "Autovetores:" << std::endl;
    eigenvectors.print();
}
int main(int argc, char const *argv[])
{

    Matrix mat(3, 3); // cria uma matriz 3x3
    mat.set_element(0, 0, 2);
    mat.set_element(0, 1, 1);
    mat.set_element(0, 2, 1);
    mat.set_element(1, 0, 1);
    mat.set_element(1, 1, 2);
    mat.set_element(1, 2, 1);
    mat.set_element(2, 0, 1);
    mat.set_element(2, 1, 1);
    mat.set_element(2, 2, 2);

    Matrix A = load_matriz_A();
    calculate_for_all_methods(mat);

}