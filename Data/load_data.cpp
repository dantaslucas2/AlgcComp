#include <iostream>
#include <fstream>
#include <vector>
#include "../src/matriz.h"
#include<string>
#include <sstream>


using namespace std;

Matrix load_matriz_A(){

    ifstream file("Data/Matriz_A.dat", ios::binary);
    
        int rows = 10;
        int cols = 10;
        Matrix matrix(rows, cols);
        if (file.is_open()) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double val;
                //file.read(reinterpret_cast<char*>(&val), sizeof(double));
                file >> val;
                matrix.set_element(i, j, val);
            }
        }

        file.close();
    } else {
        std::cout << "Unable to open file." << std::endl;
    }
    return matrix;
}

vector<double> load_vector(string name){
    string name_file = "Data/" + name + ".dat";
    ifstream infile(name_file);
    int size = 10;

    vector<double> v(size);
    if (infile.is_open()) {
        for (int i = 0; i < size; i++) {
            double val;
            infile >> val;
            v[i] = val;
        }

    }else{
        std::cout << "FAil to open vector " << name_file << endl;
    }
    return v;
}

vector<double> load_vector_1(){
    return load_vector("Vetor_B_01");
}

vector<double> load_vector_2(){
    return load_vector("Vetor_B_02");
}

vector<double> load_vector_3(){
    return load_vector("Vetor_B_03");
}