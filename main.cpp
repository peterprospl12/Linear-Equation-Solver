#include <iostream>
#include "Matrix.h"
#include <cmath>
#include </usr/include/python3.10/Python.h>

void create_plot() {
    std::string python_code = R"(
try:
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    df = pd.read_csv('jacobi_solution_error_norm.csv')

    jacobi_error_norm = df['Error norm']
    jacobi_iterations = len(jacobi_error_norm)

    df = pd.read_csv('gauss_seidel_solution_error_norm.csv')
    gauss_seidel_error_norm = df['Error norm']
    gauss_seidel_iterations = len(gauss_seidel_error_norm)

    plt.yscale("log")
    plt.plot(jacobi_error_norm)
    plt.plot(gauss_seidel_error_norm)
    plt.legend(["Jacobi solution", "Gauss-Seidel solution"])
    max_iterations = max(jacobi_iterations, gauss_seidel_iterations)
    step = max_iterations // 10
    plt.xticks(np.arange(0, max_iterations+1, step))


    plt.title("Error norms")
    plt.xlabel("Iterations")
    plt.ylabel("Error norm")

    plt.show()
except Exception as e:
    print(str(e))
)";

    PyRun_SimpleString(python_code.c_str());
}

void create_matrix(Matrix<double> &matrix, int a1, int a2, int a3) {
    for (int i = 0; i < matrix.getRowSize(); i++) {
        for (int j = 0; j < matrix.getColSize(); j++) {
            if (i == j) {
                matrix[i][j] = a1;
            } else if (i + 1 == j || i - 1 == j) {
                matrix[i][j] = a2;
            } else if (i + 2 == j || i - 2 == j) {
                matrix[i][j] = a3;
            }
        }
    }
}

void create_column_vector(Matrix<double> &matrix, int f) {
    for (int i = 0; i < matrix.getRowSize(); i++) {
        matrix[i][0] = sin(i * (f + 1));
    }
}

void test_matrix();

int main() {
    Py_Initialize();
    //test_matrix();
    int f = 2;
    int e = 5;
    int N = 963;

    Matrix<double> matrix(N, N);
    create_matrix(matrix, 5+e, -1, -1);
    auto haha = matrix.lu();


    Matrix<double> col_vector(N, 1);
    create_column_vector(col_vector, f);

    Matrix<double>::jacobi_solution(matrix, col_vector, 100, pow(10, -9));
    Matrix<double>::gauss_seidel_solution(matrix, col_vector, 100, pow(10, -9));
    //create_plot();
    Matrix<double>::LU_solution(matrix, col_vector);

    Matrix<double> matrix2(N, N);
    create_matrix(matrix, 3, -1, -1);

    Matrix<double> col_vector2(N, 1);
    create_column_vector(col_vector, f);

    Matrix<double>::jacobi_solution(matrix, col_vector, 1000, pow(10, -9));
    Matrix<double>::gauss_seidel_solution(matrix, col_vector, 1000, pow(10, -9));
    Matrix<double>::LU_solution(matrix, col_vector);
    //create_plot();
    Py_Finalize();
    return 0;
}


void test_matrix() {

    Matrix<double> matrixx(3,3);
    matrixx[0][0] = 10;
    matrixx[0][1] = -7;
    matrixx[0][2] = 0;
    matrixx[1][0] = -3;
    matrixx[1][1] = 2;
    matrixx[1][2] = 6;
    matrixx[2][0] = 5;
    matrixx[2][1] = -1;
    matrixx[2][2] = 5;

    auto mas = matrixx.lu();
    std::cout << mas.first << std::endl;
    std::cout << mas.second << std::endl;




    Matrix<double> matrix(2, 3);
    matrix[0][0] = 2;
    matrix[0][1] = 1;
    matrix[0][2] = 3;
    matrix[1][0] = -1;
    matrix[1][1] = 4;
    matrix[1][2] = 0;

    Matrix<double> matrixA(3, 3);
    matrixA[0][0] = 1;
    matrixA[0][1] = 3;
    matrixA[0][2] = 2;
    matrixA[1][0] = -2;
    matrixA[1][1] = 0;
    matrixA[1][2] = 1;
    matrixA[2][0] = 5;
    matrixA[2][1] = -3;
    matrixA[2][2] = 2;

    auto x = matrix * matrixA;
    std::cout << x << std::endl;

    Matrix<double> matrixB(2, 3);
    matrixB[0][0] = 9;
    matrixB[0][1] = 1;
    matrixB[0][2] = -4;
    matrixB[1][0] = 3;
    matrixB[1][1] = -5;
    matrixB[1][2] = 0;

    Matrix<double> matrixC(2, 3);
    matrixC[0][0] = 6;
    matrixC[0][1] = -4;
    matrixC[0][2] = 5;
    matrixC[1][0] = 7;
    matrixC[1][1] = -2;
    matrixC[1][2] = 3;

    auto y = matrixB + matrixC;
    std::cout << y << std::endl;

    Matrix<double> matrixD(3, 3);
    matrixD[0][0] = -1;
    matrixD[0][1] = 0;
    matrixD[0][2] = 2;
    matrixD[1][0] = 2;
    matrixD[1][1] = -1;
    matrixD[1][2] = 1;
    matrixD[2][0] = 3;
    matrixD[2][1] = 1;
    matrixD[2][2] = 1;
    Matrix<double> matrixE(3, 1);
    matrixE[0][0] = -1;
    matrixE[1][0] = 2;
    matrixE[2][0] = 1;
    Matrix<double> matrixF(3, 1);
    matrixF[0][0] = 10;
    matrixF[1][0] = 7;
    matrixF[2][0] = -5;
    auto z = Matrix<double>::residuum(matrixD, matrixE, matrixF);
    std::cout << z << std::endl;

    Matrix<double> matrixG(3, 1);
    matrixG[0][0] = 1;
    matrixG[1][0] = -2;
    matrixG[2][0] = 3;
    auto p = Matrix<double>::norm(matrixG);
    std::cout << p << std::endl;

    Matrix<double> matrixH(4, 4);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            matrixH[i][j] = i + j + 1; // Generowanie losowej liczby zmiennoprzecinkowej
        }
    }
    std::cout << matrixH << std::endl;

    auto a = matrixH.diag().diag();
    std::cout << a << std::endl;


    Matrix<double> matrix1(3, 3);
    matrix1[0][0] = 2;
    matrix1[0][1] = -1;
    matrix1[0][2] = 1;
    matrix1[1][0] = 1;
    matrix1[1][1] = -1;
    matrix1[1][2] = 2;
    matrix1[2][0] = 5;
    matrix1[2][1] = -2;
    matrix1[2][2] = 2;
    Matrix<double> matrix2(3, 1);
    matrix2[0][0] = 7;
    matrix2[1][0] = 6;
    matrix2[2][0] = 15;

    auto a21 = Matrix<double>::jacobi_solution(matrix1, matrix2, 10000, pow(10, -6));

}