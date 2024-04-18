//
// Created by piotr on 4/17/24.
//

#include "tasks.h"
#include "Matrix.h"
#include "plots.h"

void task_B(int f, int e, int N) {
    Matrix<double> matrix(N, N);
    int a1 = 5 + e;
    create_matrix(matrix, a1, -1, -1);

    Matrix<double> col_vector(N, 1);
    create_column_vector(col_vector, f);

    Matrix<double>::jacobi_solution(matrix, col_vector, 100, pow(10, -9), nullptr);
    Matrix<double>::gauss_seidel_solution(matrix, col_vector, 100, pow(10, -9), nullptr);

    std::string jacobi_file = "jacobi_error_norm_" + std::to_string(a1) + "_" + std::to_string(N) + ".csv";
    std::string gs_file = "gauss_seidel_error_norm_" + std::to_string(a1) + "_" + std::to_string(N) + ".csv";
    std::string plot_name = "Jacobi_GS_error_norms_" + std::to_string(a1) + "_" + std::to_string(N) + ".png";
    plot_error_norm(jacobi_file, gs_file, plot_name);
}

void task_C(int N, int f) {
    Matrix<double> matrix(N, N);
    int a1 = 3;
    create_matrix(matrix, a1, -1, -1);

    Matrix<double> col_vector(N, 1);
    create_column_vector(col_vector, f);

    Matrix<double>::jacobi_solution(matrix, col_vector, 1000, pow(10, -9), nullptr);
    Matrix<double>::gauss_seidel_solution(matrix, col_vector, 1000, pow(10, -9), nullptr);
    Matrix<double>::LU_solution(matrix, col_vector, nullptr, nullptr);

    std::string jacobi_file = "jacobi_error_norm_" + std::to_string(a1) + "_" + std::to_string(N) + ".csv";
    std::string gs_file = "gauss_seidel_error_norm_" + std::to_string(a1) + "_" + std::to_string(N) + ".csv";
    std::string plot_name = "Jacobi_GS_error_norms_" + std::to_string(a1) + "_" + std::to_string(N) + ".png";
    plot_error_norm(jacobi_file, gs_file, plot_name);
}

void task_E(int f, int e) {
    std::vector<int> N{100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
    std::vector<double> LU_error_norms(N.size());
    std::vector<double> LU_durations(N.size());
    std::vector<double> Jacobi_durations(N.size());
    std::vector<double> GS_durations(N.size());

    for (int i = 0; i < N.size(); i++) {
        int size = N[i];
        std::cout << "Size: " << size << std::endl;
        Matrix<double> matrix(size, size);
        Matrix<double> vector(size, 1);
        create_matrix(matrix, 5 + e, -1, -1);
        create_column_vector(vector, f);

        Matrix<double>::jacobi_solution(matrix, vector, size, pow(10, -9), &Jacobi_durations[i]);
        std::cout << std::endl;
        Matrix<double>::gauss_seidel_solution(matrix, vector, size, pow(10, -9), &GS_durations[i]);
        std::cout << std::endl;
        Matrix<double>::LU_solution(matrix, vector, &LU_error_norms[i], &LU_durations[i]);
        std::cout << std::endl;

    }

    std::fstream file;
    std::string name = "LU_solution_data.csv";
    file.open(name, std::fstream::out);
    file << "Size,Duration\n";
    for (int i = 0; i < N.size(); i++) {
        file << N[i] << "," << LU_durations[i] << "\n";
    }
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    file.close();

    name = "Jacobi_solution_data.csv";
    file.open(name, std::fstream::out);
    file << "Size,Duration\n";
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < N.size(); i++) {
        file << N[i] << "," << Jacobi_durations[i] << "\n";
    }
    file.close();

    name = "GS_solution_data.csv";
    file.open(name, std::fstream::out);
    file << "Size,Duration\n";
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < N.size(); i++) {
        file << N[i] << "," << GS_durations[i] << "\n";
    }
    file.close();
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
