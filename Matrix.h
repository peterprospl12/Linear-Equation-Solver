//
// Created by piotr on 4/15/24.
//

#ifndef SOLVING_LINEAR_EQUATIONS_MATRIX_H
#define SOLVING_LINEAR_EQUATIONS_MATRIX_H

#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <ctime>

template<typename T>
class Matrix {
private:
    int row_size;
    int col_size;
    std::vector<std::vector<T>> matrix;

public:
    explicit Matrix(int row_size, int col_size)
            : row_size(row_size), col_size(col_size), matrix(row_size, std::vector<T>(col_size)) {
    }

    [[nodiscard]] int getColSize() const {
        return col_size;
    }

    [[nodiscard]] int getRowSize() const {
        return row_size;
    }

    auto begin() {
        return matrix.begin();
    }

    auto end() {
        return matrix.end();
    }

    auto begin() const {
        return matrix.begin();
    }

    auto end() const {
        return matrix.end();
    }


    const std::vector<std::vector<T>> &getMatrix() const {
        return matrix;
    }

    void transpose() {
        std::vector<std::vector<T>> new_vector(this->col_size, std::vector<T>(this->row_size));
        for (int i = 0; i < this->row_size; i++) {
            for (int j = 0; j < this->col_size; j++) {
                new_vector[j][i] = this->matrix[i][j];
            }
        }
        this->matrix = new_vector;
    }

    Matrix<T> triu(int k = 0) {
        Matrix<T> new_vector(this->row_size, this->col_size);
        for (int i = 0; i < this->row_size; i++) {
            for (int j = i + k; j < this->col_size; j++) {
                new_vector[i][j] = this->matrix[i][j];
            }
        }
        return new_vector;
    }

    Matrix<T> tril(int k = 0) {
        Matrix<T> new_vector(this->row_size, this->col_size);
        for (int i = 0; i < this->row_size; i++) {
            for (int j = 0; j <= i + k; j++) {
                new_vector[i][j] = this->matrix[i][j];
            }
        }
        return new_vector;
    }

    Matrix<T> diag() {
        if (this->col_size == this->row_size) {
            Matrix<T> new_vector(this->col_size, 1);
            for (int i = 0; i < this->row_size; i++) {
                new_vector[i][0] = this->matrix[i][i];
            }
            return new_vector;
        } else if (this->col_size == 1 || this->row_size == 1) {
            int size = 0;
            if (col_size == 1) {
                size = row_size;
            } else {
                size = col_size;
            }

            Matrix<T> new_matrix(size, size);
            int counter = 0;
            for (auto row: this->matrix) {
                for (auto value: row) {
                    new_matrix[counter][counter] = value;
                    counter++;
                }
            }
            return new_matrix;
        } else {
            throw std::invalid_argument("Invalid matrix size");
        }
    }

    std::pair<Matrix<T>, Matrix<T>> lu() {
        if (this->col_size != this->row_size) {
            throw std::invalid_argument("Invalid matrix size");
        }
        int m = this->col_size;
        auto U = *this;
        auto L = eye(m);

        for (int i = 1; i < m; i++) {
            for (int j = 0; j < i; j++) {
                L[i][j] = U[i][j] / U[j][j];
                for (int z = j; z < m; z++) {
                    U[i][z] -= L[i][j] * U[j][z];
                }
            }

        }

        return std::make_pair(L, U);
    }



    /////////////////////////
    // OPERATORS
    /////////////////////////

    std::vector<T> &operator[](size_t row) {
        return matrix[row];
    }

    const std::vector<T> &operator[](size_t row) const {
        return matrix[row];
    }

    Matrix<T> operator*(const Matrix &matrixB) {
        if (this->col_size != matrixB.getRowSize()) {
            throw std::invalid_argument("Invalid matrix sizes to multiply");
        }
        Matrix<T> new_matrix(this->row_size, matrixB.getColSize());

        for (int i = 0; i < this->row_size; i++) {
            for (int j = 0; j < matrixB.getColSize(); j++) {
                for (int z = 0; z < this->col_size; z++) {
                    new_matrix[i][j] += this->matrix[i][z] * matrixB[z][j];
                }
            }
        }
        return new_matrix;
    }

    Matrix<T> operator-(const Matrix &matrixB) {
        if (this->col_size != matrixB.getColSize() || this->row_size != matrixB.getRowSize()) {
            throw std::invalid_argument("Invalid matrix sizes to multiply");
        }
        Matrix<T> new_matrix(this->row_size, this->col_size);

        for (int i = 0; i < this->row_size; i++) {
            for (int j = 0; j < this->col_size; j++) {
                new_matrix[i][j] = this->matrix[i][j] - matrixB[i][j];
            }
        }
        return new_matrix;
    }

    Matrix<T> operator+(const Matrix &matrixB) {
        if (this->col_size != matrixB.getColSize() || this->row_size != matrixB.getRowSize()) {
            throw std::invalid_argument("Invalid matrix sizes to multiply");
        }
        Matrix<T> new_matrix(this->row_size, this->col_size);

        for (int i = 0; i < this->row_size; i++) {
            for (int j = 0; j < this->col_size; j++) {
                new_matrix[i][j] = this->matrix[i][j] + matrixB[i][j];
            }
        }
        return new_matrix;
    }


    friend std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrixA) {
        for (const auto &row: matrixA.getMatrix()) {
            for (const auto &col: row) {
                os << col << " ";
            }
            os << std::endl;
        }
        return os;
    }

    /////////////////////////
    // STATIC
    /////////////////////////

    static Matrix<T> residuum(Matrix<T> &A, Matrix<T> &x, Matrix<T> &b) {
        if (b.col_size != 1) {
            throw std::invalid_argument("Invalid size of matrix b");
        } else if (x.col_size != 1) {
            throw std::invalid_argument("Invalid size of matrix x");
        }

        Matrix new_vector(x.getRowSize(), 1);

        for (int i = 0; i < A.getRowSize(); i++) {
            for (int j = 0; j < A.getColSize(); j++) {
                new_vector[i][0] += A[i][j] * x[j][0];
            }
            new_vector[i][0] -= b[i][0];
        }
        return new_vector;
    }

    static double norm(const Matrix<T> &v) {
        if (v.getColSize() > 1 && v.getRowSize() > 1) {
            throw std::invalid_argument("Invalid vector size");
        }
        double sum = 0;
        for (auto &row: v) {
            for (auto &col: row) {
                sum += pow(col, 2);
            }
        }
        return sqrt(sum);
    }

    static double norm(Matrix<T> &v) {
        return norm(const_cast<const Matrix<T> &>(v));
    }

    static Matrix<T> eye(int m) {
        Matrix<T> new_matrix(m, m);
        for (int i = 0; i < m; i++) {
            new_matrix[i][i] = 1;
        }
        return new_matrix;
    }


    static Matrix<T>
    jacobi_solution(Matrix<T> &A, Matrix<T> &b, int max_iterations, double error_norm_value, double *duration_out) {
        std::clock_t start;
        double duration;
        start = std::clock();

        std::vector<std::vector<double>> x1(A.getColSize(), std::vector<double>(A.getColSize(), 1));
        std::vector<std::vector<double>> temp_x1(A.getColSize(), std::vector<double>(A.getColSize(), 1));
        Matrix<T> x(A.getColSize(), 1);
        Matrix<T> temp_x(A.getColSize(), 1);
        x.matrix = x1;
        temp_x.matrix = temp_x1;
        std::vector<double> error_norms(max_iterations);
        int iteration;
        int size = A.getColSize();
        for (iteration = 0; iteration < max_iterations; iteration++) {
            for (int i = 0; i < size; i++) {
                double temp_sum = 0;
                for (int j = 0; j < size; j++) {
                    if (j == i) {
                        continue;
                    }
                    temp_sum += (A[i][j] * temp_x[j][0]);
                }
                x[i][0] = (b[i][0] - temp_sum) / A[i][i];


            }
            temp_x = x;

            error_norms[iteration] = norm(residuum(A, x, b));
            if (error_norms[iteration] < error_norm_value || std::isinf(error_norms[iteration])) {
                iteration++;
                break;
            }

        }
        duration = static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;
        if (duration_out != nullptr) {
            *duration_out = duration;
        }
        to_csv(error_norms, duration, A[0][0], size, iteration, "jacobi_error_norm");
        std::cout << "Jacobi solution" << std::endl;
        std::cout << "Iterations: " << iteration << std::endl;
        std::cout << "Error norm: " << error_norms[iteration - 1] << std::endl;
        std::cout << "Duration: " << duration << std::endl;
        return x;

    }

    static Matrix<T> gauss_seidel_solution(Matrix<T> &A, Matrix<T> &b, int max_iterations, double error_norm_value,
                                           double *duration_out) {
        std::clock_t start;
        double duration;
        start = std::clock();

        std::vector<std::vector<double>> x1(A.getColSize(), std::vector<double>(A.getColSize(), 1));
        Matrix<T> x(A.getColSize(), 1);
        x.matrix = x1;
        std::vector<double> error_norms(max_iterations);
        int iteration;
        int size = A.getColSize();

        for (iteration = 0; iteration < max_iterations; iteration++) {
            for (int i = 0; i < size; i++) {
                double temp_sum = 0;
                for (int j = 0; j < size; j++) {
                    if (j == i) {
                        continue;
                    }
                    temp_sum += (A[i][j] * x[j][0]);
                }
                x[i][0] = (b[i][0] - temp_sum) / A[i][i];
            }

            error_norms[iteration] = norm(residuum(A, x, b));
            if (error_norms[iteration] < error_norm_value) {
                iteration++;
                break;
            }

        }
        duration = static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;
        if (duration_out != nullptr) {
            *duration_out = duration;
        }
        to_csv(error_norms, duration, A[0][0], size, iteration, "gauss_seidel_error_norm");
        std::cout << "Gauss Seidel solution" << std::endl;
        std::cout << "Iterations: " << iteration << std::endl;
        std::cout << "Error norm: " << error_norms[iteration - 1] << std::endl;
        std::cout << "Duration: " << duration << std::endl;

        return x;

    }

    static Matrix<T> LU_solution(Matrix<T> &A, Matrix<T> &b, double *error_norm_out, double *duration_out) {
        std::clock_t start;
        double duration;
        start = std::clock();

        int m = A.getColSize();
        auto U = A;
        Matrix<T> L(m, m);
        for (int i = 0; i < m; i++) {
            L[i][i] = 1;
        }

        for (int i = 1; i < m; i++) {
            for (int j = 0; j < i; j++) {
                L[i][j] = U[i][j] / U[j][j];
                for (int z = j; z < m; z++) {
                    U[i][z] -= L[i][j] * U[j][z];
                }
            }

        }

        Matrix<T> x(A.getColSize(), 1);
        auto y = x;
        int size = A.col_size;


        // Ly = B
        for (int i = 0; i < size; i++) {
            double temp_sum = 0;
            for (int j = 0; j < i; j++) {
                temp_sum += L[i][j] * y[j][0];
            }
            y[i][0] = b[i][0] - temp_sum;
        }

        // Ux = y
        for (int i = size - 1; i >= 0; i--) {
            double temp_sum = 0;
            for (int j = i; j < size; j++) {
                temp_sum += U[i][j] * x[j][0];
            }
            x[i][0] = (y[i][0] - temp_sum) / U[i][i];
        }
        duration = static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;

        double error_norm = norm(residuum(A, x, b));

        if (error_norm_out != nullptr) {
            *error_norm_out = error_norm;
        }
        if (duration_out != nullptr) {
            *duration_out = duration;
        }
        std::cout << "LU Factorization solution" << std::endl;
        std::cout << "Error norm: " << error_norm << std::endl;
        std::cout << "Duration: " << duration << std::endl;
        return x;
    }


private:
    static void to_csv(std::vector<T> &error_norms, double duration, int a1, int size, int iterations,
                       const std::string &file_name) {
        std::fstream file;
        std::string name = file_name + "_" + std::to_string(a1) + "_" + std::to_string(size) + ".csv";
        file.open(name, std::fstream::out);
        file << "Error norm,Size,Duration\n";
        file << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (int i = 0; i < iterations; i++) {
            file << error_norms[i];
            if (i == 0) {
                file << "," << size << "," << duration << "\n";
            } else if (i < iterations - 1) {
                file << ",,\n";
            }
        }
        file << ",,";
        file.close();
    }


};


#endif //SOLVING_LINEAR_EQUATIONS_MATRIX_H