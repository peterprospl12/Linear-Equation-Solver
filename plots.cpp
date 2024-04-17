//
// Created by piotr on 4/17/24.
//

#include <string>
#include "plots.h"
#include </usr/include/python3.10/Python.h>


void plot_error_norm(const std::string& first_file, const std::string& second_file, const std::string& plot_name) {
    char python_code[1000];
    snprintf(python_code, sizeof(python_code), R"(
try:
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    df = pd.read_csv('%s')

    jacobi_error_norm = df['Error norm']
    jacobi_iterations = len(jacobi_error_norm)

    df = pd.read_csv('%s')
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

    plt.savefig('%s', dpi=300)
except Exception as e:
    print(str(e))
)", first_file.c_str(), second_file.c_str(), plot_name.c_str());

    PyRun_SimpleString(python_code);
}