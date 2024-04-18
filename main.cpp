#include "tasks.h"
#include </usr/include/python3.10/Python.h>

int main() {
    Py_Initialize();
    int f = 2;
    int e = 5;
    int N = 994;

    task_B(f, e, N);
    task_C(N, f);
    task_E(f, e);

    Py_Finalize();
    return 0;
}