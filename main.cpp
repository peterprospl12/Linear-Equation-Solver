#include "tasks.h"
#include </usr/include/python3.10/Python.h>

int main() {
    Py_Initialize();
    int f = 2;
    int e = 5;
    int N = 963;

    task_one(f, e, N);
    task_two(N, f);
    //task_three(f, e);

    Py_Finalize();
    return 0;
}