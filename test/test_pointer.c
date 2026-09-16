#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    double *a = malloc(2 * sizeof(double));
    double *b = malloc(2 * sizeof(double));
    double *c = malloc(2 * sizeof(double));
    double *pa = a; // Pointer to the first element of array a
    double *pb = b; // Pointer to the first element of array b
    for (int i = 0; i < 2; i++) {
        a[i] = i + 1.0; // Initialize array a with values 1.0 and 2.0
        b[i] = (i + 1) * 10.0; // Initialize array b with values 10.0 and 20.0
    }

    for (int i = 0; i < 2; i++) {
        printf("Value of a[%d]: %f\n", i, a[i]); // Output: 11.0 and 22.0
        printf("Value of b[%d]: %f\n", i, b[i]); // Output: 11.0 and 22.0
    }
    double *temp = pa; // temp points to the same location as pa (which is the first element of array a)
    pa = pb; // pa now points to the first element of array b
    pb = temp; // pb now points to the first element of array a
    for (int i = 0; i < 2; i++) {
        printf("Value of a[%d]: %f\n", i, pa[i]); // Output: 1.0 and 2.0
        printf("Value of b[%d]: %f\n", i, pb[i]); // Output: 10.0 and 20.0
    }
    free(a);
    free(b);
    free(c);
    return 0;
}