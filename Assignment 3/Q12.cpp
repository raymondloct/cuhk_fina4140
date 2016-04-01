#include <iostream>
#include <cstdio>
#include <cstdlib>
using namespace std;

int *LCG (int a, int b, int M, int N0){
    int *N = (int *) malloc(sizeof(int)*(M+4));
    double *U = (double *) malloc(sizeof(double)*(M+4));
    int P, start, end;
    double mean, var, sum;
    FILE *numOut;
    
    // Check initial conditions
    if (N0 == 0) {
        cout << "N_0 = 0!" << endl;
        exit(1);
    } else if (a == 0) {
        cout << "a = 0!" << endl;
        exit(1);
    }
    
    // Generate random numbers
    N[0] = N0;
    for (int i = 1; i < M+4; i++) {
        N[i] = (N[i-1] * a + b) % M;
        U[i] = (double) N[i] / (double) M;
    }
    
    // Find the period
    start = 1;
    end = 2;
    while (N[start] != N[end]) {
        if (end > M) {
            start++;
            end = start;
        }
        end++;
    }
    P = end - start;
    
    // Find mean
    sum = 0;
    for (int i = start; i < end; i++)
        sum += U[i];
    mean = sum / (double) P;
    
    // Find variance
    sum = 0;
    for (int i = start; i < end; i++)
        sum += (U[i]-mean)*(U[i]-mean);
    var = sum / (double) P;
    
    // Print out the results
    cout << "For a = " << a << ", b = " << b << \
            ", M = " << M << ", N0 = " << N0 << "," << endl;
    cout << "P = " << P << endl;
    cout << "mean = " << mean << endl;
    cout << "var = " << var << endl;
    cout << "The generated numbers are outputted in a separate file named Q12.txt.";
    
    numOut = fopen("Q12.txt","w");
    for (int i = start; i < end+3; i++)
        fprintf(numOut, "%d\n", N[i]);
    fclose(numOut);
    
    free(U);
    return N;
}

int main() {
    int *randomNumbers = LCG(1234, 2345, 1000, 11422);
    free(randomNumbers);
    
    return 0;
}
