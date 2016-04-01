#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

double *Uran (int n) {
    double *U = (double *)malloc(sizeof(double)*n);
    
    for (int i = 0; i < n; i++)
        U[i] = (double) (rand() % RAND_MAX) / (double) RAND_MAX;
    
    return U;
}

void BoxMuller1(int n){
    double *U, theta, rho, Z[2];
    
    if (n % 2 != 0) n++;
    
    for (int i = 0; i < n/2; i++) {
        U = Uran(2);
        theta = 2 * M_PI * U[1]; // M_PI is constant for pi
        rho = sqrt(-2*log(U[0]));
        Z[0] = rho * cos(theta);
        Z[1] = rho * sin(theta);
        free(U);
    }
}

void BoxMuller2(int n){
    double *U, V[2], W, Z[2];
    
    if (n % 2 != 0) n++;
    for (int i = 0; i < n/2; i++) {
        do {
            U = Uran(2);
            V[0] = U[0] * 2 - 1;
            V[1] = U[1] * 2 - 1;
            W = V[0] * V[0] + V[1] * V[1];
            free(U);
        } while (W >= 1 || W <= 0);
                
        Z[0] = V[0] * sqrt(-2*log(W)/W);
        Z[1] = V[1] * sqrt(-2*log(W)/W);
    }
}

int power(int n, int m) {
    int r = 1;
    
    while (m > 0) {
        if (m % 2 == 1)
            r *= n;
        n *= n;
        m /= 2;
    }
    return r;
}

int main() {
    double timeElapsed[2][5];
    FILE *numOut;
	
    srand(time(NULL));
    
    cout << "Box-Muller method: " << endl;
    for (int n = 2; n <= 6; n++) {
        clock_t start = clock();
        BoxMuller1(power(10,n));
        clock_t end = clock();
        timeElapsed[0][n-2] = (double)(end - start) / (double) CLOCKS_PER_SEC;
        cout << "Time elapsed for n = " << n << ", in seconds, is " << timeElapsed[0][n-2] << endl;
    }
    
    cout << "Marsaglia polar method: " << endl;
    for (int n = 2; n <= 6; n++) {
        clock_t start = clock();
        BoxMuller2(power(10,n));
        clock_t end = clock();
        timeElapsed[1][n-2] = (double)(end - start) / (double) CLOCKS_PER_SEC;
        cout << "Time elapsed for n = " << n << ", in seconds, is " << timeElapsed[1][n-2] << endl;
    }
    
	cout << "The time elapsed is also outputted to Q13.txt.";
	
	numOut = fopen("Q13.txt","w");
    for (int i = 0; i <= 4; i++)
        fprintf(numOut, "%d\t%f\t%f\n", i+2, timeElapsed[0][i], timeElapsed[1][i]);
    fclose(numOut);
	
    return 0;
}
