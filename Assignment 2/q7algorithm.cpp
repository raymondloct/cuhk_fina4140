#include <iostream>
#include <iomanip> // For setprecision, fixed
#include <cstdlib> // For malloc, free
#include <cmath> // For abs, pow, exp
using namespace std;
enum OpType {EUROCALL, EUROPUT, AMERCALL, AMERPUT}; // European or American, call or put
enum PayType {CASH, ASSET, DIFF}; // Binary cash-or-nothing, binary asset-or-nothing, ordinary option
double ** malloc2x2(int x, int y); // Allocate memory for two-dimensional array
void free2x2(double ** ptr, int x, int y); // Free memory allocated by malloc2x2
double OptionPayoff(enum OpType oType, enum PayType pType, double S, double K, double Q);

double ** malloc2x2(int x, int y) {
// Allocate memory for two-dimensional array of type double
    double **arr = (double**) malloc(x * sizeof(double*));
    
    if (arr == NULL) {
        cout << "Cannot allocate memory for " << x << "x" << y << "array!";
        exit(1);
    }
    
    for (int i = 0; i < x; i++) {
        arr[i] = (double *) malloc(y * sizeof(double));
        if (arr[i] == NULL) {
            cout << "Cannot allocate memory for " << x << "x" << y << "array!";
            exit(1);
        }
    }
    
    return arr;
}

void free2x2(double **ptr, int x, int y) {
// Free memory allocated by malloc2x2
    for (int i = 0; i < x; i++)
        free(ptr[i]);
    free(ptr);
}

double BinomialMethod(double r, double sigma, double S0, double T, double K, double bQ, \
int M, enum OpType optiontype, enum PayType payofftype) {
    double dt = T / M;
    double u, d, q, result;
    double **S = malloc2x2(M+1,M+1), **V = malloc2x2(M+1,M+1); // Allocate memory for S and V
    
    S[0][0] = S0; // Means S00
    
    // Start Parameterization
    u = exp(sigma*pow(dt,0.5));
    d = exp(-sigma*pow(dt,0.5));
    q = (exp(r*dt)-d) / (u-d);
    // End Parameterization
    
    
    for (int n = 0; n <= M; n++) { // At maturity
        S[n][M] = S[0][0] * pow(u,n) * pow(d,M-n);
        V[n][M] = OptionPayoff(optiontype,payofftype,S[n][M],K,bQ);
    }
    
    if (optiontype == EUROCALL || optiontype == EUROPUT) // European option
        for (int i = M - 1; i >= 0; i--)
            for (int n = i; n >= 0; n--)
                V[n][i] = exp(-r*dt) * (q * V[n+1][i+1] + (1 - q) * V[n][i+1]);
    else // American option
        for (int i = M - 1; i >= 0; i--)
            for (int n = i; n >= 0; n--) {
                S[n][i] = S[0][0] * pow(u,n) * pow(d,i-n);
                V[n][i] = fmax(OptionPayoff(optiontype,payofftype,S[n][i],K,bQ), \
                exp(-r * dt) * (q * V[n+1][i+1] + (1 - q) * V[n][i+1]));
            }
    
    result = V[0][0]; // Which is approximation to V(S0,0)
    free2x2(S, M+1, M+1);
    free2x2(V, M+1, M+1);
    return result;
}

// Calculate payoff at a particular point of time and given state
double OptionPayoff(enum OpType oType, enum PayType pType, double S, double K, double Q) {
    if ((oType == EUROCALL || oType == AMERCALL) && pType == DIFF)
        return fmax(S - K, 0);
    else if ((oType == EUROPUT || oType == AMERPUT) && pType == DIFF)
        return fmax(K - S, 0);
    else if ((oType == EUROCALL || oType == AMERCALL) && pType == CASH)
        return S > K ? Q : 0;
    else if ((oType == EUROPUT || oType == AMERPUT) && pType == CASH)
        return S < K ? Q : 0;
    else if ((oType == EUROCALL || oType == AMERCALL) && pType == ASSET)
        return S > K ? S : 0;
    else if ((oType == EUROPUT || oType == AMERPUT) && pType == ASSET)
        return S < K ? S : 0;
    else {
        cout << "Error in OptionPayoff!" << endl;
        return -1;
    }
}

int main(void) {
    double answerA, answerB, answerC, answerD, answerE[4];
    double r = 0.06, sigma = 0.3, S0 = 5, T = 1, K = 10, Q = 0;
    int M = 1024;
    
    cout << fixed; // For setprecision()
    
    answerA = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROPUT,DIFF);
    cout << "Part a: " << setprecision(6) << answerA << endl;
    
    S0 = 9;
    answerB = BinomialMethod(r,sigma,S0,T,K,Q,M,AMERPUT,DIFF);
    cout << "Part b: " << setprecision(6) << answerB << endl;
    
    S0 = 5;
    answerC = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROCALL,DIFF);
    cout << "Part c: " << setprecision(6) << answerC << endl;
    
    cout << "Part d: " << endl;
    for (M = 100; M <= 150; M++) {
        answerD = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROPUT,DIFF);
        cout << "M = " << M << ": " << setprecision(6) << answerD << ", ";
        cout << "error = " << setprecision(9) << abs(answerD - 4.430465) << endl;
    }
    
    cout << endl;
    cout << "Part e: " << endl;
    Q = 1;
    M = 1024;
    answerE[0] = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROCALL,CASH);
    cout << "Cash-or-nothing Call: " << setprecision(6) << answerE[0] << endl;
    answerE[1] = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROPUT,CASH);
    cout << "Cash-or-nothing Put: " << setprecision(6) << answerE[1] << endl;
    answerE[2] = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROCALL,ASSET);
    cout << "Asset-or-nothing Call: " << setprecision(6) << answerE[2] << endl;
    answerE[3] = BinomialMethod(r,sigma,S0,T,K,Q,M,EUROPUT,ASSET);
    cout << "Asset-or-nothing Put: " << setprecision(6) << answerE[3] << endl;
    
    return 0;
}
