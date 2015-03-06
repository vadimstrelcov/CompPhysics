#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "complex.h"

using namespace std;

const double a=5.0;
const double alpha=0.1;
const int f0=1.0;
const int N=1024;
const double L=16.0;

complex u_func(double x) {
    return complex(exp(-pow(x/a,10.0)),alpha*cos(2*M_PI*f0*x),0);
}

int main() {
    double dx=L/N,x;
    vector<complex> u,u_,W,W_;
    for (x=-L/2.0;x<L/2.0;x+=dx) {
        u.push_back(u_func(x));
    }
    printf("Size u = %d\n",(int)u.size());

    for (int k=0;k<N;k++) {
        complex sum=complex(0.0,0.0,1);
        for (int j=0;j<N;j++) {
            complex tmp=u[j]*complex(1.0,-2.0*M_PI*k*j/N,0);
            sum+=tmp;
        }
        W.push_back(sum/(double)N);
    }
    printf("Size of W = %d\n",(int)W.size());

    for (int k=0;k<N;k++) {
        W_.push_back(W[k]*complex(0.0,(k) ? (k<N/2) ? -1 : +1 : 0,1));
    }
    printf("Size of W_ = %d\n",(int)W_.size());

    for (int j=0;j<N;j++) {
        complex sum=complex(0.0,0.0,1);
        for (int k=0;k<N;k++) {
            complex tmp=W_[k]*complex(1.0,2*M_PI*k*j/N,0);
            sum+=tmp;
        }
        u_.push_back(sum);
    }
    printf("Size of u_ = %d\n",(int)u_.size());

    FILE* f=fopen("output.txt","wt");
    x=-L/2.0;
    for (int i=0;i<N;i++,x+=dx) {
        fprintf(f,"%lf %lf %lf\n", x, u[i].sq_mod(), u_[i].sq_mod());
    }
    fclose(f);
    return 0;
}


