#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

typedef double dbl;

struct solution {
    vector <dbl> x_1, v_1, x_2, v_2, t;
};

dbl f(dbl t, dbl v, dbl x, dbl lambda, dbl w_0, dbl A, dbl Gamma) {
    return -2*lambda*v-(w_0*w_0+A*Gamma*Gamma*cos(Gamma*t))*sin(x);
}

void solve(dbl x_0, dbl lambda, dbl w_0, dbl A, dbl *Gamma, dbl eps, dbl *T, int N, dbl *h, solution *ans ) {
    *Gamma=2*w_0+eps;
    *h=*T/N;
    dbl half_h=*h/2, t;
    ans->x_1.push_back(x_0);
    ans->v_1.push_back(0.0);
    ans->x_2.push_back(x_0);
    ans->v_2.push_back(0.0);
    ans->t.push_back(0.0);
    for (t=0.0;fabs(t-*T)>=1e-7;t+=*h) {
        //1 method
        ans->v_1.push_back(ans->v_1.back()+(*h)*f(t,ans->v_1.back(),ans->x_1.back(), lambda, w_0, A, *Gamma));
        ans->x_1.push_back(ans->x_1.back()+(*h)*ans->v_1.back());
        //2 method
        dbl x_n=ans->x_2.back(),v_n=ans->v_2.back();
        dbl v_=v_n+(*h)*f(t,v_n,x_n,lambda, w_0, A, *Gamma), x_=x_n+(*h)*v_n;
        ans->v_2.push_back(v_n+half_h*(f(t,v_n,x_n,lambda, w_0, A, *Gamma)+f(t+*h,v_,x_,lambda, w_0, A, *Gamma)));
        ans->x_2.push_back(x_n+half_h*(v_n+v_));
        ans->t.push_back(t+*h);
    }
}

int main() {
    dbl lambda, w_0, A, Gamma, x_0;
    int i,N;
    dbl T,t,h;
    x_0=0.1;
    w_0=1.0;
    lambda=0.02;

    N=100000;
    solution ans1, ans2;

    //3 point of problem
    freopen("output1.txt","wt",stdout);
    T=1000.0;
    dbl A_crit=lambda/w_0, A_1=A_crit*1.01, A_2=A_crit*0.99;
    solve(x_0, lambda, w_0, A_1, &Gamma, 0.0, &T, N, &h, &ans1);

    for (int index=0;index<(int)ans1.t.size();index+=10) {
        printf("%lf %lf %lf\n", ans1.t[index], ans1.x_1[index], ans1.x_2[index]);
    }
    /*
    for (int index=0;index<(int)ans1.t.size();index+=10) {
        printf("%lf %lf\n", ans1.x_2[index], ans1.v_2[index]);
    }
    */
    freopen("output2.txt","wt",stdout);
    solve(x_0, lambda, w_0, A_2, &Gamma, 0.0, &T, N, &h, &ans2);
    for (int index=0;index<(int)ans2.t.size();index+=) {
        printf("%lf %lf %lf\n", ans2.t[index], ans2.x_1[index], ans2.x_2[index]);
    }
    /*
    for (int index=0;index<(int)ans2.t.size();index+=10) {
        printf("%lf %lf\n", ans2.x_2[index], ans2.v_2[index]);
    }
    */

    //4 point
    solution ans3, ans4, ans5;
    dbl A_3=A_crit*1.1;
    dbl eps_crit = sqrt(pow(2*w_0*A_3,2.0)-pow(2.0*lambda,2.0));
    freopen("output3.txt","wt",stdout);
    T=700.0;
    solve(x_0,lambda, w_0, A_3, &Gamma, eps_crit*0.9, &T, N, &h, &ans3);
    for (int index=0;index<(int)ans3.t.size();index++) {
        printf("%lf %lf %lf\n",ans3.t[index], ans3.x_1[index], ans3.x_2[index]);
    }
    freopen("output4.txt","wt",stdout);
    T=1000;
    solve(x_0,lambda, w_0, A_3, &Gamma, eps_crit*1.0, &T, 10*N, &h, &ans4);
    for (int index=0;index<(int)ans4.t.size();index++) {
        printf("%lf %lf %lf\n",ans4.t[index], ans4.x_1[index], ans4.x_2[index]);
    }
    freopen("output5.txt","wt",stdout);
    T=700;
    solve(x_0,lambda, w_0, A_3, &Gamma, eps_crit*1.1, &T, 5*N, &h, &ans5);
    for (int index=0;index<(int)ans5.t.size();index++) {
        printf("%lf %lf %lf\n",ans5.t[index], ans5.x_1[index], ans5.x_2[index]);
    }

    solution ans6, ans7, ans8;
    freopen("output6.txt","wt",stdout);
    T=700.0;
    solve(x_0,lambda, w_0, A_3, &Gamma, -eps_crit*0.9, &T, N, &h, &ans6);
    for (int index=0;index<(int)ans6.t.size();index+=10) {
        printf("%lf %lf %lf\n",ans6.t[index], ans6.x_1[index], ans6.x_2[index]);
    }
    freopen("output7.txt","wt",stdout);
    T=1000;
    solve(x_0,lambda, w_0, A_3, &Gamma, -eps_crit*1.0, &T, 10*N, &h, &ans7);
    for (int index=0;index<(int)ans7.t.size();index+=10) {
        printf("%lf %lf %lf\n",ans7.t[index], ans7.x_1[index], ans7.x_2[index]);
    }
    freopen("output8.txt","wt",stdout);
    T=700;
    solve(x_0,lambda, w_0, A_3, &Gamma, -eps_crit*1.1, &T, 5*N, &h, &ans8);
    for (int index=0;index<(int)ans8.t.size();index+=10) {
        printf("%lf %lf %lf\n",ans8.t[index], ans8.x_1[index], ans8.x_2[index]);
    }

    fclose(stdout);

    return 0;
}
