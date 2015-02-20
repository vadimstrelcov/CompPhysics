#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

typedef double dbl;

dbl lambda, x_0;
int N;

dbl f(dbl t, dbl v, dbl x) {
    return -2*lambda*v-lambda*lambda*x;
}

void solve(dbl x_0, dbl lambda, int N, vector <dbl> *x_1, vector <dbl>*v_1, vector <dbl> *x_2, vector <dbl> *v_2,
           vector <dbl> *r, dbl *T, dbl *h, dbl *max_err_1, dbl *max_err_2) {
    *T=3.0/lambda;
    *h=*T/N;
    dbl half_h=*h/2.0, t;
    x_1->push_back(x_0);
    v_1->push_back(0.0);
    x_2->push_back(x_0);
    v_2->push_back(0.0);
    for (t=0.0;fabs(t-*T-(*h))>1e-7;t+=*h) {
        //1 method
        v_1->push_back(v_1->back()+(*h)*f(t,v_1->back(),x_1->back()));
        x_1->push_back(x_1->back()+(*h)*v_1->back());
        //2 method
        dbl x_n=x_2->back(),v_n=v_2->back();
        dbl v_=v_n+(*h)*f(t,v_n,x_n), x_=x_n+(*h)*v_n;
        v_2->push_back(v_n+half_h*(f(t,v_n,x_n)+f(t+(*h),v_,x_)));
        x_2->push_back(x_n+half_h*(v_n+v_));
        //right solution
        r->push_back(x_0*exp(-lambda*t)*(1.0+lambda*t));
    }

    x_1->pop_back();
    v_1->pop_back();
    x_2->pop_back();
    v_2->pop_back();

    *max_err_1=0.0;
    *max_err_2=0.0;
    int i;
    for (t=0.0, i=0;fabs(t-*T-(*h))>=1e-7;t+=(*h), i++) {
        dbl tmp_1=fabs((*r)[i]-(*x_1)[i]), tmp_2=fabs((*r)[i]-(*x_2)[i]);
        *max_err_1=(*max_err_1>tmp_1) ? *max_err_1 : tmp_1;
        *max_err_2=(*max_err_2>tmp_2) ? *max_err_2 : tmp_2;
    }
}

int main() {
    /*
    scanf("%lf %lf",&x_0,&lambda);
    getchar();
    */
    int i;
    dbl T,t,h;
    x_0=1.0;
    lambda=1.0;
    N=250;

    freopen("output1.txt","wt",stdout);
    vector <dbl> x_1,v_1,x_2,v_2,x_3,v_3,x_4,v_4,r_1,r_2;
    dbl max_err_1, max_err_2, max_err_3, max_err_4;
    solve(x_0,lambda,N,&x_1,&v_1,&x_2,&v_2,&r_1,&T,&h,&max_err_1,&max_err_2);

    for (t=0.0, i=0;fabs(t-T-h)>=1e-7;t+=h, i++) {
        printf("%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", t, fabs(r_1[i]-x_1[i]), fabs(r_1[i]-x_2[i]), x_1[i], x_2[i], r_1[i] );
    }
    //printf("error : %lf %lf\n", max_err_1, max_err_2);
    fclose(stdout);

    //(x,x')
    freopen("output3.txt","wt",stdout);
    for (t=0.0, i=0;fabs(t-T-h)>=1e-7;t+=h, i++) {
        printf("%.10lf %.10lf %.10lf %.10lf\n",x_1[i],v_1[i],x_2[i],v_2[i]);
    }
    fclose(stdout);

    freopen("output2.txt","wt",stdout);
    solve(x_0,lambda,2*N,&x_3,&v_3,&x_4,&v_4,&r_2,&T,&h,&max_err_3,&max_err_4);
    for (t=0.0, i=0;fabs(t-T-h)>=1e-7;t+=h, i++) {
        printf("%.10lf %.10lf %.10lf\n", t, fabs(r_2[i]-x_3[i]), fabs(r_2[i]-x_4[i]) );
    }
    printf("error : %lf %lf\n", max_err_3, max_err_4);

    i--;
    printf("1 net : 1 method - %.10lf %.10lf %.10lf\n",x_1.back(), r_1.back(),fabs(r_1.back()-x_1.back()));
    printf("1 net : 2 method - %.10lf %.10lf %.10lf\n",x_2.back(), r_1.back(),fabs(r_1.back()-x_2.back()));
    printf("2 net : 1 method - %.10lf %.10lf %.10lf\n",x_3.back(), r_2.back(),fabs(r_2.back()-x_3.back()));
    printf("2 net : 2 method - %.10lf %.10lf %.10lf\n",x_4.back(), r_2.back(),fabs(r_2.back()-x_4.back()));
    printf("1 method Runge's eps = %.10lf\n",fabs(x_3.back()-x_1.back()));
    printf("2 method Runge's eps = %.10lf\n",fabs(x_4.back()-x_2.back())/3.0);
    /*
    printf("1-3 : eps = %lf, eps = %lf, Runge's eps = %lf\n",fabs(r_1.back()-x_1.back()),fabs(r_1.back()-x_2.back()), 1/1.0*fabs(x_3.back()-x_1.back()) ); //1 method for h, h/2
    printf("2-4 : eps = %lf, Runge's eps = %lf\n",fabs(r_1.back()-x_4.back()), 1/3.0*fabs(x_4.back()-x_2.back()) ); //2 method for h, h/2
    */
    fclose(stdout);

    //(x,x')
    freopen("output4.txt","wt",stdout);
    for (t=0.0, i=0;fabs(t-T-h)>=1e-7;t+=h, i++) {
        printf("%.10lf %.10lf %.10lf %.10lf\n",x_3[i],v_3[i],x_4[i],v_4[i]);
    }
    fclose(stdout);

    return 0;
}
