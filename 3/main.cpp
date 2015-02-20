#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

typedef double dbl;

dbl h, N;

vector<dbl> a, b, c;


struct answer {
    vector<dbl> x, a_ans, r_ans;
    void cl() {
        x.clear();
        a_ans.clear();
        r_ans.clear();
    }
};

dbl f(dbl x, bool func) {
    return (!func) ? 1.0 : sin(M_PI*x);
}

void solve(answer *ans, int N, dbl h, bool func, bool accur) {
    /*
        func = 0 => f(x)=1;
        func = 1 => f(x)=sin(pi*x);
        accur = 0 => u_1 = u_0+u'(0)*h
        accur = 1 => u_1 = u_0+u'(0)*h+1/2*u''(0)*h^2
    */
    vector<dbl> a, b, c, F;
    a.push_back(0.0); //doesn't matter
    b.push_back(-1.0);
    c.push_back(1.0+h);
    F.push_back( (!accur) ? 0.0 : +0.5*f(0.0,func)*h*h );
    for (int i=1;i<N;i++) { //to N-1!
        a.push_back(1.0);
        c.push_back(-2.0);
        b.push_back(1.0);
        F.push_back(-f(i*h,func)*h*h);
    }
    a.push_back(1.0);
    c.push_back(-(1.0+h));
    b.push_back(0.0); //doesn't matter
    F.push_back( (!accur) ? 0.0 : -0.5*f(1.0,func)*h*h); //minus?
    vector<dbl> c_, F_;
    c_.push_back(c[0]);
    F_.push_back(F[0]);
    for (int i=1;i<=N;i++) {
        c_.push_back(c[i]-a[i]*b[i-1]/c_[i-1]);
    }
    for (int i=1;i<=N;i++) {
        F_.push_back(F[i]-a[i]*F_[i-1]/c_[i-1]);
    }

    ans->a_ans.resize(N+1);
    ans->a_ans[N]=F_[N]/c_[N];
    for (int i=N-1;i>=0;i--) {
        ans->a_ans[i]=(F_[i]-b[i]*ans->a_ans[i+1])/c_[i];
    }
    for (int i=0;i<=N;i++) {
        ans->x.push_back(i*h);
    }
    //right answer
    for (int i=0;i<=N;i++) {
        ans->r_ans.push_back( (!func) ? ( 0.5*(-pow(ans->x[i],2.0)+ans->x[i]+1.0) ) : 1.0/M_PI*(1.0+sin(M_PI*ans->x[i])/M_PI));
    }

}

int main() {
    N=1000;
    h=1.0/N;
    printf("%lf\n",h);
    answer ans1, ans2, ans3, ans4;
    solve(&ans1, N, h, 0, 0); //for f(x)=1, u_1=u_0+u'h
    freopen("output1.txt","wt",stdout);
    for (int i=0;i<=N;i++) {
        printf("%lf %.10lf %.10lf %.10lf\n",ans1.x[i],ans1.a_ans[i],ans1.r_ans[i],fabs(ans1.r_ans[i]-ans1.a_ans[i]));
    }
    answer ans_tmp;
    solve(&ans_tmp, 2*N, 0.5*h, 0, 0);
    printf("%lf\n",log( fabs(ans1.a_ans[0]-ans1.r_ans[0])/fabs(ans_tmp.a_ans[0]-ans_tmp.r_ans[0]) )/log(2.0));
    ans_tmp.cl();

    solve(&ans2, N, h, 0, 1); //for f(x)=1, u_1=u_0+u'h+0.5u''h^2
    freopen("output2.txt","wt",stdout);
    for (int i=0;i<=N;i++) {
        printf("%lf %.14lf %.14lf %.14lf\n",ans2.x[i],ans2.a_ans[i],ans2.r_ans[i],fabs(ans2.r_ans[i]-ans2.a_ans[i]));
    }
    solve(&ans_tmp, 2*N, 0.5*h, 0, 1);
    printf("%lf\n",log( fabs(ans2.a_ans[0]-ans2.r_ans[0])/fabs(ans_tmp.a_ans[0]-ans_tmp.r_ans[0]) )/log(2.0));
    ans_tmp.cl();


    solve(&ans3, N, h, 1, 0); //for f(x)=sin(pi*x), u_1=u_0+u'h
    freopen("output3.txt","wt",stdout);
    for (int i=0;i<=N;i++) {
        printf("%lf %.14lf %.14lf %.14lf\n",ans3.x[i],ans3.a_ans[i],ans3.r_ans[i],fabs(ans3.r_ans[i]-ans3.a_ans[i]));
    }
    solve(&ans_tmp, 2*N, 0.5*h, 1, 0);
    printf("%lf\n",log( fabs(ans3.a_ans[0]-ans3.r_ans[0])/fabs(ans_tmp.a_ans[0]-ans_tmp.r_ans[0]) )/log(2.0));
    ans_tmp.cl();

    solve(&ans4, N, h, 1, 1); //for f(x)=sin(pi*x), u_1=u_0+u'h+0.5u''h^2
    freopen("output4.txt","wt",stdout);
    for (int i=0;i<=N;i++) {
        printf("%lf %.14lf %.14lf %.14lf\n",ans4.x[i],ans4.a_ans[i],ans4.r_ans[i],fabs(ans4.r_ans[i]-ans4.a_ans[i]));
    }
    solve(&ans_tmp, 2*N, 0.5*h, 1, 1);
    printf("%lf\n",log( fabs(ans4.a_ans[0]-ans4.r_ans[0])/fabs(ans_tmp.a_ans[0]-ans_tmp.r_ans[0]) )/log(2.0));
    ans_tmp.cl();

    return 0;
}
