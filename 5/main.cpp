#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <map>
#include <climits>
#include <vector>
#include <cmath>

using namespace std;

const int N=100000;
const int T=1000;
double delta=0.4,a=1.0,tau=1.0,D=a*a*(1-4*delta*delta)/2.0/tau, v=delta*2.0*a/tau;

struct xt {
    int x,t;
    xt() {};
    xt(int x_, int t_) {
        x=x_;
        t=t_;
    }
};

bool operator < (const xt a, const xt b) {
    if (a.t < b.t) {
        return 1;
    } else {
        if (a.t == b.t) {
            return (a.x < b.x);
        } else {
            return 0;
        }
    }
}

bool operator == (const xt a, const xt b) {
    return (a.t == b.t)  && (a.x == b.x);
}

typedef map <xt,int> mmap;

int main(int argc, char **argv) {
    srand(time(NULL));

    mmap M; // information about count to (x;t)
    FILE* f;

    f=fopen("output.txt","wt");
    for (int n=0;n<N;n++) {
        int c=0;
        printf("%d\n",n+1);
        int x=0;
        for (int t=0;t<=T;t++) {
            mmap::iterator it=M.find(xt(x,t));
            if (it!=M.end()) {
                it->second++; //increase count
            } else {
                M.insert(pair<xt,int>(xt(x,t),1)); //if not exist, we create it
            }
            char tmp=(rand()/(double)INT_MAX<=0.5-delta) ? -1 : +1;
            c+=(tmp==1); //for statistics
            x+=tmp; //move
        }
        if ((n+1) % 100==1) fprintf(f,"%d experiment: p=%lf\n",n+1,c/(double)T);
    }
    fclose(f);

    f=fopen("output.txt","at+");
    for (mmap::iterator it=M.begin();it!=M.end();it++) {
        if (it->first.t==T) fprintf(f,"x=%+lf t=%lf count=%d\n",it->first.x*a,it->first.t*tau,it->second); //print only (x;T;count)
    }
    fclose(f);

    //for all gistogramm
    if (0) {
        f=fopen("output4.txt","wt");
        for (int t0=0;t0<=T;t0++) {
            map <int,int> mass; //(x;count)
            for (mmap::iterator it=M.begin();it!=M.end();it++) {
                if (it->first.t==t0) {
                    mass.insert(pair<int,int>(it->first.x,it->second)); //for fix t0 find all M : M.t=t0 -> mass
                }
            }
            /*
            for (map<int,int>::iterator it=mass.begin();it!=mass.end();it++) {
                fprintf(f,"x=%+lf t=%lf count=%d\n",it->first*a,t0*tau,it->second); //print (x;count) of mass
            }
            */

            map <int,double> n_mass;
            for (int pos=-t0;pos<t0;) {
                int c=0;
                for (int tmp=0;tmp<2 && pos<=t0;tmp++,pos++) {
                    map<int,int>::iterator it=mass.find(pos);
                    if (it!=mass.end()) {
                        c+=it->second;
                    }
                }
                n_mass.insert(pair<int,double>(pos-2,c/(double)N/2.0/a)); //our probablity density for [pos-2;pos]
            }
            if (t0>990*T/1000-1) {
                for (map<int,double>::iterator it=n_mass.begin();it!=n_mass.end();it++) {
                    if (it->second>1e-6) fprintf(f,"%+lf %lf %lf\n", it->first*a, t0*tau, it->second); //x time, func
                }
            }
        }
        fclose(f);
    }


    int t0=T;
    map <int,int> mass; //(x;count)
    for (mmap::iterator it=M.begin();it!=M.end();it++) {
        if (it->first.t==t0) {
            mass.insert(pair<int,int>(it->first.x,it->second)); //for fix t0 find all M : M.t=t0 -> mass
        }
    }
    f=fopen("output2.txt","wt");
    for (map<int,int>::iterator it=mass.begin();it!=mass.end();it++) {
        fprintf(f,"x=%+lf count=%d\n",it->first*a,it->second); //print (x;count) of mass
    }
    fclose(f);

    map <int,double> n_mass;
    for (int pos=-t0;pos<t0;) { //from -T to T-1?
        int c=0;
        for (int tmp=0;tmp<2 && pos<=t0;tmp++,pos++) { //to T?
            map<int,int>::iterator it=mass.find(pos);
            if (it!=mass.end()) {
                c+=it->second;
            }
        }
        n_mass.insert(pair<int,double>(pos-2,c/(double)N/2.0/a)); //our probablity density for [pos-2;pos]
    }

    f=fopen("output3.txt","wt");
    fprintf(f,"x func Delta_func tmp_p\n");
    double max_error=0.0,sq_sum=0.0;
    for (map<int,double>::iterator it=n_mass.begin();it!=n_mass.end();it++) {
        double tmp_p=1.0/sqrt(4*M_PI*D*t0*tau)*exp(-(it->first*a-v*t0*tau)*(it->first*a-v*t0*tau)/(4.0*D*t0*tau));
        if (it->second>1e-5) fprintf(f,"%lf %lf %lf %lf\n", it->first*a, it->second, it->second-tmp_p,tmp_p); //x func delta_func tmp_p
        if (fabs(tmp_p-it->second)>max_error) {
            max_error=fabs(tmp_p-it->second);
        }
    }
    for (map<int,int>::iterator it=mass.begin();it!=mass.end();it++) {
        sq_sum+=it->second*(it->first*a-v*t0*tau)*(it->first*a-v*t0*tau);
    }
    sq_sum/=(double)N;

    fprintf(f,"\t\t%lf\n",max_error);
    fprintf(f,"\t\tsq_sum=%lf 2Dt=%lf fabs(sq_sum-2Dt)=%lf eps=%lf\n",sq_sum, 2*D*t0*tau, fabs(sq_sum-2*D*t0*tau), fabs(sq_sum-2*D*t0*tau)/(2*D*t0*tau));
    fclose(f);

    return 0;
}
