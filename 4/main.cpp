#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;

typedef double dbl;

class vect {
public:
    dbl x,y,v_x,v_y;
    vect() {}
    vect(dbl X, dbl Y, dbl V_X, dbl V_Y) {
        x=X;
        y=Y;
        v_x=V_X;
        v_y=V_Y;
    }
};
vect operator +(vect P, vect Q) {
    return vect(P.x+Q.x,P.y+Q.y,P.v_x+Q.v_x,P.v_y+Q.v_y);
}
vect operator -(vect P, vect Q) {
    return vect(P.x-Q.x,P.y-Q.y,P.v_x-Q.v_x,P.v_y-Q.v_y);
}
vect operator *(vect P, dbl A) {
    return vect(P.x*A,P.y*A,P.v_x*A,P.v_y*A);
}
vect operator *(dbl A, vect P) {
    return vect(P.x*A,P.y*A,P.v_x*A,P.v_y*A);
}
vect operator /(vect P, dbl A) { //A!=0.0 !
    return vect(P.x/A,P.y/A,P.v_x/A,P.v_y/A);
}

dbl r_(vect Q) {
    return sqrt(pow(Q.x,2.0)+pow(Q.y,2.0));
}

const dbl R=6700000.0;
const dbl d=5000.0;
const dbl u=50.0;
const dbl gM=3.99533*pow(10.0,14.0);
const dbl v_c=sqrt(gM/R);
const dbl w=v_c/R;
const dbl phi_0=d/R;
const dbl TRANSF=180/M_PI;

struct answer {
    vector<vect> U;
    vector<dbl> t;
    void cl() {
        U.clear();
        t.clear();
    }
    void cp(answer *tmp) {
        for (int index=0;index<(int)tmp->U.size();index++) {
            U.push_back(tmp->U[index]);
            t.push_back(tmp->t[index]);
        }
    }
};

int is_outside(vect *Q, dbl t) {
    vect P=vect(R*cos(w*t+phi_0),R*sin(w*t+phi_0),0.0,0.0); //speed is unnessery
    dbl eps_=1E-3;
    if (r_(*Q)>R+eps_) {
        //out of orbit
        if (Q->x<P.x) return -2; //lefter
        else return -1; //righter
    }
    if (r_(*Q-P)<eps_ && Q->x>P.x) {
        return 1; //succsess
    } else {
        return 0; //in orbid, not success yet
    }
}

vect F(vect U, dbl t) {
    vect res;
    res.x=U.v_x;
    res.y=U.v_y;
    res.v_x=-gM/pow(r_(U),3.0)*U.x;
    res.v_y=-gM/pow(r_(U),3.0)*U.y;
    return res;
}

void solve(dbl *theta, dbl *h_time, answer *ans, int *result, int depth_dechotomy, int number_file=-1) {
    ans->t.push_back(0.0);
    ans->U.push_back(vect(R,0.0,-u*sin(*theta),v_c+u*cos(*theta)));
    dbl tm=0.0, half_time=*h_time*0.5;
    while (1) {
        //check
        int res = is_outside(&ans->U.back(),tm);
        if (res!=0) {
            *result=res;
            break;
        }
        //continue analysis
        vect k1, k2, k3, k4;
        k1=F(ans->U.back(),tm)*(*h_time);
        k2=F(ans->U.back()+k1*0.5,tm+half_time)*(*h_time);
        k3=F(ans->U.back()+k2*0.5,tm+half_time)*(*h_time);
        k4=F(ans->U.back()+k3,tm+*h_time)*(*h_time);
        ans->U.push_back(ans->U.back()+(k1+2*k2+2*k3+k4)/6.0);
        tm+=*h_time;
        ans->t.push_back(tm);
    }

    if (number_file!=-1 && (depth_dechotomy-1==number_file || 0)) {
        if (r_(ans->U.back())>R) {
            ans->U.pop_back();
            ans->t.pop_back();
        }
        FILE *file;
        char s[13]="output";
        if (number_file<10) {
            s[6]=number_file+'0';
            s[7]='.'; s[8]='t'; s[9]='x'; s[10]='t'; s[11]='\0'; s[12]='\0';
        } else {
            s[6]=number_file/10+'0'; s[7]=number_file%10+'0';
            s[8]='.'; s[9]='t'; s[10]='x'; s[11]='t'; s[12]='\0';
        }
        file=fopen(s,"wt");
        fprintf(file,"     %lf\n",*theta*TRANSF);
        for (int index=0;index<(int)ans->U.size();index+=(int)ans->U.size()/1000) {
            fprintf(file,"%.0lf %.0lf %.0lf %.0lf %lf\n",ans->U[index].x, ans->U[index].y,
                R*cos(w*ans->t[index]+phi_0), R*sin(w*ans->t[index]+phi_0), r_(ans->U[index]));
        }
        fclose(file);
    }

}

int main() {
    dbl epsilon_theta=1E-10, h_time=1E-4;
    dbl left_lim=0.0, right_lim=0.3, middle=(left_lim+right_lim)*0.5;
    int depth_dechotomy = (int)( log((right_lim-left_lim)*0.5/epsilon_theta)/log(2.0) ) +1;
    answer ans_l, ans_r, ans_m;
    int f_left_lim, f_right_lim, f_middle;
    solve(&left_lim,&h_time,&ans_l,&f_left_lim,depth_dechotomy);
    solve(&right_lim,&h_time,&ans_r,&f_right_lim,depth_dechotomy);
    printf("%d\n%lg\n%lg\n%lg\n%lg\n%lg\n%lg\n",depth_dechotomy,R,d,u,gM,v_c,w);
    for (int depth=0;depth<depth_dechotomy;depth++) {
        middle=(left_lim+right_lim)*0.5;
        solve(&middle,&h_time,&ans_m,&f_middle,depth_dechotomy,depth);
        //there is somewhere 1
        if (f_left_lim==1 || f_right_lim==1 || f_middle==1) {
            if (f_left_lim==1) {
                right_lim=middle;
            } else {
                if (f_right_lim==1) {
                    left_lim=middle;
                } else {
                    //nothiing :)
                }
            }
        } else {
            printf("depth=%d %02.10lf %02.10lf %02.10lf %d %d %d %lf\n",depth,left_lim*TRANSF,middle*TRANSF,
                   right_lim*TRANSF, f_left_lim, f_middle, f_right_lim, ans_m.t.back());
            if (f_left_lim==f_middle) {
                ans_l.cl();
                ans_l.cp(&ans_m);
                f_left_lim=f_middle;
                left_lim=middle;
            } else {
                ans_r.cl();
                ans_r.cp(&ans_m);
                f_right_lim=f_middle;
                right_lim=middle;
            }
        }
        if (depth!=depth_dechotomy-1) ans_m.cl();
    }
    printf("theta=%.10lf±%.10lf | %d %d %d %lf\n",middle*TRANSF, epsilon_theta*TRANSF, f_left_lim==-2,
            f_middle, f_right_lim, ans_m.t.back());
    printf("all time = %lf\n",ans_m.t.back());
    dbl min_r=R;
    for (int index=0;index<(int)ans_m.U.size();index++) {
        dbl tmp_r=r_(ans_m.U[index]);
        if (tmp_r<min_r) min_r=tmp_r;
    }
    printf("min_r=%lf\n",min_r);
    vect Ve=vect(0.0,0.0,-R*w*sin(w*ans_m.t.back()+phi_0),R*w*cos(w*ans_m.t.back())+phi_0);
    vect Va=ans_m.U.back(), Vr=Va-Ve;
    printf("Vr_x = %lf, Vr_y = %lf Vr=%lf\n",Vr.v_x,Vr.v_y,sqrt(pow(Vr.v_x,2.0)+pow(Vr.v_y,2.0)));
    dbl tmp=fabs(Ve.v_x*Va.v_y-Va.v_x*Ve.v_y), l_Va=sqrt(pow(Va.v_x,2.0)+pow(Va.v_y,2.0)), l_Ve=sqrt(pow(Ve.v_x,2.0)+pow(Ve.v_y,2.0));
    printf("%lg\n",asin(tmp/l_Va/l_Ve)*TRANSF);

    return 0;
}
