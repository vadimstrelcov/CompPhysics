#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>

using namespace std;

const int N_exp=1000000;
const int N=20;
bool show=0;

int get_count(vector<bool>* s) {
    int sum=0;
    for (int i=0;i<(int)s->size();i++) {
        sum+=(*s)[i];
    }
    return sum;
}

int f(int x, vector<int>* np) {
	if ( (*np)[x] == 0) return x;
	return f( (*np)[x] ,np);
}

void get_res(int* res, double p_) {
    for (int temp=0;temp<N_exp;temp++) {
        bool a[N+1][N+1]/*={ {0,0,0,0,0,0,0,0},
			{0,1,1,0,1,1,1,0},
			{0,1,0,0,0,1,1,1},
			{0,0,1,1,0,0,1,1},
			{0,1,1,1,1,1,1,1},
			{0,0,0,1,0,1,1,1},
			{0,1,1,0,0,1,1,1},
			{0,1,0,1,1,0,1,0}}*/;

        int i,j,c[N+1][N+1]={0},n=0,index=0;
        vector<int> np;
        np.push_back(0);
        for (i=0;i<=N;i++) {
			a[i][0]=0;
			a[0][i]=0;
		}

		for (i=1;i<=N;i++) {
			for (j=1;j<=N;j++) {
				a[i][j]=( rand()/(double)RAND_MAX < 1.0-p_ ) ? 0 : 1;
			}
		}

		if (show) {
			printf("A massive\n");
			for (i=1;i<=N;i++) {
				for (j=1;j<=N;j++) {
					printf("%d ",a[i][j]);
				}
				if (i==N) {
					int sq=0;
					for (int p=1;p<=N;p++) {
						for (int q=1;q<=N;q++) {
							sq+=a[p][q];
						}
					}
					printf("(%.3lf)",sq/(double)N/N);
				}
				printf("\n");
			}
		}
        //create c & np massive
		for (i=1;i<=N;i++) {
			for (j=1;j<=N;j++) {
				if (a[i][j]) {
					if (!a[i-1][j] && !a[i][j-1]) {
						n++;
						c[i][j]=n;
						np.push_back(0);
					}
					if (!a[i-1][j] && a[i][j-1]) {
						c[i][j]=c[i][j-1];
					}
					if (a[i-1][j] && !a[i][j-1]) {
						c[i][j]=c[i-1][j];
					}
					if (a[i-1][j] && a[i][j-1]) {
						if (c[i-1][j]==c[i][j-1]) {
							c[i][j]=c[i][j-1];
						} else {
                            c[i][j]=min(f(c[i-1][j],&np),f(c[i][j-1],&np));
                            np[max(c[i-1][j],c[i][j-1])]=c[i][j];
						}
					}
				}
			}
		}
		//correction
		for (int t=1;t<(int)np.size();t++) {
			if (np[t]!=0) {
				for (i=1;i<=N;i++) {
					for (j=1;j<=N;j++) {
						if (c[i][j]==t) {
							c[i][j]=f(np[t],&np);
						}
					}
				}
			}
		}
		if (show) {
			printf("C correction massive\n");
			for (i=1;i<=N;i++) {
				for (j=1;j<=N;j++) {
					printf("%2d ",c[i][j]);
				}
				printf("\n");
			}
		}

		//check result
		for (int j1=1;j1<=N;j1++) {
			if (c[1][j1]>0) {
				for (int j2=1;j2<=N;j2++) {
					if (c[N][j2]>0) {
						if (c[1][j1]==c[N][j2]) {
							index=c[1][j1];
							break;
						}
					}
				}
			}
			if (index>0) {
				break;
			}
		}
		(*res)+=(index>0);

		if (show) {
			printf("res = %d\n",index>0);
		}
	}
}

int main() {

    srand(time(NULL));

    double p_min=0.5,p_max=0.8,p_mid;
    int res_min=0,res_max=0,res_mid=0;
    int cou=0;

    while ( cou<10000 && fabs(p_max-p_min)>1e-4 ) {
        printf("res_min p=%lf\n",p_min);
        if (res_min==0) {
            get_res(&res_min,p_min);
        }
        printf("res_max p=%lf\n",p_max);
        if (res_max==0) {
            get_res(&res_max,p_max);
        }
        p_mid=(p_max+p_min)/2.0;
        printf("res_mid p=%lf\n",p_mid);
        if (res_mid==0) {
            get_res(&res_mid,p_mid);
        }
        printf("count=%d stat={min=%lf %d (%.3lf); max=%lf %d (%.3lf); mid=%lf %d (%.3lf)}\n",cou,p_min,res_min,
            res_min/(double)N_exp,p_max,res_max,res_max/(double)N_exp,p_mid,res_mid,res_mid/(double)N_exp);
        if (res_mid<N_exp/2.0) {
            res_min=res_mid;
            p_min=p_mid;
        }
        if (res_mid>N_exp/2.0) {
            res_max=res_mid;
            p_max=p_mid;
        }
        res_mid=0;
        cou++;
    }
    printf("stat={min=%lf %d (%.3lf); max=%lf %d (%.3lf)}\n",p_min,res_min, res_min/(double)N_exp,p_max,res_max,res_max/(double)N_exp);
    printf("mid=%lf\n",(p_max+p_min)/2.0);

    FILE* f=fopen("output.txt","at+");
    fprintf(f,"mid=%lf\n",(p_max+p_min)/2.0);
    fclose(f);

    return 0;
}
