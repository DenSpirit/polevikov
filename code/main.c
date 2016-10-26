#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int N;
int NODEC;
double EPS;
double H;

double** create() {
    int size = N;
    double** m = calloc(size, sizeof(double*));
    int i;
    for(i=0;i<size;i++) {
        m[i] = calloc(size,sizeof(double));
    }
    return m;
}

void drop(double** m, int size) {
    int i;
    for(i=0;i<size;i++) {
        free(m[i]);
    }
    free(m);
}

typedef struct state state;
struct state {
    double Gr;
    double Pr;
    double** psi;
    double** omega;
    double** T;
};

double middle(double** m, int i, int j, char pm, char dir) {
    int iindex = dir == 'x' ? 1 : 0;
    int jindex = dir == 'y' ? 1 : 0;
    double f = (m[i+iindex][j+jindex]-m[i-iindex][j-jindex])/(2*H);
    if (pm == 0 || pm == '-' && f < 0 || pm == '+' && f > 0) {
        return f;
    } else {
        return 0;
    }
}

double calc_psi(int i, int j, state s) {
    return (H*H * s.omega[i][j] + s.psi[i-1][j] + s.psi[i+1][j] 
            + s.psi[i][j-1] + s.psi[i][j+1])/4.;
}

double calc_omega(int i, int j, state s) {
    double im1= s.omega[i-1][j]*(1/(H*H) + 1/H * middle(s.psi, i, j, '+', 'y'));
    double ip1= s.omega[i+1][j]*(1/(H*H) - 1/H * middle(s.psi, i, j, '-', 'y'));
    double jm1= s.omega[i][j-1]*(1/(H*H) - 1/H * middle(s.psi, i, j, '-', 'x'));
    double jp1= s.omega[i][j+1]*(1/(H*H) + 1/H * middle(s.psi, i, j, '+', 'x'));
    double ij_coeff = (4/(H*H) + middle(s.psi,i,j,'+','y')/H
                    + middle(s.psi,i,j,'+','x')/H
                    - middle(s.psi,i,j,'-','y')/H
                    - middle(s.psi,i,j,'-','x')/H);
    return (s.Gr * middle(s.T, i, j, 0, 'x') + im1+ip1+jm1+jp1)/ij_coeff;
}

double calc_T(int i, int j, state s) {
    double im1= s.T[i-1][j]*(1/(s.Pr*H*H) + 1/H * middle(s.psi, i, j, '+', 'y'));
    double ip1= s.T[i+1][j]*(1/(s.Pr*H*H) - 1/H * middle(s.psi, i, j, '-', 'y'));
    double jm1= s.T[i][j-1]*(1/(s.Pr*H*H) - 1/H * middle(s.psi, i, j, '-', 'x'));
    double jp1= s.T[i][j+1]*(1/(s.Pr*H*H) + 1/H * middle(s.psi, i, j, '+', 'x'));
    double ij_coeff = (4/(H*H*s.Pr) + middle(s.psi,i,j,'+','y')/H
                    + middle(s.psi,i,j,'+','x')/H
                    - middle(s.psi,i,j,'-','y')/H
                    - middle(s.psi,i,j,'-','x')/H);
    return (im1+ip1+jm1+jp1)/ij_coeff;
}

void border_psi(state s) {
    int i;
    for(i=0;i<N;i++) {
        s.psi[i][0] = 0;
        s.psi[i][1] = 0;
        s.psi[0][i] = 0;
        s.psi[1][i] = 0;
        s.psi[i][NODEC] = 0;
        s.psi[i][NODEC-1] = 0;
        s.psi[NODEC][i] = 0;
        s.psi[NODEC-1][i] = 0;
    }
}

void border_omega(state s) {
    int i;
    for(i=0;i<N;i++) {
        s.omega[0][i] = - s.omega[1][i] / 2. - 3./(H*H) * s.psi[1][i];
        s.omega[NODEC][i] = - s.omega[NODEC-1][i] / 2. - 3./(H*H) * s.psi[NODEC-1][i];
        s.omega[i][0] = - s.omega[i][1] / 2. - 3./(H*H) * s.psi[i][1];
        s.omega[i][NODEC] = - s.omega[i][NODEC-1] / 2. - 3./(H*H) * s.psi[i][NODEC-1];
    }
}

void border_T(state s) {
    int i;
    for(i=0;i<N;i++) {
        s.T[i][0] = 0;
        s.T[0][i] = 0;
        s.T[NODEC][i] = 0;
        s.T[i][NODEC] = sin(M_PI*i*H);
    }
}

void output(double** m) {
    int i,j;
    for(i=0;i<=NODEC;i++) {
        for(j=0;j<=NODEC;j++) {
            //printf("%lf %lf %lf\n", (double) i*H, (double) j*H, m[i][j]);
            printf("%.3lf ", m[i][j]);
        }
        printf("\n");
    }
}


int main(int argc, char* argv[])
{
    //sscanf(argv[1],"%d",&NODEC);
    NODEC = 10;
    N = NODEC + 1;
    H = 1. / ((double) NODEC);
    EPS = H*H;
    state s;
    sscanf(argv[1],"%lf",&(s.Gr));
    sscanf(argv[2],"%lf",&(s.Pr));
    s.psi = create();
    s.omega = create();
    s.T = create();
    
    border_psi(s);
    border_T(s);
    border_omega(s);

    int i;
    int j;
    for(i=1;i<NODEC-1;i++) {
        for(j=1;j<NODEC-1;j++) {
            s.psi[i][j] = calc_psi(i, j, s);
            s.omega[i][j] = calc_omega(i, j, s);
            s.T[i][j] = calc_T(i, j, s);
        }
    }
    printf("T\n");
    output(s.T);
    printf("psi\n");
    output(s.psi);
    printf("omega\n");
    output(s.omega);

    drop(s.psi, N);
    drop(s.omega, N);
    drop(s.T, N);

    return 0;
}
