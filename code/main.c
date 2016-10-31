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
    /* Grasgoph/Prandtl and arrays */
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
    if(i==1 || j == 1 || i == NODEC-1 || j == NODEC-1) {
        return 0;
    }
    return (H*H * s.omega[i][j] + s.psi[i-1][j] 
                                + s.psi[i+1][j] 
                                + s.psi[i][j-1] 
                                + s.psi[i][j+1])/4.;
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

/**
 * psi has additional border equations
 * that will be processed in calc_psi()
 */
double border_psi(state s) {
    int i;
    for(i=0;i<N;i++) {
        s.psi[i][0] = 0;
        s.psi[0][i] = 0;
        s.psi[i][NODEC] = 0;
        s.psi[NODEC][i] = 0;
    }
    return 0;
}

double border_omega(state s) {
    int i;
    double diff = 0;
    for(i=0;i<N;i++) {
        double left = - s.omega[1][i] / 2. - 3./(H*H) * s.psi[1][i];
        double right = - s.omega[NODEC-1][i] / 2. - 3./(H*H) * s.psi[NODEC-1][i];
        double top = - s.omega[i][NODEC-1] / 2. - 3./(H*H) * s.psi[i][NODEC-1];
        double bottom = - s.omega[i][1] / 2. - 3./(H*H) * s.psi[i][1];
        diff += fabs(s.omega[0][i] - left);
        s.omega[0][i] = left;
        diff += fabs(s.omega[NODEC][i] - right);
        s.omega[NODEC][i] = right;
        diff += fabs(s.omega[i][0] - bottom);
        s.omega[i][0] = bottom;
        diff += fabs(s.omega[i][NODEC] - top);
        s.omega[i][NODEC] = top;
    }
    return diff;
}

double border_T(state s) {
    int i;
    for(i=0;i<N;i++) {
        s.T[i][0] = 0;
        s.T[0][i] = 0;
        s.T[NODEC][i] = sin(M_PI*i*H);
        s.T[i][NODEC] = 0;
    }
    return 0;
}

typedef struct outs outs;
struct outs {
    double x;
    double y;
    double z;
};
int dcomp(const double a, const double b) {
    if (fabs(a - b) < 0.000001) {
        return 0;
    } else if (a > b) {
        return 1;
    } else {
        return -1;
    }
}
int outcompare(const void* av, const void* bv) {
    outs* a = (outs*) av;
    outs* b = (outs*) bv;
    int r;
    r = dcomp(a->z, b->z);
    if(!r) {
        r = dcomp(a->x, b->x);
    }
    if(!r) {
        r = dcomp(a->y, b->y);
    }
    return r;
}
void output(const double** m) {
    int i,j,k;
    k=0;
    outs* o = calloc(N*N, sizeof(outs));
    for(i=0;i<=NODEC;i++) {
        for(j=0;j<=NODEC;j++) {
            //printf("%lf %lf %lf\n", (double) i*H, (double) j*H, m[i][j]);
            //printf("%.3lf ", m[j][i]);
            o[k].x = i*H;
            o[k].y = j*H;
            o[k++].z = m[i][j];
            printf("{%.2lf,%.2lf,%.2lf},", i*H, j*H, m[i][j]);
        }
        //printf("\n");
    }
    /*for(i=0;i<k;i++) {
        printf("(%.3lf,%.3lf,%.3lf)\n", o[k].x, o[k].y, o[k].z);
    }*/
//    qsort(o, N*N, sizeof(outs), outcompare);
/*    for(i=0;i<k;i++) {
        outs oo = o[i];
        printf("(%.2lf,%.2lf,%.2lf)", oo.x, oo.y, oo.z);
        if (i>0 && fabs(oo.z - o[i-1].z) > EPS) {
            printf("\n");
        }
    }*/
    free(o);
}

double cycle(state s) {
    double diff_psi;
    double diff_T;
    double diff_omega;
    diff_psi = border_psi(s);
    diff_T = border_T(s);
    diff_omega = border_omega(s);

    int i;
    int j;
    for(i=1;i<NODEC;i++) {
        for(j=1;j<NODEC;j++) {
            double omega = calc_omega(i, j, s);
            diff_omega += fabs(omega - s.omega[i][j]);
            s.omega[i][j] = omega;

            double psi = calc_psi(i, j, s);
            diff_psi += fabs(psi - s.psi[i][j]);
            s.psi[i][j] = psi;

            double T = calc_T(i, j, s);
            diff_T += fabs(T - s.T[i][j]);
            s.T[i][j] = T;
        }
    }
#ifdef VERBOSE
    printf("dpsi %lf domega %lf dT %lf\n", diff_psi, diff_omega, diff_T);
#endif
    return diff_psi + diff_omega + diff_T;
}

int main(int argc, char* argv[])
{
    NODEC = 20;
    N = NODEC + 1;
    H = 1. / ((double) NODEC);
    EPS = H*H;
    state s;
    sscanf(argv[1],"%lf",&(s.Gr));
    sscanf(argv[2],"%lf",&(s.Pr));
    s.psi = create();
    s.omega = create();
    s.T = create();
    
    double diff = HUGE_VAL;
    while(diff > EPS) {
        diff = cycle(s);
#ifdef VERBOSE
        printf("%lf\n", diff);
#endif
    }
#ifdef VERBOSE
    printf("T\n");
    output((const double**) s.T);
    printf("Psi\n");
    output((const double**) s.psi);
    printf("Omega\n");
    output((const double**) s.omega);
#endif

    drop(s.psi, N);
    drop(s.omega, N);
    drop(s.T, N);

    return 0;
}
