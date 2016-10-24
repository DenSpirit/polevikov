int N;
int NODEC;
double EPS;
double H;

void printv(double* v,char name[]) {
#ifdef VERBOSE
    printf("%s:",name);
    int i;
    for(i=0;i<N;i++){
        printf("%.2lf ",v[i]);
    }
    printf("\n");
#endif //VERBOSE
}
void printd(double d,char name[]){
#ifdef VERBOSE
    fprintf(stderr,"%s:%.3lf\n",name,d);
#endif //VERBOSE
}
void printm(double** A,char name[]){
#ifdef VERBOSE
    int i,j;
    printf("%s:\n",name);
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            double o;
                if(i==j)
                    o=A[1][j];
                else if(i==j-1)
                    o=A[2][j-1];
                else if(i==j+1)
                    o=A[0][j+1];
                else
                    o=0;
            printf("%.0lf\t",o);
        }
        printf("\n");
    }
#endif //VERBOSE
}

int main(int argc, char* argv[])
{
    double Gr;
    double Pr;
    //sscanf(argv[1],"%d",&NODEC);
    NODEC = 1000;
    N = NODEC + 1;
    H = 1. / ((double) NODEC);
    EPS = H*H;
    sscanf(argv[2],"%lf",&Gr);
    sscanf(argv[3],"%lf",&Pr);
    if(P<0)
    {
        P = fabs(P);
    }
    if (fabs(alpha) > 180)
    {
        return 1; // no periodical support
    }
    else
    {
        alpha = M_PI*alpha/180.; //to %pi notation
    }
    //first generation
    double* z = gen_z0();
    double* r = gen_r0();
    printv(r,"r0");
    printv(z,"z0");
    //'previous'
    double* z0 = calloc(N,sizeof(double));
    double* r0 = calloc(N,sizeof(double));
    //constant coeff tri-vectors for sweep method
    double** Ar = gen_Ar();
    double** Az = gen_Az();

    double *Fz,*Fr;
    Fz=calloc(N,sizeof(double));
    Fr=calloc(N,sizeof(double));
    int k = 1;
    double diff_r;
    double diff_z;
    double I1, I2;
    do
    {
#ifdef VERBOSE
        printf("\nk=%d\n\n",k++);
#endif
        I1 = gen_I1(z,r);
        I2 = gen_I2(alpha, r);
        printd(I1,"I1");
        double Q = gen_Q(z,r,I1,I2,P,A,alpha);
        printd(Q,"Q");
        memcpy(z0,z,sizeof(double)*N);
        memcpy(r0,r,sizeof(double)*N);
        gen_Fz(z0,r0,P,A,alpha,I1,Q,Fz);
        sweep(Az,Fz, z);
        int i;
        for(i = 0;i < N; i++) {
            z[i] = z[i] * W + (1. - W) * z0[i];
        }
        gen_Fr( z,r0,P,A,alpha,I1,Q,Fr);
        sweep(Ar,Fr, r);
        for(i = 0;i < N; i++) {
            r[i] = r[i] * W + (1. - W) * r0[i];
        }
        sanity("r", r);
        sanity("z", z);
        printv(r,"r");
        printv(z,"z");
        diff_r = vec_difference_norm(r,r0);
        diff_z = vec_difference_norm(z,z0);
#ifdef VERBOSE
        printf("diff_z: %lf\ndiff_r: %lf\n",diff_z,diff_r);
#endif
    }
    while (diff_r > EPS || diff_z > EPS);
    int i=0;
    fprintf(stderr, "I:%lf\n", I1);
    double I13 = cbrt(I1);
#ifdef LATEXOUT
    printf("r\tz\t\n");
#endif
    for(;i<N;i++){
#ifdef LATEXOUT
        printf("%lf\t%lf\n", r[i]/I13, z[i]/I13);
#else
        printf("%lf\t%lf\t%lf\n", i*H, r[i]/I13, z[i]/I13);
#endif
    }
    free(z);
    free(r);
    free(z0);
    free(r0);
    free(Fz);
    free(Fr);
    free(Ar[0]);
    free(Ar[1]);
    free(Ar[2]);
    free(Ar);
    free(Az[0]);
    free(Az[1]);
    free(Az[2]);
    free(Az);
    return 0;
}
