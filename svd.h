

double **dmatrix(int nrl, int nrh, int ncl, int nch);
double *dvector(int nl, int nh);
void free_dvector(double *v, int nl, int nh);
double** transpose(double** a);
double** matrixMul(double** a, double** b);
void svdcmp(double **a, int m, int n, double w[], double **v);
double** inverseCmp(double** input);
