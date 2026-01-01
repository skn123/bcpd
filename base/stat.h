void mvmean (double *mean, const double *X, int D, int N);
void mvsdev (double *sdev, const double *X, int D, int N);

double lnddet  (const double *sdev, int D);
double lnidet  (double sdev, int D);
double lnnormd (const double *x, const double *mu, const double *sdev, int D);
