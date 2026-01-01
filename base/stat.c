// Copyright (c) 2025-2025 Osamu Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include <math.h>
#define SQ(x) ((x)*(x))

static double mean_l(const double *x, int step, int n){
  int i; double val=0; for(i=0;i<n;i++){val+=x[step*i];} val/=n; return val;
}

static double sdev_l(const double *x, int step, int n){
  int i; double val=0,mu=mean_l(x,step,n);
  for(i=0;i<n;i++){val+=SQ(x[step*i]-mu);} val/=n; return fmax(sqrt(val),1e-6);
}

void mvmean(double *mean, const double *X, int D, int N){int d; for(d=0;d<D;d++) mean[d]=mean_l(X+d,D,N);}
void mvsdev(double *sdev, const double *X, int D, int N){int d; for(d=0;d<D;d++) sdev[d]=sdev_l(X+d,D,N);}

/* log of gaussian determinant term (diagonal) */
double lnddet(const double *sdev, int D){
  int d; double val=D*log(2*M_PI);
  for(d=0;d<D;d++) val+=2*log(sdev[d]);
  return -0.5*val;
}

/* log of gaussian determinant term (isotropic) */
double lnidet(double sdev, int D){
  return -0.5*D*(log(2*M_PI)+2*log(sdev));
}

/* log gaussian density */
double lnnormd(const double *x, const double *mu, const double *sdev, int D){
  int d; double val,sum=0;
  for(d=0;d<D;d++){val=(x[d]-mu[d])/sdev[d]; sum+=val*val;}
  return -0.5*sum+lnddet(sdev,D);
}

