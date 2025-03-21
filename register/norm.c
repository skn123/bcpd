// Copyright (c) 2018-2020 Osamu Hirose
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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<assert.h>

#define SQ(x) ((x)*(x))

static void copy(double *a, const double *b, int n){ int i; for(i=0;i<n;i++) a[i]=b[i];}

static void mean(double *mu, const double *X, int N, int D){
  int n,d; for(d=0;d<D;d++){mu[d]=0;for(n=0;n<N;n++){mu[d]+=X[d+D*n];} mu[d]/=N;}
}

static double scale(const double *X, const double *mu, int N, int D){
  int n,d; double val=0; for(d=0;d<D;d++)for(n=0;n<N;n++){val+=SQ(X[d+D*n]-mu[d]);}
  val/=N*D; return sqrt(val);
}

static void translate(double *X, const double *mu, int N, int D, int sign){
  int n,d; if(!mu) return;
  for(d=0;d<D;d++)for(n=0;n<N;n++) X[d+D*n]+=sign*mu[d];
}

static void resize(double *X, double sc, int N, int D, int inv){
  int n,d;
  if(inv) for(d=0;d<D;d++)for(n=0;n<N;n++){X[d+D*n]/=sc;}
  else    for(d=0;d<D;d++)for(n=0;n<N;n++){X[d+D*n]*=sc;}
}

static void norm_l(double *X, const double *mu, double sc, int N, int D, int revflag){
  int inv=1,plus=1,minus=-1;
  switch(revflag){
    case 0: translate(X,mu,N,D,minus); resize(X,sc,N,D,inv); break;
    case 1: resize(X,sc,N,D,0); translate(X,mu,N,D,plus);    break;
  }
}

/* alias */
void normalize  (double *X, const double *mu, double sc, int N, int D){if(!X){return;} norm_l(X,mu,sc,N,D,0);}
void denormalize(double *X, const double *mu, double sc, int N, int D){if(!X){return;} norm_l(X,mu,sc,N,D,1);}

void normalizer(double *muX, double *scX, double *muY, double *scY, const double *X, const double *Y, int N, int M, int D, const char type){
  if(type=='n'){*scY=*scX=1.0; return;} assert(type=='e'||type=='x'||type=='y');
  mean(muX,X,N,D); *scX=scale(X,muX,N,D);
  mean(muY,Y,M,D); *scY=scale(Y,muY,M,D);
  if(type=='x'){copy(muY,muX,D); *scY=*scX;}
  if(type=='y'){copy(muX,muY,D); *scX=*scY;}
}

