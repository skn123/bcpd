// Copyright (c) 2018-2019 Osamu Hirose
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
#include<string.h>
#include<math.h>
#include<assert.h>
#include"kdtree.h"
#include"misc.h"
#define SQ(x) ((x)*(x))

double genrand64_real1(void);

static void resampling(int *nums, const double *probs, int N, int M) {
  int j,i=0; double cum=probs[0], u=genrand64_real1()/(double)M;
  for(j=0;j<M;j++){
    while(u>cum&&i<N-1){i++;cum+=probs[i];}
    nums[i]++; u=((j+1)+genrand64_real1())/(double)M;
  }
}

static void downsample_a(int *U, int L, double *X, int D, int N){
  int d,l; if(L>N){goto err01;}
  randperm(U,N); return;

  err01:
  printf("\n\n  ERROR: L<=N must be satisfied in the function 'downsample_a'. Abort.\n\n");
  exit(EXIT_FAILURE);
}

static void downsample_b(int *U, int L, double *X, int D, int N, double e){
  int *S,*T,*a,*c; double *w,*v; int sd=sizeof(double),si=sizeof(int);
  int d,j,n,q,l=0,mtd=MAXTREEDEPTH; double val=0;
  /* allocation */
  T=calloc(3*N+1,si); S=calloc(N*mtd,si);
  a=calloc(6*N,  si); w=calloc(N,sd);
  v=calloc(2*N,  sd); c=calloc(N,si);
  /* build kdtree */ 
  kdtree(T,a,v,X,D,N);
  /* count #neighbors */
  #pragma omp parallel for private (j) private (q)
  for(n=0;n<N;n++){j=q=0;w[n]=0.0;
    do{eballsearch_next(&j,S+mtd*n,&q,X+D*n,e,X,T,D,N);if(j>=0){w[n]+=1.0;}} while(q);
    assert(w[n]>=1.0);
  }
  /* sampling probabilities */
  for(n=0;n<N;n++) val+=1.0/(w[n]);
  for(n=0;n<N;n++) w[n]=1.0/(w[n]*val);
  /* resampling */
  resampling(c,w,N,L);
  /* output */
  for(n=0;n<N;n++)for(j=0;j<c[n];j++) U[l++]=n;

  free(T);free(a);free(v);
  free(S);free(w);free(c);
}

static size_t voxelize(int *v, size_t *div, size_t *cum, const double *X, int D, int N, double e){
  int d,j,n; double *max,*min; double w; size_t K;
  int sd=sizeof(double),ss=sizeof(size_t); double val=0;

  /* allocation */
  max=calloc(D,sd);
  min=calloc(D,sd);
  /* bounding box */
  for(d=0;d<D;d++){min[d]=X[d];for(n=0;n<N;n++){min[d]=X[d+D*n]<min[d]?X[d+D*n]:min[d];}}
  for(d=0;d<D;d++){max[d]=X[d];for(n=0;n<N;n++){max[d]=X[d+D*n]>max[d]?X[d+D*n]:max[d];}}
  /* divide in grid & count points in a voxel */
  for(d=0;d<D;d++){w=max[d]-min[d]; div[d]=w<e?1:(int)ceil(w/e);}
  cum[0]=1; for(d=1;d<D;d++) cum[d]=cum[d-1]*div[d-1];
  for(n=0;n<N;n++){v[n]=0;for(d=0;d<D;d++){j=floor((X[d+D*n]-min[d])/e);j-=(j==div[d])?1:0;v[n]+=cum[d]*j;}}
  /* count the number of voxels (including empty voxels) */
  K=cum[D-1]*div[D-1];

  free(max);
  free(min);
  return K;
}

/* voxel grid filter */
static void downsample_c(int *U, int L, double *X, int D, int N, double e){
  int l=0,j,n,num; size_t K; int *v,*c,*np; size_t *div,*cum; double *w,val=0;
  int sd=sizeof(double),si=sizeof(int),ss=sizeof(size_t);

  /* allocation */
  v=calloc(N,si); div=calloc(D,ss);
  c=calloc(N,si); cum=calloc(D,ss);
  w=calloc(N,sd);

  /* divide into voxels */
  K=voxelize(v,div,cum,X,D,N,e); if(K>=1e8){printf("  ERROR: Voxel grid width is too small. Abort.\n\n"); exit(EXIT_FAILURE);}
  np=calloc(K,si);
  /* sampling probabilities */
  for(n=0;n<N;n++) np[v[n]]++;
  for(n=0;n<N;n++){num=np[v[n]];assert(num>0);w[n]=1.0/num;val+=w[n];}
  for(n=0;n<N;n++) w[n]/=val;
  /* resampling */
  resampling(c,w,N,L);
  /* output */
  for(n=0;n<N;n++)for(j=0;j<c[n];j++) U[l++]=n;

  free(v); free(div); free(w);
  free(c); free(cum); free(np);
}

static void geometry(double *geo, int *np, size_t K, size_t *div, size_t *cum, int D) {
  size_t k;
  #pragma omp parallel for
  for(k=0;k<K;k++) { size_t coord; int d; int empty;
    if(np[k]==0){continue;} empty=0;
    for(d=0;d<D;d++){ coord=(k/cum[d])%div[d];
      if (coord==0)       {empty++;} else if(np[k-cum[d]]==0){empty++;}
      if (coord==div[d]-1){empty++;} else if(np[k+cum[d]]==0){empty++;}
    } geo[k]=(empty==2*D)?0.0:(empty/(double)(2.0*D));
  }
}

void vgisample(int *U, int L, double *X, int D, double *fx, int Df, int N, double e, double eps, double lmd){
  int l=0,j,d,n,num; size_t k,K; int *v,*c,*np; double *w; int sd=sizeof(double),si=sizeof(int),ss=sizeof(size_t);
  double val=0,*ave,*var,*max,*sum,*wf,*wg; size_t *cum,*div; int flag=e<0?1:0; e*=e<0?-1:1; assert(e>0);

  /* allocation */
  v=calloc(N,si); div=calloc(D,ss); max=calloc(N,sd); wf=calloc(N,sd);
  c=calloc(N,si); cum=calloc(D,ss); sum=calloc(N,sd); w =calloc(N,sd);

  /* devide into voxels */
  K=voxelize(v,div,cum,X,D,N,e);
  if(K>=1e8){printf("  ERROR: Voxel grid width is too small. Abort.\n\n"); exit(EXIT_FAILURE);}

  /* allocation */
  np=calloc(K,si); ave=calloc(K,sd);
  wg=calloc(K,sd); var=calloc(K,sd);

  for(n=0;n<N; n++) np[v[n]]++;
  for(d=0;d<Df;d++){
    memset(ave,0,K*sd);
    memset(var,0,K*sd);

    for(n=0;n<N;n++) ave[v[n]]+=fx[d+Df*n];
    for(k=0;k<K;k++)if(np[k]) ave[k]/=np[k];
    for(n=0;n<N;n++) var[v[n]]+=SQ(fx[d+Df*n]-ave[v[n]]);
    for(k=0;k<K;k++)if(np[k]) var[k]/=np[k];
    for(n=0;n<N;n++) max[n] =fmax(max[n],var[v[n]]);
    for(n=0;n<N;n++) sum[n]+=var[v[n]];
  }

  /* sampling probabilities */
  if( flag)for(n=0;n<N;n++) wf[n]=sqrt(max[n]);
  if(!flag)for(n=0;n<N;n++) wf[n]=sqrt(sum[n]);

  if(lmd>0) geometry(wg,np,K,div,cum,D);
  for(n=0;n<N;n++) w[n]=fmax(wf[n],lmd*wg[v[n]])+eps;

  for(n=0;n<N;n++) val +=w[n];
  for(n=0;n<N;n++) w[n]/=val;

  /* resampling */
  resampling(c,w,N,L);
  /* output */
  for(n=0;n<N;n++)for(j=0;j<c[n];j++)if(l<L) U[l++]=n;

  free(w); free(ave); free(sum); free(np);
  free(v); free(var); free(div); free(wf);
  free(c); free(max); free(cum); free(wg);
}

void downsample(int *U, int L, double *X, int D, int N, double e){
  if     (e<0) downsample_c(U,L,X,D,N,-e);
  else if(e>0) downsample_b(U,L,X,D,N, e);
  else         downsample_a(U,L,X,D,N);
}

