// Copyright (c) 2021-2022 Osamu Hirose
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
#include<assert.h>

#define LAST (*heap)
#define ROOT 1

static void swap(int *heap, int *idx, int i, int j){
  int tmp=heap[i]; heap[i]=heap[j]; heap[j]=tmp;
  idx[heap[i]]=i;
  idx[heap[j]]=j;
}

static void shiftup(int *heap, int *idx, const double *keys, int i){
  int p; assert(i>=1&&i<=LAST);
  for(;(p=i/2)&&(keys[heap[p]]>keys[heap[i]]);i=p) swap(heap,idx,p,i);
}

static void shiftdn(int *heap, int *idx, const double *keys, int i){
  int c,cl,cr; assert(i);
  for(;(cl=2*i)<=LAST;i=c){ cr=cl+1;
    c=(LAST==cl)?cl:(keys[heap[cl]]<keys[heap[cr]]?cl:cr);
    if(keys[heap[i]]>keys[heap[c]]) swap(heap,idx,c,i);
    else break;
  }
}

void heap_init(int *heap){LAST=0;}

void heap_insert(int *heap, int *idx, double *keys, int node, double key){
  heap[++LAST]=node;
  idx [node]=LAST;
  keys[node]=key;
  shiftup(heap,idx,keys,LAST);
}

void heap_extract(int *node, int *heap, int *idx, const double *keys){
  *node=heap[ROOT];
  heap[ROOT]=heap[LAST--];
  if(LAST>=ROOT) {idx[heap[ROOT]]=ROOT; shiftdn(heap,idx,keys,ROOT);}
}

void heap_downkey(int *heap, int *idx, double *keys, int node, double key){
  int i; assert(key<keys[node]);
  keys[node]=key; i=idx[node];
  assert(i>=1&&i<=LAST);
  shiftup(heap,idx,keys,i);
}

void print_heap(const int *heap, const double *keys, int mode){
  int i,j,d=1; int L,R;
  for(i=1;i<=LAST;i++){
    if(i%d==0) { d++;
      L=1<<(d-2);
      R=(1<<(d-1))-1;
      R=LAST<R?LAST:R;
      if(R<=LAST){ 
        for(j=L;j<=R;j++){
          if(mode) printf("%d%c",   heap[j],      (j==R)?'\n':' ');
          else     printf("%.1lf%c",keys[heap[j]],(j==R)?'\n':' ');
        }
      } 
    }
  } printf("\n");
}

