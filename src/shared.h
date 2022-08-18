u32 nok(double x){
 double xx=fabs(x);
 return((xx>1e-14)&&(xx<1e300));
}

u32 *convertSEXP(struct ht *ht,u32 n,SEXP in,u32 *nout){
 if(isFactor(in)||isLogical(in)){
  //Integer-alike which needs collapsing into 1..n_levels
  u32 *out=(u32*)R_alloc(sizeof(u32),n);
  //fillHtOne includes NA check
  *nout=fillHtOne(ht,n,(u32*)INTEGER(in),out,1);
  return(out);
 }
 if(isInteger(in)){
  int *x=INTEGER(in),min=INT_MAX,max=INT_MIN;
  for(u32 e=0;e<n;e++){
   if(x[e]==NA_INTEGER) error("NAs are not accepted");
   min=min<x[e]?min:x[e];
   max=max>x[e]?max:x[e];
  }
  u32 *out=(u32*)R_alloc(sizeof(u32),n);
  if(max==min){
   //Constant stuff, won't be too useful
   *nout=1;
   for(u32 e=0;e<n;e++)
    out[e]=1;
   return(out);
  }
  if(n<6){
   *nout=2;
  }else if(n>30){
   *nout=10;
  }else{
   *nout=n/3;
  }
  double range=max-min;
  for(u32 e=0;e<n;e++){
   out[e]=(int)(((double)(x[e]-min))/range*(double)(*nout))+1;
   //Survive numerical errors
   if(out[e]>*nout) out[e]=*nout;
   if(out[e]<1) out[e]=1;
  }
  return(out);
 }
 if(isReal(in)){
  //Magically make discrete by scattering into 10-bins
  double *x=REAL(in),min=INFINITY,max=-INFINITY;
  for(u32 e=0;e<n;e++){
   if(!R_FINITE(x[e])) error("Non-finite values and NAs are not accepted");
   min=min<x[e]?min:x[e];
   max=max>x[e]?max:x[e];
  }
  u32 *out=(u32*)R_alloc(sizeof(u32),n);
  if(max==min){
   //Real value is almost constant
   *nout=1;
   for(u32 e=0;e<n;e++)
    out[e]=1;
   return(out);
  }
  if(n<6){
   *nout=2;
  }else if(n>30){
   *nout=10;
  }else{
   *nout=n/3;
  }
  //Cut, like R's cut(x,nout) would; but only if max and min are numerically ok
  if(nok(max) && nok(min)){
   min=min-(max-min)/1000.;
   max=max+(max-min)/1000.;
  }
  for(u32 e=0;e<n;e++){
   out[e]=(int)((x[e]-min)/(max-min)*(double)(*nout))+1;
   //Survive numerical errors
   if(out[e]>*nout) out[e]=*nout;
   if(out[e]<1) out[e]=1;
  }
  return(out);
 }
 //Other stuff
 return(NULL);
}

static inline void t2ij(u32 t,u32 *ip,u32 *jp){
 double k=t;
 double jf=(sqrt(8.*k+1.)-1.)/2.;
 u32 j=floor(jf)+1;
 u32 i=t-j*(j-1)/2;
 ip[0]=i;
 jp[0]=j;
}

SEXP C_convert(SEXP X){
 int N=length(X);
 struct ht *ht=R_allocHt(N);
 
 SEXP Ans; PROTECT(Ans=allocVector(INTSXP,N));
 int *ans=INTEGER(Ans);
 u32 nout;
 int *cnv=(int*)convertSEXP(ht,N,X,&nout);
 if(!cnv) error("Invalid input");
 for(int e=0;e<N;e++){
  ans[e]=cnv[e];
 }

 char buf[64];
 SEXP Levels=PROTECT(allocVector(STRSXP,nout));
 for(int e=0;e<nout;e++){
  snprintf(buf,64,"l%d",e+1);
  SET_STRING_ELT(Levels,e,mkChar(buf));
 }
 setAttrib(Ans,R_LevelsSymbol,Levels);
 setAttrib(Ans,R_ClassSymbol,mkString("factor"));

 UNPROTECT(2);
 
 return(Ans);
}
