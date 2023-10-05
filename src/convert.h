u32 *convertSEXP_mle(struct ht *ht,u32 n,SEXP in,u32 *nout){
 if(isFactor(in)||isLogical(in)){
  //Integer-alike which needs collapsing into 1..n_levels
  u32 *out=(u32*)R_alloc(sizeof(u32),n);
  //fillHtOne includes NA check
  *nout=fillHtOne(ht,n,(u32*)INTEGER(in),out,1);
  return(out);
 }
 error("Only logical and factor inputs are acccepted with the MLE estimator");
 return(NULL); //Unreachable
}

bool isBinaryFactor(SEXP x){
 if(!isFactor(x)) return(false);
 u32 l=length(getAttrib(x,R_LevelsSymbol));
 return((l==2)||(l==1)); //We accept constant factors, though they are useless
}

u32 *convertSEXP_kt(u32 n,SEXP in,u32 *nout){
 if(n>65536) error("Kendall transformation covers only up to 2^16 elements");
 if(n<2) error("Kendall transformation requires at least 2 objects");
 u32 nn=n*(n-1);
 u32 *out=(u32*)R_alloc(sizeof(u32),nn);

if(isReal(in)){
  double *v=REAL(in);
  u32 *x=out;
  for(int e=0;e<n;e++)
   if(ISNAN(v[e])) error("NAs nor NaNs are not allowed in input");

  //We need to have collapsed values of KT, hence a branch here;
  // two encodings are possible:
  //  a) = -> 1, < -> 2, > -> 3
  //  b) < -> 1, > -> 2, = -> 3
  //To have a collapsed encoding, we should have at least one state of 1..nout.
  //Thus, we have three options:
  //  1) in is constant; this requires coding a, output of 1s only and nout=1
  //  2) in is not constant and has no ties; this requires coding b, output of 1s or 2s and nout=2
  //  3) in it not constant but has ties; this works in either coding, but we'll use a, and nout=3
  if(v[0]==v[1]){
   //There is at least one tie -- coding a
   nout[0]=1;
   for(int e=0;e<n;e++)
    for(int ee=0;ee<n;ee++)
     if(e!=ee){
      u32 val=(v[e]<v[ee])+(v[e]>v[ee])*2+1;
      (x++)[0]=val;
      if(val==3) nout[0]=3;
     }
  }else{
   //There is at least one < -- coding b
   nout[0]=2;
   for(int e=0;e<n;e++)
    for(int ee=0;ee<n;ee++)
     if(e!=ee){
      u32 val=(v[e]<=v[ee])+2*(v[e]>=v[ee]);
      (x++)[0]=val;
      if(val==3) nout[0]=3;
     }
  }
  return(out);
 }
 
 if(isInteger(in)||isLogical(in)||isOrdered(in)||isBinaryFactor(in)){
  int *v=INTEGER(in);
  u32 *x=out;
  for(int e=0;e<n;e++)
   if(v[e]==NA_INTEGER) //Logical NA is NA_INTEGER
    error("NAs are not allowed in input");

  if(v[0]==v[1]){
   //There is at least one tie -- coding a
   nout[0]=1;
   for(int e=0;e<n;e++)
    for(int ee=0;ee<n;ee++)
     if(e!=ee){
      u32 val=(v[e]<v[ee])+(v[e]>v[ee])*2+1;
      (x++)[0]=val;
      if(val==3) nout[0]=3;
     }
  }else{
   //There is at least one < -- coding b
   nout[0]=2;
   for(int e=0;e<n;e++)
    for(int ee=0;ee<n;ee++)
     if(e!=ee){
      u32 val=(v[e]<=v[ee])+2*(v[e]>=v[ee]);
      (x++)[0]=val;
      if(val==3) nout[0]=3;
     }
  }
  return(out);
 }
   
 error("Only real, integer, logical, ordered and 2-level factor inputs are accepted with the KT estimator");
 return(NULL);
}

u32 *convertSEXP(struct ht *ht,u32 n,SEXP in,u32 *nout,enum estimator estimator){
 if(estimator==mle) return(convertSEXP_mle(ht,n,in,nout));
 if(estimator==kt) return(convertSEXP_kt(n,in,nout));
 return(NULL); //Unreachable due to estimator sanitisation upstream
}

SEXP C_convertTest(SEXP X,SEXP Estimator){
 enum estimator estimator=verify_estimator(asInteger(Estimator));
 u32 nout=0,
 n=length(X);
 struct ht *ht=R_allocHt(n);
 u32 *out=convertSEXP(ht,n,X,&nout,estimator);
 n=(estimator==kt)?(n*(n-1)):n;
 SEXP Ans; PROTECT(Ans=allocVector(INTSXP,n));
 int *ans=INTEGER(Ans);
 for(int e=0;e<n;e++){
  if(out[e]>nout) error("Conversion integrity error");
  ans[e]=out[e];
 }
 UNPROTECT(1);
 return(Ans);
}
