u32 *convertSEXP_mle(struct ht *ht,u32 n,SEXP in,u32 *nout){
 struct feature *f=ingestSEXP_mle(n,in);
 return(produce_mle((u32*)(f->x),ht,n,n,NULL,nout));
}

u32 *convertSEXP_kt(u32 n,SEXP in,u32 *nout){
 struct feature *f=ingestSEXP_kt(n,in);
 return(produce_kt(f->x,is_double(f),n,NULL,nout));
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
 SEXP Ans;PROTECT(Ans=allocVector(INTSXP,n));
 int *ans=INTEGER(Ans);
 for(int e=0;e<n;e++){
  if(out[e]>nout) error("Conversion integrity error");
  ans[e]=out[e];
 }
 UNPROTECT(1);
 return(Ans);
}
