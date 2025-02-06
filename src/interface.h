SEXP C_vistla(SEXP X,SEXP Y,SEXP Flow,SEXP Estimator,SEXP Threshold,SEXP Targets,SEXP Verbose,SEXP Threads){
 enum flow flow=verify_flow(asInteger(Flow));
 enum estimator estimator=verify_estimator(asInteger(Estimator));

 if(!isFrame(X)) error("X has to be a data.frame");
 u32 m=length(X);
 if(m==0) error("X has no columns to trace through");
 u32 n=length(VECTOR_ELT(X,0));
 if(n!=length(Y)) error("X and Y size mismatch");
 u32 nn=(estimator==kt)?(n*(n-1)):n;


 if(isInteger(Threads) && length(Threads)!=1) error("Invalid threads argument");
 u32 nt=asInteger(Threads);
 if(nt<0) error("Invalid threads argument");
 if(nt>omp_get_max_threads()){
  nt=omp_get_max_threads();
  warning("Thread count capped to %d",nt);
 }
 if(nt==0) nt=omp_get_max_threads();

 struct ht **ht=(struct ht**)R_alloc(sizeof(struct ht*),nt);
 for(int e=0;e<nt;e++)
  ht[e]=R_allocHt(nn);

 u32 **x=(u32**)R_alloc(sizeof(u32*),m),
     *nx=(u32*)R_alloc(sizeof(u32),m),
     *y,ny;

 double iomin=asReal(Threshold);
 if(iomin<0) error("Threshold must be at least 0!");

 int verbose=asLogical(Verbose);
 int ntargets=length(Targets);
 bool *t=(bool*)R_alloc(sizeof(bool),m);
 if(ntargets>0){
  int *vt=INTEGER(Targets);
  for(int e=0;e<m;e++) t[e]=false;
  //Convert index vector into a mask; also 1-based to 0-based
  for(int e=0;e<ntargets;e++){
   if(vt[e]<1 || vt[e]>m) error("Invalid targets -- INTERNAL PROBLEM, PLEASE REPORT");
   t[vt[e]-1]=true;
  }
 }else{
  //Synth a vector of all-true
  for(int e=0;e<m;e++) t[e]=true;
  ntargets=m;
 }

 if(verbose) Rprintf("Coercing input\n");

 for(int e=0;e<m;e++){
  SEXP Xe=PROTECT(VECTOR_ELT(X,e));
  x[e]=convertSEXP(*ht,n,Xe,nx+e,estimator);
  UNPROTECT(1);
  if(!(x[e])) error("Wrong X[,%d] type",e+1);
 }
 y=convertSEXP(*ht,n,Y,&ny,estimator);
 if(!y) error("Wrong Y type");

 SEXP Ans=PROTECT(allocVector(VECSXP,3));
 SEXP Mi=PROTECT(allocMatrix(REALSXP,m,m));
 SEXP MiY=PROTECT(allocVector(REALSXP,m));
 setAttrib(MiY,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 SEXP DN=PROTECT(allocVector(VECSXP,2));
 SEXP Xn=PROTECT(getAttrib(X,R_NamesSymbol));
 SET_VECTOR_ELT(DN,0,Xn);
 SET_VECTOR_ELT(DN,1,Xn);
 UNPROTECT(1);
 SEXP AN=PROTECT(allocVector(STRSXP,3));
 SET_STRING_ELT(AN,0,mkChar("tree"));
 SET_STRING_ELT(AN,1,mkChar("mi"));
 SET_STRING_ELT(AN,2,mkChar("miY"));
 setAttrib(Mi,R_DimNamesSymbol,DN);
 setAttrib(Ans,R_NamesSymbol,AN);
 SET_VECTOR_ELT(Ans,1,Mi);
 SET_VECTOR_ELT(Ans,2,MiY);
 UNPROTECT(4); //Mi, MiY, DN and AN

 u32 *P=(u32*)R_alloc(sizeof(u32),m*m);
 double
 *S=(double*)R_alloc(sizeof(double),m*m),
 *mi=REAL(Mi);
 double *miY=REAL(MiY);
 struct heap *queue=R_allocHeap(m*m);
 u32 *si=(u32*)R_alloc(sizeof(u32),m*m);

 //Init the rng
 struct rng rng;
 set_from_r(&rng);

 //Run the main algorithm
 u32 bc=vistla(x,nx,y,ny,m,nn,ht,&rng,nt,verbose,iomin,flow,t,ntargets,mi,miY,queue,P,S,si);


 if(verbose) Rprintf("Generating output.\n");

 //We have our result, but time to convert S, P and si tables into a tree

 //This is needed here because this is the first place we know the size (bc)
 SEXP Tree=PROTECT(allocVector(VECSXP,8));
 SEXP TreeA=PROTECT(allocVector(INTSXP,bc));
 SEXP TreeB=PROTECT(allocVector(INTSXP,bc));
 SEXP TreeC=PROTECT(allocVector(INTSXP,bc));
 SEXP TreePrv=PROTECT(allocVector(INTSXP,bc));
 SEXP TreeScore=PROTECT(allocVector(REALSXP,bc));
 SEXP TreeDepth=PROTECT(allocVector(INTSXP,bc));
 SEXP TreeLeaf=PROTECT(allocVector(LGLSXP,bc));
 SEXP TreeUsed=PROTECT(allocVector(LGLSXP,bc));
 SET_VECTOR_ELT(Tree,0,TreeA);
 SET_VECTOR_ELT(Tree,1,TreeB);
 SET_VECTOR_ELT(Tree,2,TreeC);
 SET_VECTOR_ELT(Tree,3,TreeScore);
 SET_VECTOR_ELT(Tree,4,TreeDepth);
 SET_VECTOR_ELT(Tree,5,TreeLeaf);
 SET_VECTOR_ELT(Tree,6,TreeUsed);
 SET_VECTOR_ELT(Tree,7,TreePrv);
 UNPROTECT(8); //Tree-A, B, C, Prv, Score, Depth, Leaf & Used

 u32
 *ta=(u32*)INTEGER(TreeA),
 *tb=(u32*)INTEGER(TreeB),
 *tc=(u32*)INTEGER(TreeC),
 *td=(u32*)INTEGER(TreeDepth),
 *tp=(u32*)INTEGER(TreePrv);
 Rboolean
 *tu=(Rboolean*)INTEGER(TreeUsed),
 *tl=(Rboolean*)INTEGER(TreeLeaf);
 double *ts=REAL(TreeScore);


 //Generate the (a/b/c) & depth annotations
 for(u32 e=0;e<bc;e++){
  tb[e]=si[e]%m+1;
  tc[e]=si[e]/m+1;
  tp[e]=P[si[e]];
  if(tp[e]==NA_INTEGER){
   ta[e]=NA_INTEGER;
   td[e]=1;
  }else{
   ta[e]=si[tp[e]-1]%m+1;
   td[e]=td[tp[e]-1]+1;
  }
  tu[e]=FALSE; //for now
  ts[e]=S[si[e]];
 }

 //Generate leaf and use annotations
 u32 *seen=P;//We re-use
 for(u32 e=0;e<m;e++) seen[e]=false;
 for(u32 e=0;e<bc;e++){
  if(!seen[tc[e]-1]){
   seen[tc[e]-1]=true;
   tl[e]=TRUE;//This is leaf
   int ee=e;
   for(;(tp[ee]!=NA_INTEGER) && (tu[ee]==FALSE);ee=tp[ee]-1)tu[ee]=TRUE;
   if(tp[ee]==NA_INTEGER) tu[ee]=TRUE;
  }else tl[e]=FALSE;
 }

 SET_VECTOR_ELT(Ans,0,Tree);
 UNPROTECT(1); //Tree

 UNPROTECT(1); //Ans
 return(Ans);
}

