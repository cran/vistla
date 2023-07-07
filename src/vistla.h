/*  128  64  32  16   8   4   2   1  *
 * |res|res|res|roa|hld|hlu|fwd|bkw| */
enum flow{
 backward=1, 
 forward=2, 
 hilldown=8,
 hillup=4,
 noroam=16
};

enum flow verify_flow(u32 x){
 if(x>31) error("Wrong value of the flow");
 if((x&hillup)&&(x&hilldown)) error("Cannot hill up and down at the same time");
 if(((x&hillup)|(x&hilldown))&(x&noroam)){
  warning("Force-path is redundant with up/down.");
  x&=(!noroam);
 }
 return(x);
}

SEXP C_vistla(SEXP X,SEXP Y,SEXP Flow,SEXP Threshold,SEXP Targets,SEXP Verbose,SEXP Threads){
 if(!isFrame(X)) error("X has to be a data.frame");
 u32 m=length(X);
 if(m==0) error("X has no columns to trace through");
 u32 n=length(VECTOR_ELT(X,0));
 if(n!=length(Y)) error("X and Y size mismatch");

 if(isInteger(Threads) && length(Threads)!=1) error("Invalid threads argument");
 u32 nt=asInteger(Threads);
 //TODO: Tune threads number so that it is not about m, but rather < m/<small const like 4>
 if(nt<0) error("Invalid threads argument");
 if(nt>omp_get_max_threads()){
  nt=omp_get_max_threads();
  warning("Thread count capped to %d",nt);
 }
 if(nt==0) nt=omp_get_max_threads();

 struct ht **ht=(struct ht**)R_alloc(sizeof(struct ht*),nt);
 for(int e=0;e<nt;e++)
  ht[e]=R_allocHt(n);
 
 u32 **x=(u32**)R_alloc(sizeof(u32*),m),
     *nx=(u32*)R_alloc(sizeof(u32),m),
     *y,ny;

 double iomin=asReal(Threshold);
 if(iomin<0) error("Threshold must be at least 0!");

 int verbose=asLogical(Verbose);
 int ntargets=length(Targets);
 bool *t=NULL;
 if(ntargets>0){
  int *vt=INTEGER(Targets);
  t=(bool*) R_alloc(sizeof(bool),m);
  for(int e=0;e<m;e++) t[e]=false;
  //Convert index vector into a mask; also 1-based to 0-based
  for(int e=0;e<ntargets;e++){
   if(vt[e]<1 || vt[e]>m) error("Invalid targets -- INTERNAL PROBLEM, PLEASE REPORT");
   t[vt[e]-1]=true;
  }
 }else{
  //Left as a flag
  ntargets=m;
 }

 enum flow flow=verify_flow(asInteger(Flow));
 
 if(verbose) Rprintf("Coercing input\n");

 for(int e=0;e<m;e++){
  SEXP Xe=VECTOR_ELT(X,e);
  x[e]=convertSEXP(*ht,n,Xe,nx+e);
  if(!(x[e])) error("Wrong X[,%d] type",e+1);
 }
 y=convertSEXP(*ht,n,Y,&ny);
 if(!y) error("Wrong Y type");

 SEXP Ans=PROTECT(allocVector(VECSXP,3));
 SEXP Mi=PROTECT(allocMatrix(REALSXP,m,m));
 SEXP MiY=PROTECT(allocVector(REALSXP,m));
 setAttrib(MiY,R_NamesSymbol,getAttrib(X,R_NamesSymbol));
 SEXP DN=PROTECT(allocVector(VECSXP,2));
 SEXP Xn=getAttrib(X,R_NamesSymbol);
 SET_VECTOR_ELT(DN,0,Xn);
 SET_VECTOR_ELT(DN,1,Xn);
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

 if(verbose) Rprintf("Pre-processing\n");

 //Initiate the place for MI calculations
 bool *v=(bool*)R_alloc(sizeof(bool),m*m);

 u32 *cXes=(u32*)R_alloc(sizeof(u32),n*nt),
     *cXees=(u32*)R_alloc(sizeof(u32),n*nt);
 #pragma omp parallel num_threads(nt)
 {
  int tn=omp_get_thread_num();
  u32 *cXe=cXes+n*tn,
      *cXee=cXees+n*tn;
  struct ht *htc=ht[tn];

  #pragma omp for
  for(u32 t=0;t<m*(m-1)/2;t++){
   u32 e,ee;
   t2ij(t,&e,&ee);
   fillHt(htc,n,nx[e],x[e],nx[ee],x[ee],NULL,cXe,cXee,0);
   mi[e+ee*m]=mi[ee+e*m]=miHt(htc,cXe,cXee);
  }
 }

 //Initiate the result
 for(u32 e=0;e<m;e++){
  mi[e+e*m]=NA_REAL;
  S[e+e*m]=NA_REAL;
  P[e+e*m]=NA_INTEGER;
  v[e+e*m]=true;
 }
 for(u32 e=0;e<m;e++)
  for(u32 ee=0;ee<m;ee++)
   if(e!=ee){
    v[ee+e*m]=false;
    //Unreachable; 0 is Y
    P[ee+e*m]=NA_INTEGER;
   }

 u32 *xeys=(u32*)R_alloc(sizeof(u32),n*nt),
     *cYs=(u32*)R_alloc(sizeof(u32),n*nt);
 double *miY=REAL(MiY);

 #pragma omp parallel num_threads(nt)
 {
 
  int tn=omp_get_thread_num();
  u32 *xey=xeys+n*tn,
      *cY=cYs+n*tn,
      *cXe=cXes+n*tn, //Re-used form previous loop
      *cXee=cXees+n*tn;
  u32 *cXeY=cY; //alias
  struct ht *htc=ht[tn];
  #pragma omp for
  for(u32 e=0;e<m;e++){
   u32 nxey=fillHt(htc,n,nx[e],x[e],ny,y,xey,cXe,cY,1);
   miY[e]=miHt(htc,cXe,cY);
   for(u32 ee=0;ee<e;ee++) if(ee!=e){
    fillHt(htc,n,nx[ee],x[ee],nxey,xey,NULL,cXee,cXeY,0);
    double jmiXee_XeY=miHt(htc,cXee,cXeY);
    S[ee+e*m]=mi[ee+e*m]-jmiXee_XeY;
   }
  }

 }
 for(u32 e=0;e<m;e++)
  for(u32 ee=0;ee<e;ee++) if(e!=ee) S[ee+e*m]+=miY[ee];
  
 struct heap *queue=R_allocHeap(m*m);
 
 //Set-up S for the depth-1 scores (all threes involving Y)
 for(u32 e=0;e<m;e++)
  for(u32 ee=0;ee<m;ee++) if(e!=ee){
   u32 idx=ee+e*m;
   double *cS=&(S[idx]); //cS is S[ee,e], from ee->e
   //Transpose to fill the other triangle
   if(ee>e) *cS=S[e+ee*m];
   //Apply the flow ineq for various modes -- here hill is the same as flow
   if(((flow&forward) || (flow&hilldown)) && !(miY[ee]>miY[e])) *cS=0.; 
   if(((flow&backward) || (flow&hillup)) && !(mi[e+ee*m]>miY[e])) *cS=0.; 
   if(*cS<=iomin) *cS=0.;

   if(*cS>0){
    P[ee+e*m]=NA_INTEGER; //These are reachable from Y, coded as NA
    //Push this index on the queue; we do not need it sorter so far, hence ignoring heap constraint
    addBreaking(queue,ee+e*m);
   }
  }
  
 //Restore heap for the queue after addBreaking calls
 heapify(queue,S);

 if(verbose) Rprintf("Looking for optimal paths\n");

 //Re-using some memory
 u32 *xab=xeys,//(u32*)R_alloc(sizeof(u32),n),
     *cXc=cYs,//(u32*)R_alloc(sizeof(u32),n),
     *cXab=cXes,//(u32*)R_alloc(sizeof(u32),n),
     *si=(u32*)R_alloc(sizeof(u32),m*m), //Sorted index for answer
     bc=0; //Number of branches in answer, final size of si
 while(heapLen(queue)){
  //Check if heap is tied, mix the tip if so
  if(isTied(queue,S)){
   if(verbose) Rprintf("Tie detected, breaking at random\n");
   GetRNGstate();
   breakTie(queue,S);
   PutRNGstate();
  }
  
  //Queue head is now accepted as a branch and added to the output
  u32 idx=pop(queue,S);
  v[idx]=true;
  si[bc]=idx;
  bc++;
  u32 a=idx%m,
      b=idx/m;
  if(verbose){
   const char *an=CHAR(STRING_ELT(Xn,a));
   const char *bn=CHAR(STRING_ELT(Xn,b));
   Rprintf("Established path to %s via %s of score %0.4f\n",bn,an,S[idx]);
  }
  if(t && t[b]){
   t[b]=false;
   ntargets--;
   if(!ntargets){
    if(verbose) Rprintf("Target list exhausted, breaking search\n");
    break;
   }
  }
  
  //Time to investigate its sub-branches; to this end, we
  //will need to re-create the a-b mixture. But lazily,
  //maybe it won't be necessary.
  u32 nxab=0;

  for(int c=0;c<m;c++) if((c!=a) && (c!=b)){
   //Do not re-visit
   if(v[b+c*m]) continue;

   //Check flow criterion
   if((flow&forward) && !(mi[a+b*m]>mi[a+c*m])) continue;
   if((flow&backward) && !(mi[c+b*m]>mi[c+a*m])) continue;

   //Check hill criterion
   if((flow&hilldown) && !(miY[b]>miY[c])) continue;
   if((flow&hillup) && !(mi[c+b*m]>miY[c])) continue;

   //Chek the no-roam (no sense with hill present)
   if(flow&noroam){
    u32 haveToBackOff=0;
    u32 back_bc=P[a+b*m];
    u32 cnt=10;
    while(back_bc!=NA_INTEGER){
     cnt--;
     if(cnt<1) break;
     u32 back_si=si[back_bc-1];
     u32 pa=back_si%m,
         pb=back_si/m;
     if(pa==c || pb==c){
      haveToBackOff=1;
      break;
     }
     back_bc=P[pa+pb*m];
    }
    if(haveToBackOff) continue;
   }

   //Start calculating new score
   double nS=mi[a+c*m]+mi[b+c*m];
   if((nS<=iomin) || (nS<=S[b+c*m])) continue; //Checks for immediate back-off
   if(!nxab){
    //Also generatex xab
    nxab=fillHt(*ht,n,nx[a],x[a],nx[b],x[b],xab,NULL,NULL,1);
   }
   fillHt(*ht,n,nx[c],x[c],nxab,xab,NULL,cXc,cXab,0);
   nS-=miHt(*ht,cXc,cXab);

   //Another back-off check
   if(nS<=iomin) continue;

   if(nS>S[idx]) nS=S[idx]; //Bottleneck was earlier
   //This has to be bigger than iomin still, since S element have to be

   if(nS>S[b+c*m]){
    if(verbose){
     const char *cn=CHAR(STRING_ELT(Xn,c));
     Rprintf(" widened path to %s from %0.4f to %0.4f\n",cn,S[b+c*m],nS);
    }
    //We have a relaxed path!
    S[b+c*m]=nS; //WARN: Decreasing S would violate the heap
    P[b+c*m]=bc; //Register its super-branch
    //Update the score on the queue
    update(queue,b+c*m,S);
   }
  }
 }


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
 UNPROTECT(8); //Tree-A,B,C,Prv,Score,Depth,Leaf and Used

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
 u32 *seen=(u32*)R_alloc(sizeof(u32),m); //P;//We re-use
 for(u32 e=0;e<m;e++) seen[e]=false;
 for(u32 e=0;e<bc;e++){
  if(!seen[tc[e]-1]){
   seen[tc[e]-1]=true;
   tl[e]=TRUE;//This is leaf
   int ee=e;
   for(;(tp[ee]!=NA_INTEGER)&&(tu[ee]==FALSE);ee=tp[ee]-1) tu[ee]=TRUE;
   if(tp[ee]==NA_INTEGER) tu[ee]=TRUE;
  }else tl[e]=FALSE;
 }
 
 SET_VECTOR_ELT(Ans,0,Tree);
 UNPROTECT(1); //Tree

 UNPROTECT(1); //Ans
 return(Ans);
}

