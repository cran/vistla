SEXP C_vistlaEnsemble(SEXP X,SEXP Y,SEXP Flow,SEXP Estimator,SEXP Threshold,SEXP Targets,SEXP Ensemble,SEXP Threads){
 enum flow flow=verify_flow(asInteger(Flow));
 enum estimator estimator=verify_estimator(asInteger(Estimator));

 if(!isFrame(X)) error("X has to be a data.frame");
 u32 m=length(X);
 if(m==0) error("X has no columns to trace through");
 u32 n=length(VECTOR_ELT(X,0));
 if(n!=length(Y)) error("X and Y size mismatch");

 if(length(Ensemble)!=3) error("Invalid replication options, Ens len is %d",length(Ensemble));
 int *ens=(int*)INTEGER(Ensemble);
 if(ens[0]<1) error("Replication count must be positive");
 u32 repeats=ens[0];
 if(ens[1]<0 || ens[1]>n) error("Invalid value of resample");
 u32 boot=ens[1];
 if(ens[2]<0 || ens[2]>repeats) error("Invalid value of the threshold for ensemble prune");
 u32 prune_count=ens[2];
 if((boot==0) && n<5)
  error("For bootstrap, at least five objects are required to make a practical difference");

 if(isInteger(Threads) && length(Threads)!=1) error("Invalid threads argument");
 u32 nt=asInteger(Threads);
 if(nt<0) error("Invalid threads argument");
 if(nt>omp_get_max_threads()){
  nt=omp_get_max_threads();
  warning("Thread count capped to %d",nt);
 }
 if(nt==0) nt=omp_get_max_threads();


 struct feature **aX=(struct feature**)R_alloc(sizeof(struct feature*),m),
                *aY;

 if(estimator==mle){
  for(u32 e=0;e<m;e++)
   aX[e]=ingestSEXP_mle(n,VECTOR_ELT(X,e));
  aY=ingestSEXP_mle(n,Y);
 }else if(estimator==kt){
  for(u32 e=0;e<m;e++)
   aX[e]=ingestSEXP_kt(n,VECTOR_ELT(X,e));
  aY=ingestSEXP_kt(n,Y);
 } //else impossible

 double iomin=asReal(Threshold);
 if(iomin<0) error("Threshold must be at lest 0");

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


 //Init the rng
 struct rng rng_g;
 set_from_r(&rng_g);


 //Place for a consensus tree
 struct vertex *full_tree=NULL;

 #pragma omp parallel num_threads(nt)
 {

  u32 *P=malloc(sizeof(u32)*m*m);
  double
  *S=malloc(sizeof(double)*m*m),
  *mi=malloc(sizeof(double)*m*m);
  double *miY=malloc(sizeof(double)*m);
  struct heap *queue=mallocHeap(m*m);
  u32 *si=malloc(sizeof(u32)*m*m);
  struct vertex *thread_tree=NULL;

#pragma omp for
  for(u32 i=0;i<repeats;i++){
   struct rng rng=rng_g;
   rng.stream=2*i+1;
   //Actually coerce data
   u32 nn=n;
   //For no-bootstrap, they could be made earlier and once, for all threads;
   // but given this is rarely used and gains are not huge, it is left as it is.
   u32 **x=malloc(sizeof(u32*)*m),
       *nx=malloc(sizeof(u32)*m),
       *y=NULL,ny;
   struct ht* ht;
   //If boot is on, it will fix nn to the number of selected elements
   bool *mask;
   if(boot==0) mask=boot_mask(n,&rng,&nn);
   else if(boot==n) mask=NULL;else {
    nn=boot;
    mask=count_mask(n,&rng,nn);
   }
   //ASSERT nn>=2
   if(estimator==mle){
    ht=mallocHt(nn);
    for(u32 e=0;e<m;e++){
     x[e]=produce_mle((u32*)(aX[e]->x),ht,n,nn,mask,&(nx[e]));
    }
    y=produce_mle((u32*)(aY->x),ht,n,nn,mask,&ny);
   }else if(estimator==kt){
    u32 *idx=which_mask(n,mask,nn);
    for(u32 e=0;e<m;e++)
     x[e]=produce_kt(aX[e]->x,is_double(aX[e]),nn,idx,&(nx[e]));
    y=produce_kt(aY->x,is_double(aY),nn,idx,&ny);
    free(idx);
    nn=nn*(nn-1);
    ht=mallocHt(nn);
   } //else impossible
   free(mask);

   //Initialise/reset queue object
   resetHeap(queue,m*m);

   //Run the main algorithm
   u32 bc=vistla(x,nx,y,ny,m,nn,&ht,&rng,1,false,iomin,flow,t,ntargets,mi,miY,queue,P,S,si);

   for(u32 e=0;e<m;e++) free(x[e]);
   free(x);
   free(y);
   free(nx);

   freeHt(ht);
   //We have our result, but time to convert S, P and si tables into a tree
   u32
   *ta=malloc(sizeof(u32)*bc),
   *tb=malloc(sizeof(u32)*bc),
   *tc=malloc(sizeof(u32)*bc),
   *tp=malloc(sizeof(u32)*bc);
   Rboolean
   *tu=malloc(sizeof(Rboolean)*bc),
   *tl=malloc(sizeof(Rboolean)*bc);

   //Generate the (a/b/c) & depth annotations
   for(u32 e=0;e<bc;e++){
    tb[e]=si[e]%m+1;
    tc[e]=si[e]/m+1;
    tp[e]=P[si[e]];
    if(tp[e]==NA_INTEGER){
     ta[e]=NA_INTEGER;
    }else{
     ta[e]=si[tp[e]-1]%m+1;
    }
    tu[e]=FALSE; //for now
   }

   //Generate leaf and use annotations
   u32 *seen=P;//We re-use
   for(u32 e=0;e<m;e++) seen[e]=false;
   for(u32 e=0;e<bc;e++){
    if(!seen[tc[e]-1]){
     seen[tc[e]-1]=true;
     tl[e]=TRUE;//This is leaf
     int ee=e;
     for(;(tp[ee]!=NA_INTEGER) && (tu[ee]==FALSE);ee=tp[ee]-1) tu[ee]=TRUE;
     if(tp[ee]==NA_INTEGER) tu[ee]=TRUE;
    }else tl[e]=FALSE;
   }

   struct vertex *link_tree=from_vistla_tree(bc,ta,tb,tc,tp,tu,tl);
   free(ta);free(tb);free(tc);free(tp);free(tu);free(tl);

   //This effectively frees link_tree and old thread_tree
   thread_tree=merge(thread_tree,link_tree);


  }   //FOR

  //This effectively frees thread tree & old full tree
 #pragma omp critical
  {
   full_tree=merge(full_tree,thread_tree);
  }

  free(P);free(S);free(mi);free(miY);free(si);
  freeHeap(queue);

 }

 //Trim the tree requested
 if(prune_count)
  full_tree=prune_low_count(full_tree,prune_count);

 //Return & free
 SEXP Ans=PROTECT(trie_toR(full_tree));
 free_vtx(full_tree);
 UNPROTECT(1);
 return(Ans);
}

