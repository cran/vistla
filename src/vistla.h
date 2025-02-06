static u32 vistla(
 u32 **x,u32 *nx,u32 *y,u32 ny,u32 m,u32 n, /*Dataset*/
 struct ht **ht, /*Hash tables */ struct rng *rng,
 int nt, /*Thread count*/ int verbose, double iomin, enum flow flow,
 bool *t, u32 ntargets,/*Targets*/
 double *mi,double *miY,struct heap *queue,u32 *P,double *S,u32* si /*Outputs*/
 ){
 if(verbose) Rprintf("Pre-processing\n");
 //Initiate the place for MI calculations
 bool *v=malloc(sizeof(bool)*m*m);

 u32 *cXes=malloc(sizeof(u32)*n*nt),
     *cXees=malloc(sizeof(u32)*n*nt);
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

 u32 *xeys=malloc(sizeof(u32)*n*nt),
     *cYs=malloc(sizeof(u32)*n*nt);

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
   for(u32 ee=0;ee<m;ee++) if(ee!=e){
     fillHt(htc,n,nx[ee],x[ee],nxey,xey,NULL,cXee,cXeY,0);
     double jmiXee_XeY=miHt(htc,cXee,cXeY);
     S[ee+e*m]=jmiXee_XeY;
    }
  }
  #pragma omp barrier
  #pragma omp for
  for(u32 e=0;e<m;e++)
   for(u32 ee=0;ee<m;ee++) if(e!=ee)
     S[ee+e*m]=(miY[ee]+mi[ee+e*m])-S[ee+e*m];
 }

 //Set-up S for the depth-1 scores (all threes involving Y)
 for(u32 e=0;e<m;e++)
  for(u32 ee=0;ee<m;ee++) if(e!=ee){
    u32 idx=ee+e*m;
    double *cS=&(S[idx]); //cS is S[ee,e], from ee->e
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
 u32 *xab=xeys,
     *cXc=cYs,
     *cXab=cXes,
     bc=0; //Number of branches in answer, final size of si
 while(heapLen(queue)){
  //Check if heap is tied, mix the tip if so
  if(isTied(queue,S)){
   if(verbose) Rprintf("Tie detected, breaking at random\n");
   breakTie(queue,S,rng);
  }

  //Queue head is now accepted as a branch and added to the output
  u32 idx=pop(queue,S);
  v[idx]=true;
  si[bc]=idx;
  bc++;
  u32 a=idx%m,
      b=idx/m;
  if(verbose)
   Rprintf("Established path to %d via %d of score %0.4f\n",b+1,a+1,S[idx]);

  if(t && t[b]){
   t[b]=false;
   ntargets--;
   if(verbose) Rprintf("Targets even used! Left: %d\n",ntargets);
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

    //Check the no-roam (no sense with hill present)
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
    if((nS<=iomin) || (nS<=S[b+c*m])) continue;//Checks for immediate back-off
    if(!nxab){
     //Also generatex xab
     nxab=fillHt(*ht,n,nx[a],x[a],nx[b],x[b],xab,NULL,NULL,1);
    }
    fillHt(*ht,n,nx[c],x[c],nxab,xab,NULL,cXc,cXab,0);
    nS-=miHt(*ht,cXc,cXab);

    //Another back-off check
    if(nS<=iomin) continue;

    if(nS>S[idx]) nS=S[idx];//Bottleneck was earlier
    //This has to be bigger than iomin still, since S element have to be

    if(nS>S[b+c*m]){
     if(verbose)
      Rprintf(" widened path to %d from %0.4f to %0.4f\n",c+1,S[b+c*m],nS);

     //We have a relaxed path!
     S[b+c*m]=nS; //WARN: Decreasing S would violate the heap
     P[b+c*m]=bc; //Register its super-branch
     //Update the score on the queue
     update(queue,b+c*m,S);
    }
   }
 }

 free(v);
 free(cXes);
 free(cXees);
 free(xeys);
 free(cYs);
 return(bc);
}
