#define GET_A(x) ((x)>>32)
#define GET_B(x) ((x)&0x00000000ffffffff)
#define MAKE_AB(a,b) (((uint64_t)(a)<<32)|(b))

struct ht{
 struct hte **map;
 struct hte *cnt;
 u32 N;
 u32 nAB;
};

struct hte{
 uint64_t ab;
 struct hte *nxt;
 u32 c;
};

struct ht* R_allocHt(u32 N){
 struct ht *ans=(struct ht*)R_alloc(sizeof(struct ht),1);
 ans->N=N;
 ans->map=(struct hte**)R_alloc(sizeof(struct hte*),N);
 ans->cnt=(struct hte*)R_alloc(sizeof(struct hte),N);
 return(ans);
}

uint32_t static inline fillHt(struct ht *Q,u32 N,u32 nA,u32 *a,u32 nB,u32 *b,u32 *mix,u32 *cA,u32 *cB,u32 mixOff){
 if(cA) for(u32 e=0;e<nA;e++) cA[e]=0;
 if(cB) for(u32 e=0;e<nB;e++) cB[e]=0;
 u32 Neff=nA*nB; 

 if(Neff<N && !mix){
  //HT not needed, we use it as lookup table 
  for(u32 e=0;e<Neff;e++) Q->cnt[e].c=0;
  for(u32 e=0;e<N;e++){
   uint32_t _a=a[e]-1,_b=b[e]-1,hab=_a+_b*nA;
   uint64_t _ab=MAKE_AB(_a,_b);
   Q->cnt[hab].c++;
   Q->cnt[hab].ab=_ab;
   if(cA) cA[_a]++;
   if(cB) cB[_b]++;
  }
  return(Q->nAB=Neff);
 }

 if(Neff<N && mix){
  //LT, but map used to collapse mix
  uint32_t nAB=0;
  for(u32 e=0;e<Neff;e++) Q->map[e]=NULL;
  for(u32 e=0;e<N;e++){
   uint32_t _a=a[e]-1,_b=b[e]-1,hab=_a+_b*nA;
   uint64_t _ab=MAKE_AB(_a,_b);
   struct hte **E=Q->map+hab;
   if(!*E){
    Q->cnt[nAB].ab=_ab;
    Q->cnt[nAB].c=1;
    *E=Q->cnt+nAB;
    nAB++;
   }else (*E)->c++;

   if(cA) cA[_a]++;
   if(cB) cB[_b]++;
   mix[e]=(*E-Q->cnt)+mixOff;
  }
  return(Q->nAB=nAB);
 }
 
 //Regular HT
 uint32_t nAB=0;
 for(u32 e=0;e<N;e++) Q->map[e]=NULL;
 
 for(u32 e=0;e<N;e++){
  uint32_t _a=a[e]-1,_b=b[e]-1,hab=(_a^_b)%N;
  uint64_t _ab=MAKE_AB(_a,_b);

  struct hte **E;
  for(E=Q->map+hab;(*E)&&(*E)->ab!=_ab;E=&((*E)->nxt));

  if(!*E){
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   Q->cnt[nAB].c=1;
   *E=Q->cnt+nAB;
   nAB++;
  }else (*E)->c++;

  if(cA) cA[_a]++;
  if(cB) cB[_b]++;
  if(mix) mix[e]=(*E-Q->cnt)+mixOff;
 }
 return(Q->nAB=nAB);
}

uint32_t static inline fillHtOne(struct ht *Q,u32 N,u32 *in,u32 *out,u32 mixOff){
 uint32_t nAB=0;
 for(u32 e=0;e<N;e++) Q->map[e]=NULL;
 for(u32 e=0;e<N;e++){
  if(in[e]==NA_INTEGER) error("NAs are not accepted");
  uint64_t _ab=(uint64_t)(in[e]);
  uint32_t hab=_ab%N;

  struct hte **E;
  for(E=Q->map+hab;E[0] && E[0]->ab!=_ab;E=&(E[0]->nxt));

  if(!E[0]){
   Q->cnt[nAB].ab=_ab;
   Q->cnt[nAB].nxt=NULL;
   E[0]=Q->cnt+nAB;
   nAB++;
  }
  out[e]=(E[0]-Q->cnt)+mixOff;
 }
 return(nAB);
}

double miHt(struct ht *Q,u32 *cA,u32 *cB){
 double ans=0.,N=Q->N;
 for(u32 e=0;e<Q->nAB;e++){
  if(!(Q->cnt[e].c)) continue;
  double cAB=Q->cnt[e].c,
   _cA=cA[GET_A(Q->cnt[e].ab)],
   _cB=cB[GET_B(Q->cnt[e].ab)];
  ans+=cAB*log((cAB*N)/(_cA*_cB));
 }
 return(ans/N);
}

