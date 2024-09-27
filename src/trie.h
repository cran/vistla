//Trie to merge result trees in ensemble mode

struct vertex {
 u32 count;
 u32 key;
 struct vertex *prv;
 struct vertex *nxt;
 struct vertex *dwn;
};

static struct vertex* append(struct vertex *a,struct vertex *appended,struct vertex **head){
 if(!a){
  *head=appended;
 }else{
  a->nxt=appended;
  appended->prv=a;
 }
 return(appended);
}

static struct vertex* behead(struct vertex **x){
 struct vertex* ans=*x;
 *x=ans->nxt;
 if(*x){
  (*x)->prv=NULL;
  ans->nxt=NULL;
 }
 return(ans);
}

static struct vertex* merge(struct vertex *a,struct vertex *b){
 struct vertex *ans=NULL;
 struct vertex *cur=ans;
 while(a!=NULL || b!=NULL){
  if(!a){
   append(cur,b,&ans);
   break;
  }
  if(!b){
   append(cur,a,&ans);
   break;
  }
  if(a->key==b->key){
   //We do a<-merge(a,b)
   a->count+=b->count;
   a->dwn=merge(a->dwn,b->dwn);
   //Remove b's head; its down stays, it was merged into a
   free(behead(&b));
   cur=append(cur,behead(&a),&ans);
  }else if(a->key<b->key){
   cur=append(cur,behead(&a),&ans);
  }else{
   //b key is lower than a key
   cur=append(cur,behead(&b),&ans);
  }
 }
 return(ans);
}

//There are two copies of this stuff; this one is actually useful and calculates prv
static u32 count_vtx(struct vertex *a){
 u32 c=0;
 while(a){
  c+=1+count_vtx(a->dwn);
  a=a->nxt;
 }
 return(c);
}
static void encode_vtx(struct vertex *a,int *ti,int *ts,int *td,int *tp,int depth,u32 *idx,u32 prv){
 while(a){
  u32 e=*idx;
  (*idx)++;
  ti[e]=a->key;
  td[e]=depth;
  tp[e]=prv+1;
  ts[e]=a->count;
  encode_vtx(a->dwn,ti,ts,td,tp,depth+1,idx,e);
  a=a->nxt;
 }
}

static SEXP trie_toR(struct vertex *a){
 //+1 for the root element
 u32 vn=count_vtx(a)+1;
 SEXP Tree=PROTECT(allocVector(VECSXP,5)),
      Ti=PROTECT(allocVector(INTSXP,vn)), //Name index
      Ts=PROTECT(allocVector(INTSXP,vn)), //Score (hit count)
      Td=PROTECT(allocVector(INTSXP,vn)), //Depth
      Tl=PROTECT(allocVector(LGLSXP,vn)), //Is leaf
      Tp=PROTECT(allocVector(INTSXP,vn)); //Previous link
 int *ti=INTEGER(Ti),
     *ts=INTEGER(Ts),
     *td=INTEGER(Td),
     *tp=INTEGER(Tp);
 int *tl=INTEGER(Tl);
 ti[0]=NA_INTEGER;
 ts[0]=0;
 td[0]=-1;
 tp[0]=NA_INTEGER;
 u32 idx=1;

 encode_vtx(a,INTEGER(Ti),INTEGER(Ts),INTEGER(Td),INTEGER(Tp),0,&idx,0);


 for(u32 e=0;e<vn;e++) tl[e]=ts[e]>0;

 SET_VECTOR_ELT(Tree,0,Ti);
 SET_VECTOR_ELT(Tree,1,Ts);
 SET_VECTOR_ELT(Tree,2,Td);
 SET_VECTOR_ELT(Tree,3,Tl);
 SET_VECTOR_ELT(Tree,4,Tp);
 UNPROTECT(6);
 return(Tree);
}

static void free_vtx(struct vertex *a){
 while(a){
  if(a->dwn) free_vtx(a->dwn);
  struct vertex *tmp=a;
  a=a->nxt;
  free(tmp);
 }
}

static struct vertex* array_into(u32 n,u32 *values){
 struct vertex *v=NULL,*vn;
 for(u32 e=n;e>0;e--){
  vn=malloc(sizeof(struct vertex));
  vn->key=values[e-1];
  vn->count=(e==n)?1:0;
  vn->dwn=v;
  vn->prv=NULL;
  vn->nxt=NULL;
  v=vn;
 }
 return(v);
}

static struct vertex* prune_low_count(struct vertex *a,u32 count_limit){
 struct vertex *o=a;
 while(a){
  if(a->dwn) a->dwn=prune_low_count(a->dwn,count_limit);
  if(!a->dwn && a->count<=count_limit){
   //Remove a
   if(a->prv) a->prv->nxt=a->nxt;
   if(a->nxt) a->nxt->prv=a->prv;
   struct vertex *td=a;
   a=a->nxt;
   if(td==o) o=a;
   free(td);
  }else{
   a=a->nxt;
  }
 }
 return(o);
}

static struct vertex* find_or_insert(struct vertex **ar,u32 key){
 struct vertex *ans;
 struct vertex *a=*ar;
 if(!a){
  ans=malloc(sizeof(struct vertex));
  ans->prv=ans->nxt=ans->dwn=NULL;
  ans->key=key;
  ans->count=0;
  *ar=ans;
  return(ans);
 }else{
  if(a->key>key){
   ans=malloc(sizeof(struct vertex));
   ans->dwn=ans->prv=NULL;
   ans->nxt=a;
   ans->key=key;
   ans->count=0;
   a->prv=ans;
   *ar=ans;
   return(ans);
  }
  while(a->nxt && a->key<key) a=a->nxt;
  if(a->key==key){
   return(a);
  }else if(a->key>key){
   //Insert BETWEEN a->prv (which has to exist) and a
   ans=malloc(sizeof(struct vertex));
   ans->dwn=NULL;
   ans->nxt=a;
   ans->prv=a->prv;
   a->prv->nxt=ans;
   a->prv=ans;
   ans->key=key;
   ans->count=0;
   return(ans);
  }else{
   //Append
   ans=malloc(sizeof(struct vertex));
   ans->dwn=ans->nxt=NULL;
   ans->prv=a;
   ans->key=key;
   ans->count=0;
   a->nxt=ans;
   return(ans);
  }
 }

}

static struct vertex* from_vistla_tree(u32 n,u32 *ta,u32 *tb,u32 *tc,u32 *tp,Rboolean *tu,Rboolean *tl){
 struct vertex **pvs=malloc(sizeof(struct vertex)*n);
 struct vertex *root=malloc(sizeof(struct vertex));
 root->count=0;
 root->key=777;
 root->nxt=root->prv=root->dwn=NULL;

 for(u32 e=0;e<n;e++){
  if(tu[e]){
   if(ta[e]==NA_INTEGER){
    //Root veep
    struct vertex *rv=find_or_insert(&(root->dwn),tb[e]);
    struct vertex *nv=find_or_insert(&(rv->dwn),tc[e]);
    nv->count=tl[e]?1:0;
    pvs[e]=nv;
   }else{
    //Normal vertex
    struct vertex *up=pvs[tp[e]-1];
    struct vertex *nv=find_or_insert(&(up->dwn),tc[e]);
    nv->count=tl[e]?1:0;
    pvs[e]=nv;
   }
  }
 }
 struct vertex* ans=root->dwn;
 free(root);
 free(pvs);
 return(ans);

}

SEXP C_trieTest(SEXP A){
 u32 l=length(A);
 struct vertex *ans=NULL;
 for(u32 e=0;e<l;e++){
  SEXP Z=VECTOR_ELT(A,e);
  u32 zl=length(Z);
  u32 *z=(u32*)INTEGER(Z);
  struct vertex *f=array_into(zl,z);
  ans=merge(ans,f);
 }
 SEXP Ans=trie_toR(ans);
 free_vtx(ans);
 return(Ans);
}

