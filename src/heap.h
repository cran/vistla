//Simple binary heap; used for a priority queue by vistla

struct heap {
 u32 *queue;
 u32 *map;
 u32 end;
};

struct heap* R_allocHeap(u32 N){
 struct heap *ans=(struct heap*)R_alloc(sizeof(struct heap),1);
 ans->queue=(u32*)R_alloc(sizeof(u32),N);
 ans->map=(u32*)R_alloc(sizeof(u32),N);
 for(u32 e=0;e<N;e++) ans->map[e]=NA_INTEGER;
 ans->end=0;
 return(ans);
}

//resetHeap MUST be called on mallocHeap to initialise
struct heap* mallocHeap(u32 N){
 struct heap *ans=malloc(sizeof(struct heap));
 ans->queue=malloc(sizeof(u32)*N);
 ans->map=malloc(sizeof(u32)*N);
 return(ans);
}

void resetHeap(struct heap *heap,u32 N){
 for(u32 e=0;e<N;e++) heap->map[e]=NA_INTEGER;
 heap->end=0;
}

void freeHeap(struct heap *heap){
 free(heap->queue);
 free(heap->map);
 free(heap);
}

static inline u32 parent(u32 e){
 return((e-1)/2);
}

static inline u32 childA(u32 e){
 return(2*e+1);
}

//NOTE: Only for cB>bA tie breaking
// algorithm is corrct.
static inline u32 childB(u32 e){
 return(2*e+2);
}

static inline void swap(struct heap *h,u32 a,u32 b){
 u32 ia=h->queue[a];
 u32 ib=h->queue[b];

 h->map[ia]=b;
 h->map[ib]=a;

 h->queue[b]=ia;
 h->queue[a]=ib;
}

static inline bool cmp(struct heap *h,u32 a,u32 b,double *score){
 double sa=score[h->queue[a]];
 double sb=score[h->queue[b]];
 return(sa>sb);
}

void swim(struct heap *h,u32 e,double *score){
 while(e>0 && cmp(h,e,parent(e),score)){
  swap(h,e,parent(e));
  e=parent(e);
 }
}

void sink(struct heap *h,u32 e,double *score){
 while(true){
  u32 max=e;
  u32 a=childA(e);
  u32 b=childB(e);
  if((a<h->end) && cmp(h,a,max,score)) max=a;
  if((b<h->end) && cmp(h,b,max,score)) max=b;

  if(max==e) break;
  swap(h,max,e);
  e=max;
 }
}

void update(struct heap *h,u32 e,double *score){
 //Like bump, but add if not already present
 if(h->map[e]==NA_INTEGER){
  h->map[e]=h->end;
  h->queue[h->end]=e;
  h->end++;
 }
 swim(h,h->map[e],score);
}

void addBreaking(struct heap *h,u32 e){
 //Push element breaking heap.
 //Heapify must be called after this, before
 // any pops, etc. Use for transactional push.
 h->map[e]=h->end;
 h->queue[h->end]=e;
 h->end++;
}

u32 heapLen(struct heap *h){
 return(h->end);
}

void heapify(struct heap *h,double *score){
 //Array with 0 or 1 element is already a heap
 if(h->end>1) for(u32 r=(h->end+1)/2+1;r<=h->end;r++){
   u32 e=h->end-r;
   sink(h,e,score);
  }
}

u32 pop(struct heap *h,double *score){
 u32 ans=h->queue[0];
 swap(h,0,h->end-1);
 h->end--;
 sink(h,0,score);
 h->map[ans]=NA_INTEGER;
 return(ans);
}

u32 selTied(struct heap *h,double *score,struct rng *rng){
 double rs=score[h->queue[0]];
 u32 tag=random_int(rng);
 u32 sel=0;
 u32 lasttie=0;
 for(u32 e=1;e<h->end && e<=childB(lasttie);e++)
  if(score[h->queue[e]]==rs){
   lasttie=e;
   double newtag=random_int(rng);
   if(newtag>tag){
    tag=newtag;
    sel=e;
   }
  }

 return(sel);
}

bool isTied(struct heap *h,double *score){
 double rs=score[h->queue[0]];
 return(
  ((1<h->end) && (score[h->queue[1]]==rs)) ||
  ((2<h->end) && (score[h->queue[2]]==rs))
  );
}

void breakTie(struct heap *h,double *score,struct rng *rng){
 swap(h,0,selTied(h,score,rng));
}

bool integrity_test(struct heap *h,double *score){
 bool map_integral=true;
 for(u32 e=0;e<h->end;e++) map_integral&=(h->map[h->queue[e]]==e);
 if(!map_integral) error("FATAL: Map has lost integrity!");

 bool heap_property=true;
 if(score){
  if(h->end>1){
   for(u32 e=1;e<h->end;e++){
    heap_property&=(score[h->queue[e]]<=score[h->queue[parent(e)]]);
   }
  }
  if(!heap_property) error("FATAL: Heap property violated!");
 }
 return(map_integral&heap_property);
}

SEXP C_heapTest(SEXP A,SEXP B,SEXP Torture){
 int N=length(A);
 double *a=REAL(A);
 int M=length(B);
 double *b=REAL(B);

 bool torture=(asLogical(Torture)==TRUE);
 int cap=M;
 if(M<N) error("Invalid test data, B cannot be shorter than A");

 double *x=(double*)R_alloc(sizeof(double),cap);
 for(int e=0;e<M;e++) x[e]=R_NegInf;
 for(int e=0;e<N;e++) x[e]=a[e];

 SEXP Ans;PROTECT(Ans=allocVector(REALSXP,N+M));
 double *ans=REAL(Ans);

 //Push with heapify
 struct heap *h=(struct heap*)R_allocHeap(cap);
 for(u32 e=0;e<N;e++) addBreaking(h,e);
 heapify(h,x);

 integrity_test(h,x);

 //Verify push from heapify
 for(u32 e=0;e<N;e++){
  ans[e]=a[pop(h,x)];
  if(e>0 && (ans[e-1]<ans[e])) error("FATAL: Sorting has failed (1)!");
  if(torture) integrity_test(h,x);
 }

 //Push with update
 for(u32 e=0;e<N;e++){
  update(h,e,x);
  if(torture) integrity_test(h,x);
 }

 //Verify push from update
 for(u32 e=0;e<N;e++){
  double popped=a[pop(h,x)];
  if(ans[e]!=popped) error("FATAL: Sorting has failed (2)!");
  if(torture) integrity_test(h,x);
 }

 //Get A back, with heapify since it is faster
 for(u32 e=0;e<N;e++) addBreaking(h,e);
 heapify(h,x);

 for(u32 e=0;e<M;e++){
  if(b[e]<x[e]) error("Invalid test data, cannot update to lower!");
  x[e]=b[e];
  update(h,e,x);
  if(torture) integrity_test(h,x);
 }

 //Verify update
 for(u32 e=0;e<M;e++){
  ans[e+N]=x[pop(h,x)];
  if(e>0 && (ans[N+e-1]<ans[N+e])) error("FATAL: Sorting has failed (3)!");
  if(torture) integrity_test(h,x);
 }


 //Return ans back to R so it can be compared with sort
 UNPROTECT(1);
 return(Ans);
}

SEXP C_heapTiedTest(SEXP A,SEXP B){
 int N=length(A);
 double *a=REAL(A);
 int M=length(B);
 double *b=REAL(B);

 int cap=M;
 if(M<N) error("Invalid test data, B cannot be shorter than A");

 double *x=(double*)R_alloc(sizeof(double),cap);
 for(int e=0;e<M;e++) x[e]=R_NegInf;
 for(int e=0;e<N;e++) x[e]=a[e];

 SEXP Ans;PROTECT(Ans=allocVector(INTSXP,M));
 int *ans=INTEGER(Ans);

 //Push with heapify
 struct heap *h=(struct heap*)R_allocHeap(cap);
 for(u32 e=0;e<N;e++) addBreaking(h,e);
 heapify(h,x);
 integrity_test(h,x);

 //Push new data with update
 for(u32 e=0;e<M;e++){
  if(b[e]<x[e]) error("Invalid test data, cannot update to lower");
  x[e]=b[e];
  update(h,e,x);
 }
 integrity_test(h,x);

 //Pop breaking ties
 struct rng rng;
 set_from_r(&rng);
 for(u32 e=0;e<M;e++){
  if(isTied(h,x)){
   breakTie(h,x,&rng);
   ans[e]=-(pop(h,x)+1);
  }else{
   ans[e]=pop(h,x)+1;
  }
 }

 UNPROTECT(1);
 return(Ans);
}
