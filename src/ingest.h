/*
   Internal representaiton of the input
   The basic idea goes this way:
   ingest* functions convert SEXP to struct features
   they can error, hence they only use Ralloc

   For now on, we generally don't use this stuff,
   it will be useful when NAs are on track

 */

/* 128  64  32  16   8   4   2   1 *
* ??? ??? ??? ??? INT NAs CLL MNG */
enum feature_state {
 managed=1,
 collapsed=2,
 has_missing=4,
 integer=8
};
// shuffle(permutator) *->non-collapsed+non-managed
// fromSEXP -> non-collapsed+managed
// collapse -> *->collapsed+non-managed
// kt -> *->collapsed+non-managed
// impute: with_miss *-> managed+non-collapsed
//         without_miss x->x
// for_mle -> collapse(impute(.))
// for_kt -> kt(impute(.))
// boot -> l(.) for_mle(.)|for_kt(.) (shuffle(.))
// shuffle doesn't really use re-order, just a bin mask;
//  this allows for an efficient impute index shuffle

struct feature {
 u32 l;
 void *x;
 enum feature_state state;
};

static inline bool is_double(struct feature *f){
 return(!(integer&(f->state)));
}

bool isBinaryFactor(SEXP x){
 if(!isFactor(x)) return(false);
 u32 l=length(getAttrib(x,R_LevelsSymbol));
 return((l==2) || (l==1)); //We accept constant factors, though they are useless
}

//For a future, maybe
/*
   void free_feature(struct feature *f){
        if(!(managed & f->state)){
                free(f->x);
                //Managed feature structs are R_alloced
                free(f);
        }
   }

   void free_features(u32 m,struct feature **f){
        for(u32 e=0;e<m;e++)
         if(f)
                 free_feature(f[e]);
   }
 */

//These function may error

struct feature* ingestSEXP_mle(u32 n,SEXP in){
 if(length(in)!=n) error("Incorrect feature length");
 if(isFactor(in) || isLogical(in)){
  //Integer-alike which needs collapsing into 1..n_levels
  struct feature *ans=(struct feature*)R_alloc(sizeof(struct feature),1);
  ans->state=managed|integer;
  ans->l=0; //Not collapsed
  ans->x=(void*)INTEGER(in);
  for(int e=0;e<n;e++){
   if(((u32*)(ans->x))[e]==NA_INTEGER){
    ans->state=ans->state&has_missing;
    //TODO: Till we figure out NAs
    error("NAs are not accepted");
   }
  }
  return(ans);
 }
 error("Only logical and factor inputs are acccepted with the MLE estimator");
 return(NULL); //Unreachable
}

struct feature* ingestSEXP_kt(u32 n,SEXP in){
 if(length(in)!=n) error("Incorrect feature length");
 if(n>65536) error("Kendall transformation covers only up to 2^16 elements");
 if(n<2) error("Kendall transformation requires at least 2 objects");
 struct feature *ans=(struct feature*)R_alloc(sizeof(struct feature),1);
 ans->state=managed;
 ans->l=0;

 if(isReal(in)){
  ans->x=REAL(in);
  double *v=(double*)(ans->x);

  for(int e=0;e<n;e++)
   if(ISNAN(v[e])) error("NAs nor NaNs are not allowed in input");
  return(ans);
 }

 //isBinaryFactor is from convert.h
 if(isInteger(in) || isLogical(in) || isOrdered(in) || isBinaryFactor(in)){
  ans->x=INTEGER(in);
  ans->state=managed|integer;
  u32 *v=(u32*)(ans->x);
  for(int e=0;e<n;e++)
   if(v[e]==NA_INTEGER) //Logical NA is NA_INTEGER
    error("NAs are not allowed in input");
  return(ans);
 }

 error("Only real, integer, logical, ordered and 2-level factor inputs are accepted with the KT estimator");
 return(NULL); //Unreachable
}

//TODO: Make this Fisher-Yates
bool* count_mask(u32 n,struct rng *rng,u32 count){
 bool *ans=malloc(sizeof(bool)*n);
 bool ds=false;
 if(count>n/2){
  count=n-count;
  ds=true;
 }
 for(u32 e=0;e<n;e++) ans[e]=ds;
 for(u32 e=0;e<count;){
  u32 i=random_index(rng,n);
  if(ans[i]==ds){
   ans[i]=!ds;
   e++;
  }
 }
 return(ans);
}

bool* boot_mask(u32 n,struct rng *rng,u32 *count){
 //ASSERT n>=2
 bool *ans=malloc(sizeof(bool)*n);
 *count=0;
 for(u32 e=0;e<n;e++) ans[e]=false;
 for(u32 e=0;(e<n) || ((*count)<2);e++){
  u32 i=random_index(rng,n);
  if(!ans[i]){
   ans[i]=true;
   (*count)++;
  }
 }

 return(ans);
}

u32 *which_mask(u32 n,bool *mask,u32 nn){
 u32 *ans=malloc(sizeof(u32)*nn);
 u32 ee=0;
 for(u32 e=0;e<n;e++) if(!mask || mask[e]){
   ans[ee]=e;
   ee++;
  }
 return(ans);
}

u32* produce_mle(u32 *in,struct ht *ht,u32 n,u32 nn,bool *mask,u32 *nx){//Mask & collapse
 u32 *ans=malloc(sizeof(u32)*nn);
 *nx=fillHtOneMasked(ht,n,in,mask,nn,ans,1);
 return(ans);
}

static inline u32 code(void *x,bool is_double,u32 e,u32 ee){
 u32 ans=0;
 if(is_double){
  double *xx=(double*)x;
  ans=(xx[e]<=xx[ee])+2*(xx[e]>=xx[ee]);
 }else{
  u32 *xx=(u32*)x;
  ans=(xx[e]<=xx[ee])+2*(xx[e]>=xx[ee]);
 }
 //ans is in coding "b" < -> 1, < ->2, = ->3
 //coding a (= -> 1, < -> 2, > -> 3)
 //b%3 (= -> 0 < -> 1 > -> 2) +1 (= -> 1 < -> 2 > -> 3) == a
 return(ans);
}


//Idx has to be made from mask by which_mask & re-used in the KT case;
// regular mask is good for sequential read, and here we need random
// access. Still, we don't use ht here to collapse.
u32* produce_kt(void *in,bool is_double,u32 nn,u32 *idx,u32 *nx){
 bool recode_to_a=code(in,is_double,0,1)==3;
 *nx=recode_to_a?1:2;
 //Coding a is =<>, can be = for const
 //Coding b is <>=, can be <> for no-tie
 //IF the first pair is =, we start from a and 1 level,
 // otherwise from b and 2 levels. If there is a change
 // later, we can only bump levels and left the code of
 // the former values intact.
 u32 no=nn*(nn-1);
 u32 *ans=malloc(sizeof(u32)*no);
 u32 eee=0;
 for(u32 e=0;e<nn;e++)
  for(u32 ee=0;ee<nn;ee++)
   if(e!=ee){
    u32 v=code(in,is_double,idx?idx[e]:e,idx?idx[ee]:ee);
    if(recode_to_a) v=(v%3)+1;              //Neat congruence trick
    if(v==3) *nx=3;
    ans[eee]=v;
    eee++;
   }
 return(ans);
}
