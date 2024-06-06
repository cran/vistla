enum estimator {
 mle=1,
 kt=2
};

enum estimator verify_estimator(u32 x){
 if(x<1 || x>2) error("Unknown estimator");
 return(x);
}

static inline void t2ij(u32 t,u32 *ip,u32 *jp){
 double k=t;
 double jf=(sqrt(8.*k+1.)-1.)/2.;
 u32 j=floor(jf)+1;
 u32 i=t-j*(j-1)/2;
 ip[0]=i;
 jp[0]=j;
}

/*  128  64  32  16   8   4   2   1  *
* |res|res|res|roa|hld|hlu|fwd|bkw| */
enum flow {
 backward=1,
 forward=2,
 hilldown=8,
 hillup=4,
 noroam=16
};

enum flow verify_flow(u32 x){
 if(x>31) error("Wrong value of the flow");
 if((x&hillup) && (x&hilldown))error("Cannot hill up and down at the same time");
 if(((x&hillup)|(x&hilldown))&(x&noroam)){
  warning("Force-path is redundant with up/down.");
  x&=(!noroam);
 }
 return(x);
}

/*
   For multi-thread, vistla uses its own PRNG, PCG, which is seeded from
   R's PRNG each time an rFerns model is built.
   The generator is PCG32 by M.E. O'Neill https://www.pcg-random.org
 */

struct rng {
 u64 state;
 u64 stream;
};

u32 random_int(struct rng *rng){
 rng->state=rng->state*6364136223846793005+rng->stream;
 u32 rot=(rng->state)>>59;
 u32 s=(((rng->state)>>18)^(rng->state))>>27;
 return((s<<((-rot)&31))|(s>>rot));
}

void set_rng(struct rng *rng,u64 seed,u64 stream){
 rng->state=rng->stream=stream*2+1;
 rng->state+=seed;
 random_int(rng);
}

//Sync PRNG with R; R will feel only two numbers were generated
void set_from_r(struct rng *rng){
 GetRNGstate();
 u64 a=(u32)(((double)(~((u32)0)))*unif_rand());
 u64 b=(u32)(((double)(~((u32)0)))*unif_rand());
 PutRNGstate();
 a=(a<<32)+b;
 set_rng(rng,a,0);
}

//Fast & unbiased algorithm by Daniel Lemire https://arxiv.org/pdf/1805.10941.pdf
u32 random_index(struct rng *rng,u32 upto){
 u32 x=random_int(rng);
 u64 m=((u64)x)*((u64)upto);
 u32 l=((u32)m);
 if(l<upto){
  u32 t=(-upto)%upto;
  while(l<t){
   x=random_int(rng);
   m=((u64)x)*((u64)upto);
   l=((u32)m);
  }
 }
 return(m>>32);
}

