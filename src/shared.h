enum estimator{
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

