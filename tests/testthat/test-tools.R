for(seed in 1:20)
 test_that(sprintf("subset_links works for seed %d",seed),{
  set.seed(seed)
  N<-20
  M<-15
  links<-sample(N)
  X<-data.frame(
   a=letters[1:N],
   al=letters[links],
   links=links
  )
  subset<-if(runif(1)>.5) sample(N,M,replace=TRUE) 
   else 
    sample(c(TRUE,FALSE),M,replace=TRUE)
  Xs<-X[subset,]
  Xs$nlinks<-subset_links(X$links,subset)
  Xs$alt<-Xs$a[Xs$nlinks]

  Xs$al[!is.na(Xs$alt)]->a
  Xs$alt[!is.na(Xs$alt)]->b

  expect_equal(a,b)
 })

test_that("Export to dot works",{
 data(chain)
 vistla(chain[,-6],chain$Y)->vc
 write.dot(vc)->a
 EXPECTED<-c(
  "graph {",
  "\toverlap=\"prism\";",
  "\tsplines=true;",
  "\tY [shape=egg label=\"Y\"];",
  "\tM4 [shape=ellipse label=\"M4\"];",
  "\tM3 [shape=box label=\"M3\"];",
  "\tM2 [shape=box label=\"M2\"];",
  "\tM1 [shape=box label=\"M1\"];",
  "\tX [shape=box label=\"X\"];",
  "\tY -- M4 [penwidth=3.000];",
  "\tM4 -- M3 [penwidth=3.000];",
  "\tM3 -- M2 [penwidth=2.126];",
  "\tM2 -- M1 [penwidth=2.126];",
  "\tM1 -- X [penwidth=2.101];",
  "}"
 )
 expect_equal(a,EXPECTED)
})

test_that("Prune works",{
 data(chain)
 v<-vistla(Y~.,data=chain)
 vp<-v
 vp$tree<-subset_tree(vp$tree,vp$tree$leaf|vp$tree$used)
 expect_equal(vp,prune(v))
 expect_equal(
  hierarchy(prune(v,targets="M3"))$name,
  c("Y","M4","M3")
 )
 expect_lt(min(v$tree$score),0.3)
 expect_gt(
  min(prune(v,iomin=0.3)$tree$score),
  0.3
 )
 expect_equal(
  nrow(hierarchy(prune(v,iomin=100))),
  1
 )
 expect_equal(
  nrow(hierarchy(prune(v,targets=c()))),
  1
 )
})

test_that("Prune reports errors",{
 data(chain)
 v<-vistla(Y~.,data=chain)
 expect_error(prune(v,targets="UnknownVar"),"Unknown names in targets")
 expect_error(prune(prune(v,targets="M3"),targets="M4"),"Prune can only remove targets")
 vp<-prune(v,iomin=0.3)
 expect_error(prune(vp,iomin=0.2),"Prune can only increase")
 expect_error(prune(iris))
})

test_that("Agreement works",{
 data(chain)
 va<-vistla(Y~.,data=chain)
 vb<-vistla(Y~.,data=chain[,sample(6)])
 expect_gt(agreement(va,vb),.99)
 expect_error(agreement(list(va)),"Comparison requires")
 agreement(list(v1=va,v2=vb,v3=va)->s)->r
 expect_equal(r,matrix(rep(1,9),3,dimnames=list(names(s),names(s))))
 expect_equal(names(agreement(s,raw=TRUE)),c("a","b","c",names(s)))
})

synth_tree<-data.frame(
 name=sprintf("v%d",1:15),
 prv=c(NA,1,2,3,3,2,6,6,1,9,10,10,9,13,13),
 depth=c(-1,0,1,2,2,1,2,2,0,1,2,2,1,2,2),
 score=15:1,
 leaf=  c(F,T,T,T,T,T,T,T,T,T,T,T,T,T,T)
)
class(synth_tree)<-c("vistla_hierarchy",class(synth_tree))
 
test_that("TreeLayout works",{
 organise_plot(synth_tree)->st
 expect_equal(st$y,c(0,2,3,3.5,2.5,1,1.5,0.5,-2,-1,
  -0.5,-1.5,-3,-2.5,-3.5))
})

test_that("Plot doesn't error",{
 plot(synth_tree,circular=TRUE)->L
 expect_true(inherits(L,"gTree"))
})
