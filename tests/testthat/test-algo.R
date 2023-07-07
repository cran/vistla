verify_tree<-function(tree){
 diff(tree$depth)->dd
 expect_true(all(dd<=1))
 expect_true(is.na(tree$prv[1]))
 p<-tree$prv[-1]-1
 expect_true(all(p<1:length(p)))
}

verify_leaf<-function(x){
 x$tree->tr
 wrongc<-c()
 for(x in split(tr,tr$c)){
  if(sum(x$leaf)!=1){
   wrongc<-c(wrongc,x$c[1])
  }else if(x$score[x$leaf]!=max(x$score)) wrongc<-c(wrongc,x$c[1])
 }
 unique(wrongc)
}

rebuild_used<-function(x){
 x$tree->tr
 u<-rep(FALSE,nrow(tr))
 #Start with leaves
 to_use<-which(tr$leaf)
 while(length(to_use)>0){
  u[to_use]<-TRUE
  tr$prv[to_use]->to_use
  unique(to_use[!is.na(to_use)])->to_use
 }
 u
}

test_that("vistla works on chain data",{
 data(chain)
 vistla(Y~.,data=chain)->ans
 expect_error(
  vistla(Y~.,data=chain,targets="J"),
  "Unknown variables specified as targets"
 )
 expect_named(ans,c("tree","mi","miY","yn","iomin"))
 expect_s3_class(ans,"vistla")
 expect_named(ans$tree,c("a","b","c","score","depth","leaf","used","prv"))
 expect_equal(verify_leaf(ans),c())
 expect_equal(rebuild_used(ans),ans$tree$used)
 verify_tree(hierarchy(ans))
})

test_that("vistla works on junction data",{
 data(junction)
 expect_equal(
  paths(vistla(A2~.,data=junction)),
  list(A1=c("A1","X","A2"))
 )
 expect_equal(
  paths(vistla(B2~.,data=junction)),
  list(B1=c("B1","X","B2"))
 )
 expect_equal(
  paths(vistla(A1~.,data=junction)),
  list(A2=c("A2","X","A1"))
 )
 expect_equal(
  paths(vistla(B1~.,data=junction)),
  list(B2=c("B2","X","B1"))
 )
})

test_that("path extraction works",{
 vistla(Y~.,data=chain,iomin=0.2)->v_nox
 vistla(Y~.,data=chain)->v
 expect_equal(
  path_to(v_nox,"X"),
  character(0)
 )
 expect_equal(
  path_to(v_nox,"X",detailed=TRUE),
  data.frame(
   a=character(0),
   b=character(0),
   c=character(0),
   score=numeric(0)
  )
 )
 path_to(v,"X",detailed=TRUE)->dp
 expect_equal(
  dp[,-4],
  data.frame(
   a=c("M2","M3","M4","Y"),
   b=c("M1","M2","M3","M4"),
   c=c("X","M1","M2","M3")
  )
 )
 expect_equal(
  round(dp$score,3),
  c(0.197,0.201,0.201,0.308)
 )
 expect_error(
  path_to(v,"J"),
  "Target feature J was not"
 )
})

test_that("vistla works on iris data",{
 vistla(iris[,-5],iris[,5])->ans
 expect_named(ans,c("tree","mi","miY","yn","iomin"))
 expect_s3_class(ans,"vistla")
 expect_named(ans$tree,c("a","b","c","score","depth","leaf","used","prv"))
 expect_equal(verify_leaf(ans),c())
 expect_equal(rebuild_used(ans),ans$tree$used)
 verify_tree(hierarchy(ans))
 organise_plot(ans)->xp
 expect_equal(xp$prv,c(NA,1,2,3,1,5))
 expect_equal(xp$y,c(0,.5,.5,.5,-.5,-.5))
})

test_that("vistla works properly with prune",{
 data(chain)
 vistla(Y~.,data=chain,targets="M3")->ans
 expect_named(ans,c("tree","mi","miY","yn","iomin","targets"))
 expect_s3_class(ans,"vistla")
 expect_named(ans$tree,c("a","b","c","score","depth","leaf","used","prv"))
 expect_equal(verify_leaf(ans),c())
 expect_equal(rebuild_used(ans),ans$tree$used)
 verify_tree(hierarchy(ans))

 vistla(Y~.,data=chain)->ref
 expect_equal(prune(ref,targets="M3"),ans)
})

test_that("vistla verbose prints something",{
 data(chain)
 expect_output(vistla(Y~.,data=chain,verbose=TRUE))
})

test_that("empty trees do not explode",{
 data(chain)
 vistla(Y~.,data=chain,iomin=1)->empty
 he<-hierarchy(empty)
 expect_equal(he$name,empty$yn)
 expect_equal(names(he),c("name","score","depth","leaf","prv"))
 expect_output(print(empty))
 expect_equal(nrow(as.data.frame(empty)),0)
})

test_that("mi matrix export works",{
 data(chain)
 vistla(Y~.,data=chain)->v
 mi_scores(v)->mv
 expect_equal(rownames(mv),colnames(mv))
 expect_true(all(is.na(diag(mv))))
 expect_equal(mv,t(mv))
})

test_that("leaf score matrix export works",{
 data(chain)
 vistla(Y~.,data=chain)->v
 leaf_scores(v)->mv
 expect_equal(rownames(mv),colnames(mv))
 expect_true(all(is.na(diag(mv))))
 expect_true(all(mv<v$mi,na.rm=TRUE))
})

test_that("printing works on chain",{
 data(chain)
 vistla(Y~.,data=chain)->v

 expect_output(print(v),"Vistla tree")
 expect_output(print(v),"M2 \\(0\\.2\\) ~ M3")
 expect_output(print(v,1),"3 more")
 
 hv<-hierarchy(v)
 expect_output(print(hv),"Vistla hierarchy")
 expect_output(print(hv),"| \\+ M2 \\(0\\.2\\)")
 
 pv<-prune(v,targets="M3")
 expect_output(print(pv),"rooted in Y, pruned")
 expect_output(print(pv),"specified target")
 
 pv<-prune(v,targets=c("M3","M2"))
 expect_output(print(pv),"specified targets")
})

test_that("vistla rejects bad input class",{
 expect_error(vistla(list()),"Expecting a formula or a data.frame as an input, got list")
 expect_error(vistla(matrix(1:12,4)),"Expecting a formula or a data.frame as an input, got matrix/array")
})

test_that("flow works",{
 expect_error(flow("eee"))
 expect_error(flow(pi))
 expect_output(print(flow(from=TRUE,into=TRUE,down=FALSE,up=FALSE,forcepath=TRUE)),"both!")
 expect_output(print(flow("from",from=FALSE)),"from") #code overrides
 expect_output(print(flow(6L)),"intoup")
})

test_that("targets are robust",{
 a<-vistla(A1~.,data=junction,targets="B1")
 expect_output(print(a),"No paths")
})
