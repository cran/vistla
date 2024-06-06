set.seed(1)

library(vistla)

data(chain)

vistla(Y~.,data=chain,threads=2)->V

expected_paths<-list(
 M2=c("M2","M1","Y"),
 M3=c("M3","M2","M1","Y"),
 M4=c("M4","M3","M2","M1","Y"),
 T=c("T","M4","M3","M2","M1","Y")
)

stopifnot(identical(
 paths(V),
 expected_paths
))

stopifnot(identical(
 hierarchy(V)$name,
 names(chain)
))

