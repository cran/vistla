set.seed(1)

library(vistla)

data(chain)

vistla(Y~.,data=chain,threads=2)->V

expected_paths<-list(
 M3=c("M3","M4","Y"),
 M2=c("M2","M3","M4","Y"),
 M1=c("M1","M2","M3","M4","Y"),
 X=c("X","M1","M2","M3","M4","Y")
)

stopifnot(identical(
 paths(V),
 expected_paths
))

stopifnot(identical(
 hierarchy(V)$name,
 rev(names(chain))
))

