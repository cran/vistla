subset_links<-function(links,subset)
 match(
  links,
  (1:length(links))[subset],
 )[subset]

subset_tree<-function(tree,subset){
 st<-tree[subset,]
 st$prv<-subset_links(tree$prv,subset)
 st
}

