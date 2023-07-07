#' Overview plot of the vistla tree
#'
#' Plots a vistla tree, using layout derived by a Buchheim et al. extension of the standard Reingold-Tilford method.
#' The tree root is placed on the left, while the paths extend to the right, with all branches of the same depth at the same horizontal coordinate.
#' The path are sorted vertically, from strongest on top to weakest on the bottom.
#' Link weight indicates, by default, the link's score.
#' A feature name in parentheses indicates that is is only a way-point in a path to some other feature.
#' @note The graph is rendered using the grid graphics system, in a manner similar to \code{ggplot2}; the output of the \code{plot.vistla} function is only a grid graphical object, while the actual plotting is done when this object is printed or plotted.
#' Yet, said object can be used with other functions in the grid ecosystem for rendering into files, being edited, combined with other plots, etc.
#' @method plot vistla
#' @param x vistla, vistla hierarchy or vistla plot object.
#' @param slant arrange vertices in a slanted way. 
#'  Can be given as a number, possibly negative, indicating the amount of slant, or as \code{TRUE}, for an auto value.
#'  No slant is applied when set to 0 or omitted.
#' @param circular if given \code{TRUE}, switches to circular layout; alternatively, can be given two numbers, then the first one will be interpreted as an angle to fit the whole graph in (\code{2*pi} when using \code{TRUE}), and the second one as an initial angle offset (0 when using \code{TRUE}), which can be used to rotate the whole graph around the root.
#' Both angles are expected to be in radians.
#' It is recommended to add \code{asp=TRUE} parameter to make this layout truly circular, otherwise lines of equal depth are going to be elliptical.
#'  When \code{FALSE}, linear layout is enforced.
#' @param edge_col edge colour; can be given as vector, then mapping order adheres to the one in hierarchy object; please note that the edge towards first feature, the root, is not drawn, so the first element is effectively ignored.
#' If given as a function, it is called on the internally generated extended hierarchy object, and the result is used as an aesthetic.
#' @param edge_lwd edge width; behaves similarly to \code{edge_col}, yet also accepts special value \code{'scale'}, which triggers default scaling of edge width to be proportional to score.
#' @param edge_lty edge line-type; behaves similarly to \code{edge_col}.
#' @param label_text vertex label text, feature name by default.
#'  Behaves similarly to \code{edge_col}.
#' @param label_border_col vertex label border colour; behaves similarly to \code{edge_col}, can be set to 0 for no border.
#' @param label_border_lty vertex label border line-type; behaves similarly to \code{edge_col}, can be set to 0 for no border.
#' @param label_fill vertex label fill colour; behaves similarly to \code{edge_col}, can be set to 0 for no fill.
#' @param ... ignored.
#' @param asp1 if \code{TRUE}, scales on both axes are the same, like with \code{asp=1} in base graphics.
#' @param pmar Specifies margins as a fraction of graph size; expects a 4-element vector, in standard R bottom-left-top-right order.
#' @return Grid object with the graph.
#' @export
plot.vistla<-function(
 x,...,
 slant,
 circular,
 asp1=FALSE,
 pmar=c(0.05,0.05,0.05,0.05),
 edge_col=1,
 edge_lwd='scale',
 edge_lty=1,
 label_text=function(x) x$name,
 label_border_col=1,
 label_border_lty=function(x) ifelse(x$leaf,1,2),
 label_fill='white'
 ){

 #Plot layout planning 
 organise_plot(x)->trg

 #Vertex and edge generators
 wmar<-grid::unit(1,"strwidth","m")
 hmar<-grid::unit(1,"strheight","m")
 
 uu<-ifelse(asp1,'snpc','npc')
  
 vertex<-function(x,y,lab,lty,fill,col){
  if(asp1){
   x<-grid::unit(.5,'npc')+grid::unit(x-0.5,'snpc')
   y<-grid::unit(.5,'npc')+grid::unit(y-0.5,'snpc')
  }else{
   x<-grid::unit(x,'npc')
   y<-grid::unit(y,'npc')
  }
  lb<-grid::textGrob(label=lab,x=x,y=y)
  bkg<-grid::rectGrob(
   x=x,y=y,
   width=grid::unit(1,"grobwidth",lb)+wmar,
   height=grid::unit(1,"grobheight",lb)+hmar,
   gp=grid::gpar(fill=fill,col=col,lty=lty))
  grid::grobTree(bkg,lb)
 }

 edge<-function(x0,y0,x1,y1,lwd,lty,col){
  if(asp1){
   hh<-grid::unit(.5,'npc')+grid::unit(c(x0,x1)-0.5,'snpc')
   vv<-grid::unit(.5,'npc')+grid::unit(c(y0,y1)-0.5,'snpc')
  }else{
   hh<-grid::unit(c(x0,x1),'npc')
   vv<-grid::unit(c(y0,y1),'npc')
  }
  grid::linesGrob(
   hh,
   vv,
   gp=grid::gpar(lwd=lwd,col=col))
 }

 #Slant layout modifier
 if(missing(slant)) slant<-0
 if(identical(slant,TRUE)) slant<-0.5/max(trg$cidx)
 if(slant!=0) trg$x<-trg$x-slant*trg$cidx

 #Circular layout modifier
 if(missing(circular)) circular<-FALSE
 if(identical(circular,TRUE)) circular<-c(2*pi,0)
 if(is.numeric(circular) && length(circular)==2){
  trg$r<-trg$x+1
  #Add 0.5 units of margin on both sides on the tree
  trg$phi<-(trg$y-min(trg$y)+.5)/(diff(range(trg$y))+1)*circular[1]+circular[2]
  trg$x<-trg$r*cos(trg$phi)
  trg$y<-trg$r*sin(trg$phi)
 }

 stopifnot(is.numeric(pmar))
 stopifnot(length(pmar)==4)
 trg$x<-trg$x-min(trg$x)
 trg$x<-if(max(trg$x)>0) trg$x/max(trg$x) else 0.5
 trg$x<-trg$x*(1-pmar[2]-pmar[4])+pmar[2]
 trg$y<-trg$y-min(trg$y)
 trg$y<-if(max(trg$y)>0) trg$y/max(trg$y) else 0.5
 trg$y<-trg$y*(1-pmar[1]-pmar[3])+pmar[1]

 #Vertex styling
 trg$lab<-if(is.function(label_text)) label_text(trg) else label_text
 trg$label_border_col<-if(is.function(label_border_col)) label_border_col(trg) else label_border_col
 trg$label_border_lty<-if(is.function(label_border_lty)) label_border_lty(trg) else label_border_lty
 trg$label_fill<-if(is.function(label_fill)) label_fill(trg) else label_fill

  #Special case for root only
 if(nrow(trg)==1)
  return(
   as_vistla_plot(grid::grobTree(
    vertex(.5,.5,trg$lab,trg$lty,trg$fill,trg$col)
   ))
  )

 #Edge styling
 trg$wd<-(trg$score-min(trg$score,na.rm=TRUE))/(max(trg$score,na.rm=TRUE)-min(trg$score,na.rm=TRUE))*3+1
 trg$wd[!is.finite(trg$wd)]<-1

 trg$edge_col<-if(is.function(edge_col)) edge_col(trg) else edge_col
 trg$edge_lty<-if(is.function(edge_lty)) edge_lty(trg) else edge_lty
 trg$edge_lwd<-if(is.function(edge_lwd))
  edge_lwd(trg) else if(identical(edge_lwd,'scale')) trg$wd else edge_lwd

 #Edge ends
 trg$xp<-trg$x[trg$prv]
 trg$yp<-trg$y[trg$prv]

 #Order by score to have important branches on top
 trg<-trg[order(trg$score,na.last=FALSE),]
 trg$prv<-NULL
 
 #Make a grid object
 V<-do.call(grid::gList,lapply(1:nrow(trg),function(e){
  vertex(
   trg$x[e],trg$y[e],
   lab=trg$lab[e],
   col=trg$label_border_col[e],
   lty=trg$label_border_lty[e],
   fill=trg$label_fill[e])
 }))
 
 E<-do.call(grid::gList,lapply(2:nrow(trg),function(e){
  edge(trg$x[e],trg$y[e],trg$xp[e],trg$yp[e],col=trg$edge_col[e],lty=trg$edge_lty[e],lwd=grid::unit(trg$edge_lwd[e],'mm'))
 }))
 
 as_vistla_plot(grid::grobTree(E,V))
}

as_vistla_plot<-function(x){
 class(x)<-c("vistla_plot",class(x))
 x
}

#' @rdname plot.vistla
#' @method plot vistla_plot
#' @export
plot.vistla_plot<-function(x,...){
 grid::grid.newpage()
 grid::grid.draw(x)
}

#' @rdname plot.vistla
#' @method print vistla_plot
#' @export
print.vistla_plot<-function(x,...)
 plot.vistla_plot(x,...)


#' @method plot vistla_hierarchy
#' @export
plot.vistla_hierarchy<-function(x,...)
 plot.vistla(x,...)

organise_plot<-function(x){
 if(inherits(x,"vistla")) x<-hierarchy(x)
 stopifnot(inherits(x,"vistla_hierarchy"))

 x$prv[is.na(x$prv)]<-0
 x$idx<-1:nrow(x)
 split(x,x$prv)->xs
 xs<-lapply(xs,function(x){
   x$cidx<-1:nrow(x)
   x$nxtu<-c(0,x$idx[-length(x$idx)]) #Neighbour up
   x$nxtd<-c(x$idx[-1],0) #Neighbour down
   x$nxtf<-0 #First sub
   x$nxtfd<-0 #Last sub
   x
 })
 unsplit(xs,x$prv)->x
 for(e in xs) #Try xs here
   if(e$prv[1]!=0){
     x$nxtf[e$prv[1]]<-e$idx[1]
     x$nxtfd[e$prv[1]]<-e$idx[nrow(e)]
   }

 x<-TreeLayout(x,distance=1)

 x<-x[,c("name","score","depth","leaf","prv","cidx","y")]
 x$prv[1]<-NA

 #Flip vertically
 x$y<--x$y+x$y[1]
 x$x<-x$depth

 x
}

#R implementation of the pseudo-code from
# Buchheim, C., Jünger, M. and Leipert, S. (2006), Drawing rooted trees in linear time.
#  Softw: Pract. Exper., 36: 651-665. https://doi.org/10.1002/spe.713

TreeLayout<-function(x,distance=1){
 x$mod<-0
 x$thread<-0
 x$ancestor<-x$idx
 x$shift<-0
 x$change<-0
 x$prelim<-0
 defaultAncestor<-0

 x$y<-0

 NextUp<-function(v)
  if(x$nxtf[v]>0) x$nxtf[v] else x$thread[v]
   
 NextDown<-function(v)
  if(x$nxtfd[v]) x$nxtfd[v] else x$thread[v]

 Ancestor<-function(vim,v)
  if(x$prv[x$ancestor[vim]]==x$prv[v]) x$ancestor[vim] else defaultAncestor

 MoveSubtree<-function(wm,wp,shift){
  subtrees<-(x$cidx[wp]-x$cidx[wm])
  x$change[wp]<<-x$change[wp]-shift/subtrees
  x$shift[wp]<<-x$shift[wp]+shift
  x$change[wm]<<-x$change[wm]+shift/subtrees
  x$prelim[wp]<<-x$prelim[wp]+shift
  x$mod[wp]<<-x$mod[wp]+shift
 }

 Apportion<-function(v){
  w<-x$nxtu[v]
  if(w>0){
   vip<-v0p<-v
   vim<-w
   v0m<-x$nxtf[x$prv[vip]]
   sip<-x$mod[vip]
   s0p<-x$mod[v0p]
   sim<-x$mod[vim]
   s0m<-x$mod[v0m]
   while(NextDown(vim)!=0 && NextUp(vip)!=0){
    vim<-NextDown(vim)
    vip<-NextUp(vip)
    v0m<-NextUp(v0m)
    v0p<-NextDown(v0p)
    x$ancestor[v0p]<<-v
    shift<-(x$prelim[vim]+sim)-(x$prelim[vip]+sip)+distance
    if(shift>0){
     MoveSubtree(Ancestor(vim,v),v,shift)
     sip<-sip+shift
     s0p<-s0p+shift
    }
    sim<-sim+x$mod[vim]
    sip<-sip+x$mod[vip]
    s0m<-s0m+x$mod[v0m]
    s0p<-s0p+x$mod[v0p]
   } 
   if(NextDown(vim)!=0 && NextDown(v0p)==0){
    x$thread[v0p]<<-NextDown(vim)
    x$mod[v0p]<<-x$mod[v0p]+sim-s0p
   }
   if(NextUp(vip)!=0 && NextUp(v0m)==0){
    x$thread[v0m]<<-NextUp(vip)
    x$mod[v0m]<<-x$mod[v0m]+sip-s0m
    defaultAncestor<<-v
   }
  }
 }

 ExecuteShifts<-function(v){
  shift<-0
  change<-0
  w<-x$nxtfd[v]
  while(w>0){
   x$prelim[w]<<-x$prelim[w]+shift
   x$mod[w]<<-x$mod[w]+shift
   change<-change+x$change[w]
   shift<-shift+x$shift[w]+change
   w<-x$nxtd[w]
  }
 }

 FirstWalk<-function(v){
  if(x$nxtf[v]==0){
   x$prelim[v]<<-0
   w<-x$nxtu[v]
   if(w>0)
    x$prelim[v]<<-x$prelim[w]+distance
  }else{
   w<-defaultAncestor<<-x$nxtf[v]
   while(w!=0){
    FirstWalk(w)
    Apportion(w) 
    w<-x$nxtd[w]
   }
   ExecuteShifts(v)
   midpoint<-.5*(x$prelim[x$nxtf[v]]+x$prelim[x$nxtfd[v]])
   w<-x$nxtu[v]
   if(w>0){
    x$prelim[v]<<-x$prelim[w]+distance
    x$mod[v]<<-x$prelim[v]-midpoint
   }else{
    x$prelim[v]<<-midpoint
   }
  }
 }

 SecondWalk<-function(v,m){
  x$y[v]<<-x$prelim[v]+m
  w<-x$nxtf[v]
  while(w!=0){
   SecondWalk(w,m+x$mod[v])
   w<-x$nxtd[w]
  }
 }
 
 FirstWalk(1)
 SecondWalk(1,x$prelim[1])
 x
}
