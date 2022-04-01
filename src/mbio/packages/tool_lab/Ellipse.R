# zouguaniqng 20181213  
library(MASS)
library(getopt)


ellipse <- function(center, shape, radius, log="", center.pch=19, center.cex=1.5, segments=51, draw=FALSE, add=draw, 
		xlab="", ylab="",  lwd=2, fill=FALSE, fill.alpha=0.3,
		grid=TRUE, ...) {
	
	logged <- function(axis=c("x", "y")){
		axis <- match.arg(axis)
		0 != length(grep(axis, log))
	}
	
	if (! (is.vector(center) && 2==length(center))) stop("center must be a vector of length 2")
	if (! (is.matrix(shape) && all(2==dim(shape)))) stop("shape must be a 2 by 2 matrix")
	if (max(abs(shape - t(shape)))/max(abs(shape)) > 1e-10) stop("shape must be a symmetric matrix")
	angles <- (0:segments)*2*pi/segments 
	unit.circle <- cbind(cos(angles), sin(angles)) 
	Q <- chol(shape, pivot=TRUE)
	order <- order(attr(Q, "pivot"))
	ellipse <- t( center + radius*t( unit.circle %*% Q[,order]))
	colnames(ellipse) <- c("x", "y")
	if (logged("x")) ellipse[, "x"] <- exp(ellipse[, "x"])  #做对数处理
	if (logged("y")) ellipse[, "y"] <- exp(ellipse[, "y"])
	
	if (draw) {
		if (add) {
			lines(ellipse, lwd=lwd, ...)  
			if (fill) polygon(ellipse ) 
		}
		else {
			plot(ellipse, type="n", xlab = xlab, ylab = ylab, ...) 
			if(grid){
				grid(lty=1, equilogs=FALSE)
				box()}
			lines(ellipse, lwd=lwd, ... )  
			if (fill) polygon(ellipse ) 
		} 	
		if ((center.pch != FALSE) && (!is.null(center.pch))) points(center[1], center[2], pch=center.pch, cex=center.cex ) 
	}
	v<-radius*Q[,order]
    #ret <- list(c=center,v=v, u=unit.circle, r=radius)
	#ret <- list(m1=center[1],m2=center[2],c11=v[1],c12=v[2],c21=v[3],c22=v[4])
	#ret <- list(d1=center[1],d2=v[1],d3=v[2],d4=center[2],d5=v[3],d6=v[4])
	cent <-as.matrix(center)
	ret <-c(cent[1],v[1],v[2],cent[2],v[3],v[4])
    invisible(ret)
}


dataEllipse <- function(x, y, groups,
    group.labels=group.levels, weights, log="", levels=c(0.5, 0.95), center.pch=19, 
    center.cex=1.5, draw=FALSE,plot.points=draw, add=!plot.points, segments=51, robust=FALSE, 
    xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), 
    pch=if (missing(groups)) 1 else seq(group.levels),
    lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE, ...) {
	
    if(missing(y)){
        if (is.matrix(x) && ncol(x) == 2) {
            if (missing(xlab)) xlab <- colnames(x)[1]
            if (missing(ylab)) ylab <- colnames(x)[2]
            y <- x[,2]
            x <- x[,1]
        }
        else stop("x and y must be vectors, or x must be a 2 column matrix")
    }
    else if(!(is.vector(x) && is.vector(y) && length(x) == length(y))){
        stop("x and y must be vectors of the same length")
    }
    if (missing(weights)) weights <- rep(1, length(x))
    if (length(weights) != length(x)) stop("weights must be of the same length as x and y")
    if (!missing(groups)){
        xlab
        ylab
        if (!is.factor(groups)) stop ("groups must be a factor")
        if (!(length(groups) == length(x))) stop ("groups, x, and y must all be of the same length")
        labels <- seq(length(x))
        valid <- complete.cases(x, y, groups)
        x <- x[valid]
        y <- y[valid]
        weights <- weights[valid]
        groups <- groups[valid]
        labels <- labels[valid]
        group.levels <- levels(groups)
        result <- vector(length(group.levels), mode="list")
        names(result) <- group.levels
        if(draw) {
            if (!add) {
                plot(x, y, type="n", xlab=xlab, ylab=ylab,  ...) 
                if(grid){
                    grid(lty=1, equilogs=FALSE)
                    box()
                }
            }
        }
       
        for (lev in 1:length(group.levels)){
            level <- group.levels[lev]
            sel <- groups == level
            #id_lev_labels <- labels[sel]
          
            result[[lev]] <- dataEllipse(x[sel], y[sel],
                weights=weights[sel], log=log, levels=levels, center.pch=center.pch,
                center.cex=center.cex, draw=draw, plot.points=plot.points, add=TRUE, segments=segments,
                robust=robust, pch=pch[lev], lwd=lwd, fill=fill, fill.alpha=fill.alpha)
        }
        return(invisible(result))
    }
  
    if(draw) {
        if (!add) {
            plot(x, y, type="n", xlab=xlab, ylab=ylab,  ...) 
            if(grid){
                grid(lty=1, equilogs=FALSE)
                box()}
        }
        if (plot.points)  points(x, y,  pch=pch[1], ...)
    }
    dfn <- 2
    dfd <- length(x) - 1
    if (robust) {
        use <- weights > 0
        v <- MASS::cov.trob(cbind(x[use], y[use]), wt=weights[use])
        shape <- v$cov
        center <- v$center
    }
    else {
        v <- cov.wt(cbind(x, y), wt=weights)
        shape <- v$cov
        center <- v$center
    }
    result <- vector("list", length=length(levels))
    names(result) <- levels
    for (i in seq(along=levels)) {
        level <- levels[i]
        radius <- sqrt(dfn * qf(level, dfn, dfd ))
        result[[i]] <- ellipse(center, shape, radius, log=log,
            center.pch=center.pch, center.cex=center.cex, segments=segments, lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
    }
    invisible(if (length(levels) == 1) result[[1]] else result)
}
    

	
	
command <- matrix(c("group", "g", 1, "character",
					"infile", "f", 1, "character",
					"level", "l", 1, "character",
					"out", "o", 1, "character",
					"meta", "m", 2, "character"), byrow = T, ncol = 4)  # add meta by houshuang 20190924
args <- getopt(command)

#args <- commandArgs(T)
#pfile <- args[1]
#gfile <- args[2]
#level <- args[3]
#draw.file <- args[4]
pfile <- args$infile
gfile <- args$group
level <- as.numeric(args$level)
draw.file <- args$out

data <- read.table(pfile, header=TRUE, sep='\t', comment.char = "")

# by houshuang 20190924 calculate two or three dimensions for meta >>>
if (!is.null(args$meta)){
    if (length(colnames(data)) > 3) {
        data <- data[,1:4]
    } else {
        data <- data[,1:3]
    }
}
# <<<

len <- length(names(data))
new_names <- paste0(replicate(len-1, "pc"),as.character(1:(len-1)))
names(data) <- c('Sample_ID',new_names)
if (!is.null(gfile)) {
    group <- read.table(gfile, header=TRUE, sep='\t', comment.char = "")
    names(group) <- c('Sample_ID','group')
    tmp <- merge(data, group, by='Sample_ID', sort=FALSE)
    group <- tmp$group
}

head <- names(data)
len <- length(head)
m1 <- c()
c11 <- c()
c12 <- c()
m2 <- c()
c21 <- c()
c22 <- c()
pc <- c()
g <- c()
for (id1 in 2:(len-1)){
	h1 <- head[id1]
	for (id2 in (id1+1):len){
		h2 <- head[id2]
		#result <- dataEllipse(as.vector(data[h1]), as.vector(data[h2]), group, levels = level)
        x <- as.vector(unlist(data[h1]))
        y <- as.vector(unlist(data[h2]))
		if (is.null(gfile)){
			result <- dataEllipse(x, y, levels = level)
			m1 <-c(m1, result[1])
			c11 <- c(c11, result[2])
			c12 <- c(c12, result[3])
			m2 <- c(m2, result[4])
			c21 <- c(c21, result[5])
			c22 <- c(c22, result[6])
			pc <- c(pc, paste(h1,h2,sep='_'))
            g <- c(g, 'All')
		}else{
			result <- dataEllipse(x, y, group, levels = level)
            name <- names(result)
            for (r in name){
                rr <- as.vector(unlist(result[r]))
                m1 <-c(m1, rr[1])
                c11 <- c(c11, rr[2])
                c12 <- c(c12, rr[3])
                m2 <- c(m2, rr[4])
                c21 <- c(c21, rr[5])
                c22 <- c(c22, rr[6])
                pc <- c(pc, paste(h1, h2, sep='_'))
                g <- c(g, r)
            }		
		}
    }
}



out <- data.frame(PC=pc, group=g, m1=m1, c11=c11, c12=c12, m2=m2, c21=c21, c22=c22)
write.table(out, draw.file, sep='\t', quote=FALSE, row.names=FALSE)
#print(warnings())
