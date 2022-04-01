# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# __version__ = 'v1.0'
# __last_modified__ = '20151111'
import sys
import os


def pan_core(otutable, dowhat, groupfile='none', work_dir=None):
    """
    计算pan或core OTU
    """
    if not work_dir:
        output = os.path.dirname(otutable)
    else:
        output = work_dir
    if dowhat not in ("pan", "core"):
        raise Exception("只能是计算core或pan!")
        sys.exit()
    cmd = "setwd('" + output + "')"
    cmd += '''

    library(vegan)
    pan_core <- function (comm, dowhat = "pan" , method = "random", permutations = 100, conditioned = TRUE, gamma = "jack1", w = NULL, subset, ...)
    {
      METHODS <- c("collector", "random", "exact", "rarefaction", "coleman")
      method <- match.arg(method, METHODS)
      if (!is.null(w) && !(method %in% c("random", "collector")))
        stop(gettextf("weights 'w' can be only used with methods 'random' and 'collector'"))
      if (!missing(subset)) {
        comm <- subset(comm, subset)
        w <- subset(w, subset)
      }
      x <- comm
      x <- as.matrix(x)
      x <- x[, colSums(x) > 0, drop = FALSE]
      n <- nrow(x)
      p <- ncol(x)
      if (p == 1) {
        x <- t(x)
        n <- nrow(x)
        p <- ncol(x)
      }
      accumulator <- function(x, ind) {
        rowSums(apply(x[ind, ], 2, switch(dowhat, pan = cumsum, core = cumprod)) > 0)
      }
      specaccum <- sdaccum <- sites <- perm <- NULL
      if (n == 1 && method != "rarefaction")
        message("No actual accumulation since only 1 site provided")
      switch(method, collector = {
        sites <- 1:n
        xout <- weights <- cumsum(w)
        specaccum <- accumulator(x, sites)
      }, random = {
        permat <- vegan:::getPermuteMatrix(permutations, n)
        perm <- apply(permat, 1, accumulator, x = x)
        if (!is.null(w)) weights <- apply(permat, 1, function(i) cumsum(w[i]))
        sites <- 1:n
        if (is.null(w)) {
          specaccum <- apply(perm, 1, mean)
          sdaccum <- apply(perm, 1, sd)
        } else {
          sumw <- sum(w)
          xout <- seq(sumw/n, sumw, length.out = n)
          intx <- sapply(seq_len(n), function(i) approx(weights[, i], perm[, i], xout = xout)$y)
          specaccum <- apply(intx, 1, mean)
          sdaccum <- apply(intx, 1, sd)
        }
      }, exact = {
        freq <- colSums(x > 0)
        freq <- freq[freq > 0]
        f <- length(freq)
        ldiv <- lchoose(n, 1:n)
        result <- array(dim = c(n, f))
        for (i in 1:n) {
          result[i, ] <- ifelse(n - freq < i, 0, exp(lchoose(n - freq, i) - ldiv[i]))
        }
        sites <- 1:n
        specaccum <- rowSums(1 - result)
        if (conditioned) {
          V <- result * (1 - result)
          tmp1 <- cor(x > 0)
          ind <- lower.tri(tmp1)
          tmp1 <- tmp1[ind]
          tmp1[is.na(tmp1)] <- 0
          cv <- numeric(n)
          for (i in 1:n) {
            tmp2 <- outer(sqrt(V[i, ]), sqrt(V[i, ]))[ind]
            cv[i] <- 2 * sum(tmp1 * tmp2)
          }
          V <- rowSums(V)
          sdaccum <- sqrt(V + cv)
        } else {
          Stot <- specpool(x)[, gamma]
          sdaccum1 <- rowSums((1 - result)^2)
          sdaccum2 <- specaccum^2/Stot
          sdaccum <- sqrt(sdaccum1 - sdaccum2)
        }
      }, rarefaction = {
        freq <- colSums(x)
        freq <- freq[freq > 0]
        tot <- sum(freq)
        ind <- round(seq(tot/n, tot, length = n))
        result <- matrix(NA, nrow = 2, ncol = n)
        for (i in 1:n) {
          result[, i] <- rarefy(t(freq), ind[i], se = TRUE)
        }
        specaccum <- result[1, ]
        sdaccum <- result[2, ]
        sites <- ind/tot * n
      }, coleman = {
        freq <- colSums(x > 0)
        result <- array(dim = c(n, p))
        for (i in 1:n) {
          result[i, ] <- (1 - i/n)^freq
        }
        result <- 1 - result
        sites <- 1:n
        specaccum <- apply(result, 1, sum)
        sdaccum <- sqrt(apply(result * (1 - result), 1, sum))
      })
      out <- list(call = match.call(), method = method, sites = sites, richness = specaccum, sd = sdaccum, perm = perm)
      if (!is.null(w)) {
        out$weights <- weights
        out$effort <- xout
      }
      if (method == "rarefaction")
        out$individuals <- ind
      if (method == "random")
        attr(out, "control") <- attr(permat, "control")
      class(out) <- "specaccum"
      out
    }

    data = read.table("''' + otutable + '''",sep="\\t",check.names=F,head=T,comment.char="")
    rownames(data) <-data[,1]
    data <- data[,-1]
    gdata <-list()
    gfile = "''' + groupfile + '''"
    dowhat = "''' + dowhat + '''"
    g = 1
    groups = "none"
    if (gfile !="none"){
        map <- read.table(gfile,sep="\\t",head=F,check.names=F, colClasses = c("character"))
        #map <- as.matrix(map)
        groups <- as.character(unique(map[,2]))
        for(i in 1:length(groups))
        {
            samples <- as.character(map[which(map[,2] %in% groups[i]),1])
            data_pick <- data[,which(colnames(data) %in% samples)]
            gdata[[i]] <- data_pick
        }
        g = length(groups)
    }else{
        gdata[[1]] <- data
    }
    richness_out<-matrix(nrow=g,ncol=ncol(data)+1)
    if (length(gdata) >0){
        for (i in 1:length(gdata)){
            sp <- pan_core(t(gdata[[i]]), dowhat)
            groupname = "all"
            if (gfile !="none"){
              groupname = groups[i]
            }
            richness_out[i,1:(length(sp$richness)+1)] <- c(groupname,sp$richness)
            write.table(sp$perm,paste(groupname,".perm.txt",sep=""),sep="\\t",col.names=F,row.names=F,quote=F)
        }
    }else{
        stop("no proper data to run.")
    }
    colnames(richness_out) <- c("group",seq(1,ncol(data)))
    write.table(richness_out,paste(dowhat,".richness.xls",sep=""),sep="\\t",col.names=T,row.names=F,quote=F)
    '''
    output = os.path.join(output, dowhat + ".r")
    with open(output, "wb") as w:
        w.write(cmd)
    return output
    # os.system("Rscript tmp.r")


# test
if __name__ == "__main__":
    pan_core("otu_table.xls", dowhat="core")
    pan_core("otu_table.xls", dowhat="pan")
