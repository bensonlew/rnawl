#!/usr/bin/env Rscript

# load package
library(getopt)
library(maSigPro)
library(mclust)

# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing normalized gene expression data",
  "design", "d", 1, "character", "input file describing experimental design",
  "cluster", "k", 1, "integer", "number of clusters for data partioning",
  "method", "m", 1, "character", "clustering method for data partioning",
  "output", "o", 1, "character", "output directory for results",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# function



self.PlotGroups <-function (data, edesign = NULL, time = edesign[, 1], groups = edesign[,
    c(3:ncol(edesign))], repvect = edesign[, 2], show.fit = FALSE,
    dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05,
    nvar.correction = FALSE, summary.mode = "median", show.lines = TRUE,
    groups.vector = NULL, xlab = "Time", ylab = "Expression value",
    cex.xaxis = 1, ylim = NULL, main = NULL, cexlab = 0.8, legend = TRUE,
    sub = NULL, item = NULL) {
    if (!is.vector(data)) {
        if (summary.mode == "representative") {
            distances <- apply(as.matrix(dist(data, diag = TRUE,
                upper = TRUE)), 1, sum)
            representative <- names(distances)[distances == min(distances)]
            yy <- as.numeric(data[rownames(data) == representative,
                ])
            sub <- paste("Representative:", representative)
        }
        else if (summary.mode == "median") {
            yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
            if (is.null(sub)) {
                sub <- paste("Median profile of", nrow(data),
                  item, sep = " ")
            }
        }
        else stop("not valid summary.mode")
        if (dim(data)[1] == 1) {
            sub <- rownames(data)
        }
    }
    else if (length(data) != 0) {
        yy <- as.numeric(data)
        sub <- rownames(data)
    }
    else stop("empty data")
    if (is.null(ncol(groups))) {
        ncol = 1
        legend = FALSE
        codeg = "group"
    }
    else {
        ncol = ncol(groups)
        codeg <- as.character(colnames(groups))
    }
    reps <- i.rank(repvect)
    y <- vector(mode = "numeric", length = length(unique(reps)))
    x <- vector(mode = "numeric", length = length(unique(reps)))
    g <- matrix(nrow = length(unique(reps)), ncol = ncol)
    for (k in 1:length(y)) {
        y[k] <- mean(yy[reps == k], na.rm = TRUE)
        x[k] <- mean(time[reps == k])
        for (j in 1:ncol) {
            g[k, j] <- mean(groups[reps == k, j])
        }
    }
    if (is.null(ylim))
        ylim = c(min(as.numeric(yy), na.rm = TRUE), max(as.numeric(yy),
            na.rm = TRUE))
    abcissa <- x
    xlim = c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE) *
        1.3)
    color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
    color2 <- groups
    for (j in 1:ncol) {
        color2[, j] <- color2[, j] * j
    }
    color2 <- as.vector(apply(color2, 1, sum) + 1)
    plot(x = time, y = yy, pch = 21, xlab = xlab, ylab = ylab,
        xaxt = "n", main = main, sub = sub, ylim = ylim, xlim = xlim,
        cex = cexlab, col = color2)
    axis(1, at = unique(abcissa), labels = unique(abcissa), cex.axis = cex.xaxis)
    if (show.fit) {
        rm <- matrix(yy, nrow = 1, ncol = length(yy))
        rownames(rm) <- c("ratio medio")
        colnames(rm) <- rownames(dis)
        fit.y <- T.fit(rm, design = dis, step.method = step.method,
            min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction)
        betas <- fit.y$coefficients
    }
    if ( ncol(groups) <= 5 ){
        for (i in 1:ncol(groups)) {
            group <- g[, i]
            if ((show.fit) && !is.null(betas)) {
                li <- c(2:6)
                a <- reg.coeffs(coefficients = betas, groups.vector = groups.vector,
                    group = colnames(groups)[i])
                a <- c(a, rep(0, (7 - length(a))))
                curve(a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) +
                    a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5), from = min(time),
                    to = max(time), col = color1[i], add = TRUE,
                    lty = li[i])
            }
            if (show.lines) {
                lx <- abcissa[group != 0]
                ly <- y[group != 0]
                ord <- order(lx)
                lxo <- lx[ord]
                lyo <- ly[ord]
                lines(lxo, lyo, col = color1[i])
            }
        }
    }
    else {
        for (i in 1:ncol(groups)) {
            group <- g[, i]
            if ((show.fit) && !is.null(betas)) {
                li <- c(2:6)
                a <- reg.coeffs(coefficients = betas, groups.vector = groups.vector,
                    group = colnames(groups)[i])
                a <- c(a, rep(0, (7 - length(a))))
                curve(a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) +
                    a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5), from = min(time),
                    to = max(time), col = color1[i], add = TRUE,
                    lty = 2)
            }
            if (show.lines) {
                lx <- abcissa[group != 0]
                ly <- y[group != 0]
                ord <- order(lx)
                lxo <- lx[ord]
                lyo <- ly[ord]
                lines(lxo, lyo, col = color1[i])
            }
        }
    }
    op <- par(bg = "white")
    if (legend)
        legend(max(abcissa, na.rm = TRUE) * 1.02, ylim[1], legend = codeg,
            text.col = color1, col = color1, cex = cexlab, lty = 1,
            yjust = 0)
    par(op)
}

self.clust <- function(object.sig.genes, char.method, num.cluster) {
  dat.cluster <- object.sig.genes[["sig.profiles"]]
  num.cluster <- min(num.cluster, nrow(dat.cluster), na.rm = TRUE)
  if (char.method == "hclust") {
    mat.dcorrel <- matrix(rep(1, nrow(dat.cluster)^2), nrow(dat.cluster), nrow(dat.cluster)) -
      cor(t(dat.cluster), use = "pairwise.complete.obs")
    obj.hclust <- hclust(as.dist(mat.dcorrel), method = "ward.D")
    vec.cluster <- cutree(obj.hclust, k = num.cluster)
  } else if (char.method == "kmeans") {
    obj.kmeans <- kmeans(dat.cluster, centers = num.cluster, iter.max = 500)
    vec.cluster <- obj.kmeans$cluster
  } else if (char.method == "Mclust") {
    mc <- Mclust(dat.cluster, G = num.cluster)
    vec.cluster <- mc$classification
  }
  vec.cluster
}

self.plot.groups <- function(df.data, lst.design) {
  df.design <- lst.design$edesign
  vec.time <- df.design[, 1]
  vec.reps <- i.rank(df.design[, 2])
  codeg <- as.character(colnames(df.design)[3])
  df.groups <- df.design[codeg]
  vec.sample.mean <- apply(as.matrix(df.data), 2, median, na.rm = TRUE)
  vec.y <- vector(mode = "numeric", length = length(unique(vec.reps)))
  vec.x <- vector(mode = "numeric", length = length(unique(vec.reps)))
  for (i in unique(vec.reps)) {
    vec.x[i] <- mean(vec.time[vec.reps == i])
    vec.y[i] <- mean(vec.sample.mean[vec.reps == i])
  }
  xlim <- c(min(vec.x, na.rm = TRUE), max(vec.x, na.rm = TRUE) * 1.3)
  ylim <- c(min(as.numeric(vec.sample.mean), na.rm = TRUE), max(as.numeric(vec.sample.mean), na.rm = TRUE))
  lst.glm <- glm(vec.sample.mean ~ ., family = gaussian, data = as.data.frame(lst.design$dis))
  lst.result <- summary(lst.glm)
  vec.coeff <- as.vector(lst.result$coefficients[, 1])
  vec.coeff <- c(vec.coeff, rep(0, (7 - length(vec.coeff))))
  plot(x = vec.time, y = vec.sample.mean, xlim = xlim, ylim = ylim, cex = 0.8, pch = 21, xaxt = "n",
       main = NULL, sub = paste("Median profile of", nrow(df.data)),
       xlab = "Time", ylab = "Expression value")
  axis(1, at = unique(vec.x), labels = unique(vec.x), cex.axis = 1)
  lines(vec.x, vec.y)
  curve(vec.coeff[1] + vec.coeff[2]*x + vec.coeff[3]*(x^2) + vec.coeff[4]*(x^3) +
        vec.coeff[5]*(x^4) + vec.coeff[6]*(x^5) + vec.coeff[7]*(x^6),
        from = min(vec.time), to = max(vec.time), lty = 2, add = TRUE)
  legend(max(vec.time, na.rm = TRUE) * 1.02, ylim[1], legend = codeg, cex = 0.8, lty = 1, yjust = 0)
}

# check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$data) || is.null(opts$design)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}
if (is.null(opts$cluster) || is.null(opts$method)) {
  cat(getopt(command, usage = TRUE))
  q(status=-3)
}

# prepare data
print("start reading data table")
df.count <- read.delim(opts$data, header = TRUE, stringsAsFactors = FALSE)
mat.count <- as.matrix(df.count)

# prepare design
print("start reading design table")
df.design <- read.delim(opts$design, header = TRUE, stringsAsFactors = FALSE)
mat.design <- as.matrix(df.design)
lst.design <- make.design.matrix(mat.design)


# calculate minimum observations
print("start calculating requried arguments")
#degree <- sort(as.vector(table(df.design$Time)))[1]
degree <- length(unique(df.design$Time)) - 1
groups <- ncol(df.design) - 2
min.obs <- (degree + 1) * groups + 1
sample_num <- length(df.design$Time)
if (min.obs >=  sample_num) {
  min.obs <- sample_num
} else if (min.obs < sample_num) {
  min.obs <- min.obs
}

print(degree)
print(groups)
print(min.obs)

# p.vector performs a regression fit for each gene taking all variables present in the model given by
# a regression matrix and returns a list of FDR corrected significant genes.
print("start making regression fit for time series gene expression experiments")
NBp <- p.vector(mat.count, lst.design, counts = TRUE, min.obs = min.obs)

# reshape NBp result
df.pvt <- data.frame(seq.id = rownames(NBp$p.vector), p.value = NBp$p.vector, row.names = NULL)
df.pad <- data.frame(seq.id = rownames(NBp$p.vector), p.adjust = NBp$p.adjusted, row.names = NULL)
df.p <- dplyr::inner_join(df.pvt, df.pad, by = "seq.id")

# T.fit selects the best regression model for each gene using stepwise regression.
print("start making a stepwise regression fit for time series gene expression experiments")
NBt <- T.fit(NBp, min.obs = min.obs)

# get.siggenes creates lists of significant genes for a set of variables
# whose significance value has been computed with the T.fit function.
print("start extracting significant genes for sets of variables in time series gene expression experiments")
NBg <- get.siggenes(NBt, vars = "groups")
NBa <- get.siggenes(NBt, vars = "all")

# determine group situation and get R squared values of the sigificant gene
if (length(colnames(mat.design)) == 3) {
  comparison <- colnames(mat.design)[3]
  flag = 1
} else if (length(colnames(mat.design)) >= 4) {
  comparison <- paste(colnames(mat.design)[4], colnames(mat.design)[3], sep = "vs")
  flag = 0
}


if (length(colnames(mat.design)) <= 4){sig.rsq <- NBg[["sig.genes"]][[comparison]][["sig.pvalues"]]["R-squared"]
    df.r <- data.frame(seq.id = rownames(sig.rsq), r.squared = sig.rsq[["R-squared"]], row.names = NULL)
    result_dir <- paste(opts$output, "/", comparison, sep = "")
    dir.create(result_dir)
    # see.genes provides visualisation tools for gene expression values in a time course experiment (deprecated)
    # NBs <- see.genes(NBg[["sig.genes"]][[comparison]], k = opts$cluster, cluster.method = opts$method, min.obs = min.obs)
    print("start making cluster")
    if (nrow( NBg[["sig.genes"]][[comparison]][["sig.profiles"]])< 2){
      print("number of sig.genes greater than or equal to 2")
    }
    if (nrow( NBg[["sig.genes"]][[comparison]][["sig.profiles"]]) <  opts$cluster){
      print("number of sig.genes greater than or equal to num of cluster")
    }

    vec.cluster <- self.clust(NBg[["sig.genes"]][[comparison]], opts$method, opts$cluster)
    df.c <- data.frame(seq.id = names(vec.cluster), cluster = as.vector(vec.cluster))

    # plot profiles
    print("start exporting figures")
    df.sig.profiles <- data.frame(seq.id = rownames(NBg[["sig.genes"]][[comparison]][["sig.profiles"]]),
                                  NBg[["sig.genes"]][[comparison]][["sig.profiles"]])
    df.used <- dplyr::inner_join(df.c, df.sig.profiles, by = "seq.id")
    for (cluster in unique(vec.cluster)) {
      df.select <- df.used[df.used["cluster"] == cluster,]
      rownames(df.select) <- as.vector(df.select[["seq.id"]])
      df.plot <- df.select[c(-1, -2)]
      vec.ylim <- c(0, max(as.matrix(df.plot), na.rm = TRUE) * 1.1)
      pdf(paste(result_dir, paste("profile", cluster, "pdf", sep = "."), sep = "/"))
      PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
      dev.off()
      png(paste(result_dir, paste("profile", cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
      PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
      dev.off()
      svg(paste(result_dir, paste("profile", cluster, "svg", sep = "."), sep = "/"))
      PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
      dev.off()
      if (flag) {
        pdf(paste(result_dir, paste("groups", cluster, "pdf", sep = "."), sep = "/"))
        self.plot.groups(df.plot, lst.design)
        dev.off()
        png(paste(result_dir, paste("groups", cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
        self.plot.groups(df.plot, lst.design)
        dev.off()
        svg(paste(result_dir, paste("groups", cluster, "svg", sep = "."), sep = "/"))
        self.plot.groups(df.plot, lst.design)
        dev.off()
      } else {
        pdf(paste(result_dir, paste("groups", cluster, "pdf", sep = "."), sep = "/"))
        self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
        dev.off()
        png(paste(result_dir, paste("groups", cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
        self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
        dev.off()
        svg(paste(result_dir, paste("groups", cluster, "svg", sep = "."), sep = "/"))
        self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
        dev.off()
      }
    }

    # merge tables and export result
    print("start merging table")
    df.pr <- dplyr::inner_join(df.p, df.r, by = "seq.id")
    df.out <- dplyr::inner_join(df.pr, df.c, by = "seq.id")
    write.table(df.out, paste(result_dir, 'result.tsv', sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE)
} else if (length(colnames(mat.design)) > 4) {
    comparison <- colnames(mat.design)[3:length(colnames(mat.design))]
    control <-  colnames(mat.design)[3]
    tt <- opts$output
    for (i in comparison[2:length(comparison)]){
        comparison_u <- paste(i, control, sep = "vs")
        dir.create(paste(tt, "/", comparison_u, sep = ""))
        result_dir <- paste(tt, "/", comparison_u, sep = "")
        sig.rsq <- NBg[["sig.genes"]][[comparison_u]][["sig.pvalues"]]["R-squared"]
        df.r <- data.frame(seq.id = rownames(sig.rsq), r.squared = sig.rsq[["R-squared"]], row.names = NULL)

        # see.genes provides visualisation tools for gene expression values in a time course experiment (deprecated)
        # NBs <- see.genes(NBg[["sig.genes"]][[comparison_u]], k = opts$cluster, cluster.method = opts$method, min.obs = min.obs)
        print("start making cluster")

        vec.cluster <- self.clust(NBg[["sig.genes"]][[comparison_u]], opts$method, opts$cluster)

        #vec.cluster <- self.clust(NBg[["sig.genes"]][[comparison_u]], "hclust", 2)
        df.c <- data.frame(seq.id = names(vec.cluster), cluster = as.vector(vec.cluster))

        # plot profiles
        print("start exporting figures")
        df.sig.profiles <- data.frame(seq.id = rownames(NBg[["sig.genes"]][[comparison_u]][["sig.profiles"]]),
                                      NBg[["sig.genes"]][[comparison_u]][["sig.profiles"]])
        df.used <- dplyr::inner_join(df.c, df.sig.profiles, by = "seq.id")

        for (cluster in unique(vec.cluster)) {
          df.select <- df.used[df.used["cluster"] == cluster,]
          rownames(df.select) <- as.vector(df.select[["seq.id"]])
          df.plot <- df.select[c(-1, -2)]
          vec.ylim <- c(0, max(as.matrix(df.plot), na.rm = TRUE) * 1.1)
          pdf(paste(result_dir, paste("profile",comparison_u, cluster, "pdf", sep = "."), sep = "/"))
          PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
          dev.off()
          png(paste(result_dir, paste("profile", comparison_u,cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
          PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
          dev.off()
          svg(paste(result_dir, paste("profile", comparison_u,cluster, "svg", sep = "."), sep = "/"))
          PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
          dev.off()
          if (flag) {
            pdf(paste(result_dir, paste("groups",comparison_u, cluster, "pdf", sep = "."), sep = "/"))
            self.plot.groups(df.plot, lst.design)
            dev.off()
            png(paste(result_dir, paste("groups",comparison_u, cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
            self.plot.groups(df.plot, lst.design)
            dev.off()
            svg(paste(result_dir, paste("groups",comparison_u, cluster, "svg", sep = "."), sep = "/"))
            self.plot.groups(df.plot, lst.design)
            dev.off()
          } else {
            pdf(paste(result_dir, paste("groups", comparison_u,cluster, "pdf", sep = "."), sep = "/"))
            self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
            dev.off()
            png(paste(result_dir, paste("groups",comparison_u, cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
            self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
            dev.off()
            svg(paste(result_dir, paste("groups",comparison_u, cluster, "svg", sep = "."), sep = "/"))
            self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
            dev.off()
          }
        }
        # merge tables and export result
        print("start merging table")
        df.pr <- dplyr::inner_join(df.p, df.r, by = "seq.id")
        df.out <- dplyr::inner_join(df.pr, df.c, by = "seq.id")
        write.table(df.out, paste(result_dir,'result.tsv', sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE)
    }


    if (length(colnames(mat.design)) > 4) {
        sig.rsq <- NBa[["sig.genes"]][["sig.pvalues"]]["R-squared"]
        df.r <- data.frame(seq.id = rownames(sig.rsq), r.squared = sig.rsq[["R-squared"]], row.names = NULL)

        # see.genes provides visualisation tools for gene expression values in a time course experiment (deprecated)
        # NBs <- see.genes(NBg[["sig.genes"]][[comparison]], k = opts$cluster, cluster.method = opts$method, min.obs = min.obs)
        print("start making cluster")
        vec.cluster <- self.clust(NBa[["sig.genes"]], opts$method, opts$cluster)
        df.c <- data.frame(seq.id = names(vec.cluster), cluster = as.vector(vec.cluster))

        # plot profiles
        print("start exporting figures")
        df.sig.profiles <- data.frame(seq.id = rownames(NBa[["sig.genes"]][["sig.profiles"]]),
                                      NBa[["sig.genes"]][["sig.profiles"]])
        df.used <- dplyr::inner_join(df.c, df.sig.profiles, by = "seq.id")
        result_dir <- paste(opts$output, "/", "all", sep = "")
        dir.create(paste(opts$output, "/", "all", sep = ""))
        for (cluster in unique(vec.cluster)) {
          df.select <- df.used[df.used["cluster"] == cluster,]
          rownames(df.select) <- as.vector(df.select[["seq.id"]])
          df.plot <- df.select[c(-1, -2)]
          vec.ylim <- c(0, max(as.matrix(df.plot), na.rm = TRUE) * 1.1)
          pdf(paste(result_dir, paste("profile", "all",cluster, "pdf", sep = "."), sep = "/"))
          PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
          dev.off()
          png(paste(result_dir, paste("profile", "all",cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
          PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
          dev.off()
          svg(paste(result_dir, paste("profile","all", cluster, "svg", sep = "."), sep = "/"))
          PlotProfiles(df.plot, cond = rownames(df.design), main = cluster, ylim = vec.ylim, repvect = df.design[["Replicate"]])
          dev.off()
          if (flag) {
            pdf(paste(result_dir, paste("groups","all", cluster, "pdf", sep = "."), sep = "/"))
            self.plot.groups(df.plot, lst.design)
            dev.off()
            png(paste(result_dir, paste("groups","all", cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
            self.plot.groups(df.plot, lst.design)
            dev.off()
            svg(paste(result_dir, paste("groups","all", cluster, "svg", sep = "."), sep = "/"))
            self.plot.groups(df.plot, lst.design)
            dev.off()
          } else {
            pdf(paste(result_dir, paste("groups", "all",cluster, "pdf", sep = "."), sep = "/"))
            self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
            dev.off()
            png(paste(result_dir, paste("groups", "all",cluster, "png", sep = "."), sep = "/"), width = 4200, height = 4200, res = 600)
            self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
            dev.off()
            svg(paste(result_dir, paste("groups","all", cluster, "svg", sep = "."), sep = "/"))
            self.PlotGroups(df.plot, edesign = df.design, show.fit = TRUE, dis = lst.design$dis, groups.vector = lst.design$groups.vector, min.obs = min.obs)
            dev.off()
          }
        }

        # merge tables and export result
        print("start merging table")
        df.pr <- dplyr::inner_join(df.p, df.r, by = "seq.id")
        df.out <- dplyr::inner_join(df.pr, df.c, by = "seq.id")
        write.table(df.out, paste(result_dir, 'result.tsv', sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE)
    }
}