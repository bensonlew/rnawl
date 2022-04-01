add_ellipse <- function(pls.result,classFc,confidence,parCompVi=c(1,2)){
    x <- pls.result
    if(x@summaryDF[, "ort"] > 0) {
        if(parCompVi[2] > x@summaryDF[, "ort"] + 1){
            stop("Selected orthogonal component for plotting (ordinate) exceeds the total number of orthogonal components of the model", call. = FALSE)
        }
        tCompMN <- cbind(x@scoreMN[, 1], x@orthoScoreMN[, parCompVi[2] - 1])
        pCompMN <- cbind(x@loadingMN[, 1], x@orthoLoadingMN[, parCompVi[2] - 1])
        colnames(pCompMN) <- colnames(tCompMN) <- c("h1", paste("o", parCompVi[2] - 1, sep = ""))
    }else{
        if(max(parCompVi) > x@summaryDF[, "pre"]){
            stop("Selected component for plotting as ordinate exceeds the total number of predictive components of the model", call. = FALSE)
        }
        tCompMN <- x@scoreMN[, parCompVi, drop = FALSE]
        pCompMN <- x@loadingMN[, parCompVi, drop = FALSE]
    }
    cxtCompMN <- cor(x@suppLs[["xModelMN"]], tCompMN,use = "pairwise.complete.obs")

    if(!is.null(x@suppLs[["yModelMN"]])){
        cytCompMN <- cor(x@suppLs[["yModelMN"]], tCompMN, use = "pairwise.complete.obs")
    }

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]])) {
        pexVi <- integer(x@suppLs[["topLoadI"]] * ncol(pCompMN) * 2) ## 'ex'treme values
        for(k in 1:ncol(pCompMN)) {
            pkVn <-  pCompMN[, k]
            pexVi[1:(2 * x@suppLs[["topLoadI"]]) + 2 * x@suppLs[["topLoadI"]] * (k - 1)] <- c(order(pkVn)[1:x@suppLs[["topLoadI"]]],
                                                                         rev(order(pkVn, decreasing = TRUE)[1:x@suppLs[["topLoadI"]]]))
        }
    }else{ pexVi <- 1:ncol(x@suppLs[["xModelMN"]]) }

    pxtCompMN <- cbind(pCompMN, cxtCompMN)

    if(ncol(pCompMN) == 1) {
       colnames(pxtCompMN)[2] <- paste0("cor_", colnames(pxtCompMN)[2])
    } else{
        colnames(pxtCompMN)[3:4] <- paste0("cor_", colnames(pxtCompMN)[3:4])
    }

    topLoadMN <- pxtCompMN
    topLoadMN <- topLoadMN[pexVi, , drop = FALSE]

    if(x@suppLs[["topLoadI"]] * 4 < ncol(x@suppLs[["xModelMN"]]) && ncol(pCompMN) > 1) {
        topLoadMN[(2 * x@suppLs[["topLoadI"]] + 1):(4 * x@suppLs[["topLoadI"]]), c(1, 3)] <- NA
        topLoadMN[1:(2 * x@suppLs[["topLoadI"]]), c(2, 4)] <- NA
    }

    parAsColFcVn = classFc

    if(!any(is.na(parAsColFcVn))) {
        obsColVc <- .plotColorF(as.vector(parAsColFcVn))[["colVc"]]
    } else if(!is.null(x@suppLs[["yMCN"]]) && ncol(x@suppLs[["yMCN"]]) == 1) { ## (O)PLS of single response
        obsColVc <- .plotColorF(c(x@suppLs[["yMCN"]]))[["colVc"]]
    } else { ## PCA
        obsColVc <- rep("black", nrow(tCompMN))
    }

    ploMN <- tCompMN
    colC <- unique(obsColVc)
    ploColVc <- obsColVc
    all_ellipse <- c("group","m1","m2","c11","c12","c21","c22")

    for(each in colC){
        if(table(ploColVc)[each]<2){next}  ##zouguanqing add
        each_ellipse <- .each_ellipse(ploMN,ploColVc,each)
        group <- names(ploColVc[ploColVc==each])[1]
        each_ellipse <- c(group,each_ellipse)
        all_ellipse <- rbind(all_ellipse,each_ellipse)
    }
    all_sample_ellipse <- .each_ellipse(ploMN,ploColVc,"all",all="T")
    all_sample_ellipse <- c("All",all_sample_ellipse)
    all_ellipse <- rbind(all_ellipse,all_sample_ellipse)
    return (all_ellipse)
}

.each_ellipse <- function(ploMN,ploColVc,each_colC,all="F"){
    if(all=="T"){
        xMN <- ploMN
    }else{
        xMN <- ploMN[ploColVc == each_colC, , drop = FALSE]
    }

    ## Hotteling's T2 (Tenenhaus98, p86)
    ##----------------------------------
    radVn <- seq(0, 2 * pi, length.out = 100)

    if(ncol(xMN) != 2){
        stop("Matrix must have two columns", call. = FALSE) }
    csqN <- qchisq(confidence, 2) ## ncol(xMN) == 2

    xMeaVn <- colMeans(xMN)
    xCovMN <- cov(xMN)
    xCovSvdLs <- svd(xCovMN, nv = 0)

    mahMN <- matrix(1, nrow = length(radVn), ncol = 1) %*% xMeaVn +
             cbind(cos(radVn), sin(radVn)) %*% diag(sqrt(xCovSvdLs[["d"]] * csqN)) %*% t(xCovSvdLs[["u"]])

    m1 =  as.matrix(xMeaVn[1])
    m2 = as.matrix(xMeaVn[2])
    c11 = (sqrt(xCovSvdLs[["d"]] * csqN)[1]) * (t(xCovSvdLs[["u"]])[1,1])
    c12 = (sqrt(xCovSvdLs[["d"]] * csqN)[2]) * (t(xCovSvdLs[["u"]])[1,2])
    c21 = (sqrt(xCovSvdLs[["d"]] * csqN)[1]) * (t(xCovSvdLs[["u"]])[2,1])
    c22 = (sqrt(xCovSvdLs[["d"]] * csqN)[2]) * (t(xCovSvdLs[["u"]])[2,2])
    return (c(m1,m2,c11,c12,c21,c22))
}

.plotColorF <- function(namVcn) {
    ## 16 color palette without 'gray'
    #palVc <- c("blue", "red", "green3", "cyan", "magenta", "#FF7F00", "#6A3D9A", "#B15928", "aquamarine4", "yellow4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#FFFF99")
    # miodify in 20181126, change 16 groups limit
    palVc <- paste("mycolor", 1:500, sep = "")
    if(is.null(namVcn) || all(is.na(namVcn))) {
        if(!is.null(namVcn)) {
            palNamVc <- paste0(1:length(palVc),"_",palVc)
            #print(matrix(palNamVc, ncol = 1))
        }
        return(palVc)
}
else {
        if(is.character(namVcn)) {
            namFcn <- factor(namVcn)
            if(length(levels(namFcn)) <= length(palVc)) {
                scaVc <- palVc[1:length(levels(namFcn))]
            } else{
                scaVc <- c(palVc,rep("gray",length(levels(namFcn)) - length(palVc)))
            }
            names(scaVc) <- levels(namFcn)
            colVc <- scaVc[unlist(sapply(namVcn,
                                         function(scaleC) {
                                             if(is.na(scaleC)){
                                             return(NA)}else{which(levels(namFcn) == scaleC)}
                                         }))]
        }else if(is.numeric(namVcn)) {
            scaVc <- rev(rainbow(100, end = 4/6))
            if(length(namVcn) > 1) {
                colVc <- scaVc[round((namVcn - min(namVcn, na.rm = TRUE)) / diff(range(namVcn, na.rm = TRUE)) * 99) + 1]
            }else{
                colVc <- rep("black", length(namVcn))
            }
        }else{
            stop("'namVcn' argument must be a vector of either character or numeric mode", call. = FALSE)
        }
        colVc[is.na(colVc)] <- "grey"
        names(colVc) <- namVcn
    }
    return(list(colVc = colVc, scaVc = scaVc))
}





