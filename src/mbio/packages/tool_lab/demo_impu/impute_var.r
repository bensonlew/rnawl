# Title     : function (ImputeVar) from MetaboAnalystR package
# Objective : DEMO
# Created by: jincheng.qin
# Created on: 2019/7/23

function (mSetObj = NA, method = "min")
{
  mSetObj <- .get.mSet(mSetObj)
  int.mat <- mSetObj$dataSet$preproc
  new.mat <- NULL
  msg <- mSetObj$msgSet$replace.msg
  if (method == "exclude") {
    good.inx <- apply(is.na(int.mat), 2, sum) == 0
    new.mat <- int.mat[, good.inx]
    msg <- c(msg, "Variables with missing values were excluded.")
  }
  else if (method == "min") {
    minConc <- min(int.mat[int.mat > 0], na.rm = T)/2
    int.mat[int.mat == 0 | is.na(int.mat)] <- minConc
    new.mat <- int.mat
    msg <- c(msg, "Variables with missing values were replaced with a small value.")
  }
  else if (method == "colmin") {
    new.mat <- apply(int.mat, 2, function(x) {
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- min(x, na.rm = T)/2
      }
      x
    })
    msg <- c(msg, "Missing variables were replaced with the half of minimum values for each feature column.")
  }
  else if (method == "mean") {
    new.mat <- apply(int.mat, 2, function(x) {
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- mean(x, na.rm = T)
      }
      x
    })
    msg <- c(msg, "Missing variables were replaced with the mean value for each feature column.")
  }
  else if (method == "median") {
    new.mat <- apply(int.mat, 2, function(x) {
      if (sum(is.na(x)) > 0) {
        x[is.na(x)] <- median(x, na.rm = T)
      }
      x
    })
    msg <- c(msg, "Missing variables were replaced with the median for each feature column.")
  }
  else {
    if (method == "knn") {
      new.mat <- t(impute::impute.knn(t(int.mat))$data)
    }
    else {
      if (method == "bpca") {
        new.mat <- pcaMethods::pca(int.mat, nPcs = 5,
          method = "bpca", center = T)@completeObs
      }
      else if (method == "ppca") {
        new.mat <- pcaMethods::pca(int.mat, nPcs = 5,
          method = "ppca", center = T)@completeObs
      }
      else if (method == "svdImpute") {
        new.mat <- pcaMethods::pca(int.mat, nPcs = 5,
          method = "svdImpute", center = T)@completeObs
      }
    }
    msg <- c(msg, paste("Missing variables were imputated using",
      toupper(method)))
  }
  mSetObj$dataSet$proc <- as.data.frame(new.mat)
  mSetObj$msgSet$replace.msg <- msg
  return(.set.mSet(mSetObj))
}

