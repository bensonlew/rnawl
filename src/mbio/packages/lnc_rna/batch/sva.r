# load package

library(getopt)

# set options
command <- matrix(c(
    "json", "i", 1, "character", "setting file in JSON format",
    "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (! is.null(opts$help)) {
    cat(getopt(command, usage = TRUE))
    q(status = 2)
}
if (is.null(opts$json)) {
    cat(getopt(command, usage = TRUE))
    q(status = - 1)
}

library(rjson)
library(sva)
svaBatchCor2 <- function(dat, mmi, mm0,n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  if (is.null(n.sv))   {n.sv <- num.sv(dat,mmi,method="leek")}
o <- sva(dat,mmi,mm0,n.sv=n.sv)
W <- cbind(mmi,o$sv)
alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
o$corrected <- t(Y - W %*% alpha)
return(o)
}

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
output.path <- json.lst$output
counts <- read.delim(count.path, row.names = 1)
counts = as.matrix(counts)

keep_lst=list(which(apply(counts, 1, function(x){!all(x==0)})))
keep <- Reduce(intersect, keep_lst)
rm <- setdiff(1:nrow(counts), keep)
countsOri <- counts
counts <- counts[keep, ]

print("start creating group.factor")
group = as.factor(group.vec)

print("start batch effects analysis")
mod1 = model.matrix(~group)
print('step1')
mod0 = cbind(rep(1,length(group)))
print('step2')
n.sv <- num.sv(counts,mod1,method="leek")
print(n.sv)
print('step3')
count <- 0
repeat {
  if (count > n.sv) break
  count <- count + 1
  x.inv <- try(svaBatchCor2(counts,mod1,mod0,n.sv=count))
  if ('try-error' %in% class(x.inv)){
  next
  } else {
  lnj.corr <- svaBatchCor2(counts,mod1,mod0,n.sv=count)
  break
  }
}

#lnj.corr <- svaBatchCor2(counts,mod1,mod0,n.sv=1)
print('step4')
co <- lnj.corr$corrected
co[co<0]=0
adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
dimnames(adjust_counts_whole) <- dimnames(countsOri)
adjust_counts_whole[keep, ] <- co
adjust_counts_whole[rm, ] <- countsOri[rm, ]
col_name = colnames(counts)
write.table(t(c('seq_id',col_name)), file = file.path(output.path, 'count_batch.txt'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(adjust_counts_whole, file = file.path(output.path, 'count_batch.txt'), sep = '\t', quote = FALSE, row.names = TRUE, col.name = FALSE, append = TRUE)


