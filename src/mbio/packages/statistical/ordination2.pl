#!/usr/bin/perl
# 计算分组椭圆
# 只计算前三个维度

use strict;
use warnings;
use Getopt::Long;

my $VERSION = "1.0.1";
my $DATE = "2019-09-23";
my $AUTHOR = "shuang.hou\@majorbio.com";

my ( $type, $community, $sites, $group, $outdir);
GetOptions(
    "type=s"          => \$type,
    "community|com=s"     => \$community,
    "sites=s"     => \$sites,
    "group|g=s"         => \$group,
    "outdir|o=s"        => \$outdir
);


my $usage = <<"USAGE";
    Program : $0
    Discription :
    Version : $VERSION
    Contact : $AUTHOR
    Lastest : 2019-09-23
    Usage : perl $0 -type rda-plsda-dbrda [options]
        [REQUIRED options]
        -type           [rda-plsda-dbrda] which type of the ordination analysis you want to do
        -community      the input community table
        -sites      the input sites table
        -group          the group design file
        -outdir         [ordination] the output file dir

USAGE
#die $usage if( !defined $type );
#die $usage if( !defined $community );
#die $usage if( !defined $sites );
#die $usage if( !defined $group );
#
#die "You must have a otu file !\n" if ($community eq "");
#die "You must have a sites file !\n" if ($sites eq "");
#die "You must have a group file !\n" if ($group eq "");


####

my $run_dbrda = "F";
my $run_plsda = "F";
my $run_rda = "F";

my @types = split(/-/, $type);
foreach my $type_tmp (@types) {
    if ($type_tmp =~ /dbrda/i) {
        $run_dbrda = "T";
        `mkdir -p $outdir/dbrda`;
    } elsif ($type_tmp =~ /plsda/i) {
        $run_plsda = "T";
        `mkdir -p $outdir/plsda`;
    } elsif ($type_tmp =~ /rda/i) {
        $run_rda = "T";
        `mkdir -p $outdir/rda`;
    }
}

open RCMD, ">cmd2.r";

print RCMD "

##
r_dbrda <- \"$run_dbrda\"
r_plsda <- \"$run_plsda\"
r_rda <- \"$run_rda\"

##
g_design    <- \"$group\"
i_sites <- \"$sites\"
o_dir <- \"$outdir\"
i_community <- \"$community\"
type <- \"$type\"

####    read group_design file
if ( g_design != \"\" ) {
    group_design <- read.table(g_design, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"!\",colClasses = c(\"character\"))
    colnames(group_design)[1] <- \"Sample_ID\"
    group_design\$Sample_ID <- as.character(group_design\$Sample_ID)
    rownames(group_design) <- group_design\$Sample_ID
    sample_sort <- sort(rownames(group_design))
    group_design <- group_design[sample_sort, ]
    data_group <- data.frame(Sample_ID = group_design\$Sample_ID)
    rownames(data_group) <- as.character(data_group\$Sample_ID)
}

####
if ( i_community != \"\" ) {   # otu_id   sam1    sam2    sam3    ...
    data_community <- read.table(i_community, row.names = 1, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"\")
    data_community_temp <- read.table(i_community, row.names = 1, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"\", colClasses = c(\"character\"))
    rownames(data_community) <- row.names(data_community_temp)
    data_community <- t(data_community)

    if ( g_design != \"\" ) {
        inter_samples <- sort(intersect(rownames(data_community), rownames(data_group)))
        data_community <- data_community[inter_samples, ]
        group_real <- group_design[inter_samples, ]
    }
}

library(\"vegan\")
library(\"maptools\")
library(\"ade4\")

fac2disj<- function(fac, drop = FALSE) {
    ## Returns the disjunctive table corrseponding to a factor
    n <- length(fac)
    fac <- as.factor(fac)
    if(drop)
        fac <- factor(fac)
    x <- matrix(0, n, nlevels(fac))
    x[(1:n) + n * (unclass(fac) - 1)] <- 1
    dimnames(x) <- list(names(fac), as.character(unique(fac)))
    return(data.frame(x, check.names = FALSE))
}

add_ellipse <- function(dfxy, fac, wt = rep(1, length(fac)), pc, xax = 1, yax = 2, col = rep(1, length(levels(fac))), cellipse = 1.5){
    dfxy <- data.frame(dfxy)
    dfdistri <- fac2disj(fac) * wt
    coul <- col
    w1 <- unlist(lapply(dfdistri, sum))
    dfdistri <- t(t(dfdistri)/w1)
    if (nrow(dfxy) != nrow(dfdistri))
        stop(paste(\"Non equal row numbers\", nrow(dfxy), nrow(dfdistri)))
    results <- as.data.frame(matrix(0, 1, ncol(dfdistri)+1))
    colnames(results) <- c(\"PC\", colnames(dfdistri))
    results[, 1] <- pc
    for (i in 1:ncol(dfdistri)) {
        z<- dfdistri[, i]
	if (length(which( !z == 0.0))<3){
            results[, i+1]  <- 0
            next
        }
        z <- z/sum(z)
        x <- dfxy[, xax]
        y <- dfxy[, yax]
        m1 <- sum(x * z)
        m2 <- sum(y * z)
        v1 <- sum((x - m1) * (x - m1) * z)
        v2 <- sum((y - m2) * (y - m2) * z)
        cxy <- sum((x - m1) * (y - m2) * z)

        lig <- 100
        epsi <- 1e-10
        x <- 0
        y <- 0
        if (v1 < 0)
            v1 <- 0
        if (v2 < 0)
            v2 <- 0
        if (v1 == 0 && v2 == 0){
            results[, i+1]  <- 0
            next
	}
        delta <- (v1 - v2) * (v1 - v2) + 4 * cxy * cxy
        delta <- sqrt(delta)
        l1 <- (v1 + v2 + delta)/2
        l2 <- v1 + v2 - l1
        if (l1 < 0)
            l1 <- 0
        if (l2 < 0)
            l2 <- 0
        l1 <- sqrt(l1)
        l2 <- sqrt(l2)
        test <- 0
        if (v1 == 0) {
            a0 <- 0
            b0 <- 1
            test <- 1
        }
        if ((v2 == 0) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (((abs(cxy)) < epsi) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (test == 0) {
            a0 <- 1
            b0 <- (l1 * l1 - v1)/cxy
            norm <- sqrt(a0 * a0 + b0 * b0)
            a0 <- a0/norm
            b0 <- b0/norm
        }
        a1 <- 2 * pi/lig
        c11 <- cellipse * a0 * l1
        c12 <- (-cellipse) * b0 * l2
        c21 <- cellipse * b0 * l1
        c22 <- cellipse * a0 * l2
        tmp1 <- paste(m1, c11, c12, m2, c21, c22, sep=\",\")
        results[, i+1]  <- tmp1
	if ((c12==0) && (c22==0))
		results[, i+1]  <- 0
    }
    return(results)
}


if (r_plsda == \"T\" || r_dbrda == \"T\" || r_rda == \"T\") {
    if ( g_design != \"\" ) {
        pca_sites <- read.table(i_sites, sep=\"\\t\", header=TRUE, row.names=1)
        fac <- as.factor(group_real[ ,2])
	    fac <- factor(fac,levels = unique(group_real[, 2]))
	    if (r_plsda == \"T\") {
	        colnames(pca_sites) <- gsub(\"\\\\.\", \"\", colnames(pca_sites))
	    }
	    if (length(colnames(pca_sites)) > 2) {
            ellipse12 <- add_ellipse(pca_sites[, c(1,2)], fac, pc=paste(colnames(pca_sites)[1], colnames(pca_sites)[2], sep=\"\"))
            ellipse13 <- add_ellipse(pca_sites[, c(1,3)], fac, pc=paste(colnames(pca_sites)[1], colnames(pca_sites)[3], sep=\"\"))
            ellipse23 <- add_ellipse(pca_sites[, c(2,3)], fac, pc=paste(colnames(pca_sites)[2], colnames(pca_sites)[3], sep=\"\"))
            ellipse_all <- rbind(colnames(ellipse12), ellipse12, ellipse13, ellipse23)
	    } else {
            ellipse12 <- add_ellipse(pca_sites[, c(1,2)], fac, pc=paste(colnames(pca_sites)[1], colnames(pca_sites)[2], sep=\"\"))
            ellipse_all <- rbind(colnames(ellipse12), ellipse12)
	    }
	}
        write.table(ellipse_all, paste(o_dir, \"/\", type, \"/\", \"ellipse.xls\", sep = \"\"), sep = \"\\t\", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
    ";

    close RCMD;
