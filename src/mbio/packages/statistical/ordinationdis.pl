#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my ( $type, $community, $dist, $environment, $outdir, $group, $color, $shape, $point_lab, $priorSpe_pick, $pick_number, $rda_threshold, $width, $height, $label_font, $label_cex, $pca_env, $env_labs, $envfit,$scale);
GetOptions(
    "type=s"          => \$type,

    "community|com=s"     => \$community,
    "dist=s"          => \$dist,
    "environment|env=s"   => \$environment,
    "outdir|o=s"        => \$outdir,

    "group|g=s"         => \$group,
    "color|col=s"         => \$color,
    "shape|pch=s"         => \$shape,

    "point_lab|pl=s"     => \$point_lab,
    "width|pw=f"           => \$width,
    "height|ph=f"           => \$height,
    "label_font|lf=s"    => \$label_font,
    "label_cex|size=f"     => \$label_cex,

    "priorSpe_pick=s" => \$priorSpe_pick,
    "pick_number=i"   => \$pick_number,
    "rda_threshold=f" => \$rda_threshold,
    "pca_env=s" => \$pca_env,
    "env_labs=s" => \$env_labs,
    "envfit=s" => \$envfit,
    "scale=s" => \$scale,
);


my $usage = <<"USAGE";
    Program : $0
    Usage : perl $0 -type pca-pcoa-nmds-rda [options]
        [REQUIRED options]
        -type           [pca-pcoa-nmds-rda] which type of the ordination analysis you want to do
        [Optional options]
        -outdir         [ordination] the output file dir

        -community      the input community table
        -dist           the input dist matrix
        -pca_env        [T/F] add environment data to pca
        -environment    the env data

        -group          the group design file
        -color          the column title to defined the sample color in group file
        -shape          the column title to defined the sample shape in group file

        -point_lab      [T] whether to plot point lables
        -label_font     [1/2/3/4] the font of the point label in the graph
        -label_cex      [0.5] the cex of the point lable in the graph
        -width          [8] the width of the graph
        -height         [6] the height of the graph

        -priorSpe_pick  [T/F] in rda analysis, whether to pick the main otus of the data
        -pick_number    [10] the tip n otus of this table in graph
        -rda_threshold  [3.5] the threshold to decide rda or cca
        -env_labs       choose the columns of env_labs
        -envfit         [T/F] whether to get envfit_table, in rda/cca
        -scale         [T/F] whether to scale

USAGE
die $usage if( !defined $type );

if (!defined $community) {
    $community = "";
}

if (!defined $dist) {
    $dist = "";
}

if (!defined $pca_env) {
    $pca_env = "";
}else{
    die "You must have a environment file !\n" if ($environment eq "");
}

if (!defined $environment) {
    $environment = "";
}

if (!defined $outdir) {
    $outdir = "ordination";
}
`mkdir -p $outdir`;

if (!defined $group) {
    $group = "";
    $color = "";
    $shape = "";
}

if (!defined $point_lab) {
    $point_lab = "T";
}
if (!defined $width) {
    $width = 8;
}
if (!defined $height) {
    $height = 6;
}
if (!defined $label_font) {
    $label_font = 1;
}
if (!defined $label_cex) {
    $label_cex = 0.5;
}

if (!defined $priorSpe_pick) {
    $priorSpe_pick = "T";
}
if (!defined $pick_number) {
    $pick_number = 10;
}
if (!defined $rda_threshold) {
    $rda_threshold = 3.5;
}
if (!defined $env_labs) {
    $env_labs = "nothing";
}
if (!defined $envfit) {
    $envfit = "F"
}
if (!defined $scale) {
    $scale = "F";
}


####

my $run_pca = "F";
my $run_pcoa = "F";
my $run_nmds = "F";
my $run_rda = "F";

my @types = split(/-/, $type);
foreach my $type_tmp (@types) {
    if ($type_tmp =~ /pca/i) {
        $run_pca = "T";
        `mkdir -p $outdir/pca`;
    } elsif ($type_tmp =~ /pcoa/i) {
        $run_pcoa = "T";
        `mkdir -p $outdir/pcoa`;
    } elsif ($type_tmp =~ /nmds/i) {
        $run_nmds = "T";
        `mkdir -p $outdir/nmds`;
    } elsif ($type_tmp =~ /rda/i) {
        $run_rda = "T";
        `mkdir -p $outdir/rda`;
    }
}

if ($run_pca eq "T") {
    die "You must have a community file !\n" if ($community eq "");
}
if ($run_rda eq "T") {
    die "You must have a community file !\n" if ($community  eq "");
    die "You must have an environment file !\n" if ($environment  eq "");
}
if (($run_pcoa eq "T") || ($run_nmds eq "T")) {
    die "You must have a dist matrix file !\n" if ($dist  eq "");
}

if ($group ne "") {
    #die "You must set which column to specific the color or shape, at least one of them !" if ((!defined $color) && (!defined $shape ));

    if (!defined $color) {
        $color = "";
    }

    if (!defined $shape) {
        $shape = "";
    }
}


open RCMD, ">cmd.r";

print RCMD "

####    ²ÎÊýÊäÈë Óë ´¦Àí
##
r_pca <- \"$run_pca\"
r_pcoa <- \"$run_pcoa\"
r_nmds <- \"$run_nmds\"
r_rda <- \"$run_rda\"
envfit <- \"$envfit\"


# Èç¹û´æÔÚ g_design£¬ Ôò g_col Óë g_pch ±ØÐëÖÁÉÙ´æÔÚÒ»¸ö
g_design    <- \"$group\"
g_col <- \"$color\"
g_pch <- \"$shape\"

##
i_community <- \"$community\"
i_dist      <- \"$dist\"
i_env       <- \"$environment\"
o_dir <- \"$outdir\"
pca_env <- \"$pca_env\"
# env_labs <- \"$env_labs\"

##
point_lab <- \"$point_lab\"

priorSpe_pick <- \"$priorSpe_pick\"
pick_number <- $pick_number
rda_threshold <- $rda_threshold

p_w <- $width
p_h <- $height
label_font <- $label_font
label_cex <- $label_cex


####    read group_design file £ºSample_ID   g1    g2    g3  ...
if ( g_design != \"\" ) {
    # 12 ÖÖÑÕÉ«
    mycol <- c(\"#E41A1C\", \"#377EB8\", \"#4DAF4A\", \"#984EA3\", \"#FF7F00\", \"#E7298A\", \"#66C2A5\", \"#FC8D62\", \"#8DA0CB\", \"#FFFF33\", \"#A65628\", \"#F781BF\")
    # 12 ÖÖÐÎ×´
    mypch <- c(15:18, 7:14)

    group_design <- read.table(g_design, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"!\",colClasses = c(\"character\"))
    colnames(group_design)[1] <- \"Sample_ID\"
    group_design\$Sample_ID <- as.character(group_design\$Sample_ID)
    rownames(group_design) <- group_design\$Sample_ID
    sample_sort <- sort(rownames(group_design))
    group_design <- group_design[sample_sort, ]

    gd_titles <- colnames(group_design)
    data_group <- data.frame(Sample_ID = group_design\$Sample_ID)

    if ( g_col != \"\" ) {
        data_group\$g_col <- as.character(group_design[, which(gd_titles %in% g_col)])
        data_group\$color <- factor(as.factor(data_group\$g_col), labels = mycol[1:length(levels(as.factor(data_group\$g_col)))])

        data_group\$group <- data_group\$g_col
    }

    if ( g_pch  != \"\" ) {
        data_group\$g_pch <- group_design[, which(gd_titles %in% g_pch)]
        data_group\$shape <- factor(as.factor(data_group\$g_pch), labels = mypch[1:length(levels(as.factor(data_group\$g_pch)))])

        if ( g_col != \"\" ) {
            data_group\$group <- paste(data_group\$g_col, data_group\$g_pch, sep = \"_\")
        } else {
            data_group\$group <- data_group\$g_pch
        }
    }

    rownames(data_group) <- as.character(data_group\$Sample_ID)
    # Sample_ID    g_col   color   group    g_pch   shape

}

####    ¶ÁÊý¾ÝÎÄ¼þ, ²ÎÊý´¦Àí
if ( i_community != \"\" ) { # otu_id   sam1    sam2    sam3    ...
    i_splitc <- unlist(strsplit(x = i_community, split = \"/\", fixed = FALSE, perl = TRUE))
    i_community_clean <- i_splitc[length(i_splitc)]
    i_split1 <- unlist(strsplit(x = i_community_clean, split = \"[-_.]\", fixed = FALSE, perl = TRUE))
    all_levels <- c(\"phylum\", \"class\", \"order\", \"family\", \"genus\", \"species\", \"otu\", \"Pathway\", \"pathway\", \"KO\", \"ko\", \"Module\", \"module\", \"Enzyme\", \"enzyme\", \"EC\", \"ec\", \"COG\", \"cog\", \"CAZy\", \"cazy\", \"ARDB\", \"ardb\")

    i_lab_cpos <- which(i_split1 %in% all_levels)
    if(length(i_lab_cpos) == 0){
        i_lab_comm <- paste(i_split1, collapse = \"_\")
    } else {
        if(length(i_lab_cpos) > 1){
            i_lab_comm <- i_split1[i_lab_cpos[length(i_lab_cpos)]]
        }else{
            i_lab_comm <- i_split1[i_lab_cpos]
        }
    }


    data_community <- read.table(i_community, row.names = 1, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"\")
    data_community_temp <- read.table(i_community, row.names = 1, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"\", colClasses = c(\"character\"))
    rownames(data_community) <- row.names(data_community_temp)
    data_community <- t(data_community)

    if ( g_design != \"\" ) {
        inter_samples <- sort(intersect(rownames(data_community), rownames(data_group)))
        data_community <- data_community[inter_samples, ]
        group_community <- data_group[inter_samples, ]
        group_real <- group_design[inter_samples, ]
    }

    filename_community <- i_lab_comm

}

if ( i_dist != \"\" ) { # sample_id sam1    sam2    sam3    ...

    i_splitd <- unlist(strsplit(x = i_dist, split = \"/\", fixed = FALSE, perl = TRUE))
    i_dist_clean <- i_splitd[length(i_splitd)]
    i_split2 <- unlist(strsplit(x = i_dist_clean, split = \"[-_.]\", fixed = FALSE, perl = TRUE))
    all_levels <- c(\"phylum\", \"class\", \"order\", \"family\", \"genus\", \"species\", \"otu\", \"Pathway\", \"pathway\", \"KO\", \"ko\", \"Module\", \"module\", \"Enzyme\", \"enzyme\", \"EC\", \"ec\", \"COG\", \"cog\", \"CAZy\", \"cazy\", \"ARDB\", \"ardb\")

    i_lab_cpod <- which(i_split2 %in% all_levels)
    if (length(i_lab_cpod) == 0){
        i_lab_dist <- paste(i_split2, collapse = \"_\")
    }else{
        if(i_lab_cpod > 1) {
            i_lab_dist <- i_split2[i_lab_cpod[length(i_lab_cpod)]]
        } else {
            i_lab_dist <- i_split2[i_lab_cpod]
        }
    }    

    data_dist <- read.table(i_dist,row.names = 1, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"\")
    data_dist_temp <- read.table(i_dist,row.names = 1, sep = \"\\t\", header = TRUE, check.names = FALSE, comment.char = \"\", colClasses = c(\"character\"))
    rownames(data_dist) <- row.names(data_dist_temp)
    # rownames(data_dist) <- as.character(data_dist[, 1])
    # data_dist <- data_dist[, -1]

    if ( g_design != \"\" ) {
        inter_samples <- sort(intersect(rownames(data_dist), rownames(data_group)))
        data_dist <- data_dist[inter_samples, inter_samples]
        group_dist <- data_group[inter_samples, ]
        group_real <- group_design[inter_samples, ]
    }

    filename_dist <- i_lab_dist

}

if ( i_env != \"\" ) { # sample_id  env1    env2    env3    ...
    data_env <- read.table(i_env, sep = \"\\t\", header = TRUE, check.names = FALSE, row.names=1,comment.char = \"\")
    data_env_temp <- read.table(i_env, sep = \"\\t\", header = TRUE, check.names = FALSE, row.names=1,comment.char = \"\",colClasses = c(\"character\"))
    rownames(data_env) <- row.names(data_env_temp)
    # rownames(data_env) <- as.character(data_env[, 1])
    # data_env <- data_env[, -1]

    inter_samples <- sort(intersect(rownames(data_env), rownames(data_community)))
    data_env_c <- data_env[inter_samples, ]
    data_env_c <- data.frame(data_env_c)
    rownames(data_env_c) <- inter_samples
    colnames(data_env_c) <- colnames(data_env)
    data_env <- data_env_c
    data_community <- data_community[inter_samples, ]
    if ( g_design != \"\" ) {
        group_community <- data_group[inter_samples, ]
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

if (r_pca == \"T\") {
    data_pca <- data_community

    filename_pca <- paste(o_dir, \"/pca/\", filename_community, \"_pca\", sep = \"\")
    if (\"$scale\" == \"T\"){
        scale = TRUE
    }else{
        scale = FALSE
    }

    pca_analysis <- prcomp(data_pca, scale = scale)
    pca_sites <- pca_analysis\$x
    pca_summary <- summary(pca_analysis)
    pca_importance <- pca_summary\$importance[2, ]
    pca_rotation <- pca_summary\$rotation
    pca_rotation_all <- pca_rotation
    if(nrow(pca_rotation) > 100) {
        pca_rotation <- pca_rotation[1:100, ]
    }

    if ( g_design != \"\" ) {
        group_pca <- group_community
        fac <- as.factor(group_real[ ,2])
	fac<-factor(fac,levels = unique(group_real[, 2]))
        ellipse12 <- add_ellipse(pca_sites[, c(1,2)], fac, pc=\"12\")
        ellipse13 <- add_ellipse(pca_sites[, c(1,3)], fac, pc=\"13\")
        ellipse23 <- add_ellipse(pca_sites[, c(2,3)], fac, pc=\"23\")
        ellipse_all <- rbind(colnames(ellipse12), ellipse12, ellipse13, ellipse23)
        write.table(ellipse_all, paste(o_dir, \"/pca/\", \"ellipse.xls\", sep = \"\"), sep = \"\\t\", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    ####    write table
    pca_sitesname <- paste(filename_pca, \"sites.xls\", sep = \"_\")
    pca_imporname <- paste(filename_pca, \"importance.xls\", sep = \"_\")
    pca_rotatname <- paste(filename_pca, \"rotation.xls\", sep = \"_\")
    pca_rotatname_all <- paste(filename_pca, \"rotation_all.xls\", sep = \"_\")

    if (pca_env == \"T\") {
        ef <- envfit(pca_analysis, data_env, permu = 999)
        ef_vector <- scores(ef, \"vectors\")
        ef_factor <- scores(ef, \"factors\")
        if (is.matrix(ef_vector)){
            ef_vector_scores <- paste(filename_pca, \"envfit_vector_scores.xls\", sep = \"_\")
            ef_vector_pr <- paste(filename_pca, \"envfit_vector.xls\", sep = \"_\")
            write.table(cbind(rownames(ef\$vectors\$arrows),ef\$vectors\$arrows,ef\$vectors\$r,ef\$vectors\$pvals),ef_vector_pr, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", \"PC1\", \"PC2\", \"r2\", \"Pr(>r)\"))
            write.table(cbind(rownames(ef_vector),ef_vector),ef_vector_scores, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", \"PC1\", \"PC2\"))
        }
        if (is.matrix(ef_factor)){
            ef_factor_scores <- paste(filename_pca, \"envfit_factor_scores.xls\", sep = \"_\")
            ef_factor_pr <- paste(filename_pca, \"envfit_factor.xls\", sep = \"_\")
            factor_r <- data.frame(ef\$facotrs\$r)
            write.table(matrix(c(ef\$factors\$r,ef\$factors\$pvals),ncol=2,dimnames = list(names(ef\$factors\$pvals),c(\"Sample_ID\\tr2\",\"Pr(>r)\"))),ef_factor_pr, sep = \"\\t\", quote = FALSE, row.names = TRUE, col.names=TRUE)
            write.table(cbind(rownames(ef_factor),ef_factor),ef_factor_scores, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", \"PC1\", \"PC2\"))
        }
    }
    write.table(cbind(rownames(pca_sites), pca_sites), pca_sitesname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", colnames(pca_sites)))
    write.table(cbind(names(pca_importance), pca_importance), pca_imporname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", \"Proportion of Variance\"))
    write.table(cbind(rownames(pca_rotation), pca_rotation), pca_rotatname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", colnames(pca_rotation)))
    write.table(cbind(rownames(pca_rotation_all), pca_rotation_all), pca_rotatname_all,sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", colnames(pca_rotation_all)))  #add by zouxuan
}


if (r_pcoa == \"T\") {
    data_pcoapre <- as.dist(data_dist, diag = TRUE)
    filename_pcoa <- paste(o_dir, \"/pcoa/\", filename_dist, \"_pcoa\", sep = \"\")
    if (\"$scale\" == \"T\"){
        data_pcoa <- scale(data_pcoapre)
    }else{
        data_pcoa <- data_pcoapre
    }
    #data_pcoa <- vegdist(data_pcoap, method = d_method)
    pcoa_analysis <- cmdscale(data_pcoa, k=length(data_dist[1,]) - 1, eig=T)

    pcoa_point <- pcoa_analysis\$points
    eigenvaluespre <- abs(pcoa_analysis\$eig)\/sum(abs(pcoa_analysis\$eig))
    eigenvalues <- pcoa_analysis\$eig
    pcoa_sitesname <- paste(filename_pcoa, \"sites.xls\", sep = \"_\")
    pcoa_eigname <- paste(filename_pcoa, \"eigenvalues.xls\", sep = \"_\")#zhangpeng-xiugai-baifenlv-2017/1/3
    pcoa_eignamepre <- paste(filename_pcoa, \"eigenvaluespre.xls\", sep = \"_\")
    write.table(pcoa_point, pcoa_sitesname, sep=\"\\t\", col.names = NA, quote = FALSE)
    write.table(as.data.frame(eigenvalues), pcoa_eigname, sep=\"\\t\", col.names = NA, quote = FALSE)
    write.table(as.data.frame(eigenvaluespre), pcoa_eignamepre, sep=\"\\t\", col.names = NA, quote = FALSE)

    if ( g_design != \"\" ) {
        ### 添加分组椭圆
        fac <- as.factor(group_real[,2])
	fac<-factor(fac,levels = unique(group_real[, 2]))
        ellipse12 <- add_ellipse(pcoa_point[, c(1,2)], fac, pc=\"12\")
        ellipse13 <- add_ellipse(pcoa_point[, c(1,3)], fac, pc=\"13\")
        ellipse23 <- add_ellipse(pcoa_point[, c(2,3)], fac, pc=\"23\")
        ellipse_all <- rbind(colnames(ellipse12), ellipse12, ellipse13, ellipse23)
        write.table(ellipse_all, paste(o_dir, \"/pcoa/\", \"ellipse.xls\", sep = \"\"), sep = \"\\t\", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}


if (r_nmds == \"T\") {
    data_nmdspre <- as.dist(data_dist, diag = TRUE)
    filename_nmds <- paste(o_dir, \"/nmds/\", filename_dist, \"_nmds\", sep = \"\")

    if(\"$scale\" == \"T\"){
        data_nmds <- scale(data_nmdspre)
    }else{
        data_nmds <- data_nmdspre
    }
    nmds_analysis <- metaMDS(data_nmds, k = 2)
    nmds_sites <- nmds_analysis\$points
    nmds_stress <- nmds_analysis\$stress

    if ( g_design != \"\" ) {
        group_nmds <- group_dist
        fac <- as.factor(group_real[,2])
	fac<-factor(fac,levels = unique(group_real[, 2]))
        ellipse12 <- add_ellipse(nmds_sites[, c(1,2)], fac, pc=\"12\")
        ellipse12 <- rbind(colnames(ellipse12), ellipse12)
        write.table(ellipse12, paste(o_dir, \"/nmds/\", \"ellipse.xls\", sep = \"\"), sep = \"\\t\", row.names = FALSE, col.names = FALSE, quote = FALSE)
    } 

    nmds_sitesname <- paste(filename_nmds, \"sites.xls\", sep = \"_\")
    nmds_stressfile <- paste(filename_nmds, \"stress.xls\", sep = \"_\")
    write.table(cbind(rownames(nmds_sites), nmds_sites), nmds_sitesname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", colnames(nmds_sites)))
    write.table( nmds_stress, nmds_stressfile, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = TRUE )
}


if (r_rda == \"T\") {
    data_rda <- data_community
    data_rda <- data_rda[, (apply(data_rda, 2 , function(y) sum(y) != 0 ))]
    env_rda <- data_env
    if ( g_design != \"\") {
        group_rda <- group_community
    }
    filename_rda <- paste(o_dir, \"/rda/\", filename_community, sep = \"\")
    dca_normal_run <- function(){dca_analysis <<- decorana(data_rda)}
    dca_error_run <- function(){
        scale_data <<- apply(data_rda,2,function(x) as.numeric(x)/sum(as.numeric(x)))
        dca_analysis <<- decorana(scale_data)
    }
    tryCatch(dca_normal_run(), error=function(e){print(\"error info\");print(e);dca_error_run()})
    # dca_analysis <- decorana(data_rda)
    dca_scores <- scores(dca_analysis)
    axis_lengths <- max(dca_scores[, 1]) - min(dca_scores[, 1])

    dca_name <- paste(filename_rda, \"dca.xls\",sep = \"_\")
    scores_dca <- scores(dca_analysis)
    Axis_lengths <- c()
    for (i in 1:length(scores_dca[1,])) Axis_lengths=c(Axis_lengths, max(scores_dca[,i]) - min(scores_dca[,i]))
    Eigenvalues <- dca_analysis\$evals
    Decorana_values <- dca_analysis\$evals.decorana
    dca_table <- rbind(Eigenvalues, Decorana_values,Axis_lengths)
    write.table(dca_table, dca_name, sep = '\\t', col.names = NA, quote = FALSE)
    # sink(file = dca_name, append = FALSE)
    # print(dca_analysis)
    # print(\"axis_lengths:\")
    # print(axis_lengths)
    # sink(file = NULL)

    if(axis_lengths < rda_threshold) {
        i_method <- \"rda\"
        cons_analysis <- rda(data_rda~$env_labs, data=data_env)
    }

    if(axis_lengths >= rda_threshold) {
        i_method <- \"cca\"
        cons_analysis <- cca(data_rda~$env_labs, data=data_env)
    }
    filename_rda <- paste(filename_rda, i_method, sep = \"_\")

    ####    rda/cca analysis and scores
    cons_summary <- summary(cons_analysis)
    if (is.matrix(cons_summary\$biplot)) {
        cons_biplot <- cons_summary\$biplot
        cons_biplotname <- paste(filename_rda, \"biplot.xls\",sep = \"_\")
        write.table(cbind(rownames(cons_biplot), cons_biplot), cons_biplotname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Environment\", colnames(cons_biplot)))
    }
    if (is.matrix(cons_summary\$centroids)) {
        cons_centroids <- cons_summary\$centroids
        cons_centroidsname <- paste(filename_rda, \"centroids.xls\",sep = \"_\")
        write.table(cbind(rownames(cons_centroids), cons_centroids), cons_centroidsname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Environment\", colnames(cons_centroids)))
    }
    if (is.matrix(cons_summary\$biplot)){
        cons_environment <- cons_summary\$biplot
    }else{
        cons_environment <- cons_summary\$centroids
    }
    cons_species <- cons_summary\$species
    cons_sites <- cons_summary\$sites
    cons_importance <- cons_summary\$cont\$importance[2, ]

    cons_envname <- paste(filename_rda, \"environment.xls\",sep = \"_\")
    cons_spename <- paste(filename_rda, \"species.xls\",sep = \"_\")
    cons_sitesname <- paste(filename_rda, \"sites.xls\",sep = \"_\")
    cons_importname <- paste(filename_rda, \"importance.xls\",sep = \"_\")

    ####envfit.txt##add by zhouxuan 20170401
    if (envfit == \"T\"){
        cons_envfit <- paste(filename_rda, \"envfit.xls\",sep = \"_\")
        ef<-envfit(cons_analysis,data_env,permu=999)
        ef<-data.frame(ef\$vectors\$arrows,ef\$vectors\$r,ef\$vectors\$pvals)
        write.table(ef,file = cons_envfit)
    }
    #cons_envfit <- paste(filename_rda, \"envfit.xls\",sep = \"_\")
    #ef<-envfit(cons_analysis,data_env,permu=999)
    #ef<-data.frame(ef\$vectors\$arrows,ef\$vectors\$r,ef\$vectors\$pvals)
    #write.table(ef,file = cons_envfit)

    # write.table(cbind(rownames(cons_environment), cons_environment), cons_envname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Environment\", colnames(cons_environment)))
    write.table(cbind(rownames(cons_species), cons_species), cons_spename, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Species\", colnames(cons_species)))
    write.table(cbind(rownames(cons_sites), cons_sites), cons_sitesname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Sample_ID\", colnames(cons_sites)))
    write.table(cbind(names(cons_importance), cons_importance), cons_importname, sep = \"\\t\", quote = FALSE, row.names = FALSE, col.names = c(\"Axes\", \"Proportion of Variance\"))
    }

    ";

    close RCMD;
    # `Rscript cmd.r`;
    # `rm cmd.r`;
