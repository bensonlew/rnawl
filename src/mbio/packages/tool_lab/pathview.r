library(methods)
library(pathview)
library(gage)
library(clusterProfiler)
library(getopt)
library(stringr)

command <- matrix(c(
  'type', 't', 1, 'character', 'data type. One of continuous, genelist, discrete',
  'gene', 'g', 2, 'character', 'gene file',
  'cpd', 'c', 2, 'character', 'compound file',
  'gene_id', 'i', 2, 'character', 'gene id type, entrezid as default',
  'cpd_id', 'd', 2, 'character', 'compound id type, kegg as default',
  'species', 's', 2, 'character', 'species. if missing, Homo Sapiens will be used',
  'pathway_id', 'p', 2, 'character', 'kegg pathway ids, separated by comma',
  'ko_dir', 'k', 2, 'character', 'path for species-specific ko file',
  'kegg_dir', 'e', 2, 'character', 'path for kegg png repository',
  'if_ko', 'y', 0, 'logical', 'state whether input pathways starts with ko',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$gene) && is.null(opts$cpd)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}



# params for pathview
in_species <- ifelse(is.null(opts$species), 'hsa', opts$species)
pathview_species<-ifelse(is.null(opts$if_ko), in_species, 'ko')
# pathview_species<- in_species
gene_id<-ifelse(is.null(opts$gene_id), 'entrezid', opts$gene_id)
gene_id<-ifelse(gene_id=='entrez', 'entrezid', gene_id)
cpd_id<-ifelse(is.null(opts$cmp_id), 'kegg', opts$cmp_id)
both_dir<-list(gene=F, cpd=F)
col_bin<-list(gene=10, cpd=10)
if_discrete<-list(gene=ifelse(opts$type=='continuous',F,T), cpd=F)
colour_key<-ifelse(opts$type=='genelist', F, T)
kegg_dir<-ifelse(is.null(opts$kegg_dir), '.', opts$kegg_dir)



# read data
if(!is.null(opts$gene)){
  if(opts$type=='genelist'){
    gene_raw<-read.delim(opts$gene)
    colnames(gene_raw)[1]<-'name'
    gene_raw$data<-1
    gene_df<-subset(gene_raw, select=c(data))
    col_bin$gene<-1
  }else if(opts$type=='discrete'){
    gene_raw<-read.delim(opts$gene, row.names = 1)
    colnames(gene_raw)[1]<-'raw_data'
    classes<-unique(gene_df[,1])
    gene_raw$data<-' '
    for (i in c(1:length(classes))){
      gene_raw[gene_raw$raw_data == classes[i],]$data<-i
    }
    gene_df<-subset(gene_raw, select=c(data))
    col_bin$gene<-length(classes)
  }else{
    gene_df<-read.delim(opts$gene, row.names = 1)
  }
}

if(!is.null(opts$cpd)){
  cpd_df<-read.delim(opts$cpd, row.names = 1)
}


# Convert genename/entrezid to K_id
convert_for_ko<-function(id_type, df){
    ko_path <- paste0(opts$ko_dir, '_allGene_ko.xls')
    ko_df<-read.delim(ko_path)
    if (id_type == 'gene_entry'){
       ko_df$name<-sapply(ko_df[,1], function(x) strsplit(as.character(x), ':')[[1]][2])
    }else{
       ko_df$name<-ko$gene_name
    }
    df$name<-rownames(df)
    new_df<-merge(df, ko_df[, c('name', 'KO_id')], by='name')
    new_df<-aggregate(new_df[ ,2:ncol(new_df)-1],by=list(ID=new_df$KO_id),FUN=mean)
    rownames(new_df)<-new_df$ID
    new_df<-new_df[,3:ncol(new_df)]
    return(new_df)
}

# Obtain KO list for GAGE
ko_list<-function(file_path, df, df_type){
    kolist<-list()
    ko_path<-ifelse(df_type =='gene',paste0(file_path,'_allGene_ko.xls'), paste0(file_path,'_compound.xls'))
    ko_df<-read.delim(ko_path)
    if(df_type == 'gene'){
      ko_df$name<-' '
      for (i in c(1:nrow(ko_df))){
        pathway<-strsplit(as.character(ko_df[i,4]), ';')[[1]]
        ko_df$name[i]<-gene_id<-strsplit(as.character(ko_df[i,1]), ':')[[1]][2]
        for (j in pathway){
          kolist[[j]]<-c(kolist[[j]], gene_id)
        }
      }
    }else if(df_type=='cpd'){
      for (i in c(1:nrow(ko_df))){
        pathway<-as.character(ko_df[i,1])
        cpd_id<-strsplit(as.character(ko_df[i,2]), '  ')[[1]]
        for (j in cpd_id){
          id<-strsplit(j, ':')[[1]]
          if (length(id)==1 && which(cpd_id==j) != 1){
            next
          }else{
            kolist[[pathway]]<-c(kolist[[pathway]], ifelse(length(id)>1, id[2:length(id)], id[1]))
          }
        }
      }
    }
    return(kolist)
}

# Enrichment Analysis
run_enrich<-function(df, df_type){
  if (opts$type == 'continuous'){
    if (df_type == 'gene' && tolower(gene_id) != 'entrezid'){
      s<-bitr(rownames(df), fromType = toupper(gene_id), toType ="ENTREZID" ,OrgDb = "org.Hs.eg.db")
      col_num<-dim(df)[2]
      df$name<-rownames(df)
      lncRNA_target<-merge(df,s,by.x="name",by.y=toupper(gene_id))
      rownames(df)<-df$ENTREZID
      df<-df[,1:col_num]
    }
    gs_data<-ko_list(opts$ko_dir, df, df_type)
    raw<-gage(df, gsets = gs_data)
    result<-sigGeneSet(raw, cutoff = 0.05, qpval = 'p.val')
    result<-as.data.frame(rbind(result$greater, result$less))
    result<-result[order(abs(result$p.val)),]
  }else{
    if (tolower(gene_id) != 'entrezid'){
        s<-bitr(rownames(df), fromType = toupper(gene_id), toType ="ENTREZID" ,OrgDb = "org.Hs.eg.db")
    }
    kegg <- enrichKEGG(gene = ifelse(gene_id !='entrezid', s$ENTREZID, rownames(df)),
                       organism = in_species,
                       pvalueCutoff = 0.05,
                       pAdjustMethod="none",
                       qvalueCutoff = 1,
                       minGSSize = 2)
    result<-kegg@result
  }
  write.table(result, paste0(df_type, '_GAGE_result.txt'), sep = '\t', quote=F)
  return(result)
}



# Obtain pathway ids
if(!is.null(opts$pathway_id)){
  pathways<-strsplit(opts$pathway_id, ',')
  if (pathview_species == 'ko' && !is.null(opts$gene)){
    if (tolower(gene_id) %in% c('entrez', 'entrezid')){
        gene_df<-convert_for_ko('gene_entry', gene_df)
        gene_id<-'kegg'
    }else if (tolower(gene_id) %in% c('genename', 'symbol')){
        gene_df<-convert_for_ko('gene_name', gene_df)
        gene_id<-'kegg'
    }
  }
}else if (is.null(opts$pathway_id) && opts$type=='continuous'){
  gene_enrich<-cpd_enrich<-NULL
  if (!is.null(opts$gene)){
    gene_enrich<-run_enrich(gene_df, 'gene')
  }
  if (!is.null(opts$cpd)){
    cpd_enrich<-run_enrich(cpd_df, 'cpd')
  }
  if (!is.null(cpd_enrich) && !is.null(gene_enrich)){
    inter_path<-intersect(rownames(gene_enrich), rownames(cpd_enrich))
    pathways<-inter_path[1:ifelse(length(inter_path)>5, 5, length(inter_path))]
  }else if (!is.null(gene_enrich)){
    pathways<-rownames(gene_enrich)[1:ifelse(length(rownames(gene_enrich))>5, 5, length(rownames(gene_enrich)))]
  }else{
    pathways<-rownames(cpd_enrich)[1:ifelse(length(rownames(cpd_enrich))>5, 5, length(rownames(cpd_enrich)))]
  }
}else{
  stop('请输入关注pathwayID')
}



# run pathview
for (each in pathways){
  each_id<-strsplit(each, ' ')[[1]][1]
  each_id<-str_extract_all(each_id,"[0-9]+")[[1]]
  print(each_id)
  if(!is.null(opts$gene) && is.null(opts$cpd)){
    both_dir$gene<-ifelse(min(gene_df)<0, T, F)
    pv.out <- pathview(gene.data = gene_df,
                       pathway.id = each_id,
                       species = pathview_species,
                       gene.idtype = gene_id,
                       both.dirs = both_dir,
                       discrete=if_discrete,
                       bins = col_bin,
                       plot.col.key= colour_key,
                       kegg.dir = kegg_dir,
                       same.layer=T)
  }else if(!is.null(opts$cpd) && is.null(opts$gene)){
    both_dir$cpd<-ifelse(min(cpd_df)<0, T, F)
    pv.out <- pathview(cpd.data = cpd_df,
                       pathway.id = each_id,
                       species = pathview_species,
                       cpd.idtype = cpd_id,
                       both.dirs = both_dir,
                       discrete=if_discrete,
                       bins = col_bin,
                       plot.col.key= colour_key,
                       kegg.dir = kegg_dir,
                       same.layer=T)
  }else{
    both_dir$gene<-ifelse(min(gene_df)<0, T, F)
    both_dir$cpd<-ifelse(min(cpd_df)<0, T, F)
    pv.out <- pathview(gene.data = gene_df,
                       cpd.data = cpd_df,
                       pathway.id = each_id,
                       species = pathview_species,
                       gene.idtype = gene_id,
                       both.dirs = both_dir,
                       discrete=if_discrete,
                       bins = col_bin,
                       cpd.idtype = cpd_id,
                       plot.col.key= colour_key,
                       kegg.dir = kegg_dir,
                       same.layer=T)
  }
}