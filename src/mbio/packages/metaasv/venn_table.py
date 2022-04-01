# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import argparse

"""
功能：根据输入的otu表和group表生成cmd.r，花瓣图
"""

def venn_table(otu_table, group_table,r_file):
    """
    输入OTU表和group表，生成cmd.r
    :param otu_table: OTU表的全路径
    :param group_table: group表的全路径
    """
    workpath = os.path.dirname(group_table)
    with open(r_file, 'w') as f:
        # 生成输入表格
        my_dict = {"otu_table": otu_table, "group_table": group_table}
        r_str = """
        lst <-list()
        table <-read.table(file="%(otu_table)s",sep="\\t",head=T,check.names = FALSE, comment.char = "")
        rownames(table) <- as.character(table[,1])
        table <-table[,-1]
        group <- read.table("%(group_table)s", colClasses=c('character'))
        glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] %%in%% unique(group[,2])[x]),1])
        names(glst) <-unique(group[,2])
        tab <-sapply(1:length(glst),function(x) apply(table[as.character(as.vector(glst[[x]]))],1,sum))
        table <-tab[apply(tab,1,function(x)any(x>0)),]
        colnames(table) <-unique(group[,2])
        for(i in 1:length(colnames(table))){
            samp <-colnames(table)[i]
            lst[[samp]] <- rownames(table)[which(table[,i] !=0)]
        }
        """ % my_dict
        f.write(r_str + "\n")
        # 生成输出表格
        r_str = """
        intersect_n <-function(x) {
            n <-length(x)
            x2 <- x[[1]]
            for(i in 2:n){
                    x2 <- intersect(x[[i]],x2)
            }
            return(x2)
        }

        getset <- function(lst){   # lst is a list
            sn <-length(names(lst))
            sets <- list()
            sls <- list()
            maxl <-0
            # get all intersect
            inter <- lapply(sn,function(x) intersect_n(lst))
            new_name <- names(lst)
            names(inter) <- paste(new_name,collapse =\" & \")
            sets <- c(sets,inter)
            # sls <- c(sls,sl)
            # get all unique
            for(i in 1:sn){
                    uniq <- list(setdiff(unlist(lst[names(lst)[i]]),unlist(lst[names(lst)[-i]])))
                    names(uniq) <- paste(names(lst)[i],\" only\")
                    sets <- c(sets,uniq)
            }
            return(sets)
        }
        print(names(lst))
        sets <- getset(lst)
        lst_len <- length(sets)
        venn_matrix <- matrix(, lst_len, 3)
        for(i in 1:lst_len){
            venn_matrix[i,1] <- names(sets[i])
            venn_matrix[i,2] <- length(sets[[i]])
            venn_matrix[i,3] <- paste(sets[[i]],collapse =",")
        }
        # write sets to file
        otuset <- '"""+ workpath +"""/venn_table.xls'
        write.table(venn_matrix,otuset,sep = "\\t",eol="\\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
        
        # add by khl 20170601        
        new_table <- '"""+ workpath +"""/new_venn_graph.xls'
        new_matrix = matrix(,length(lst),2)
        for(i in 1:length(lst)){
            new_matrix[i,1] = names(lst)[i]
            new_matrix[i,2] = paste(lst[[i]],collapse=",")
        }
        colnames(new_matrix)= c("#group_name","species_name")
        write.table(new_matrix,new_table,sep = "\\t",eol = "\\n",row.names=FALSE,col.names=TRUE,quote=FALSE)
        
        """ % my_dict
        f.write(r_str + "\n")

def main():
    _i = ''
    _g = ''
    _o = ''
    parse = argparse.ArgumentParser()
    parse.add_argument('-i','--inputfile_otu',help='input otu_table')
    parse.add_argument('-g','--inputfile_group',help='input group_table')
    parse.add_argument('-o','--outputfile',help='outputfile is cmd.r')
    args = parse.parse_args()
    _i = args.inputfile_otu
    _g = args.inputfile_group
    _o = args.outputfile
    venn_table(_i,_g,_o)

if __name__ == '__main__':
    main()
