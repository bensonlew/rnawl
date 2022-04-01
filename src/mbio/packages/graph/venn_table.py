# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import os
import subprocess


def venn_table(otu_table, group_table, R_path):
    """
    输入OTU表和group表，生成venn表

    :param otu_table: OTU表的全路径
    :param group_table: group表的全路径
    :R_path: 软件R的安装路径
    :retuen: venn表格的全路径
    """
    workpath = os.path.dirname(group_table)
    r_file = os.path.join(workpath, 'cmd.r')
    with open(r_file, 'w') as f:
        # 生成输入表格
        my_dict = {"otu_table": otu_table, "group_table": group_table}
        r_str = """
        # pdf="venn.%(otu_table)s.tiff"
        lst <-list()
        table <-read.table(file="%(otu_table)s",sep="\\t",head=T,check.names = FALSE)
        rownames(table) <- as.character(table[,1])
        table <-table[,-1]
        # pdf <-"venn.%(otu_table)s.%(group_table)s.tiff"
        group <- read.table("%(group_table)s")
        glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] %%in%% unique(group[,2])[x]),1])
        names(glst) <-unique(group[,2])
        tab <-sapply(1:length(glst),function(x) apply(table[as.character(as.vector(glst[[x]]))],1,sum))
        table <-tab[apply(tab,1,function(x)any(x>0)),]
        colnames(table) <-unique(group[,2])
        for(i in 1:length(colnames(table))){
            samp <-colnames(table)[i]
            lst[[samp]] <- rownames(table)[which(table[,i] !=0)]
        }
        # tiff(file=pdf,width=600,height=800,pointsize=16)
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
            for(i in sn:2){
                    sl <- combn (names(lst),i,simplify=FALSE)
                    inter <-lapply(sl,function(x) intersect_n(lst[x]))
                    names(inter) <-lapply(sl,function(x) paste(x,collapse =\" & \"))
                    sets <- c(sets,inter)
                    sls <- c(sls,sl)
            }
            # get all unique
            for(i in 1:sn){
                    uniq <- list(setdiff(unlist(lst[names(lst)[i]]),unlist(lst[names(lst)[-i]])))
                    names(uniq) <- paste(names(lst)[i],\" only\")
                    sets <- c(sets,uniq)
            }
            return(sets)
        }
        sets <- getset(lst)
        lst_len <- length(sets)
        venn_matrix <- matrix(, lst_len, 3)
        for(i in 1:lst_len){
            venn_matrix[i,1] <- names(sets[i])
            venn_matrix[i,2] <- length(sets[[i]])
            #venn_matrix[i,3] <- paste(sets[[i]],collapse =",")
            venn_matrix[i,3] <- paste(sets[[i]],collapse =";")
        }

        # write sets to file
        otuset <- '"'workpath'"'+"venn_table.%(otu_table)s.%(group_table)s.xls"
        write.table(sets_table,otuset,sep = "\\t",eol="\\n",row.names=FALSE,col.names=FALSE,quote=FALSE)

        # add new_table by khl
        new_table <- "new_venn_table.%(otu_table)s.xls"
        new_matrix = matrix(,length(lst),2)
        for(i in 1:length(lst)){
                new_matrix[i,1] = names(lst)[i]
                #new_matrix[i,2] = paste(lst[[i]],collapse=",")
                new_matrix[i,2] = paste(lst[[i]],collapse=";")
        }
        colnames(new_matrix)= c("#group_name","species_name")
        write.table(new_matrix,new_table,sep = "\t",eol = "\n",row.names=FALSE,col.names=TRUE,quote=FALSE)
        
        """ % my_dict
        f.write(r_str + "\n")

    try:
        my_cmd = R_path + " --restore --no-save < " + r_file
        subprocess.check_call(my_cmd, shell=True)
        venn_table = os.path.join(workpath, "venn_table." + otu_table + '.' + group_table + '.xls')
        return venn_table
    except subprocess.CalledProcessError:
        raise Exception(u"生成Venn表出错！")


def venn_graph(otu_table, group_table, output):
    group = {}
    sets = {}
    with open(otu_table, "r") as o, open(group_table, "r") as g, open(output, "w") as w:
        samples = o.readline().strip().split("\t")[1:]
        g.readline()
        for line in g:
            line = line.strip().split("\t")
            if line[1] in group:
                group[line[1]].append(line[0])
            else:
                group[line[1]] = [line[0]]
        # print group
        for gp in group:
            sets[gp] = set()
        for line in o:
            line = line.strip().split("\t")
            for n, s in enumerate(samples):
                for gp in group:
                    if float(line[n+1]) == 0:
                        continue
                    elif s in group[gp]:
                        sets[gp].add(line[0].split("; ")[-1])
        w.write("#group_name\tspecies_name\n")
        for s in sets:
            w.write("{}\t{}\n".format(s, ",".join(list(sets[s]))))

def venn_graph_mg(otu_table, group_table, output):
    group = {}
    sets = {}
    with open(otu_table, "r") as o, open(group_table, "r") as g, open(output, "w") as w:
        samples = o.readline().strip().split("\t")[1:]
        g.readline()
        for line in g:
            line = line.strip().split("\t")
            if line[1] in group:
                group[line[1]].append(line[0])
            else:
                group[line[1]] = [line[0]]
        # print group
        for gp in group:
            sets[gp] = set()
        for line in o:
            line = line.strip().split("\t")
            for n, s in enumerate(samples):
                for gp in group:
                    if float(line[n+1]) == 0:
                        continue
                    elif s in group[gp]:
                        sets[gp].add(line[0].split("; ")[-1])
        w.write("#group_name\tspecies_name\n")
        for s in sets:
            #w.write("{}\t{}\n".format(s, ",".join(list(sets[s]))))
            w.write("{}\t{}\n".format(s, ";".join(list(sets[s]))))  # modify by zhujuan 20171107 一些分类名中含有","会导致结果错误
