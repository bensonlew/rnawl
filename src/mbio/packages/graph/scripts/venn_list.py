# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
import os
import argparse


def venn_list(list_file, r_file,skip_num=-1):
    """
    输入list表，生成cmd.r
    list文件格式：无表头，两列\t分隔；第一列为集合名，第二列为元素名（元素之间逗号分隔）
    :param list_file: list文件的全路径

    """
    workpath = os.path.dirname(r_file)
    #r_file = os.path.join(workpath, 'cmd.r')
    if skip_num >= 6:
        skip_str = '''
            if(sn>=%s){
                end=sn
            }else{
                end=2
            }
        '''%(skip_num)
    else:
        skip_str = '''
            end = 2
        '''

    with open(r_file, 'w') as f:
        # 生成输入表格
        my_dict = {"list_file": list_file}
        r_str = """
        lst <-list()
        listf <- read.table(file = "%(list_file)s",sep = "\\t", head=F, check.names = FALSE, comment.char="")
        rownames(listf) <- as.character(listf[,1])

        for(i in 1:length(rownames(listf))){
            samp <- rownames(listf)[i]
            elements <- strsplit(as.character(listf[i,2]), ",")
            lst[[samp]] <- as.vector(unlist(elements))
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
            maxl <-0 """
        f.write(r_str + "\n")
        f.write(skip_str)


        r_str = """
            # get all intersect
            for(i in sn:end){
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
        outset <- '"""+ workpath +"""/venn_table.xls'
        write.table(venn_matrix,outset,sep = "\\t",eol="\\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
        
        # add by khl 20170601
        new_table <- '"""+ workpath +"""/new_venn_graph.xls'
        new_matrix = matrix(,length(lst),2)
        for(i in 1:length(lst)){
            new_matrix[i,1] = names(lst)[i]
            new_matrix[i,2] = paste(lst[[i]],collapse=",")
        }
        colnames(new_matrix)= c("#set_name","attributes_name")
        write.table(new_matrix,new_table,sep = "\\t",eol = "\\n",row.names=FALSE,col.names=TRUE,quote=FALSE)
        
        """
        f.write(r_str + "\n")

def main():
    _i = ''
    _o = ''
    parse = argparse.ArgumentParser()
    parse.add_argument('-i','--inputfile_list',help='input list_file')
    parse.add_argument('-o','--outputfile',help='outputfile is cmd.r')
    parse.add_argument('-n','--skip_num',default=-1,help='if not -1 , then More than skip_num groups, only common and unique are counted')
    args = parse.parse_args()
    _i = args.inputfile_list
    _o = args.outputfile
    venn_list(_i,_o, args.skip_num)

if __name__ == '__main__':
    main()
