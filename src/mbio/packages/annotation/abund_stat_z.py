#!/usr/bin/python
#guanqing.zou 20181105
import pandas as pd
import sys

class Table2stat:
    def __init__(self,num_file,anno_file,nlink, alink,modify,anno_rm):
        self.reads_file =  pd.read_table(num_file,sep='\t').drop_duplicates(nlink)
        self.rlink = nlink
        self.anno_file = pd.read_table(anno_file, sep='\t').drop_duplicates(alink)
        self.alink = alink
        self.anno_file_rm_item = anno_rm
        if modify == True:
            self.modify_alink()
        if self.anno_file_rm_item !=[]:
            self.rm_item(self.anno_file,self.anno_file_rm_item)
        self.merge2table()

    def modify(self,k):
        self.anno_file[k] = self.anno_file[k].str[:-2]

    def modify_alink(self):
        self.modify(self.alink)

    def merge2table(self):
        self.mdata = pd.merge(self.reads_file, self.anno_file, left_on=self.rlink, right_on=self.alink, how='inner')

    def level_sum(self,level):
        level_stat = self.mdata.groupby(by=[level]).sum()
        return level_stat

    def add_desc(self,sub_item,id,level_stat):
        id_desc =self.anno_file.ix[:,sub_item]
        id_desc.drop_duplicates(id, keep='first', inplace=True)
        id_desc.set_index(id,inplace=True)
        out_data = level_stat.join(id_desc)
        return out_data

    def rm_item(self,data,rm_list):
        for i in rm_list:
            if i in data.columns:
                del data[i]

    def write_file(self,data,name):    
        data.to_csv( name +'_out.xls',sep='\t')

    def main_func(self,level_list):  #level_list =[{''id':ID','out_name':'ID_out','desc':['ID','ID Description'],'rm_item':[]},]
        for i in level_list:
            id = i['id']
            out_name = i['out_name']
            level_stat = self.level_sum(id)
            if 'desc' in i.keys():
                sub_item = i['desc']
                out_data = self.add_desc(sub_item,id,level_stat)
            else:
                out_data = level_stat
            if 'rm_item' in i.keys():
                self.rm_item(out_data, i['rm_item'])
            self.write_file(out_data,out_name)

    def level_contain_gene(self,level_gene):   # level_gene = [gene ,level,outname]
        gene = level_gene[0]
        level = level_gene[1]
        out_name = level_gene[2]
        data = pd.DataFrame()
        data[gene] = self.anno_file[gene]
        data[level] = self.anno_file[level]
        data_dic = dict(list(data.groupby(by=[level])))
        with open(out_name + '.xls','w') as f:
            f.write(level+'\tGene_counts\tGene_list\n')
            for k in data_dic.keys():
                gene_num = len(data_dic[k][gene])
                gene_detail = ';'.join(data_dic[k][gene])
                f.write(str(k)+'\t'+str(gene_num)+'\t'+gene_detail+'\n')

    def batch_level_contain_gene(self,level_gene_list):
        for i in level_gene_list:
            self.level_contain_gene(i)



def TCDB_stat(read_xls,anno_xls):   
    cls = Table2stat(read_xls,anno_xls,'GeneID','#Query',False, ['Identity(%)','Evalue','Score','Align_len'])
    level_list = [{'id':'TCDB ID','out_name':'tcdb_abund','desc':['TCDB ID','TCDB Description'],'rm_item':['class_id']}]
    level_list.append({'id':'class_id','out_name':'class_abund','desc':['class_id','TCDB Class']})
    level_list.append({'id':'subclass_id','out_name':'subclass_abund','desc':['subclass_id','TCDB Subclass'],'rm_item':['class_id']})
    level_list.append({'id':'family_id','out_name':'family_abund','desc':['family_id','TCDB Family'],'rm_item':['class_id']})
    cls.main_func(level_list)

    level_stat_list = [['#Query','TCDB ID','TCDB_gene_stat']]
    level_stat_list.append(['#Query','TCDB Class','Class_gene_stat'])
    level_stat_list.append(['#Query','TCDB Subclass','SubClass_gene_stat'])
    level_stat_list.append(['#Query','TCDB Family','Family_gene_stat'])
    cls.batch_level_contain_gene(level_stat_list)

def Mvirdb_stat(read_xls,anno_xls):
    cls = Table2stat(read_xls,anno_xls,'GeneID','#Query',False, ['Identity(%)','Evalue','Score','Align_len'])
    level_list = []
    level_list.append({'id':'Virulence Factor ID','out_name':'Factor_abund','desc':['Virulence Factor ID','Short Description']})
    level_list.append({'id':'Virulence Factor Type','out_name':'Type_abund' ,'rm_item':['Virulence Factor ID']})
    cls.main_func(level_list)

    level_stat_list = [['#Query','Virulence Factor ID','Factor_gene_stat']]
    level_stat_list.append(["#Query",'Virulence Factor Type','Type_gene_stat'])
    cls.batch_level_contain_gene(level_stat_list)



if __name__ == '__main__':
    read_xls = sys.argv[1]
    anno_xls = sys.argv[2]
    anno_type = sys.argv[3]
    if anno_type == 'tcdb':
        TCDB_stat(read_xls,anno_xls)
    elif anno_type == 'mvirdb':
        Mvirdb_stat(read_xls,anno_xls)
