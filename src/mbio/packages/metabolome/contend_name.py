# -*- coding:utf-8 -*-
# /usr/bin/python
__author__ = 'guanqing.zou'
import argparse
import pandas as pd
import os

"""
作用：
阳离子表和阴离子表的Metabolite列有相同的代谢名称，按照规则保留一个，另一个名称改成‘-‘。
规则如下：
A：如果上传的数据表格有：Fragmentation Score（实际匹配得分：第一优先级）//Theoretical Fragmentation Score（理论匹配得分）
那么根据得分进行判断清除：
如果两个均为Fragmentation Score 则清除得分低的；
如果两个均为Theoretical Fragmentation Score 则清除得分低的；
如果只有一个有Fragmentation Score，另一个为Theoretical Fragmentation Score或没有，清除有Theoretical Fragmentation Score或没有得分的那个代谢物；
如果只有一个有Theoretical Fragmentation Score，另一个没有得分值，清除没有得分值的那个代谢物；
B：如果两个均没有这两个得分值（得分值为零，等价于没有得分值），那么根据QC样本的RSD（峰度值的标准差除以平均数据） 进行判断清除：
清除RSD 大的那个代谢物；
C：如果没有得分值也没有QC 样本，那么根据两个峰在所有样本中的强度和进行判断清除：
清除强度和小的那个代谢物；
D:如果根据峰强度无法判断，那么清除阴离子模式下的代谢物；

"""
class Contend_Name():
    def __init__(self,pos_file,neg_file,target='Metabolite',not_dispose=None):
        if not  os.path.exists(pos_file):
            raise("not exists: "+ pos_file)
        if not os.path.exists(neg_file):
            raise("not exists: "+ neg_file)

        self.not_dispose = None
        new_neg_file = self.check_change_name(pos_file,neg_file, not_dispose)
        self.pos_data = pd.read_table(pos_file,sep='\t')
        self.neg_data = pd.read_table(new_neg_file,sep='\t')
        self.allot = {'pos':[],'neg':[]}
        self.target = target
        self.get_same_name_data()



    def get_same_name_data(self):
        target = self.target
        pos_meta = self.pos_data[target]
        neg_meta = self.neg_data[target]
        self.same_names = set(pos_meta.tolist()) & set(neg_meta.tolist())


        self.pos_sub = self.pos_data[self.pos_data[target].apply(lambda x: x in self.same_names)]
        self.neg_sub = self.neg_data[self.neg_data[target].apply(lambda x: x in self.same_names)]
        if '-' in self.same_names:
            self.same_names.remove('-')


    # 可继承重写
    '''
    group:
    #sample 	group
    IH207_2A	IH_2A
    IH201_3A	IH_3A
    IH205_3A	IH_3A
    QC_1	QC
    QC_2	QC

    '''
	#sample_info 样本分组文件，从中获取样本的名称和QC的样本名称
    def contend_pip(self,project_name='metabolome',sample_info=None):
        if project_name=='metabolome':
            priority = ['Fragmentation Score','Theoretical Fragmentation Score','rsd','sum']
            if sample_info:
                if os.path.exists(sample_info):
                    group_data = pd.read_table(sample_info,sep='\t',header=0)
                    group_data.columns = ['sample','group']
                    all_samples = group_data['sample'].tolist()
                    qc_samples = group_data[group_data['group'] == 'QC']['sample'].tolist()
                    if len(qc_samples)==0:
                        priority.remove('rsd')
                    if len(all_samples) ==0:
                        priority.remove('sum')
                else:
                    raise('not exisits: %s'% sample_info)
            else:
                priority.remove('rsd')
                priority.remove('sum')

            for sifter in priority:
                if len(self.same_names) == 0 :
                    break
                if sifter in ['Fragmentation Score','Theoretical Fragmentation Score']:
                    self.by_one_column(sifter)
                elif sifter == 'rsd':
                    def rsd_fun(v):
                        return v.std()/v.mean()
                    self.pick_fun(qc_samples,rsd_fun,pick_big=False)

                elif sifter == 'sum':
                    def sum_fun(v):
                        return sum(v)
                    self.pick_fun(all_samples,sum_fun,pick_big=True)
            if len(self.same_names) > 0:   #剩下的归于pos
                self.update(self.same_names,pos=True)

            ## 名称改成’-‘，并把它的相关信息改成’-‘
            clear_k = set(['Metabolite','Adducts','Formula','Fragmentation Score','Theoretical Fragmentation Score',
                       'Library ID','CAS ID','KEGG Compound ID'])

            pos_clear_k = clear_k - set(self.pos_data.columns.tolist())
            neg_clear_k = clear_k - set(self.neg_data.columns.tolist())
            current_clear = {'pos':pos_clear_k,'neg':neg_clear_k}

            for lab, df in [('neg',self.pos_data),('pos',self.neg_data)]:
                ori_columns = df.columns
                df = df.set_index(self.target)
                for name in self.allot[lab]:
                    for ck in current_clear[lab]:
                        df.loc[name,ck] = '-'
                df = df.reset_index()
                if self.allot[lab] != []:
                    df[self.target].replace(self.allot[lab],'-',inplace = True)

                if lab == 'pos':
                    out_file = './neg.rm_rep.csv'
                else:
                    out_file = './pos.rm_rep.csv'

                # print('write out file')
                # if lab == 'neg':
                #     if self.not_dispose:
                #         print(self.not_dispose)
                #         #df[self.target].replace(self.not_dispose,self.ori_name,inplace=True)
                #         for i in range(len(df)):
                #             if df['Metabolite'][i] == self.not_dispose:
                #                 print('get %s'% self.not_dispose)
                #                 df['Metabolite'][i] = self.ori_name
                #                 print('change to %s'%self.ori_name)
                #         print(self.target)
                #         print('change name %s to %s'%(self.not_dispose, self.ori_name))
                df.to_csv(out_file,sep='\t',index=False,columns=ori_columns, quoting=3)

            if self.not_dispose:
                self.change_back_not_dispose('neg.rm_rep.csv')


    def by_one_column(self,column,pick_big=True):
        target = self.target
        if column in self.pos_sub.columns and column not in self.neg_sub.columns:
            pick_out = self.pos_sub[self.pos_sub[column].apply(lambda x: x  not in ['-',0])]
            tmp = pick_out[target].tolist()
            self.update(tmp,pos=True)

        elif column not in self.pos_sub.columns and column  in self.neg_sub.columns:
            pick_out = self.neg_sub[self.neg_sub[column].apply(lambda x: x not in ['-',0])]
            tmp = pick_out[target].tolist()
            self.update(tmp,pos=False)

        elif column in self.pos_sub.columns and column  in self.neg_sub.columns:
            def change_v(v):
                if  v == '-':
                    return 0
                else:
                    return v
            self.pick_fun(column,change_v,pick_big)


    def update(self,in_list,pos=True):
        target = self.target
        if pos:
            self.allot['pos'].extend(in_list)
        else:
            self.allot['neg'].extend(in_list)
        self.same_names = self.same_names - set(in_list)
        if len(self.same_names) > 0:
            self.pos_sub = self.pos_sub[self.pos_sub[target].apply(lambda x: x in self.same_names)]
            self.neg_sub = self.neg_sub[self.neg_sub[target].apply(lambda x: x in self.same_names)]

    def _pick_extremum(self,obj, max=True):  #obj [(22,'a'),(33,'b'),(44,'c')], return lab
        if max:
            v =sorted(obj, key=lambda x:x[0],reverse=True)
            if v[0][0] > v[1][0]:
                return v[0][1]
        else:
            v =sorted(obj, key=lambda x:x[0],reverse=False)
            if v[0][0] < v[1][0]:
                return v[0][1]
        return None

    def pick_fun(self,column,change_fun,pick_big):
        pos_tmp = []
        neg_tmp = []
        self.pos_sub = self.pos_sub.set_index(self.target)
        self.neg_sub = self.neg_sub.set_index(self.target)
        for name in self.same_names:
            sort_ = []
            for sub,lab in [(self.pos_sub,'pos'),(self.neg_sub,'neg')]:
                v = sub.loc[name,column]
                v = change_fun(v)
                sort_.append((v,lab))
            pick_ret = self._pick_extremum(sort_, max=pick_big)
            if not pick_ret:
                if pick_ret == 'pos':
                    pos_tmp.append(name)
                else:
                    neg_tmp.append(name)
        self.pos_sub = self.pos_sub.reset_index()
        self.neg_sub = self.neg_sub.reset_index()
        self.update(pos_tmp,pos=True)
        self.update(neg_tmp,pos=False)

    def check_change_name(self,pos_file,neg_file, not_dispose=None):   # pos和neg代谢物名称 不区分大小写和末尾的空格。neg保持和pos一致

        pos = pd.read_table(pos_file,sep='\t')
        pos_meta = pos['Metabolite'].tolist()
        new_map = {}
        for p in pos_meta:
            new_p = p.strip().lower()
            if new_p in ['-','']:
                continue
            if new_p not in new_map:
                new_map[new_p] = p

        neg = pd.read_table(neg_file,sep='\t')
        for i in range(len(neg)):
            new_n = neg['Metabolite'][i].strip().lower()
            if new_n in new_map:
                neg['Metabolite'][i] = new_map[new_n]
                if not_dispose:
                    if neg['Metabolite'][i] == not_dispose:
                        self.ori_name = not_dispose
                        neg['Metabolite'][i] = not_dispose +'_NegDispose'
                        self.not_dispose = not_dispose +'_NegDispose'


        neg.to_csv('neg_file_check_1',sep='\t',index=False, quoting=3)
        return 'neg_file_check_1'

    def change_back_not_dispose(self, infile):
        data = pd.read_table(infile,sep='\t',header=0)
        data[self.target].replace(self.not_dispose, self.ori_name,inplace=True)
        print('change_back_not_dispose')
        data.to_csv(infile,sep='\t',index=False, quoting=3)


if __name__ ==  '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-p','--pos',required=True, help='pos file')
    parse.add_argument('-n','--neg',help='neg file',required=True)
    parse.add_argument('-g','--group', help='group file')
    args = parse.parse_args()
    CN = Contend_Name(args.pos,args.neg)
    if args.group:
        CN.contend_pip(sample_info=args.group)
    else:
        CN.contend_pip()




