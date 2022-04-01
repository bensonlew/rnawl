#! -*- coding:utf-8 -*-
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
import json
import pandas as pd
import datetime
import glob
import re
import os


class ProjectStatistics(Base):
    def __init__(self, bind_object):
        super(ProjectStatistics, self).__init__(bind_object)
        self._project_type = "metabolome"


    def get_group_info(self):
        self.group_info = dict()
        with open(self.group_file) as f:
            f.readline()
            for line in f:
                if line.strip() == '':
                    continue
                s,g = line.strip().split('\t')
                if g not in self.group_info:
                    self.group_info[g] = []
                self.group_info[g].append(s)

        specimen_count = 0
        for g in self.group_info:
            specimen_count += len(self.group_info[g])

        self.overview_data['specimen_count'] = specimen_count
        if 'QC' in self.group_info:
            self.overview_data['qc_specimen_count'] = len(self.group_info['QC'])
        else:
            self.overview_data['qc_specimen_count'] = 0


    def deal_file(self,abund):
        #statistic = self.output_dir + '/Preprocess/statistic.xls'
        #raw_statistic = self.output_dir + '/Preprocess/raw_statistic.xls'
        #stat_r = pd.read_table(raw_statistic,sep='\t', header=0)


        exp_data = pd.read_table(abund,sep='\t',header=0)
        has_exp_list = exp_data['metab_id'].tolist()


        anno_view = self.output_dir + '/AnnoOverview/anno.xls'

        # stat_o = pd.read_table(statistic,sep='\t',header=0)
        # hmdb_n = 0
        # keggc_n = 0
        # all_n = 0
        # if len(stat_o) == 3:
        #     e = 2
        # else:
        #     e= 1
        # for i in stat_o.index[0:e]:
        #     all_n += stat_o.loc[i, 'all']
        #     hmdb_n += stat_o.loc[i, 'hmdb']
        #     keggc_n += stat_o.loc[i, 'kegg']

        anno_view_data = pd.read_table(anno_view,sep='\t',header=0)
        #统计有表达量的
        anno_view_data = anno_view_data[anno_view_data['metab_id'].apply(lambda x: True if x in has_exp_list else False)]

        all_n = len(anno_view_data)
        keggc_n = len(anno_view_data[anno_view_data['compound_id']!='-'])
        hmdb_n = len(anno_view_data[anno_view_data['hmdb_id'].apply(lambda x: 'HMDB' in x )])

        self.overview_data['kegg_c_count'] = keggc_n
        self.overview_data['hmdb_h_count'] = hmdb_n
        self.overview_data['kegg_c_percent'] = round(float(keggc_n)/all_n,4)
        self.overview_data['hmdb_h_percent'] = round(float(hmdb_n)/all_n,4)

        keggp_n = len(anno_view_data[anno_view_data['pathway_id']!='-'])
        self.overview_data['kegg_map_count'] = keggp_n
        self.overview_data['kegg_map_percent'] = round(float(keggp_n)/all_n,4)

        ##差异代谢物分析信息  最少，最多，中位数
        metab_set_num = dict()
        #metabset_data =  pd.read_table(self.metabset_file,sep='\t',header=0)
        with open(self.metabset_file) as f:
            for line in f:
                g,metabs = line.strip().split('\t')
                metab_set_num[g] = len(metabs.split(','))

        sort_metab_set = sorted(metab_set_num.items(),key=lambda x:x[1],reverse=False)
        mid_k = len(sort_metab_set)/2

        self.overview_data.update({
            'set_min' : sort_metab_set[0][1],
            "set_middle": sort_metab_set[mid_k][1],
            'set_max' : sort_metab_set[-1][1],
            "set_min_name" : sort_metab_set[0][0],
            "set_middle_name" : sort_metab_set[mid_k][0],
            "set_max_name" : sort_metab_set[-1][0]
        })


    def deal_nead_qc(self):
        ## RSD  output/Preprocess/cv_percent.xls  30%
        cv_percent = self.output_dir + '/Preprocess/cv_percent.xls'
        cv_data = pd.read_table(cv_percent,sep='\t',header=0,index_col=0)
        if self.project == 'GC':
            pos_rsd = cv_data.loc['per3','pos']
            self.overview_data['rsd'] = pos_rsd
        else:
            if self.mix == 'T':
                self.overview_data['rsd'] = cv_data.loc['per3','mix']
            else:
                self.overview_data['rsd'] = cv_data.loc['per3','pos']
                self.overview_data['neg_rsd'] = cv_data.loc['per3','neg']


        #
        #pca 离散度和跨度
        #output/ExpPCA/pos/PCA.sites.xls
        #remote_input/group_table/group.txt

        self.expcorr = list()

        if self.project == 'GC':
            type_list = ['pos']
        else:
            if self.mix == 'T':
                type_list = ['mix']
            else:
                type_list = ['pos','neg']
        for type in type_list:
            tmp = {
                'task_id' : self.task_id,
                'ion_mode' : type,
            }
            pca_site = self.output_dir + '/ExpPCA/'+type+ '/PCA.sites.xls'
            pca_data = pd.read_table(pca_site,sep='\t',header=0)
            pca_data.rename(columns={pca_data.columns[0]:'sample',pca_data.columns[1]:'pc1',pca_data.columns[2]:'pc2'},inplace=True)
            group_data = pd.read_table(self.group_file,sep='\t',header=0)
            group_data.columns = ['sample','group']
            data_add_group = pd.merge(group_data,pca_data,on='sample')
            qc_df = data_add_group[data_add_group.group == 'QC']
            extremum_pca_p1 = data_add_group['pc1'].max() - data_add_group['pc1'].min()
            extremum_pca_p2 = data_add_group['pc2'].max() - data_add_group['pc2'].min()
            extremum_qc_p1 = qc_df['pc1'].max() - qc_df['pc1'].min()
            extremum_qc_p2 = qc_df['pc2'].max() - qc_df['pc2'].min()
            extremum_p1 = round(extremum_qc_p1 / extremum_pca_p1,2)
            extremum_p2 = round(extremum_qc_p2 / extremum_pca_p2,2)
            tmp['qc_span'] = str(extremum_p1)+',' + str(extremum_p2)

            pca_g = data_add_group.groupby('group')[['pc1', 'pc2']].std()
            tmp['qc_seperate'] = str(round(pca_g.loc['QC', 'pc1'],2)) + ',' + str(round(pca_g.loc['QC', 'pc2'],2))
            self.expcorr.append(tmp)


    def deal_oplsda(self):
        # oplsda
        #
        # pos : pls_dir , neg: pls_dir
        # DiffPls/DiffMulStat__3/output/lic_vs_Cd/OPLS-DA.permMN.xls  第一行  R2Y(cum)        Q2(cum) sim， Q2(cum) 大于下面所有的值
        # DiffPls/DiffMulStat__3/output/lic_vs_Cd/OPLS-DA.intercept.xls  第一和第二行 的第2列 分别对应 r,q 的截距, q 的截距
        self.diff_metab= list()

        if self.oplsda_dir_neg:
            type_dir_map = [(' pos',self.oplsda_dir),(' neg',self.oplsda_dir_neg)]
        else:
            if self.project == 'LC':
                type_dir_map = [('mix',self.oplsda_dir)]
            else:
                type_dir_map = [('pos',self.oplsda_dir)]

        for t, oplsda_dir in type_dir_map:
            for d in os.listdir(oplsda_dir):
                diff_g = d
                tmp = {
                    "task_id" : self.task_id,
                    "name" : diff_g,
                    "type" : t
                }
                sub_d = os.path.join(oplsda_dir, d)
                perm_file = sub_d + '/OPLS-DA.permMN.xls'
                intercept_file = sub_d  + '/OPLS-DA.intercept.xls'
                with open(intercept_file) as f:
                    f.readline()
                    f.readline()
                    l3 = f.readline()
                    x0y = round(float(l3.strip().split('\t')[1]),4)

                perm = pd.read_table(perm_file,sep='\t',header=0)
                x1y = perm['Q2(cum)'][0]
                sy = sorted(perm['Q2(cum)'].tolist(), reverse=True)
                if sy[0] == x1y:
                    if sy[1] != x1y:
                        q2_left_high = 'yes'
                    else:
                        q2_left_high = 'no'
                else:
                    q2_left_high = 'no'

                slope = round((x1y-x0y)/1,2)
                q2_y_intercept = x0y
                tmp.update({
                    "q2_left_high" : q2_left_high,
                    "slope" : slope,
                    "q2_y_intercept": q2_y_intercept
                })
                self.diff_metab.append(tmp)


    ##其他数据完整性
    def check_result(self):
        self.error = list()
        no_workflow = ['diff_sample', 'metabset_kegg_heatmap', 'metabset_roc']
        main_map_detail = dict()
        main_map_key = dict()
        for d in self.db['table_relation'].find({"is_detail" : "y"}):
            if self.project == 'GC':
                if d['table'] in ['exp_neg','exp_mix']:
                    continue
            else:
                if self.mix == 'T':
                    if d['table'] in ['exp_neg','exp_pos']:
                        continue
                else:
                    if d['table'] in ['exp_mix']:
                        continue

            main = d['analysis']
            detail = d['table']
            if main in no_workflow:
                continue
            if main not in main_map_detail:
                main_map_detail[main] = []
                main_map_key[main] = d['asso_id']
            main_map_detail[main].append(detail)

        for main in main_map_detail:
            main_ret = self.db[main].find({'task_id':self.task_id,"status" : "end"})
            if main_ret:
                for mr in main_ret:
                    try:
                        submit_loc = json.loads(mr['params'])['submit_location']
                        main_id = mr['main_id']

                    except Exception as e:
                        continue
                    dts = main_map_detail[main]
                    for dt in dts:
                        d_ret = self.db[dt].find_one({main_map_key[main]: main_id})
                        if not d_ret:
                            tmp = {
                                'task_id' : self.task_id,
                                'sub_location' : submit_loc,
                                'main_id' : str(main_id),
                                'loss_table': dt,
                                'key' : main_map_key[main]
                            }
                            self.error.append(tmp)

    def project_check_pip(self, out_dir, group_file, metabset_dir,oplsda_dir,oplsda_dir_neg=None,abund=None):
        self.output_dir = out_dir
        self.group_file = group_file
        self.metabset_file = metabset_dir + '/merge_mul.metabset.xls'
        self.oplsda_dir = oplsda_dir
        if oplsda_dir_neg:
            self.oplsda_dir_neg = oplsda_dir_neg
        else:
            self.oplsda_dir_neg = None

        self.project = self.bind_object.sheet.option('project_type')
        self.mix = 'F'
        if self.project =='LC':
            self.mix = self.bind_object.sheet.option('mix_table')


        self.task_id = self.bind_object.sheet.id
        #样本数，项目类型，注释物种， pvalue/fdr值， vip ，fc
        try:
            vip_value = self.bind_object.sheet.option('vip_value')
        except Exception as e:
            vip_value = '0'

        self.overview_data = {
            'task_id' : self.task_id,
            'kegg_class' : self.bind_object.sheet.option('organism_type'),
            'vip' :  self.bind_object.sheet.option('vip_condition') + vip_value,
            'p_value_fdr' : self.bind_object.sheet.option('p_fdr_condition') + self.bind_object.sheet.option('p_fdr_value')
        }

        if self.project =='LC':
            self.overview_data['lc_merge'] = self.mix

        try:
            self.overview_data['fc_1'] =  self.bind_object.sheet.option('up_condition') +' '+ str(self.bind_object.sheet.option('up_value'))
        except Exception as e:
            pass

        try:
            self.overview_data['fc_2'] = self.bind_object.sheet.option('down_condition')+ ' ' + str(self.bind_object.sheet.option('down_value'))
        except Exception as e:
            pass

        try:
            self.overview_data['kegg_name'] = self.bind_object.sheet.option('organism')
        except Exception as e:
            pass


        self.get_group_info()
        self.deal_file(abund)
        if 'QC' in self.group_info:
            self.deal_nead_qc()

        if self.oplsda_dir:
            self.deal_oplsda()

        self.check_result()
        self.db['project_overview'].insert_one(self.overview_data)

        if 'QC' in self.group_info:
            if self.expcorr:
                self.db['project_expcorr'].insert_many(self.expcorr)

        if self.oplsda_dir:
            if self.diff_metab:
                self.db['project_diff_metab'].insert_many(self.diff_metab)

        if self.error:
            self.db['project_error'].insert_many(self.error)

















