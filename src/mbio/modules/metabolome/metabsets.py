# -*- coding: utf-8 -*-
# __author__ = zouguanqing
# last_modifiy = modified 201908

import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest
import pandas as pd
import numpy as np
import re

class MetabsetsModule(Module):
    """
    单个代谢集分析
    """
    
    def __init__(self,work_id):
        super(MetabsetsModule, self).__init__(work_id)
        options = [
            #公用
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "mul_metabset", "type": "infile", "format": "metabolome.mul_metabset"},
            ##keggp
            {"name": "ko_overview","type": "infile", "format":"sequence.profile_table" },
            #hmdb
            {"name": "hmdb_overview","type": "infile", "format":"sequence.profile_table" },
            # vip
            {'name': 'diff_dir', 'type': 'infile', 'format': 'annotation.mg_anno_dir'},
            ## 公用
            {"name": "metab_trans", "type": "infile", "format": "sequence.profile_table"},
            {"name": "metab_table", "type": "infile", "format": "sequence.profile_table"}, # 标准化时用，用原始丰度表
            ##cluster
            # {"name": "scale_exp","type": "infile", "format": "sequence.profile_table"},
            ###{"name": "skip_analysis","type": "string", "default":""}
            {"name": "species", "type": "string", "default": "all"},
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老工作流
            {"name": "database_version", "type": "string", "default": ""}, # 用于区分新老工作流

        ]
        self.add_option(options)



    def check_options(self):
        if not self.option("anno_overview").is_set:
            raise OptionError("请传入anno_overview！")
        if not self.option("mul_metabset").is_set:
            raise OptionError("请传入mul_metabset！")
        if not self.option("ko_overview").is_set:
            raise OptionError("请传入ko_overview！")
        # if not self.option("hmdb_overview").is_set:
        #     raise OptionError("请传入hmdb_overview！")
        return True

    def scale_data(self, file_path):
        table = pd.read_table(file_path, sep="\t", index_col=0)
        scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
        exp_profile = os.path.join(self.work_dir, "scale_data.xls")
        scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
        return exp_profile

    def divide_metabset(self):
        mul_metabset = self.option('mul_metabset').path
        dir = self.work_dir
        self.map_file = {}
        with  open(mul_metabset) as f:
            for line in f:
                line = line.strip()
                if line =="" :
                    continue
                spline = line.split('\t')
                set_name = spline[0]
                metab_num = len(spline[1].split(','))
                self.map_file[set_name] = {
                    'file':dir+'/set_'+set_name,
                    'file_row': dir+'/set_'+set_name+'.row',
                    'metab_num': metab_num,
                    'metab_list' : spline[1].split(',')
                }
                with open(dir+'/set_'+set_name,'w') as fw , open(dir+'/set_'+set_name+'.row','w') as fw2:
                    fw2.write(line+'\n')
                    fw.write('\n'.join(spline[1].split(',')))

    def rm_qc_sample(self):
        metab_table = self.option('metab_table').path
        data = pd.read_table(metab_table,sep='\t',header=0)
        cols = data.columns
        rm_qc_cols = []
        for c in cols:
            if not re.match('^QC\d*$', c):
                rm_qc_cols.append(c)
        new_data = data[rm_qc_cols]
        self.rm_qc_abund = self.work_dir+'/rm_QC_abund.xls'
        new_data.to_csv(self.rm_qc_abund,sep='\t',index=False)

    def _check_skip_hmdb(self,db,metab_set):
        sub_db = db[db['Metab_id'].apply( lambda  x:  x in metab_set)]
        if len(sub_db) > 0:
            return False
        else:
            return  True



    def run_modules(self):
        self.divide_metabset()
        self.rm_qc_sample()
        self.scale_exp = self.scale_data( self.rm_qc_abund)
        if self.option('hmdb_overview').is_set:
            self.hmdb_all_data = pd.read_table(self.option('hmdb_overview').path,sep='\t',header=0)
        else:
            self.hmdb_all_data = ''


        for set_name in  self.map_file.keys():

            if self.map_file[set_name]['metab_num'] == 0:
                self.logger.info('%s metab_num is 0, skip metabset_analysis module'%set_name)
                skip = 'all'
                continue
            elif self.map_file[set_name]['metab_num'] < 4:

                skip = 'cluster,vip,corr_cluster'
                if isinstance(self.hmdb_all_data, str):
                    skip += ',hmdb'
                else:
                    if self._check_skip_hmdb(self.hmdb_all_data,self.map_file[set_name]['metab_list']):
                        skip += ',hmdb'
            else:
                skip = ''
                if isinstance(self.hmdb_all_data ,str):
                    skip = 'hmdb'
                else:
                    if self._check_skip_hmdb(self.hmdb_all_data,self.map_file[set_name]['metab_list']):
                        skip = 'hmdb'

            metab_module = self.add_module("metabolome.metabset_analysis")
            opts = {
                'metabset' : self.map_file[set_name]['file'],
                'metabset_row' : self.map_file[set_name]['file_row'],
                'anno_overview': self.option('anno_overview'),
                'ko_overview': self.option('ko_overview'),
                'metab_trans': self.option('metab_trans'),
                'metab_table': self.rm_qc_abund,       #self.option('metab_table'),
                'scale_exp' : self.scale_exp,
                'species' : self.option("species")
            }
            if 'hmdb' not in skip:
                opts['hmdb_overview'] = self.option('hmdb_overview')

            #如果self.option("diff_dir") 路径下不存在diff.exp.xls ，则不作vip分析。

            if not self.option("diff_dir").is_set or not os.path.exists(self.option("diff_dir").path+'/{}.diff.exp.xls'.format(set_name)):
                skip = skip+',vip'
                opts['skip_analysis'] = skip
            else:
                opts['diff_dir'] = self.option('diff_dir')
                opts['skip_analysis'] = skip
            opts['database_version'] = self.option("database_version") ## add by qingchen.zhang用于区分新老工作流
            self.logger.info('%s skip %s'%(set_name,skip))
            metab_module.set_options(opts)
            self.map_file[set_name]['module'] = metab_module



    def run(self):
        super(MetabsetsModule, self).run()
        self.run_modules()
        module_list = [ self.map_file[k]['module'] for k in self.map_file]
        if len(module_list) == 1:
            module_list[0].on('end',self.set_output)
        elif len(module_list) > 1:
            self.on_rely(module_list,self.set_output)
        for m in module_list:
            m.run()

    def set_output(self):
        dirs = [
            'MetabsetCluster',
            'MetabsetVip',
            'MetabsetCorr',
            'MetabsetKeggc',
            'MetabsetKeggp',
            'MetabsetEnrich',
            'MetabsetHmdb' ,
            'MetabsetIpath'
        ]
        for d in dirs:
            out_dir = self.output_dir +'/' +d
            if os.path.exists(out_dir):
                shutil.rmtree(out_dir)
            os.mkdir(out_dir)

        for set_name in self.map_file:
            if 'module' not in self.map_file[set_name].keys():
                 continue
            for d in dirs:
                ori_dir = self.map_file[set_name]['module'].output_dir + '/' + d
                if not os.path.exists(ori_dir):
                    continue
                tar_dir = self.output_dir +'/'+ d + '/' + set_name
                shutil.copytree(ori_dir,tar_dir)

        self.end()

    def end(self):
        super(MetabsetsModule, self).end()



class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "metabset2",
            "type": "module",
            "name": "metabolome.metabsets",
            "instant": False,
            "options": dict(
                anno_overview = '/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/Anno/AnnoOverview/output/anno.xls',
                mul_metabset = '/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/CreatTable/output/merge_mul.metabset.xls',
                ko_overview = '/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/Anno/AnnoOverview/output/ko.xls',
                hmdb_overview = "/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/AnnoHmdb/output/HmdbLevel_Origin.xls",
                diff_dir = "/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/DiffPls/MergeDiff/output/tmp_DiffStat",
                metab_trans = "/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/Preprocess/output/mix/metab_desc.txt",
                metab_table="/mnt/ilustre/users/sanger-dev/workspace/20190816/Metabolome_tsg_35242/Preprocess/output/mix/metab_abund.txt",
               # scale_exp=
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

