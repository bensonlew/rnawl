# -*- coding: utf-8 -*-


from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import json
import pandas as pd


class RocBatchModule(Module):
    def __init__(self, work_id):
        super(RocBatchModule, self).__init__(work_id)
        options = [
            {"name": "exp_table", "type": "infile", "format": "sequence.profile_table"},  #表达量表
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组表，只接受两个分组。接口用to_file 生成
            {"name": "metab_set_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "confidence_interval", "type": "float", "default": 0.95}
        ]
        self.add_option(options)
        self.map_tool = {}

    def check_options(self):

        if not self.option("exp_table").is_set:
            raise OptionError("请传入丰度矩阵！", code="24700201")
        if not self.option("metab_desc").is_set:
            raise OptionError("请传入metab_desc文件！", code="24700202")
        if not self.option("metab_set_table").is_set:
            raise OptionError("请传入metab_set_table文件！", code="24700202")

        return True

    def run_roc(self):
        table = pd.read_table(self.option('exp_table').path, sep='\t',header=0)
        metab_table = pd.read_table(self.option('metab_set_table').path, sep='\t',index_col=0)
        for metab_id in metab_table.index:
            new_table = table[table['metab_id'] == metab_id]
            if len(new_table) !=0:
                metab_id_table = self.work_dir+'/%s_exp_table.xls'%metab_id
                new_table.to_csv(metab_id_table, sep='\t',index=False)
            else:
                self.logger.info('%s not exist in exp_table'%metab_id)
                continue

            roc_tool = self.add_tool("statistical.roc")
            options = {
                'abu_table': metab_id_table,
                'group_table': self.option('group_table'),
                'method': 'sum',
                'top_n': 1,
                'confidence_interval': str(self.option('confidence_interval')),
            }
            roc_tool.set_options(options)
            self.map_tool[metab_id] = roc_tool

        self.on_rely(self.map_tool.values(), self.set_output)

        with open(self.work_dir+'/metab_id_map_tool.xls','w') as fw:
            for metab_id in self.map_tool:
                fw.write('%s\t%s\n'%(metab_id, self.map_tool[metab_id]._id))

        for tool in self.map_tool.values():
            tool.run()

    def set_output(self):

        files = ['best_loc.xls','roc_auc_smooth.xls','roc_auc.xls','roc_curve_smooth.xls','roc_curve.xls','roc_interval.xls']
        self.data = dict()
        for file in files:
            self.data[file] = []
        #files 第一列添加metab_id ,添加到 self.data
        for metab in self.map_tool:
            tool = self.map_tool[metab]
            for file in self.data.keys():
                file_path = tool.output_dir + '/' + file
                if not os.path.exists(file_path):
                    continue
                data = pd.read_table(file_path,sep='\t')
                cols = data.columns.tolist()
                cols.insert(0,'metab_id')
                data['metab_id'] = metab
                data = data.reindex(columns=cols)
                self.data[file].append(data)

        for file in self.data:
            if len(self.data[file]) == 0:
                continue
            cat_data = pd.concat(self.data[file],axis=0)
            out_file = self.output_dir +'/' + file
            cat_data.to_csv(out_file, sep='\t', index=False)

        self.end()

    def run(self):
        super(RocBatchModule, self).run()
        self.run_roc()

    def end(self):
        super(RocBatchModule, self).end()

