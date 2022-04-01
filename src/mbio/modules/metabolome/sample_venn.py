# -*- coding: utf-8 -*-


from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import json
import pandas as pd


class SampleVennModule(Module):
    def __init__(self, work_id):
        super(SampleVennModule, self).__init__(work_id)
        options = [
            {"name": "metab_table", "type": "infile", "format": "metabolome.express"},
            {"name": "group_detail", "type": "string"},
            {"name": "threshold", "type": "float", "default":50},   #指 50%
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)
        self.venn = self.add_tool("metabolome.metabset.venn")

    def check_options(self):
        if not self.option("metab_table").is_set:
            raise OptionError("请传入丰度矩阵！", code="24700201")

        if self.option('group_detail'):
            try:
                self.group_detail = json.loads(self.option('group_detail'))
            except Exception as e:
                raise OptionError("group_detail  Format Is Wrong")
        else:
            if not self.option("group_table").is_set:
                 raise OptionError("请传入分组文件！", code="24700201")
            self.group_detail = self.get_group_detail(self.option('group_table').path)

        return True

    def get_group_detail(self,infile):
        data = pd.read_table(infile,sep='\t')
        data.columns = ['sample','group']
        j = dict(list(data.groupby(by=['group'])))
        for g in j.keys():
            j[g] = j[g]['sample'].tolist()
        return j

    def create_list_file(self):
        group_detail = self.group_detail
        data = pd.read_table(self.option('metab_table').path, sep='\t')
        threshold = float(self.option("threshold"))/100
        group_metab = []
        group_names = []
        def fill_fun(one,g):
            samples = group_detail[g]
            s_num = len(samples)
            count = 0
            for s in samples:
                if one[s]  in ['-',0]:
                    count +=1
            per = float(count)/s_num
            if per >= threshold:
                return False   #
            else:
                return True

        for g  in group_detail:
            group_names.append(g)
            sub = data[data.apply(lambda x: fill_fun(x,g), axis=1)]
            group_metab.append(','.join(sub['metab_id']))

        out_data = pd.DataFrame({"group": group_names, "metab_id": group_metab})
        out_data.to_csv(self.work_dir+'/group_metab.list',sep='\t', header=None,index=False)

    def run(self):
        super(SampleVennModule, self).run()
        self.create_list_file()
        options = {
            "list_file": self.work_dir+'/group_metab.list'
        }
        self.venn.set_options(options)
        self.venn.on('end', self.set_output)
        self.venn.run()

    def set_output(self):
        # if os.path.exists(self.output_dir):
        #     shutil.rmtree(self.output_dir)
        if os.path.exists(self.output_dir+'/metabset_detail.xls'):
            os.remove(self.output_dir+'/metabset_detail.xls')
        os.link(self.venn.output_dir+'/metabset_detail.xls' , self.output_dir+'/metabset_detail.xls')
        self.rm_other()

        self.end()

    def rm_other(self):
        ori_venn = self.venn.output_dir+'/venn_table.xls'
        group_ele = self.work_dir+'/group_metab.list'
        group_map = {}
        with open(group_ele) as f:
            for line in f:
                spline = line.strip().split('\t')
                group_map[spline[0]] = set(spline[1].split(','))
        groups = set(group_map.keys())

        new_venn = self.output_dir + '/venn_table.xls'
        with open(new_venn,'w') as fw, open(ori_venn) as fr:
            for line in fr:
                spline = line.strip().split('\t')
                if '&' in spline[0]:
                    cg = set([i.strip() for i in spline[0].split('&')])
                    other_g = groups - cg
                    cur_ele = set(spline[2].split(','))
                    for og in other_g:
                        cur_ele  =cur_ele - group_map[og]

                    s_cur_ele = sorted(cur_ele,key=lambda x: int(x.split('_')[1]))
                    fw.write(spline[0]+'\t'+str(len(s_cur_ele))+'\t'+','.join(s_cur_ele)+'\n')

                else:
                    fw.write(line)

    def end(self):
        super(SampleVennModule, self).end()

