# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modifies 20180409

"""circos 图"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.config import Config


class CircosWorkflow(Workflow):
    """
    报告中使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(CircosWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pre_file", "type": "infile", "format": "bacgenome.circos_pre_dir"},
            {"name": "task_id", "type": "string"},
            {"name": "specimen_id", "type": "string"},
            {"name": "location", "type": "string", "default": "Scaffold"},
            {"name": "para1", "type": "infile", "format": "sequence.profile_table"},
            {"name": "para2", "type": "infile", "format": "sequence.profile_table"},
            {"name": "para3", "type": "infile", "format": "sequence.profile_table"},
            {"name": "para4", "type": "infile", "format": "sequence.profile_table"},
            {"name": "para5", "type": "infile", "format": "sequence.profile_table"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "para_def", "type": "infile", "format": "sequence.profile_table"}, # 自定义的区段
            {"name": "labs", "type": "string", "default": ""}, #分号分割
            {"name": "color_version", "type": "string","default":""},
            {"name": "seq_type", "type": "string", "default": "Circular"},  # 开闭环

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.plot = self.add_tool("bacgenome.plot_circos")
        self.file_path = self._sheet.output

    def change_seq_name(self, infile):
        with open(infile) as def_p, open(infile+'_new','w') as fw:
            for line in def_p:
                spline = line.split(' ')
                spline[0] = self.seq_name
                fw.write(' '.join(spline))


    def run_plot(self):
        list = []
        self.antisense_cog = self.option('pre_file').prop['path'] + "/antisense_strand_cog.txt"
        self.sense_cog = self.option('pre_file').prop['path'] + "/sense_strand_cog.txt"

        karyotype =self.option('pre_file').prop['path'] + "/karyotype.txt"
        with open(karyotype) as f:
            spline = f.readline().split(' ')
            self.seq_name = spline[2]
            self.seq_all_length = float(spline[5])

        if self.option('color_version') in ['1.0','1'] :
            self.change_cog_v1_color()
        if self.option('para_def').is_set:  #zouguanqing 20190410
            para_def = self.option('para_def').prop['path']
            self.change_seq_name(para_def)
            list.append(para_def+'_new')

        if self.option('para1').is_set:
            para1 = self.option('para1').prop['path']
            self.change_seq_name(para1)
            list.append(para1+"_new")
        if self.option('para2').is_set:
            para2 = self.option('para2').prop['path']
            self.change_seq_name(para2)
            list.append(para2+"_new")
        if self.option('para3').is_set:
            para3 = self.option('para3').prop['path']
            self.change_seq_name(para3)
            list.append(para3+'_new')
        if self.option('para4').is_set:
            para4 = self.option('para4').prop['path']
            self.change_seq_name(para4)
            list.append(para4+'_new')
        if self.option('para5').is_set:
            para5 = self.option('para5').prop['path']
            self.change_seq_name(para5)
            list.append(para5+'_new')
        f = ','.join(list)
        options = {
            "k": self.option('pre_file').prop['path'] + "/karyotype.txt",
            "c":  self.sense_cog,
            "t": self.option('pre_file').prop['path'] + "/temp.txt",
            "ac": self.antisense_cog ,
            "pgc": self.option('pre_file').prop['path'] + "/positive_gc_count.txt",
            "ngc": self.option('pre_file').prop['path'] + "/negative_gc_count.txt",
            "pgs": self.option('pre_file').prop['path'] + "/positive_gc_skew.txt",
            "ngs": self.option('pre_file').prop['path'] + "/negative_gc_skew.txt",
            "f": f,
            "labs": 'cog;'+self.option('labs'),
            "seq_type": self.option('seq_type')
        }
        #if self.option('location').startswith('p'):
        if self.seq_all_length < 500000:
            options['type'] = 2
        else:
            options['type'] = 1
        self.plot.set_options(options)
        self.plot.on('end', self.set_db)
        self.plot.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_plot()
        super(CircosWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        result = self.plot.output_dir
        api_circos = self.api.api('bacgenome.circos')

        api_circos.add_circos(result, main=False, main_id=self.option('main_id'),
                              location=self.option('location'), specimen_id=self.option('specimen_id'),
                              task_id=self.option('task_id'), update=False, link_path=self.file_path)
        self.end()

    def end(self):
        repaths = [
            [".", "", ""],
        ]
        regexps = [
        ]
        sdir = self.add_upload_dir(self.plot.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CircosWorkflow, self).end()


    def change_cog_v1_color(self,new_color_version=2):
        db_name = 'bacgenome'
        db = Config().get_mongo_client(mtype=db_name, ref=True)[Config().get_mongo_dbname(db_name, ref=True)]
        pre_map_new_color = {}
        for i in db.COG_color.find({"version":new_color_version}):
            new_color = i['rgb']
            cog = i['fid']
            pre_color_info = db.COG_color.find_one({"_id":cog})
            if pre_color_info:
                pre_color = pre_color_info['rgb']
                pre_map_new_color[pre_color] = new_color

        for f in [self.sense_cog,self.antisense_cog]:
            with open(f) as fr, open(f+'_new','w') as fw:
                for line in fr:
                    line = line.rstrip()
                    spline = line.split(' ')
                    pre_c = spline[3].split('=')[1]
                    if pre_c in pre_map_new_color:
                        spline[3] = 'fill_color=%s'%pre_map_new_color[pre_c]
                    fw.write(' '.join(spline)+"\n")
        self.sense_cog = self.sense_cog +"_new"
        self.antisense_cog = self.antisense_cog + "_new"