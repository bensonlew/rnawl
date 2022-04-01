# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'


import os
import re
import types
import codecs
import chardet
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class RandomForestWorkflow(Workflow):
    """
    RandomForestWorkflow
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RandomForestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"input_table","type":"infile", "format":"meta.otu.otu_table"},
            {"name":"group_table","type":"infile", "format":"toolapps.group_table,meta.otu.group_table"},
            {"name":"norm_method","type": "string", "default": ""},
            {"name":"method","type": "string", "default": "CV"},
            {"name":"main_id","type":"string"},
            {"name":"update_info","type":"string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.random_forest = self.add_tool("tool_lab.randomforest")
    
    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_table"):
            raise OptionError("必须输入表格数据")
        if not self.option("group_table"):
            raise OptionError("必须输入分组表格")
        if self.option("norm_method") not in ["none","z","minmax","log","relative"]:
            raise OptionError("输入的数据标准化方法不在可选范围内")


    def change_otuname(self, tablepath, file_name):
        newtable = self.work_dir + "/" + file_name + "_input_abund.xls"
        content = open(tablepath,'r').read()
        encoding = chardet.detect(content)['encoding']
        oldFile = codecs.open(tablepath,'r',encoding=encoding)
        content = oldFile.read()
        oldFile.close()
        newFile = codecs.open(tablepath,'w','ascii')
        newFile.write(content)
        newFile.close()
        with open(tablepath, "r") as f, open(newtable, "w") as g:
            head = f.readline()
            g.write(head)
            for line in f:
                lines = line.split("\t", 1)
                specimen = re.subn("^.*; ", "", lines[0])[0]
                g.write(specimen + "\t" + lines[1])
        return newtable

    def run_randomforest(self):
        newtable = self.change_otuname(self.option('input_table').prop['path'], "otutable")
        if self.option('method') == "none":
            method = "CV"
        else:
            method = self.option('method')
        if self.option("norm_method") == "none":
            norm_method = ""
        else:
            norm_method = self.option("norm_method")
        self.group = self.check_group()
        options = {
            'otutable': newtable,
            'method': method,
            'grouptable': self.group,
            'ntree': 500,
            'norm_method': norm_method
        }
        self.random_forest.set_options(options)
        self.random_forest.on('end', self.set_output)
        self.random_forest.run()

    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        os.link(os.path.join(self.random_forest.output_dir,"randomForest_imptance_table.xls")
                ,os.path.join(self.output_dir,"randomForest_imptance_table.xls"))
        if self.option('method') == "CV":
            os.link(os.path.join(self.random_forest.output_dir,"randomForest_10-fold_CV.xls")
                    ,os.path.join(self.output_dir,"randomForest_10-fold_CV.xls"))
        elif self.option('method') == "AUC":
            os.link(os.path.join(self.random_forest.output_dir,"randomForest_AUC.xls")
                    ,os.path.join(self.output_dir,"randomForest_AUC.xls"))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_random_forest = self.api.api("tool_lab.random_forest")
        api_random_forest.add_column_detail(self.option("main_id"),os.path.join(self.output_dir,"randomForest_imptance_table.xls"))
        api_random_forest.add_table_detail(self.option("main_id"),os.path.join(self.output_dir,"randomForest_imptance_table.xls"))
        if self.option('method') == "CV":
            api_random_forest.add_line_detail(self.option("main_id"),os.path.join(self.output_dir,"randomForest_10-fold_CV.xls"))
            api_random_forest.add_tooltip(self.option("main_id"),"CV")
        elif self.option('method') == "AUC":
            api_random_forest.add_line_detail(self.option("main_id"),os.path.join(self.output_dir,"randomForest_AUC.xls"))
            api_random_forest.add_tooltip(self.option("main_id"),"AUC")
        else:
            api_random_forest.add_tooltip(self.option("main_id"),"none")
        self.logger.info("导表结束")
        self.end()

    def check_group(self):
        group_path = os.path.join(self.work_dir,"group_tem.xls")
        old_path = self.option("group_table").prop['path']
        content = open(old_path,'r').read()
        encoding = chardet.detect(content)['encoding']
        oldFile = codecs.open(old_path,'r',encoding=encoding)
        content = oldFile.read()
        oldFile.close()
        newFile = codecs.open(old_path,'w','ascii')
        newFile.write(content)
        newFile.close()
        with open(group_path,"w") as new_group:
            with open(self.option("group_table").prop['path'],"r") as group:
                line = group.readline()
                if line[1] == "#":
                    new_group.write(line)
                else:
                    new_group.write("#")
                    new_group.write(line)
                title = line.rstrip().split('\t')[1:]
                group_num = {}
                group_list = {}
                for i in title:
                    group_num[i] = []
                    group_list[i] = []
                while 1:
                    line = group.readline()
                    if not line:
                        break
                    new_group.write(line)
                    fd = line.rstrip().split('\t')
                    for i in range(len(title)):
                        group_list[title[i]].append(fd[i+1])
                        if fd[i+1] not in group_num[title[i]]:
                            group_num[title[i]].append(fd[i+1])
                for i in title:
                    if len(group_num[i]) < 2:
                        self.set_error("分组列表的组{}的分组数小于2".format(i))
                    if self.option("method") == "AUC" and len(group_num[i]) > 2:
                        self.set_error("当验证方式为AUC时，分组数只能等于2")
                    for j in group_num[i]:
                        if group_list[i].count(j) < 5:
                            self.set_error("分组列表的组{}的{}组小于5".format(i,j))
        return group_path


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(RandomForestWorkflow,self).end()

    def run(self):
        self.run_randomforest()
        super(RandomForestWorkflow, self).run()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    data = {
        'name': 'test_randomfrest',
        'id': 'randomfrest_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "input_table":"/mnt/ilustre/users/sanger-dev/workspace/20210430/RandomForest_tuqs_ogn4f8vk6guect81n28tkr_0430085604115228_2529/remote_input/input_table/otu_new.txt",
        "group_table":"/mnt/ilustre/users/sanger-dev/workspace/20210430/RandomForest_tuqs_ogn4f8vk6guect81n28tkr_0430085604115228_2529/remote_input/group_table/group_30.txt",
        "norm_method":"z",
        "method":"none",
        "main_id" : "5e9e6a6017b2bf2049a81b26"
        }
    }
    wsheet = Sheet(data=data)
    wf = RandomForestWorkflow(wsheet)
    wf.run()