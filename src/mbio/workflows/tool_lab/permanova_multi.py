# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

import os
import types
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class PermanovaMultiWorkflow(Workflow):
    """
    PERMANOVA分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PermanovaMultiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"input_table","type":"infile","format":"meta.otu.otu_table"},
            {"name":"factor_type","type":"string"},
            {"name": "env_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string", "default": 'bray_curtis'},
            {"name":"main_id","type":"string"},
            {"name":"update_info","type":"string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.permanova = self.add_tool("metaasv.permanova")

    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_table"):
            raise OptionError("必须输入表格数据")
        if not self.option("factor_type"):
            raise OptionError("请输入因素类型")
        elif self.option("factor_type") == "both":
            if not self.option("env_table") or not self.option("group_table"):
                raise OptionError("必须输入环境因子表格和分组表格")
        elif self.option("factor_type") == "env":
            if not self.option("env_table") :
                raise OptionError("必须输入环境因子表格")
        elif self.option("factor_type") == "group":
            if not self.option("group_table"):
                raise OptionError("必须输入分组表格")
        else:
            raise OptionError("因素类型格式不正确")

    def run_permanova(self):
        """
        运行permanova
        """
        env_table= self.combine_env()
        [binary, distance_method] = self.convert_dis_method(self.option('method'))
        options = {
            'dis_method': distance_method,
            'permutations': 999,
            'envtable': env_table,
            'binary': binary
        }
        options['abu_table'] = self.option('input_table').prop['path']
        self.permanova.set_options(options)
        self.permanova.on('end', self.set_output)
        self.permanova.run()

    def combine_env(self):
        if self.option("factor_type") =="both":
            env_info = {}
            env_title = []
            env_list =[]
            env_path = os.path.join(self.work_dir,"env_tem.xls")
            with open(self.option("env_table").prop['path'],"r") as env:
                line = env.readline()
                title = line.rstrip().split('\t')
                env_title.extend(title)
                while 1:
                    line = env.readline()
                    if not line:
                        break
                    fd = line.rstrip().split('\t')
                    for i in fd[1:]:
                        try:
                            float(i)
                        except:
                            self.set_error("环境因子表中有非数值型的环境因子")
                    env_info[fd[0]] = fd
                    env_list.append(fd[0])
            with open(self.option("group_table").prop['path'],"r") as group:
                line = group.readline()
                title = line.rstrip().split('\t')[1:]
                env_title.extend(title)
                group_num = {}
                group_list = {}
                for i in title:
                    group_num[i] = []
                    group_list[i] = []
                while 1:
                    line = group.readline()
                    if not line:
                        break
                    fd = line.rstrip().split('\t')
                    env_info[fd[0]].extend(fd[1:])
                    for i in range(len(title)):
                        group_list[title[i]].append(fd[i+1])
                        if fd[i+1] not in group_num[title[i]]:
                            group_num[title[i]].append(fd[i+1])
                for i in title:
                    if len(group_num[i]) < 2:
                        self.set_error("分组列表的组{}的分组数小于2".format(i))
                    for j in group_num[i]:
                        if group_list[i].count(j) < 3:
                            self.set_error("分组列表的组{}的{}组小于3".format(i,j))
            with open(env_path,"w") as new_env:
                new_env.write("\t".join(env_title))
                new_env.write('\n')
                for i in env_list:
                    new_env.write("\t".join(env_info[i]))
                    new_env.write('\n')
            return env_path
        elif self.option("factor_type") == "env":
            with open(self.option("env_table").prop['path'],"r") as env:
                env.readline()
                while 1:
                    line = env.readline()
                    if not line:
                        break
                    fd = line.rstrip().split('\t')
                    for i in fd[1:]:
                        try:
                            float(i)
                        except:
                            self.set_error("环境因子表中有非数值型的环境因子")
            return self.option("env_table").prop['path']
        elif self.option("factor_type") == "group":
            with open(self.option("group_table").prop['path'],"r") as group:
                line = group.readline()
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
                    fd = line.rstrip().split('\t')
                    for i in range(len(title)):
                        group_list[title[i]].append(fd[i+1])
                        if fd[i+1] not in group_num[title[i]]:
                            group_num[title[i]].append(fd[i+1])
                for i in title:
                    if len(group_num[i]) < 2:
                        self.set_error("分组列表的组{}的分组数小于2".format(i))
                    for j in group_num[i]:
                        if group_list[i].count(j) < 3:
                            self.set_error("分组列表的组{}的{}组小于3".format(i,j))
            return self.option("group_table").prop['path']

    def run(self):
        self.run_permanova()
        super(PermanovaMultiWorkflow, self).run()

    def _get_samplenames(self, groupfile):
        try:
            with open(groupfile, 'rb') as f:
                alllines = f.readlines()
                all_names = [i.split('\t')[0] for i in alllines]
            return all_names[1:]
        except IOError:
            self.set_error('无法打开分组文件或者文件不存在')

    def convert_dis_method(self, dis_method):
        """
        将距离算法进行转换
        :param dis_method:
        :return:
        """
        distance_name = ['euclidean', 'binary_euclidean', 'manhattan', 'binary_manhattan', 'gowerM', 'binary_gowerM',
                        'altGower', 'binary_altGower', 'canberraNZ', 'binary_canberraNZ', 'bray_curtis', 'binary_bray_curtis',
                        'kulczynski', 'binary_kulczynski', 'morisita_horn', 'binary_morisita_horn', 'morisita',
                        'binomial', 'binary_binomial', 'cao', 'binary_cao', 'chao', 'jaccard', 'binary_jaccard',
                        'raup_crick', 'mountford', 'mahalanobis']
        if dis_method not in distance_name:
            self.set_error('不支持所选的距离算法')
        splite = dis_method.split("_")
        if splite[0] == "binary":
            binary = "true"
            distance_method = dis_method[7:]
        else:
            binary = "false"
            distance_method = dis_method
        if distance_method == 'gowerM':
            distance_method = 'gower'
        elif distance_method == 'bray_curtis':
            distance_method = 'bray'
        elif distance_method == 'morisita_horn':
            distance_method = 'horn'
        elif distance_method == 'raup_crick':
            distance_method = 'raup'
        elif distance_method == 'canberraNZ':
            distance_method = 'canberra'
        return binary, distance_method

    def set_output(self):
        self.logger.info('strat set_output as {}'.format(
            self.__class__.__name__))
        # os.link(os.path.join(self.permanova.output_dir,"Multi-factor_PERMANOVA.xls")
        #         ,os.path.join(self.output_dir,"Multi-factor_PERMANOVA.xls"))
        with open(os.path.join(self.permanova.output_dir,"Multi-factor_PERMANOVA.xls"), "r") as inputFile,open(os.path.join(self.output_dir,"Multi-factor_PERMANOVA.xls"),"w") as output:
            while 1:
                line = inputFile.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                for i in range(len(fd)):
                    if fd[i] == "NA":
                        fd[i] = "-"
                output.write("\t".join(fd))
                output.write("\n")
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_permanmova = self.api.api("tool_lab.permanova_multi")
        api_permanmova.add_detail(self.option("main_id"),os.path.join(self.output_dir,"Multi-factor_PERMANOVA.xls"))
        self.logger.info("导表结束")
        self.end()
    
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(PermanovaMultiWorkflow, self).end()

    def filter_otu_sample(self, otu_path, filter_samples, newfile):
        """
        根据otu表和group表过滤掉样本
        :param otu_path: otu表
        :param filter_samples: 要过滤的样本list
        :param newfile: 生成的结果表
        :return:
        """
        if not isinstance(filter_samples, types.ListType):
            self.logger.error('过滤otu表样本的样本名称应为列表')
            self.set_error("filter_otu_sample错误")
        try:
            with open(otu_path, 'rb') as f, open(newfile, 'wb') as w:
                one_line = f.readline()
                all_samples = one_line.strip().split('\t')[1:]
                if not ((set(all_samples) & set(filter_samples)) == set(filter_samples)):
                    self.logger.error('提供的过滤样本存在otu表中不存在的样本all:%s,filter_samples:%s' % (all_samples, filter_samples))
                    self.set_error("otu表中不存在过滤样本")
                if len(all_samples) == len(filter_samples):
                    return otu_path
                samples_index = [all_samples.index(i) + 1 for i in filter_samples]
                w.write('#OTU\t' + '\t'.join(filter_samples) + '\n')
                for line in f:
                    all_values = line.rstrip().split('\t')
                    new_values = [all_values[0]] + [all_values[i] for i in samples_index]
                    w.write('\t'.join(new_values) + '\n')
                return newfile
        except IOError:
            self.set_error('无法打开OTU相关文件或者文件不存在')

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    data = {
        'name': 'test_permanova',
        'id': 'permanova_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "input_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_permanova/input_data.txt",
        "env_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_permanova/env.txt",
        "group_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_permanova/group.txt",
        "factor_type":"both",
        "method":"binary_euclidean",
        "main_id" : "5e9e6a6017b2bf2049a81be1"
        }
    }
    wsheet = Sheet(data=data)
    wf = PermanovaMultiWorkflow(wsheet)
    wf.run()