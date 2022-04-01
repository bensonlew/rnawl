# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
# import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.graph.venn_table import venn_graph
from mbio.packages.graph.venn_table import venn_graph_mg
from mbio.packages.taxon.mask_taxon import mask_taxon
import re

class VennTableAgent(Agent):
    """
    metaasv分析 功能：花瓣图不再计算两两组之间的结果，Venn图正常呈图
    """
    def __init__(self, parent):
        super(VennTableAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir,denovo_rna.express.express_matrix,toolapps.table"},
            {"name": "group_table", "type": "infile", "format": " meta.otu.group_table, toolapps.group_table"},  # 输入的group表格
            {"name": "level", "type": "string", "default": "otu"},  # 物种水平
            {"name": "analysis_model", "type": "string", "default": ""}  # venn结果中物种/功能的连接符，因宏基因组的功能名称中含有逗号与默认连接符冲突，特设mg
        ]
        self.add_option(options)
        self.step.add_steps('venn_table')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.venn_table.start()
        self.step.update()

    def step_end(self):
        self.step.venn_table.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("otu_table"):
            raise OptionError("参数otu_table不能为空")
        if self.option("level") not in ['otu', 'domain', 'kindom', 'phylum', 'class',
                                        'order', 'family', 'genus', 'species']:
            raise OptionError("请选择正确的分类水平")
        if not self.option("group_table"):
            raise OptionError("参数group_table不能为空")
        if self.option("group_table").format == 'toolapps.group_table':
            for i in self.option('group_table').prop['sample_name']:
                if i not in self.option('otu_table').prop['col_sample']:
                    raise OptionError('分组文件中的样本不存在于表格中，查看是否是数据取值选择错误')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        """
        结束和上传结果文件
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "venn图结果目录"],
        ])
        super(VennTableAgent, self).end()


class VennTableTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(VennTableTool, self).__init__(config)
        self.R_path = '/program/R-3.3.1/bin/'
        self.venn_path = self.config.PACKAGE_DIR + '/graph/scripts/venn_table.py'
        self.petal_path =  self.config.PACKAGE_DIR + '/metaasv/venn_table.py'
        self.venn_path_mg = self.config.PACKAGE_DIR + '/graph/scripts/venn_table_mg.py'
        self.python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.python_path2 = '/program/Python/bin/'  # 用于小工具
        self.software = 'program/parafly-r2013-01-21/bin/bin/ParaFly'
        self._version = 1.0

    def _create_venn_table(self):
        """
        调用脚本venn_table.py,输出venn表格
        """
        if self.option("group_table").format == 'toolapps.group_table':  # add by wzy 2017.6.23
            group_file = self.option("group_table").prop['new_table']
        else:
            group_file = self.option("group_table").prop['path']
        if self.option("otu_table").format == 'toolapps.table':
            otu_table = self.option("otu_table").prop['new_table']
        elif self.option("otu_table").format is "meta.otu.tax_summary_dir":
            otu_table = self.option("otu_table").get_table(self.option("level"))
        else:
            otu_table = self.option("otu_table").prop['path']
        get_cmd_list = []
        cmd_list = []
        venn_path = ""
        self.name_to_name = mask_taxon(otu_table, self.work_dir + "/tmp_mask_otu.xls")  # 将分类名称替换成特定的name名,防止R计算出错 2017.11.07 by zhujuan
        otu_table = self.work_dir + '/tmp_mask_otu.xls'
        self.logger.info(">>>>>>>>>>>>>>>>>>>>>")

        ##分组方案的名称的list：
        group_scheme_list = self.option("group_table").prop['group_scheme']
        if len(group_scheme_list) == 1:   # 判断分组方案的个数
            group_number = self.option("group_table").group_num(group_scheme_list[0])
            ##判断分组个数，如果大于等于2小于等于5做venn图，大于等于6做花瓣图
            if self.option("analysis_model") == "":
                if group_number >= 2 and group_number <= 5:
                    venn_path = self.venn_path
                elif group_number >= 6:
                    venn_path = self.petal_path
            else:
                if group_number >= 2 and group_number <= 5:
                    venn_path = self.venn_path_mg
                elif group_number >= 6:
                    venn_path = self.petal_path
            self.logger.info(venn_path)
            os.link(group_file, self.work_dir + '/group_table')  # venn_table的结果与分组文件的目录一致，所以需要将分组文件放在工作目录下
            if self.option("otu_table").format == 'toolapps.table':
                self.option("otu_table").get_table_of_main_table(otu_table, self.work_dir + '/new_input.xls',group_file)
                venn_cmd = '%spython %s -i %s -g %s -o cmd.r' % (
                self.python_path, venn_path, self.work_dir + '/new_input.xls', self.work_dir + '/group_table')
            else:
                venn_cmd = '%spython %s -i %s -g %s -o cmd.r' % (
                self.python_path, venn_path, otu_table, self.work_dir + '/group_table')
            self.logger.info(venn_cmd)
            os.system(venn_cmd)
            self.logger.info('运行venn_cmd')
            cmd = self.R_path + 'Rscript cmd.r'
            self.logger.info("开始运行venn_table")
            command = self.add_command("venn_table", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行venn_table完成")
            else:
                self.set_error("运行venn_table运行出错!")
                raise Exception("运行venn_table运行出错，请检查输入的表格是否正确")
            # 统计各组所有otu/物种名 add by qindanhua
            if self.option("analysis_model") == "":
                venn_graph(otu_table, group_file, "venn_graph.xls")
            else:
                venn_graph_mg(otu_table, group_file, "venn_graph.xls")


    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171031
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        add func by guhaidong 20171031
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                #line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def remove_other_label(self,infile, graph_infile):
        """
        去掉之前共有的部分
        """
        ori_venn = infile
        group_ele = graph_infile
        group_map = {}
        with open(group_ele, 'r') as f:
            for line in f:
                spline = line.strip().split('\t')
                group_map[spline[0]] = set(spline[1].split(','))
        groups = set(group_map.keys())

        new_venn = self.work_dir + '/venn_table2.xls'
        with open(new_venn,'w') as fw, open(ori_venn, 'r') as fr:
            for line in fr:
                spline = line.strip().split('\t')
                if '&' in spline[0]:
                    cg = [i.strip() for i in spline[0].split('&')]
                    cg = set([i.split(" only")[0] for i in cg])
                    other_g = groups - cg
                    if len(spline) == 3:
                        cur_ele = set(spline[2].split(','))
                    else:
                        cur_ele = set([])
                    for og in other_g:
                        cur_ele = cur_ele - group_map[og]
                    # s_cur_ele = sorted(cur_ele,key=lambda x: int(x.split('_')[1]))
                    s_cur_ele = sorted(cur_ele)
                    fw.write(spline[0]+'\t'+str(len(s_cur_ele))+'\t'+','.join(s_cur_ele)+'\n')
                else:
                    fw.write(line)

    def set_output(self):
        """
        将结果文件链接至output
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        if len(self.option("group_table").prop['group_scheme']) == 1:
            if os.path.exists(self.work_dir + "/venn_graph.xls") and self.option("analysis_model") == "":  # 去除多样性物种名前的多余水平(如d__Archaea; k__norank; p__Thaumarchaeota) 20180112
                self.remove_other_label(self.work_dir + '/venn_table.xls', self.work_dir + '/venn_graph.xls')
                venn_graph(self.option("otu_table").prop['path'], self.option("group_table").prop['path'], "venn_graph.xls")
                self.add_taxon(self.work_dir + '/venn_table2.xls', self.work_dir + '/venn_table_new.xls')  # add_taxon换回原名称
                os.link(self.work_dir + '/venn_table_new.xls', self.output_dir + '/venn_table.xls')
                os.link(self.work_dir + '/venn_graph.xls', self.output_dir + '/venn_graph.xls')
            elif os.path.exists(self.work_dir + "/venn_graph.xls") and self.option("analysis_model") != "":
                self.add_taxon(self.work_dir + '/venn_graph.xls', self.work_dir + '/venn_graph_new.xls')
                self.remove_other_label(self.work_dir + '/venn_table.xls', self.work_dir + '/venn_graph.xls')
                self.add_taxon(self.work_dir + '/venn_table2.xls', self.work_dir + '/venn_table_new.xls')  # add_taxon换回原名称
                os.link(self.work_dir + '/venn_table_new.xls', self.output_dir + '/venn_table.xls')
                os.link(self.work_dir + '/venn_graph_new.xls', self.output_dir + '/venn_graph.xls')
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(VennTableTool, self).run()
        self.logger.info("开始运行！")
        self._create_venn_table()
        self.logger.info("set out put")
        self.set_output()
        self.end()
