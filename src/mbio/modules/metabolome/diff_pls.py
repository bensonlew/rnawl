# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class DiffPlsModule(Module):
    def __init__(self, work_id):
        super(DiffPlsModule, self).__init__(work_id)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'metabolome.express,metabolome.metab_abun'},  # 表达矩阵文件
            {"name": "metab_desc", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {'name': 'group_name', 'type': 'string', 'default': ''},  # 差异组名，eg："A|B";"A|C" ,两两比较
            {'name': 'mul_type', 'type': 'string', 'default': 'pca;plsda;oplsda'},  # 多元统计类型，pca，plsda, oplsda
            {'name': 'confidence', 'type': 'string', 'default': '0.95;0.95;0.95'},  # 置信度，与mul_type对应
            {'name': 'perm', 'type': 'string', 'default': ';200;200'},  # 置换次数，与mul_type对应
            {'name': 'data_trans', 'type': 'string', 'default': 'UV;Par;Par'},  # 数据转化方法："UV","Ctr","Par"，"", 与mul_type对应个数
            {'name': 'test_method', 'type': 'string', 'default': 't-test'},  # 差异检验方法
            {'name': 'side_type', 'type': 'string', 'default': 'two-tailed'},  # 单尾或双尾检验 two-tailed,left-tailed,right-tailed
            {'name': 'creat_metabset', 'type': 'bool', 'default': False},  # 是否创建基因集
            {'name': 'metabset_dir', 'type': 'outfile', 'format': 'annotation.mg_anno_dir'},  # 
            {'name': 'pls_dir', 'type': 'outfile', 'format': 'annotation.mg_anno_dir'}, #
            {'name': 'diffStat_dir', 'type': 'outfile', 'format': 'annotation.mg_anno_dir'}, #
            {'name': 'id_diffStat_dir', 'type': 'outfile', 'format': 'annotation.mg_anno_dir'}, #
            {'name': 'diff_plot_dir', 'type': 'outfile', 'format': 'annotation.mg_anno_dir'},
            {'name': 'filter_k', 'type':'string','default':''},
            {"name": "filter_t","type":"string"},
            {"name": "filter_v", "type": "string"}
        ]
        self.add_option(options)
        self.test_tool = self.add_tool("metabolome.diff.diff_test")
        self.pls_tool = self.add_tool("metabolome.diff.diff_mul_stat")
        self.creat_metab_tool = self.add_tool("metabolome.diff.merge_diff")

    def check_options(self):
        if not self.option("exp_file").is_set:
            raise OptionError("请传入丰度矩阵！", code="24700201")
        if not self.option("group_file").is_set:
            raise OptionError("请传入group文件！", code="24700202")
        return True

    def set_output(self):
        '''
        '''
        try:
            DiffStat_dir = os.path.join(self.creat_metab_tool.output_dir, "DiffStat")
            id_diffStat_dir = os.path.join(self.creat_metab_tool.output_dir, "tmp_DiffStat")
            pls_dir = self.pls_tool.output_dir
            if self.option("creat_metabset"):
                metabset_dir = os.path.join(self.creat_metab_tool.output_dir, "Metabset")
                self.option("metabset_dir", metabset_dir)
            self.option("diffStat_dir", DiffStat_dir)
            self.option("pls_dir", pls_dir)
            self.option("id_diffStat_dir", id_diffStat_dir)
            self.option('diff_plot_dir',self.test_tool.work_dir + '/plot_dir' )
        except Exception as e:
            self.set_error("输出结果文件异常——%s", variables=(e), code="24700201")
        self.end()

    def run(self):
        super(DiffPlsModule, self).run()
        self.sub_groups = self.option("group_name").split(";")
        self.logger.info(self.option("group_file"))
        self.group_file = self.option("group_file")
        self.exp_file = self.option("exp_file")
        self.run_test()
        self.run_pls()
        self.diff_tools = [self.test_tool, self.pls_tool]
        self.on_rely(self.diff_tools, self.run_creat_metab)
        self.creat_metab_tool.on("end", self.set_output)

    def run_test(self):
        """
        两两差异检验
        """
        self.logger.info("run diff_pls module t-test")
        self.logger.info(self.group_file)
        options = {
            "exp_file": self.exp_file,
            "test_method": self.option("test_method"),
            "group_file": self.group_file,
            "group_name": self.option("group_name"),
            "side_type": self.option("side_type"),
            "mul_test" : "fdr"
        }
        self.test_tool.set_options(options)
        self.test_tool.run()

    def run_pls(self):
        """
        pca、plsda、oplsda分析
        """
        self.logger.info("module pls analysis start")
        self.logger.info(self.exp_file.prop["path"])
        self.logger.info(self.group_file.prop["path"])
        self.logger.info(self.option("mul_type"))
        options = {
            "exp_file": self.exp_file,
            "metab_desc": self.option("metab_desc"),
            "mul_type": self.option("mul_type"),
            "group_file": self.group_file,
            "group_name": self.option("group_name"),
            "confidence": self.option("confidence"),
            "perm": self.option("perm"),
            "data_trans": self.option("data_trans"),
        }
        self.pls_tool.set_options(options)
        self.pls_tool.run()

    def run_creat_metab(self):
        if self.option("creat_metabset"):
            metabset = True
        else:
            metabset = False
        metab_trans = self.option("metab_desc").prop["path"]
        options = {
            "test_dir": self.test_tool.output_dir,
            "pls_dir": self.pls_tool.output_dir,
            "sub_group": self.option("group_name"),
            "metabset": metabset,
            "metab_trans": metab_trans
        }
        if self.option('filter_k'):
            options['filter_k'] = self.option('filter_k')
            options['filter_t'] = self.option('filter_t')
            options['filter_v'] = self.option('filter_v')
        self.creat_metab_tool.set_options(options)
        self.creat_metab_tool.run()

    def end(self):
        super(DiffPlsModule, self).end()

    def mk_dir(self,dir_names):
        for each in dir_names:
            if os.path.exists(each):
                pass
            else:
                os.mkdir(each)
