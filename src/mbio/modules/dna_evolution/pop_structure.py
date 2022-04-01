# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last_modify:20180910

import os
import re
import gevent
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class PopStructureModule(Module):
    """
    structure 结果分析, 包含了admistructure 与cverror
    """
    def __init__(self, work_id):
        super(PopStructureModule, self).__init__(work_id)
        options = [
            {"name": "pop_fam", "type": "infile", "format": "dna_evolution.bed", "required": True},
            {"name": "pop_bed", "type": "infile", "format": "dna_evolution.bed", "required": True},
            {"name": "k_min", "type": "int", "default": 2},
            {"name": "k_max", "type": "int", "default": 20}
        ]
        self.add_option(options)
        self.cverror = self.add_tool("dna_evolution.cverror")
        self.structure_tools = []
        self.sample_num = 0
        self.k_max_ = 21

    def check_options(self):
        if type(self.option("k_min")) != int:
            raise OptionError("k_value一定要是整型数据", code="123345")
        else:
            if self.option("k_min") > 100 or self.option("k_min") < 1:
                raise OptionError("k值要大于%s且小于%s", variables=(1, 100), code='12345')
        if type(self.option("k_max")) != int:
            raise OptionError("k_value一定要是整型数据", code="123345")
        else:
            if self.option("k_max") > 100 or self.option("k_max") < 1:
                raise OptionError("k值要大于%s且小于%s", variables=(1, 100), code='12345')

    def admi_structure_run(self):
        self.set_k_max()
        for i in range(self.option("k_min"), self.k_max_):
            admi_structure = self.add_tool("dna_evolution.admi_structure")
            admi_structure.set_options({
                "pop_bed": self.option("pop_bed").prop['path'],
                "k_value": i,
                "sample_list": self.work_dir + "/sample_list"
            })
            self.structure_tools.append(admi_structure)
        for j in range(len(self.structure_tools)):
            self.structure_tools[j].on('end', self.set_output, 'structure')
        if len(self.structure_tools) > 1:
            self.on_rely(self.structure_tools, self.set_structure)
        elif len(self.structure_tools) == 1:
            self.structure_tools[0].on('end', self.set_structure)
        else:
            self.set_error("structure_tools列表为空！")
        for t in self.structure_tools:
            gevent.sleep(1)
            t.run()

    def set_structure(self):
        """
        设置structure列表
        2\tpop.2.log\tpop.2.xls
        :return:
        """
        if os.path.exists(os.path.join(self.output_dir, "structure.list")):
            os.remove(os.path.join(self.output_dir, "structure.list"))
        if len(self.structure_tools) > 1:
            path_ = self.output_dir + "/structure"
        else:
            path_ = self.structure_tools[0].output_dir
        with open(os.path.join(self.output_dir, "structure.list"), 'w') as w:
            for i in range(self.option("k_min"), self.option('k_max') + 1):
                log = os.path.join(path_, "pop.{}.log".format(i))
                xls = os.path.join(path_, "pop.{}.xls".format(i))
                self.is_exists(log)
                self.is_exists(xls)
                w.write("{}\t{}\t{}\n".format(i, log, xls))
        self.cverror_run()

    def is_exists(self, file_path):
        if not os.path.exists(file_path):
            self.set_error("文件{}不存在！".format(file_path))

    def cverror_run(self):
        self.cverror.set_options({
            "structure_list": os.path.join(self.output_dir, "structure.list")
        })
        self.cverror.on("end", self.set_output, "cverror")
        self.cverror.on("end", self.end)
        self.cverror.run()

    def get_sample_list(self):
        """
        获取样本列表
        cut -f 1 -d " " pop.fam  > sample.list
        :return:
        """
        cmd = "cut -f 1 -d ' ' {} > {}".format(self.option('pop_fam').prop['path'], self.work_dir + "/sample_list")
        code = os.system(cmd)
        if code == 0:
            self.logger.info("命令{}执行成功！".format(cmd))
        else:
            self.set_error("命令{}执行失败！".format(cmd))

    def get_file_len(self, file_path):
        """
        获取文件行数
        :return:
        """
        return len(open(file_path, 'rU').readlines())

    def set_k_max(self):
        if self.sample_num < self.k_max_:
            self.k_max_ = self.sample_num + 1

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'structure':
            self.linkdir(obj.output_dir, 'structure')
        elif event['data'] == 'cverror':
            self.linkdir(obj.output_dir, 'cverror')
        else:
            pass

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(PopStructureModule, self).run()
        self.get_sample_list()
        self.sample_num = self.get_file_len(self.work_dir + "/sample_list")
        self.logger.info("样本个数:{}".format(self.sample_num))
        self.admi_structure_run()

    def end(self):
        super(PopStructureModule, self).end()
