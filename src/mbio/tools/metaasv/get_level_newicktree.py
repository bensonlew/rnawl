# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file
from mbio.packages.metaasv.filter_newick import get_level_newicktree


class GetLevelNewicktreeAgent(Agent):
    """
    功能： 将不同类型的文件进行合并
    """

    def __init__(self, parent):
        super(GetLevelNewicktreeAgent, self).__init__(parent)
        options = [
            {"name": "asv_id", "type": "string"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "otu_file", "type": "string"},
            {"name": "temp_otu_file", "type": "outfile", "format": "meta.otu.otu_table"},
            {"name": "temp_tree_file", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("asv_id"):
            raise OptionError("请传入asv_id!")
        if not self.option("level"):
            raise OptionError("必须设置level!")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(GetLevelNewicktreeAgent, self).end()

class GetLevelNewicktreeTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(GetLevelNewicktreeTool, self).__init__(config)

    def get_newicktree(self):
        """
        获取进化树
        :return:
        """
        if self.option('level') != 9:
            newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),
                                              tempdir=self.work_dir, return_file=False, bind_obj=self)
            all_find = re.findall(r'\'.+?\'', newicktree)  # 找到所有带引号的进化树中复杂的名称
            for n, m in enumerate(all_find):
                all_find[n] = m.strip('\'')
            all_find = dict((i[1], i[0]) for i in enumerate(all_find))  # 用名称做键，找到的位置数字做值

            def match_newname(matchname):
                """
                随着自身被调用，自身的属性count随调用次数增加，返回OTU加次数，用于重命名进化树复杂的名称
                """
                if hasattr(match_newname, 'count'):
                    match_newname.count = match_newname.count + 1
                else:
                    match_newname.count = 1
                return 'ASV' + str(match_newname.count)  # 后面替换OTU中名称用同样的命名规则

            newline = re.sub(r'\'.+?\'', match_newname,
                             newicktree)  # 替换树种的复杂名称用 OTU 加数字代替 , 选哟注意的是这里的sub查找与findall查到方式是一致的
            temp_tree_file = self.work_dir + '/temp.tree'
            tempfile = open(temp_tree_file, 'w')
            tempfile.write(newline)
            tempfile.close()
            self.logger.info('get_newick:' + temp_tree_file)
            otu_table = self.option('otu_file')
            temp_otu_file = self.option('otu_file') + '.temp'
            all_lines = open(otu_table, 'r').readlines()
            if len(all_lines) < 3:
                self.logger.error('分类水平：%s,otu表数据少于2行：%s' % (self.option('level'), len(all_lines)))
                self.set_error("otu表数据少于2行")
            new_all = []
            new_all.append(all_lines[0])
            for line in all_lines[1:]:  # 遍历OTU表，将OTU表的复杂OTU名称改为之前find到的复杂名称对应的字典
                name = line.split('\t')
                if name[0] not in all_find:
                    self.set_error('ASV表中有原始表不存在的ASV名：%s', variables=(name[0]))
                name[0] = 'ASV' + str(all_find[name[0]] + 1)
                new_all.append('\t'.join(name))
            otu_file_temp = open(temp_otu_file, 'w')
            otu_file_temp.writelines(new_all)
            otu_file_temp.close()
            self.option("temp_otu_file", temp_otu_file)
            self.option("temp_tree_file", temp_tree_file)
        else:
            newicktree = get_level_newicktree(self.option('asv_id'), level=self.option('level'),
                                              tempdir=self.work_dir, return_file=False, bind_obj=self)
            temp_tree_file = self.work_dir + '/temp.tree'
            tempfile = open(temp_tree_file, 'w')
            tempfile.write(newicktree)
            tempfile.close()
            otu_table = self.option('otu_file')
            temp_otu_file = self.option('otu_file') + '.temp'
            all_lines = open(otu_table, 'r').readlines()
            new_all = []
            new_all.append(all_lines[0])
            for line in all_lines[1:]:  # OTU表中有复杂的名称OTU名称，包含进化物种类型，进化树种只有OTU名称
                name = line.split('\t')
                name[0] = name[0].split(';')[-1].strip()
                new_all.append('\t'.join(name))
            otu_file_temp = open(temp_otu_file, 'w')
            otu_file_temp.writelines(new_all)
            otu_file_temp.close()
            self.option("temp_otu_file", temp_otu_file)
            self.option("temp_tree_file", temp_tree_file)
        self.end()

    def run(self):
        super(GetLevelNewicktreeTool, self).run()
        self.get_newicktree()
