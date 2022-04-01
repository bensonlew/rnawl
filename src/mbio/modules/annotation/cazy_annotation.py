# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last_modify:2017.09.29
# last_modify by: linmeng.liu 20200303 add best/add_score option

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CazyAnnotationModule(Module):
    def __init__(self, work_id):
        super(CazyAnnotationModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集的输出
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            # gene_profile.reads_number.txt
            {"name": "cazy_family_profile", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "cazy_class_profile", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "best", "type": "bool", "default": True},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "add_score", "type": "string", "default": "True"},  #解析结果是否增加score和identity
        ]
        self.add_option(options)
        self.hmmscan = self.add_module("align.hmmscan")  # 使用统一的质控模块，用参数控制使用该模块
        self.anno = self.add_tool("annotation.cazy_anno")  # 质控前后序列信息统计
        self.step.add_steps('hmmscan', 'anno')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21200501")
        return True

    def run_hmmscan(self):
        os.system("sed -n '$=' {} > {}".format(self.option('query').prop['path'], self.work_dir + "/list.txt"))
        with open(self.work_dir + "/list.txt", 'r') as l:
            line = l.readline().strip('\n')
        self.logger.info(line)
        n = int(line)/20
        if n >= 400000:
            n = 400000
        self.hmmscan.set_options({
            'query': self.option('query'),
            'lines': n,
            "database": "cazy_v8",## add by qingchen.zhang@20200922
        })
        self.hmmscan.on('start', self.set_step, {'start': self.step.hmmscan})
        self.hmmscan.on('end', self.set_step, {'end': self.step.hmmscan})
        self.hmmscan.on('end', self.set_step, {'start': self.step.anno})
        self.hmmscan.on('end', self.set_output, "hmmscan")
        self.hmmscan.on('end', self.run_anno)
        self.hmmscan.run()

    def run_anno(self):
        self.anno.set_options({
            'hmmscan_result': self.hmmscan.option('align_result'),
            'reads_profile_table': self.option('reads_profile_table'),
            'best': self.option('best'),
            'add_score': self.option('add_score'),
            'evalue': self.option('evalue'),
            'version': "cazy_v8"## add by qingchen.zhang@20200922
        })
        self.anno.on('start', self.set_step, {'start': self.step.anno})
        self.anno.on('end', self.set_step, {'end': self.step.anno})
        self.anno.on('end', self.set_output, "anno")
        self.anno.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'hmmscan':
            self.linkdir(obj.output_dir, "hmmscan")
        if event['data'] == 'anno':
            self.linkdir(obj.output_dir, 'anno_result')
            self.option('cazy_family_profile', obj.output_dir + "/cazy_family_profile.xls")
            self.option('cazy_class_profile', obj.output_dir + "/cazy_class_profile.xls")
            self.end()

    def run(self):
        super(CazyAnnotationModule, self).run()
        self.run_hmmscan()

    def end(self):
        super(CazyAnnotationModule, self).end()

    def linkdir(self, dirpath, dirname):  # 暂时无用
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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))
