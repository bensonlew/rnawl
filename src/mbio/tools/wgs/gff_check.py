# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last modify 20190603
import re
import json
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class GffCheckAgent(Agent):
    """
    用于将基因组配置中对gff进行检查，因为gff文件比较大，所以将检查放到tool中进行
    """
    def __init__(self, parent):
        super(GffCheckAgent, self).__init__(parent)
        options = [
            {"name": "gff_file", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps('GOanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.GOanno.start()
        self.step.update()

    def step_end(self):
        self.step.GOanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("gff_file"):
            raise OptionError("请设置gff_file参数")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(GffCheckAgent, self).end()


class GffCheckTool(Tool):
    def __init__(self, config):
        super(GffCheckTool, self).__init__(config)

    def check_gff(self):
        """
        检查gff的最后一列注释信息中，是否含有冒号和空格，如果有，用下划线替换。暂时只是检查是否有问题，如果有问题，反馈给产品线
        :return:''
        """
        self.logger.info('refgff check is start!')
        with open(self.option("gff_file"), 'r') as r, open(self.output_dir + "/ref.gff", 'w') as w:
            for line in r:
                if re.match('#.*', line):
                    pass
                else:
                    new_anno = []
                    temp = line.strip().split('\t')
                    for m in temp[-1].split(';'):
                        if re.match(r'ID=(.*) (.*)', m) or re.match(r'Parent=(.*) (.*)', m) or \
                                re.match(r'Name=(.*) (.*)', m):
                            temp_ = '_'.join(re.split(r'\s+', m))
                            if re.match('(.*):(.*)', temp_):
                                new_anno.append('_'.join(temp_.split(':')))
                            else:
                                new_anno.append(temp_)
                        elif re.match(r'ID=(.*):(.*)', m) or re.match(r'Parent=(.*):(.*)', m) or \
                                re.match(r'Name=(.*):(.*)', m):
                            new_anno.append('_'.join(m.split(':')))
                        else:
                            new_anno.append(m)
                    w.write('{}\t{}\n'.format('\t'.join(temp[:-1]), ';'.join(new_anno)))
        self.logger.info('refgff check is end!')

    def run(self):
        super(GffCheckTool, self).run()
        self.check_gff()
        self.end()
