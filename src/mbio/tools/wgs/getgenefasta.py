# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180916

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class GetgenefastaAgent(Agent):
    """
    提取ref.fa里面的gene.fa文件
    lasted modified by hongdong@20190315 添加对geneid的长度过滤
    """
    def __init__(self, parent):
        super(GetgenefastaAgent, self).__init__(parent)
        options = [
            {"name": "reffa", "type": "string"},  # 更名后的reffa
            {"name": "refgff", "type": "string"},  # 更名后的refgff
        ]
        self.add_option(options)
        self.step.add_steps('Getgenefasta')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Getgenefasta.start()
        self.step.update()

    def step_end(self):
        self.step.Getgenefasta.finish()
        self.step.update()

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置reffa参数") 
        if not self.option("refgff"):
            raise OptionError("请设置refgff参数")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '50G'

    def end(self):
        super(GetgenefastaAgent, self).end()


class GetgenefastaTool(Tool):
    def __init__(self, config):
        super(GetgenefastaTool, self).__init__(config)
        self.gene_path = self.config.PACKAGE_DIR + "/wgs/getGeneFasta.pl"
        self.perl_path = 'miniconda2/bin/perl '

    def getgenefasta(self):
        """
        要重新写下！！！
        :return:
        """
        cmd = "{}{} -i {} -g {} -o {}"\
            .format(self.perl_path, self.gene_path, self.option("reffa"), self.option("refgff"),
                    self.work_dir + "/ref.gene.fa")
        self.logger.info(cmd)
        self.logger.info("开始进行Getgenefasta")
        command = self.add_command("getgenefasta", cmd).run() 
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Getgenefasta完成！")
        else:
            self.set_error("Getgenefasta出错！")

    def filter_id(self):
        """
        '>rna4552:NM_001182509.1|855026|NM_001182509.1;NR_132237.1;NR_132238.1;NR_132239.1;NR_132240.1;NR_132241.1
        ;NR_132242.1;NR_132243.1|NM_001182509.1;NR_132237.1;NR_132238.1;NR_132239.1;NR_132240.1:chr13:295179:296738'
        当id的长度特别长的时候会导致interpro运行失败， 所以在对注释的序列id进行过滤处理
        ERROR - StepExecution with errors - stepName: stepLoadOrfFromFasta
        :return:
        """
        with open(self.work_dir + "/ref.gene.fa", 'r') as r, open(self.output_dir + "/ref.gene.fa", 'w') as w:
            for line in r:
                if re.match('>:.*--:--:$', line) or line == '\n':
                    continue
                if re.match('>.*', line):
                    if len(line) >= 212:
                        temp = line.split(':')
                        new_id = ':'.join([temp[0], self.get_shorts_name(temp[1]), ':'.join(temp[2:])])
                        w.write('{}'.format(new_id))
                        self.logger.info('异常基因id：{}, 改为：{}'.format(line, new_id))
                    else:
                        w.write(line)
                else:
                    w.write(line)

    def get_shorts_name(self, ids):
        if ids.split('|'):
            new_id = ids.split('|')[0]
            for id_ in ids.split('|'):
                if len(id_) < len(new_id):
                    new_id = id_
        else:
            new_id = 'None'
        return new_id

    def run(self):
        super(GetgenefastaTool, self).run()
        self.getgenefasta()
        self.filter_id()
        self.end()
