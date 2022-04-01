# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
from biocluster.core.exceptions import OptionError
import subprocess


class CatTableAgent(Agent):
    """
    cat多个文件为一个文件
    1. 如果文件存在header，先删掉header，如果没有则直接进行合并，最后合并再加上header
    2. 可以合并一个文件夹下相同字段的文件通过table_name进行判断

    """

    def __init__(self, parent):
        super(CatTableAgent, self).__init__(parent)
        options = [
            {"name": "merge_dir", "type": "infile", "format": "paternity_test.data_dir"},  # 输入合并文件夹
            {"name": "prefix", "type": "string", "default": "table"}, #合并后的文件前缀名
            {"name": "header", "type": "string", "default": 'False'}, #默认为没有header
            {"name": "table_name", "type": "string"}, #合并前的文件的名称或者部分字段
        ]
        self.add_option(options)
        self.step.add_steps('cat_table')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cat_table.start()
        self.step.update()

    def step_end(self):
        self.step.cat_table.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("merge_dir").is_set:
            raise OptionError("请传入merge_dir文件夹路径")
        if not self.option('table_name'):
            raise OptionError("请输入需要合并的文件名称或部分字段")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(CatTableAgent, self).end()

class CatTableTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CatTableTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'

    def cat_table(self):
        """
        合并文件
        :return:
        """
        file_list = os.listdir(self.option('merge_dir').prop['path'])
        header = self.work_dir + '/head.txt'
        cmd = self.sh_path + 'cat_seq.sh'
        if self.option('header') == "True":
            for file in file_list:
                if re.search(r'%s'%(self.option('table_name')), file):
                    file_path = self.option('merge_dir').prop['path'] + '/' + file
                    os.system('head -n 1 %s > %s'%(file_path, header))
                    break
        for file in file_list:
            if re.search(r'%s'%(self.option('table_name')), file):
                file_path = self.option('merge_dir').prop['path'] + '/' + file
                if self.option('header') == "True":
                    os.system("sed -i '1d' %s"%(file_path))
                cmd += ' ' + file_path
        cmd += ' ' + self.work_dir + '/{}'.format(self.option("prefix") + self.option('table_name'))
        self.logger.info('运行cat_seq，将table进行合并')
        self.logger.info(cmd)
        command = self.add_command("cat_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            pre_path = self.work_dir + '/{}'.format(self.option("prefix") + self.option('table_name'))
            if self.option('header') == "True":
                with open(header, 'r') as h:
                    head = h.readline()
                with open(pre_path, 'r+') as f:
                    content = f.read()
                    f.seek(0, 0)
                    f.write('{}{}'.format(head, content))
            self.logger.info("合并文件完成啦")
        else:
            self.set_error("合并文件失败！")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if os.path.exists(self.output_dir + '/{}'.format(self.option("prefix") + self.option('table_name'))):
            os.remove(self.output_dir + '/{}'.format(self.option("prefix") + self.option('table_name')))
        os.link(self.work_dir + '/{}'.format(self.option("prefix")+ self.option('table_name')),
                self.output_dir + '/{}'.format(self.option("prefix")+ self.option('table_name')))

    def run(self):
        super(CatTableTool, self).run()
        self.cat_table()
        self.set_output()
        self.end()