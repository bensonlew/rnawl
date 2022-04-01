# -*- coding: utf-8 -*-
# __author__ = 'JieYap'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import types
import subprocess
from biocluster.core.exceptions import OptionError


class OtunetworkAgent(Agent):
    """
    需要calc_otu_network.py
    version 1.0
    author: JieYao
    last_modified:2016.8.1
    """
    
    def __init__(self, parent):
        super(OtunetworkAgent, self).__init__(parent)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "string", "default": "otu"},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)
        self.step.add_steps('OtunetworkAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.OtunetworkAnalysis.start()
        self.step.update()
        
    def step_end(self):
        self.step.OtunetworkAnalysis.finish()
        self.step.update()
        
    def gettable(self):
        """
        根据输入的otu表和分类水平计算新的otu表
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            return self.option('otutable').get_table(self.option('level'))
        else:
            return self.option('otutable').prop['path']
        
    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('otutable').is_set:
            raise OptionError('必须提供otu表', code="32705101")
        self.option('otutable').get_info()
        if self.option('grouptable').is_set:
            self.option('grouptable').get_info()
        if self.option('otutable').prop['sample_num'] < 2:
            raise OptionError('otu表的样本数目少于2，不可进行网络分析', code="32705102")
        if self.option('grouptable').is_set:
            if self.option('grouptable').prop['sample_number'] < 2:
                raise OptionError('分组表的样本数目少于2，不可进行网络分析', code="32705103")
        samplelist = open(self.gettable()).readline().strip().split('\t')[1:]
        if self.option('grouptable').is_set:
            if len(self.option('grouptable').prop['sample']) > len(samplelist):
                raise OptionError('OTU表中的样本数量:%s少于分组表中的样本数量:%s', variables=(len(samplelist),
                                  len(self.option('grouptable').prop['sample'])), code="32705104")
            for sample in self.option('grouptable').prop['sample']:
                if sample not in samplelist:
                    raise OptionError('分组表的样本中存在OTU表中未知的样本%s', variables=(sample), code="32705105")
        table = open(self.gettable())
        if len(table.readlines()) < 4:
            raise OptionError('数据表信息少于3行', code="32705106")
        table.close()
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        self._memory = '15G'
        
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "网络分析结果输出目录"],
            ["./real_node_table.txt", "txt", "网络节点属性表"],
            ["./real_edge_table.txt", "txt", "网络边的属性表"],
            ["./real_dc_otu_degree.txt", "txt", "网络物种节点度分布表"],
            ["./real_dc_sample_degree.txt", "txt", "网络sample节点度分布表"],
            ["./real_dc_sample_otu_degree.txt", "txt", "网络所有节点度分布表"],
            ["./network_centrality.txt", "txt", "网络中心系数表"],
            ["./network_attributes.txt", "txt", "网络单值属性表"],
            ["./network_degree.txt", "txt", "网络度统计总表"]
        ])
        # print self.get_upload_files()
        super(OtunetworkAgent, self).end()


class OtunetworkTool(Tool):
    def __init__(self, config):
        super(OtunetworkTool, self).__init__(config)
        self._version = "1.0.1"
        self.cmd_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/calc_otu_network.py'
        if self.option('grouptable').is_set:
            self.group_table = self.option('grouptable').prop['path']
        self.otu_table = self.get_otu_table()
        self.out_files = ['real_node_table.txt', 'real_edge_table.txt', 'real_dc_otu_degree.txt',
                          'real_dc_sample_degree.txt', 'real_dc_sample_otu_degree.txt', 'network_centrality.txt',
                          'network_attributes.txt', 'network_degree.txt']

        
    def get_otu_table(self):
        """
        根据调用的level参数重构otu表
        :return:
        """
        if self.option('otutable').format == "meta.otu.tax_summary_dir":
            otu_path = self.option('otutable').get_table(self.option('level'))
        else:
            otu_path = self.option('otutable').prop['path']
        return otu_path

    def run(self):
        """
        运行
        """
        super(OtunetworkTool, self).run()
        self.run_otu_network_py()

    def formattable(self, tablepath):
        with open(tablepath) as table:
            if table.read(1) == '#':
                newtable = os.path.join(self.work_dir, 'temp_format.table')
                with open(newtable, 'w') as w:
                    w.write(table.read())
                return newtable
        return tablepath

    def run_otu_network_py(self):
        """
        运行calc_otu_network.py
        """
        real_otu_path = self.formattable(self.otu_table)
        cmd = self.config.SOFTWARE_DIR + '/program/Python/bin/python '
        cmd += self.cmd_path
        cmd += ' -i %s -o %s' % (real_otu_path, self.work_dir + '/otu_network')
        if self.option('grouptable').is_set:
            cmd += ' -m %s' % (self.group_table)
        print cmd
        self.logger.info('开始运行calc_otu_network生成OTU网络并进行计算')
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('OTU_Network计算完成')
        except subprocess.CalledProcessError:
            self.logger.info('OTU_Network计算失败')
            self.set_error('运行calc_otu_network.py失败', code="32705101")
        allfiles = self.get_filesname()
        for i in range(len(self.out_files)):
            self.linkfile(allfiles[i], self.out_files[i])
        self.end()

    def linkfile(self, oldfile, newname):
        """
        link文件到output文件夹
        :param oldfile: 资源文件路径
        :param newname: 新的文件名
        :return:
        """
        newpath = os.path.join(self.output_dir, newname)
        if os.path.exists(newpath):
            os.remove(newpath)
        os.link(oldfile, newpath)

    def get_filesname(self):
        filelist = os.listdir(self.work_dir + '/otu_network')
        files_status = [None, None, None, None, None, None, None, None]
        for paths, d, filelist in os.walk(self.work_dir + '/otu_network'):
            for filename in filelist:
                name = os.path.join(paths, filename)
                print name
                for i in range(len(self.out_files)):
                    if self.out_files[i] in name:
                        files_status[i] = name
        for i in range(len(self.out_files)):
            if not files_status[i]:
                self.set_error('未知原因，结果文件生成出错或丢失', code="32705102")
        return files_status
