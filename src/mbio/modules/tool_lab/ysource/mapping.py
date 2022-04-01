# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import random
import datetime
import subprocess
import re
import glob
import os
import sys
import shutil

class MappingModule(Module):
    '''
    Yfull流程mapping
    '''
    def __init__(self, work_id):
        super(MappingModule, self).__init__(work_id)
        options = [
            {'name':'bam_list','type':'infile','format':'denovo_rna_v2.bamlist'},
        ]
        self.add_option(options)
        # self.split_list()
        self.fq_file = {}
        self.samples = list()
        self.mkbam_tools = list()
        self.mapping_tools = list()
        self.mapping_option = {}
        self.mkchrybam_tools = list()
        self.mkchrybam_option = {}
        self.fastqstat_tools = list()

    def check_options(self):
        if not self.option('bam_list').is_set:
            raise OptionError("没有导入list文件，模块mapping启动失败")
        with open(self.option("bam_list").prop["path"],'r') as bl:
            while 1:
                line = bl.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                if not os.path.exists(fd[1]) and not os.path.exists(fd[2]):
                    raise OptionError("mapping需要的fq1和fq2文件缺失，请检查bam_list")
                self.fq_file[fd[0]] = [fd[1],fd[2]]
                self.samples.append(fd[0])
                self.logger.info(len(self.samples))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(MappingModule, self).run()
        self.run_tool()
        
    def run_tool(self):
        self.run_mkbam()
        self.logger.info(len(self.mkbam_tools))
        if len(self.mkbam_tools) == 1 :
            self.mkbam_tools[0].on('end',self.run_mapping)
        else:
            self.on_rely(self.mkbam_tools, self.run_mapping)
        for tool in self.mkbam_tools:
            tool.run()

    def run_mkbam(self):
        for i in self.samples:
            sample_name = i
            fq1_path = self.fq_file[i][0]
            fq2_path = self.fq_file[i][1]
            mkbam = self.add_tool('tool_lab.ysource.mkbam')
            options ={
                "fq1":fq1_path,
                "fq2":fq2_path,
                "sample_name":sample_name
            }
            self.mapping_option[i] = {
                "sample_name" : i,
                "bam_file": mkbam.output_dir + "/aln-pe.bam"
            }
            mkbam.set_options(options)
            self.mkbam_tools.append(mkbam)

    def run_mapping(self):
        for i in self.samples:
            index = self.samples.index(i)
            sample_name = i 
            mapping = self.add_tool('tool_lab.ysource.mapping')
            # options = {
            #     "sample_name" : i,
            #     "bam_file": self.mkbam_tools[index].output_dir + "/aln-pe.bam"
            # }
            options = self.mapping_option[i]
            mapping.set_options(options)
            self.mapping_tools.append(mapping)
            self.mkchrybam_option[i] =  {
                "sample_name" : i,
                "bam_file" : mapping.output_dir + "/align_final_sort.bam"
            }
        if len(self.mapping_tools) == 1:
            self.mapping_tools[0].on('end', self.run_mkchrybam)
        else:
            self.on_rely(self.mapping_tools,self.run_mkchrybam)
        for tool in self.mapping_tools:
                tool.run()            

    def run_mkchrybam(self):
        for i in self.samples:
            sample = i 
            index = self.samples.index(i)
            mkchrybam = self.add_tool('tool_lab.ysource.mkchrybam')
            # options = {
            #     "sample_name" : sample,
            #     "bam_file" : self.mapping_tools[index].output_dir + "/align_final_sort.bam"
            # }
            options = self.mkchrybam_option[i]
            mkchrybam.set_options(options)
            self.mkchrybam_tools.append(mkchrybam)
        if len(self.mkchrybam_tools) == 1 :
            self.mkchrybam_tools[0].on('end', self.run_fastqstat)
        else:
            self.on_rely(self.mkchrybam_tools, self.run_fastqstat)
        for tool in self.mkchrybam_tools:
            tool.run() 

    def run_fastqstat(self):
        for i in self.samples:
            sample = i 
            fastqstat = self.add_tool("tool_lab.ysource.fastqstat")
            input_dir = os.path.dirname(os.path.abspath(self.fq_file[sample][0]))
            sn_raw_input = os.path.join(input_dir,"{}_raw_input.txt".format(sample))
            sn_clean_input = os.path.join(input_dir,"{}_clean_input.txt".format(sample))
            options = {
                "raw_input": sn_raw_input,
                "clean_input": sn_clean_input,
                "sample_name": sample
            }
            fastqstat.set_options(options)
            self.fastqstat_tools.append(fastqstat)
        if len(self.fastqstat_tools) == 1 :
            self.fastqstat_tools[0].on('end',self.set_output)
        else:
            self.on_rely(self.fastqstat_tools, self.set_output)
        for tool in self.fastqstat_tools:
            tool.run()

    def set_output(self):
        self.logger.info('start set_output at {}'.format(
            self.__class__.__name__))
        with open(os.path.join(self.output_dir, "bam_list"), 'w') as bl:
            for tool in self.mapping_tools:
                sample_name = tool.option('sample_name')
                text1 = ""
                for source in glob.glob(os.path.join(tool.output_dir, '*')):
                    link_name = os.path.join(
                        self.output_dir, "{}_{}".format(sample_name,os.path.basename(source)))
                    if os.path.exists(link_name):
                        os.remove(link_name)
                    os.link(source, link_name)
                    if link_name[-3:] != "bai":
                        text1 = link_name
                    # sample_name = os.path.basename(source)[:-8]
                    # if os.path.basename(source)[-7] == "1":
                    #     text1 = link_name
                    # elif os.path.basename(source)[-7] == "2":
                    #     text2 = link_name

                    # text1 = text1 + "\t" + link_name
                    # bl.write(link_name)
                    # bl.write('\t')
                    self.logger.info(
                        'succeed in linking {} to {}'.format(source, link_name))
                text = '\t' + text1
                bl.write(sample_name)
                bl.write(text)
                bl.write('\n')
            self.logger.info('finish set_output at {}'.format(
                self.__class__.__name__))
        for tool in self.fastqstat_tools:
            sample_name = tool.option('sample_name')
            for source in glob.glob(os.path.join(tool.output_dir, '*')):
                link_name = os.path.join(
                    self.output_dir, os.path.basename(source)
                )
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
        source_size = {}
        now_date =  datetime.date.today().strftime("%Y%m%d")
        now_time = datetime.datetime.now().strftime("%H%M%S") + '_' + str(random.randint(1, 10000))
        for tool in self.mkchrybam_tools:
                sample_name = tool.option('sample_name')
                # text1 = ""
                for source in glob.glob(os.path.join(tool.output_dir, '*')):
                    link_name = os.path.join(
                        self.output_dir, os.path.basename(source))
                    if os.path.exists(link_name):
                        os.remove(link_name)
                    if sample_name[3] == "B":
                        if os.path.exists("/mnt/ilustre/users/sanger-dev/yoogene/kehu/{}".format(now_date)):
                            ybam_dir = "/mnt/ilustre/users/sanger-dev/yoogene/kehu/{}/{}_bam".format(now_date,now_time)
                            if not os.path.exists(ybam_dir):
                                os.mkdir(ybam_dir)
                            os.link(source, os.path.join(ybam_dir,os.path.basename(source)))
                        else:
                            ybam_dir = "/mnt/ilustre/users/sanger-dev/yoogene/kehu/{}/{}_bam".format(now_date,now_time)
                            os.mkdir("/mnt/ilustre/users/sanger-dev/yoogene/kehu/{}".format(now_date))
                            os.mkdir(ybam_dir)
                            os.link(source, os.path.join(ybam_dir,os.path.basename(source)))
                        os.link(source, link_name)
                    else:
                        if os.path.getsize(source) > 2097152:
                            if os.path.exists("/mnt/ilustre/users/sanger-dev/yoogene/fastq/{}".format(now_date)):
                                ybam_dir = "/mnt/ilustre/users/sanger-dev/yoogene/fastq/{}/{}_bam".format(now_date,now_time)
                                if not os.path.exists(ybam_dir):
                                    os.mkdir(ybam_dir)
                                os.link(source, os.path.join(ybam_dir,os.path.basename(source)))
                            else:
                                ybam_dir = "/mnt/ilustre/users/sanger-dev/yoogene/fastq/{}/{}_bam".format(now_date,now_time)
                                os.mkdir("/mnt/ilustre/users/sanger-dev/yoogene/fastq/{}".format(now_date))
                                os.mkdir(ybam_dir)
                                os.link(source, os.path.join(ybam_dir,os.path.basename(source)))
                            os.link(source, link_name)
                    source_size[os.path.basename(source)] = (os.path.getsize(source)/1024)/1024

                    # text1 = link_name
                    # sample_name = os.path.basename(source)[:-8]
                    # if os.path.basename(source)[-7] == "1":
                    #     text1 = link_name
                    # elif os.path.basename(source)[-7] == "2":
                    #     text2 = link_name
                    # text1 = text1 + "\t" + link_name
                    # bl.write(link_name)
                    # bl.write('\t')
                    self.logger.info(
                        'succeed in linking {} to {}'.format(source, link_name))
        with open(os.path.join(self.output_dir,"size.xls"),"w") as sz:
            for i in source_size.keys():
                sz.write("{}M\t{}\n".format(source_size[i],i))
        self.end()

    def end(self):
        super(MappingModule,self).end()

class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Mapping_' + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.ysource.mapping",
            "options": {
                # "bam_list": "/mnt/ilustre/users/sanger-dev/workspace/20201224/Single_Qc_7579/Qc/output/bam_list"
                "bam_list": "/mnt/ilustre/users/sanger-dev/workspace/20201227/Single_Qc_6186/Qc/output/bam_list"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
