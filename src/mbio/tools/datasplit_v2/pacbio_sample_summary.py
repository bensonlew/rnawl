# -*_ coding: utf-8 -*-
# __author__ = 'XueQinwen'
# last_modified: 20210907
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os
import datetime
import random
import unittest
import time
import gzip

class PacbioSampleSummaryAgent(Agent):
    def __init__(self,parent):
        super(PacbioSampleSummaryAgent,self).__init__(parent)
        options = [
            {"name":"input_bam","type":"infile","format":"align.bwa.bam"},
            {"name":"raw_fastq","type":"string"},
            {"name":"project_type","type":"string"},
            {"name":"clean_fastq","type":"string"},
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_option(self):
        """
        参数检查
        """
        if not self.option("input_bam"):
            raise OptionError("没有输入bam文件")
        if not self.option("project_type"):
            raise OptionError("没有输入项目类型")
        elif self.option("project_type") == "diversity" and not self.option("raw_fastq"):
            raise OptionError("多样性项目需要输入fastq文件")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 8
        self._memory = "20G"

    def end(self):
        super(PacbioSampleSummaryAgent,self).end()

class PacbioSampleSummaryTool(Tool):
    def __init__(self,config):
        super(PacbioSampleSummaryTool,self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/library/glibc-2.14/lib")
        self.samtools = self.config.SOFTWARE_DIR+"/bioinfo/align/samtools-1.8/samtools"

    def samtool_stat(self, inputbam):
        output = os.path.join(self.work_dir, "stat.txt")
        cmd = "{} stats {} > {}".format(self.samtools,inputbam, output)
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/datasplit/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "samstat_{}.sh".format(now_time)
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/datasplit/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始统计")
        self.logger.info(shell)
        command = self.add_command("stat",shell).run()
        self.wait(command)
        total_num = 0
        avg_num = 0
        max_num = 0
        if command.return_code == 0 :
            self.logger.info("运行samtools stat统计完成")
        else:
            self.set_error("统计出现错误")
        with open(output,'r') as rs:
            while 1:
                line = rs.readline()
                if not line:
                    break
                if line[:3] == "FFQ":
                    break
                if line[3:17] == "average length":
                    avg_num = int(line.split(":")[1])
                if line[3:17] == "maximum length":
                    max_num = int(line.split(":")[1])
                if line[3:15] == "total length":
                    total_num = int(line.split(":")[1].split('\t')[1])
        return total_num,avg_num,max_num

    def getFastqReadsNum(self, input_fastq):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        self.readfq = os.path.join( self.config.PACKAGE_DIR,"datasplit/readfq")
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/datasplit/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "fastqstat_{}.sh".format(now_time)
        stat_file = os.path.join(self.work_dir,"{}.txt".format(os.path.basename(input_fastq)))
        cmd = "{} {} >{}".format(self.readfq,input_fastq,stat_file)
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/datasplit/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始统计{}的reads数".format(os.path.basename(input_fastq)))
        self.logger.info(shell)
        command1 = self.add_command("stat{}".format(os.path.basename(input_fastq).split('.')[0].lower()), shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行统计完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        reads_num = 0 
        st_file = open(stat_file,"r")
        if st_file.read() == '\n':
            self.logger.info("结果没有获取到，等待1分钟后重新获取")
            time.sleep(60)
            command2 = self.add_command("stat{}1".format(os.path.basename(input_fastq).split('.')[0].lower()), shell).run()
            self.wait(command2)
            if command2.return_code == 0:
                self.logger.info("运行统计完成")
            else:
                self.set_error("运行错误，请重新检查输入参数")
        with open(stat_file,'r') as sf:
            reads_num = int(sf.readline().split('\t')[0].split(':')[1])
        # os.system('rm {}'.format(file_path))
        return reads_num
    
    # def getFastqReadsNum(self, input_fastq):
    #     read_num = 0 
    #     # fastq_file = 
    #     self.logger.info(input_fastq)
        
    #     return read_num

    def run(self):
        super(PacbioSampleSummaryTool,self).run()
        if self.option('project_type') == "diversity":
            self.run_diversity_statitic()
        else:
            self.run_non_diversity_statistic()
        self.end()
           
    def run_diversity_statitic(self):
        sample_name = os.path.basename(self.option('clean_fastq')).split('.')[0]
        total_len,ave_len,_ = self.samtool_stat(self.option("input_bam").prop['path'])
        reads = self.getFastqReadsNum(self.option('raw_fastq'))
        qc_reads= self.getFastqReadsNum(self.option('clean_fastq'))
        reads_base = 0
        self.output_file = os.path.join(self.output_dir,"{}.txt".format(sample_name))
        with open(self.output_file,"w") as of:
            of.write('total_len\tave_len\treads\tqc_reads\treads_base\tatype\n')
            of.write('{}\t{}\t{}\t{}\t{}\t{}'.format(total_len,ave_len,reads,qc_reads,reads_base,"diversity"))

    def run_non_diversity_statistic(self):
        sample_name = os.path.basename(self.option('input_bam').prop['path']).split('.')[1]
        total_len,ave_len,_ = self.samtool_stat(self.option("input_bam").prop['path'])
        reads = 0
        qc_reads = 0
        reads_base = self.friendly_size(total_len)
        self.output_file = os.path.join(self.output_dir,"{}.txt".format(sample_name))
        with open(self.output_file,"w") as of:
            of.write('total_len\tave_len\treads\tqc_reads\treads_base\tatype\n')
            of.write('{}\t{}\t{}\t{}\t{}\t{}'.format(total_len,ave_len,reads,qc_reads,reads_base,"non_diversity"))

    def friendly_size(self, size):
        gb = 1000 * 1000* 1000.0
        mb = 1000 * 1000.0
        kb = 1000.0
        if size > gb:
            new_size = round(float(size) / gb, 4)
            return str(new_size) + "G"
        else:
            new_size = round(float(size) / mb,4)
            return str(new_size) + "M"

            
class TestFunctionf(unittest.TestCase):
    '''
    测试脚本
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "pacsam_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "datasplit_v2.pacbio_sample_summary",
            "options": {
                "input_bam":"/mnt/ilustre/users/sanger-dev/workspace/20210915/Single_lima_9691/Lima/output/split/split.test11111.bam",
                "raw_fastq":"/mnt/ilustre/users/sanger-dev/workspace/20210918/PacbioSplit_SP20210915-1631687156_20210918_173310/PacbioQcStat/PacbioQc1/output/EFH061965.ccs.fastq.gz",
                "clean_fastq":"/mnt/ilustre/users/sanger-dev/workspace/20210918/PacbioSplit_SP20210915-1631687156_20210918_173310/PacbioQcStat/PacbioQc1/output/EFH061965_value.fastq.gz",
                "project_type":"diversity"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

    