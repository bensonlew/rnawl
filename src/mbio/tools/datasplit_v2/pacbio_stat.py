# -*_ coding: utf-8 -*-
# __author__ = 'XueQinwen'
# last_modified: 20210907
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os
import datetime
import unittest
import random

class PacbioStatAgent(Agent):
    def __init__(self, parent):
        super(PacbioStatAgent, self).__init__(parent)
        options = [
            {"name":"ccsbam","type":"infile","format":"align.bwa.bam"},
            {"name":"split_dir","type":"infile","format":"denovo_rna_v2.common_dir"},
            {"name":"origin_bam","type":"infile","format":"align.bwa.bam"},
            {"name":"index_path","type":"string"},
            # {"name":"clean_dir","type":"infile"}
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_option(self):
        """
        参数检查
        """
        if not self.option("ccsbam"):
            raise OptionError("缺少ccs.new.bam")
        if not self.option("split_dir"):
            raise OptionError("没有找到拆分数据后的文件夹")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(PacbioStatAgent,self).end()

class PacbioStatTool(Tool):
    def __init__(self,config):
        super(PacbioStatTool, self).__init__(config)
        self.samtools = self.config.SOFTWARE_DIR+"/miniconda2/bin/samtools"
        self.ccs_bam  = self.option("ccsbam").prop['path']
        self.split_dir = self.option('split_dir').prop['path']
        self.seqkit = os.path.join(self.config.SOFTWARE_DIR,"bioinfo/meta/seqkit/seqkit")
        self.bam2fastq = self.config.SOFTWARE_DIR+ "/program/SmrtLink/smrtlink/smrtcmds/bin/bam2fastq"
        if self.option('origin_bam'):
            self.origin_bam = self.option("origin_bam").prop['path']

    def samtool_getreads(self, inputbam):
        output = os.path.join(self.work_dir, "statreads{}.txt".format(os.path.basename(inputbam)))
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/datasplit/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "getreads_{}.sh".format(now_time)
        cmd = "{} view -c {} > {}".format(self.samtools,inputbam, output)
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
        self.logger.info("开始统计{}的reads数".format(os.path.basename(inputbam)))
        self.logger.info(shell)
        command = self.add_command("stat_getreads{}{}".format(os.path.basename(inputbam).split(".")[-2].lower(),now_time),shell).run()
        self.wait(command)
        num = 0
        if command.return_code == 0 :
            self.logger.info("统计成功")
        else:
            self.set_error("统计出现错误")
        with open(output,'r') as sf:
            num = int(sf.read())
        return num

    def samtool_stat(self, inputbam):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/datasplit/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "samtool_stat_{}.sh".format(now_time)
        output = os.path.join(self.work_dir, "samtool_stat.txt")
        cmd = "{} stats {} > {}".format(self.samtools,inputbam, output)
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
        self.logger.info("开始统计{}的bam".format(os.path.basename(inputbam)))
        self.logger.info(shell)
        command = self.add_command("stat",shell).run()
        self.wait(command)
        total_num = 0
        avg_num = 0
        max_num = 0
        if command.return_code == 0 :
            self.logger.info("samtools statistic统计成功")
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

    def run_statistic_ccs(self):
        subreads_num = self.samtool_getreads(self.origin_bam)
        ccsreads_num = self.samtool_getreads(self.ccs_bam)
        if subreads_num != 0:
            get_ratio = float(ccsreads_num)/subreads_num
        recon_reads = 0
        num_ccs= 0
        for i in self.index.keys():
            file = os.path.join(self.split_dir,"split.{}.bam".format(i))
            if not os.path.exists(file):
                continue
            num = self.samtool_getreads(file)
            recon_reads += num
            num_ccs+=1
        if ccsreads_num != 0:
            recon_retio = float(recon_reads)/ccsreads_num
        else:
            recon_retio = 0
        avg_len_ccs,_,_ = self.samtool_stat(self.ccs_bam)
        self.output_name = os.path.basename(self.origin_bam).split('.')[0]
        self.output_sta = os.path.join(self.output_dir,self.output_name)
        with open(self.output_sta,'w') as ots:
            line = "subreads_num\tccsreads_num\tget_ratio\treconize_reads\treconize_retio\tavg_len_ccs\n"
            ots.write(line)
            ots.write("{}\t{}\t{}\t{}\t{}\t{}".format(subreads_num,ccsreads_num,get_ratio,recon_reads,recon_retio,avg_len_ccs))
    
    def run_statistic_reads(self):
        """
        bam2fastq -u -o name bam
        seqkit fx2tab -l -n -i -H fastq|cut -f 2|sort|uniq -c|sort -n -k2> txt
        """
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/datasplit/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "fastqstat_{}.sh".format(now_time)
        stat_file = os.path.join(self.work_dir,"{}.txt".format("seqLength.txt"))
        cmd = "{} -u -o temp {}\n".format(self.bam2fastq,self.option("ccsbam").prop['path'])
        cmd += "{} fx2tab -l -n -i -H temp.fastq|cut -f 2|sort|uniq -c|sort -n -k2> {}".format(self.seqkit,stat_file)
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
        self.logger.info("开始统计reads数")
        self.logger.info(shell)
        command1 = self.add_command("statreads", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行统计完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))
        self.stato_file = os.path.join(self.output_dir, "statReads.new.txt")
        with open(stat_file,'r') as instat,open(self.stato_file,'w') as ostat:
            readDistribution = {}
            line = instat.readline()
            for i in range(500,8001,500):
                readDistribution[i] = 0
            while 1:
                line = instat.readline()
                if not line:
                    break
                fd = line.rstrip().split(' ')
                read_num = int(fd[-1])
                if read_num <=500:
                    readDistribution[500] += int(fd[-2])
                elif read_num <=1000:
                    readDistribution[1000] += int(fd[-2])
                elif read_num <=1500:
                    readDistribution[1500] += int(fd[-2])
                elif read_num <=2000:
                    readDistribution[2000] += int(fd[-2])
                elif read_num <=2500:
                    readDistribution[2500] += int(fd[-2])
                elif read_num <=3000:
                    readDistribution[3000] += int(fd[-2])
                elif read_num <=3500:
                    readDistribution[3500] += int(fd[-2])
                elif read_num <=4000:
                    readDistribution[4000] += int(fd[-2])
                elif read_num <=4500:
                    readDistribution[4500] += int(fd[-2])
                elif read_num <=5000:
                    readDistribution[5000] += int(fd[-2])
                elif read_num <=5500:
                    readDistribution[5500] += int(fd[-2])
                elif read_num <=6000:
                    readDistribution[6000] += int(fd[-2])
                elif read_num <=6500:
                    readDistribution[6500] += int(fd[-2])
                elif read_num <=7000:
                    readDistribution[7000] += int(fd[-2])
                elif read_num <=7500:
                    readDistribution[7500] += int(fd[-2])
                else:
                    readDistribution[8000] += int(fd[-2])
            for key in sorted(readDistribution.keys()):
                ostat.write("{}\t{}\n".format(key,readDistribution[key]))

    def get_index(self):
        self.index = {}
        with open(self.option("index_path"),'r') as ip:
            while 1:
                line = ip.readline()
                if not line:
                    break 
                fd = line.rstrip().split('\t')
                self.index[fd[0]]={}
                self.index[fd[0]][1]= fd[1]
                self.index[fd[0]][2]= fd[2]
                self.index[fd[0]]['type'] = fd[3]
                self.index[fd[0]]['primer_type'] = fd[4]

    def run(self):
        super(PacbioStatTool,self).run()
        self.get_index()
        self.run_statistic_ccs()
        self.run_statistic_reads()
        self.end()


class TestFunctionf(unittest.TestCase):
    '''
    测试脚本
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "pastat_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "datasplit_v2.pacbio_stat",
            "options": {
                "ccsbam":"/mnt/ilustre/users/sanger-dev/workspace/20210915/Single_pbmerge_6488/Pbmerge/ccs.new.bam",
                "split_dir":"/mnt/ilustre/users/sanger-dev/workspace/20210915/Single_lima_9691/Lima/output/split",
                "origin_bam":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/checkSmrtLink/split.F7--R10.bam",
                "index_path":'/mnt/ilustre/users/sanger-dev/workspace/20210915/PacbioSplit_20181016PE300_20210915_105314/sample_list.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

    
        

        
