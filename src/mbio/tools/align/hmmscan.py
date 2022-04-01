# -*- coding: utf-8 -*-
# __author__ :zhouxuan

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess
import math
from mbio.packages.denovo_rna.gene_structure.pfam_domtblout import pfam_out


class HmmscanAgent(Agent):
    """
    为cazy注释写的比对软件 hummscan
    version 1.0
    author: zhouxuan
    last_modify: 20170524
    """

    def __init__(self, parent):
        super(HmmscanAgent, self).__init__(parent)
        options = [
            {"name": "faa_file", "type": "infile", "format": "sequence.fasta"},
            {"name": "database", "type": "string", "default": "cazy"},  # 数据库：cazy, pfam  add by zhujuan10180102
            {"name": "hmmscan_out_dm", "type": "outfile", 'format': "meta_genomic.hmmscan_table"},  # 结果文件的设置
        ]
        self.add_option(options)
        self.step.add_steps("hmmscan")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.hmmscan.start()
        self.step.update()

    def stepfinish(self):
        self.step.hmmscan.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('faa_file'):
            raise OptionError("必须输入fastq文件", code="31101401")
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        """
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""],
        ])
        super(HmmscanAgent, self).end()


class HmmscanTool(Tool):
    def __init__(self, config):
        super(HmmscanTool, self).__init__(config)
        self._version = "v1.0"
        self.soft_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan"
        # self.data_path = self.config.SOFTWARE_DIR + "/database/CAZyDB/dbCAN-fam-HMMs.txt.v5"
        self.data_path = self.config.SOFTWARE_DIR + "/database/CAZyDB/dbCAN-fam-HMMs.txt.v6"  # last_modify  by zhujuan10180102 更新数据库至V6版本
        if self.option("database") in ['cazy_v8']:
            self.data_path = self.config.SOFTWARE_DIR + "/database/CAZyDB/dbCAN-HMMdb-V8.txt" ## qingchen.zhang @ 20201009
        self.pfam_db = self.config.SOFTWARE_DIR + "/database/pfam_31/Pfam-A.hmm"
        if self.option("database") in ['pfam_v33.1']:
            self.pfam_db = self.config.SOFTWARE_DIR + "/database/pfam_33.1/Pfam-A.hmm"

    def run(self):
        """
        运行
        :return:
        """
        super(HmmscanTool, self).run()
        if self.option("database") in ["cazy", 'cazy_v8']:## fix by qingchen.zhang @20201009
            self.run_hmm()
        elif self.option("database") in ["pfam", "pfam_v33.1"]:
            self.run_pfam_hmmscan()
        else:
            self.set_error("不存在该数据库名称有误！", code="31101401")
        self.set_output()
        self.end()

    def run_hmm(self):
        hmm_out_dm = os.path.join(self.output_dir, "dbCAN.hmmscan.out.dm")
        hmm_out = os.path.join(self.output_dir, "dbCAN.hmmscan.out")
        for i in os.listdir(self.output_dir):
            i_path = os.path.join(self.output_dir, i)
            os.remove(i_path)
        cmd = "{} --domtblout {} -o {} --cpu 6 {} {}" \
            .format(self.soft_path, hmm_out_dm, hmm_out, self.data_path, self.option("faa_file").prop['path'])
        self.logger.info("start run_hmm")
        command = self.add_command("hmmscan", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_hmm done")
        else:
            command.rerun()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_hmm done")
            else:
                self.set_error("run_hmm error", code="31101402")
                # raise Exception("run_hmm error")

    def run_pfam_hmmscan(self):
        cmd = "{} --cpu 10 --noali --cut_nc --acc --notextw --domtblout {} {} {}".format(
            self.soft_path, 'pfam.domtblout', self.pfam_db, self.option("faa_file").prop['path'])
        print(cmd)
        self.logger.info("开始运行hmmscan")
        command = self.add_command("hmmscan", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行hmmscan结束")
        else:
            self.set_error("运行hmmscan出错", code="31101403")

        pfam_out("pfam.domtblout",score=True)  # 增加score


    def set_output(self):
        if self.option("database") in ["cazy", "cazy_v8"] and len(os.listdir(self.output_dir)) == 2:
            self.logger.info("hmm比对的结果文件正确生成")
            number = os.path.basename(self.option('faa_file').prop['path']).split('_')[-1]
            if os.path.exists(os.path.join(self.output_dir, "dbCAN.hmmscan.out")):
                os.remove(os.path.join(self.output_dir, "dbCAN.hmmscan.out"))
            os.link(os.path.join(self.output_dir, "dbCAN.hmmscan.out.dm"),
                    os.path.join(self.output_dir, "dbCAN.hmmscan.out.dm_{}".format(number)))
            os.remove(os.path.join(self.output_dir, "dbCAN.hmmscan.out.dm"))
            self.logger.info('dbCAN.hmmscan.out原始比对结果删除')
            self.option("hmmscan_out_dm").set_path(self.output_dir + '/dbCAN.hmmscan.out.dm_{}'.format(number))
        # self.option('hmmscan_out_dm', self.output_dir + '/dbCAN.hmmscan.out.dm_{}'.format(number))
        elif self.option("database") in ["pfam", "pfam_v33.1"]:
            if os.path.exists(self.output_dir + "/" + "pfam_domain"):
                os.remove(self.output_dir + "/" + "pfam_domain")
            os.link(self.work_dir + "/" + "pfam_domain", self.output_dir + "/" + "pfam_domain")
        else:
            self.set_error("hmm比对结果出错", code="31101404")
