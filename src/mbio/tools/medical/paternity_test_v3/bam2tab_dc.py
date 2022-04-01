# -*- coding: utf-8 -*-
# __author__ = 'hongyu.chen'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import subprocess
import datetime


class Bam2tabDcAgent(Agent):
    def __init__(self, parent):
        super(Bam2tabDcAgent, self).__init__(parent)
        options = [
            {"name": "sample_id", "type": "string"},                 # 输入F/M/S的fastq文件的样本名
            {"name": "bam_path", "type": "string"},                  # bam文件的路径
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # hg38.chromosomal_assembly/ref.fa
            {"name": "targets_bedfile", "type": "infile", "format": "paternity_test.rda"},
            {"name": "batch_id", "type": "string"},                    # 拆分主表id
            {"name": "board_batch", "type": "string"},                  # 版号170928_TPNB500180_0112_AHMWJKAFXX
            {"name": "cpu_number", "type": "int", "default": 4},
            {"name": "memory", "type": "int", "default": 10},
            {"name": "is_update", "type": "string", "default": "true"}  # add by hd 20180204
        ]
        self.add_option(options)
        self.step.add_steps("Bam2tabDc")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = "gada"

    def stepstart(self):
        self.step.Bam2tabDc.start()
        self.step.update()

    def stepfinish(self):
        self.step.Bam2tabDc.finish()
        self.step.update()

    def check_options(self):
        if not self.option("sample_id"):
            raise OptionError("必须输入样本编号")
        if not self.option("bam_path"):
            raise OptionError("必须输入bam文件的所在路径")
        if not self.option("ref_fasta").is_set:
            raise OptionError("必须输入参考基因组序列fasta文件")
        if not self.option('targets_bedfile'):
            raise OptionError('必须提供target_bedfile文件')
        if not self.option("batch_id"):
            raise OptionError("必须输入batch_id")
        if not self.option("board_batch"):
            raise OptionError("必须输入board_batch")
        return True

    def set_resource(self):
        self._cpu = self.option("cpu_number")
        self._memory = "{}G".format(self.option("memory"))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [r".mem.sort.hit.vcf.tab", "tab", "所有位点的信息"],
            [r".qc", "qc", "质控文件"]
        ])
        super(Bam2tabDcAgent, self).end()


class Bam2tabDcTool(Tool):
    def __init__(self, config):
        super(Bam2tabDcTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ruby-2.4.1/bin')  # 正式机
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/ruby-2.3.1')  # 测试机
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bioruby-vcf-master/bin')  # 测试机
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bioawk')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/seqtk-master')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/bwa-0.7.15')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/samblaster-0.1.24')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/medical/bedtools-2.24.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bcftools-1.4/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/vt-0.577')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/vcflib-master/bin')
        self.cmd_path = "bioinfo/medical/scripts"
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.bcftools = "bioinfo/seq/bcftools-1.4/bin/bcftools"
        self.vt = "bioinfo/seq/vt-0.577/vt"
        self.script_path = "bioinfo/medical/scripts/pt_scripts"
        self.script_absolute_path = self.config.SOFTWARE_DIR + '/' + self.script_path
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_data"
        self.ref_all = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_all"

    def bam2vcf(self, rerun=False):
        sample_id = self.option("sample_id")
        bam_path = self.option("bam_path")
        ref = self.option("ref_fasta").prop["path"]
        bed = self.option("targets_bedfile").prop["path"]
        if not os.path.exists(self.script_absolute_path):
            os.mkdir(self.script_absolute_path)
        now_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        full_path = self.script_absolute_path + '/' + sample_id + '_bam2vcf_dcpt.sh' + now_time
        path = self.script_path + '/' + sample_id + '_bam2vcf_dcpt.sh' + now_time
        cmd_pool = list()
        with open(full_path, 'w') as out:
            header_line = '#!/bin/bash'
            if not rerun:
                cmd1 = 'samtools mpileup -uvf {} -l {} {} -o {}.bcf\n'.format(ref, bed, bam_path, sample_id)
                cmd1 += "[ $(less {sample_id}.bcf|grep -v '^#'|wc -l) -eq 0 ] && samtools mpileup -Auvf {ref} -l " \
                        "{bed} {bam_path} -o {sample_id}.bcf".format(sample_id=sample_id, ref=ref, bed=bed,
                                                                     bam_path=bam_path)  # bcf文件位点信息为空，用-A参数重新运行。
                # 解决正式机上因Ruby版本不同产生的异常
            else:
                cmd1 = 'samtools mpileup -Auvf {} -l {} {} -o {}.bcf   # 单端运行'.format(ref, bed, bam_path, sample_id)
            cmd_pool.append(cmd1)

            cmd2 = 'bcftools call --multiallelic-caller --keep-alts --targets-file {} -Oz {}.bcf -o ' \
                   '{}.mem.sort.hit.vcf.gz'.format(bed, sample_id, sample_id)
            cmd_pool.append(cmd2)

            cmd3 = "bcftools view -i 'INDEL=1' {}.mem.sort.hit.vcf.gz -o {}.mem.sort.hit.indel.vcf".format(
                sample_id, sample_id)
            cmd_pool.append(cmd3)

            cmd4 = "bioawk -c vcf '{{print $chrom\"\\t\"$pos}}' {}.mem.sort.hit.indel.vcf > {}.indel.reg.txt".format(
                sample_id, sample_id
            )
            cmd_pool.append(cmd4)

            out.write('{}\n'.format(header_line))
            for i in range(0, len(cmd_pool)):
                out.write('{}\n'.format(cmd_pool[i]))
        self.run_single_command('chmod +x {}'.format(full_path))
        return path, full_path

    def check_indel_output(self):
        sample_id = self.option("sample_id")
        number = 0
        with open(sample_id + ".mem.sort.hit.indel.vcf", 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    number += 1
                if number != 0:
                    break
        if number == 0:
            return True
        else:
            return False

    def deal_with_vcf(self):
        sample_id = self.option("sample_id")
        ref = self.option("ref_fasta").prop["path"]
        if not os.path.exists(self.script_absolute_path):
            os.mkdir(self.script_absolute_path)
        now_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        full_path = self.script_absolute_path + '/' + sample_id + '_vcf2tab_dcpt.sh' + now_time
        path = self.script_path + '/' + sample_id + '_vcf2tab_dcpt.sh' + now_time
        cmd_pool = list()
        with open(full_path, 'w') as out:
            header_line = '#!/bin/bash'
            if self.check_indel_output():
                cmd1 = "bcftools view -e 'INDEL=1' {}.mem.sort.hit.vcf.gz -o {}.mem.sort.hit.snp.vcf".format(
                    sample_id, sample_id)
            else:
                step1 = "bcftools view -e 'INDEL=1' {}.mem.sort.hit.vcf.gz -o {}.temp.vcf\n".format(sample_id, sample_id
                                                                                                    )
                step2 = "bcftools view -T ^{}.indel.reg.txt {}.temp.vcf -o {}.mem.sort.hit.snp.vcf".format(
                    sample_id, sample_id, sample_id)
                cmd1 = step1 + step2
            cmd_pool.append(cmd1)

            cmd2 = 'cat {}.mem.sort.hit.snp.vcf <(bcftools view -H {}.mem.sort.hit.indel.vcf)' \
                   '>{}.mem.sort.hit.snp.indel.vcf'.format(sample_id, sample_id, sample_id, sample_id)
            cmd_pool.append(cmd2)

            cmd3 = 'vcfstreamsort < {}.mem.sort.hit.snp.indel.vcf > {}.mem.sort.hit.vcf '.format(sample_id, sample_id)
            cmd_pool.append(cmd3)

            cmd4 = 'vt normalize -r {} {}.mem.sort.hit.vcf -o {}.mem.sort.hit.dedup.vcf'.format(
                ref, sample_id, sample_id)
            cmd_pool.append(cmd4)

            cmd5 = "bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(\",\"),r.info.dp," \
                   "r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]' " \
                "< {}.mem.sort.hit.dedup.vcf > {}.mem.sort.hit.dedup.reformed.vcf".format(sample_id, sample_id)
            cmd_pool.append(cmd5)

            out.write('{}\n'.format(header_line))
            for i in range(0, len(cmd_pool)):
                out.write('{}\n'.format(cmd_pool[i]))
        self.run_single_command('chmod +x {}'.format(full_path))
        return path, full_path

    def vcf2tab(self):
        sample_id = self.option("sample_id")
        out = open(sample_id + '.mem.sort.hit.vcf.tab1', 'w')
        out1 = open(sample_id + '.chrY.tab', 'w')
        with open(sample_id + '.mem.sort.hit.dedup.reformed.vcf', 'r') as f:
            for line in f:
                line_split = line.split('\t')
                out.write(sample_id + '\t' + line)
                if line_split[0] == 'chrY' and float(line_split[4]) > 3:   # 修复深度过滤错位 add by hd 20180129
                    out1.write(sample_id + '\t' + line)
        out.close()
        out1.close()
        dcpt_site_pri = self.config.SOFTWARE_DIR + "/database/human/pt_ref/dcpt_site_pri.txt"  # modified by HD 20180112
        self.logger.info(dcpt_site_pri)
        # dcpt_site_pri = '/mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/dcpt_site_pri.txt'
        # dcpt_site_pri = '/mnt/ilustre/users/sanger/app/database/human/pt_ref/dcpt_site_pri.txt'
        sites = list()
        with open(dcpt_site_pri, 'r') as f1:
            for line in f1:
                line_split = line.rstrip('\n').split('\t')
                sites.append(line_split[0])

        out2 = open(sample_id + '.mem.sort.hit.vcf.tab', 'w')
        with open(sample_id + '.mem.sort.hit.vcf.tab1', 'r') as f2:
            for line in f2:
                line_split = line.split('\t')
                if line_split[2] in sites:
                    out2.write(line)
        out2.close()

    def run_bam2tab(self, rerun=False):
        cmd = self.bam2vcf(rerun)
        cmd_name = 'bam2vcf_rerun' if rerun else 'bam2vcf'
        self.run_command(cmd_name, cmd[0])
        os.system('rm {}'.format(cmd[1]))                 # 清理程序生成的临时sh脚本

        cmd1 = self.deal_with_vcf()
        cmd_name1 = 'deal_with_vcf_rerun' if rerun else 'deal_with_vcf'
        self.run_command(cmd_name1, cmd1[0])
        os.system('rm {}'.format(cmd1[1]))                # 清理程序生成的临时sh脚本

        self.vcf2tab()

    def qc(self):
        sample_id = self.option("sample_id")
        bam_path = self.option("bam_path")
        bam_dir = os.path.dirname(bam_path)
        mem_sort_bam = os.path.join(bam_dir, '{}.mem.sort.bam'.format(sample_id))
        mem_sort_hit_bam = os.path.join(bam_dir, '{}.mem.sort.hit.bam'.format(sample_id))
        cpu = self.option("cpu_number") - 1
        cmd = '{}/dcpt_qc.sh {} {} {} {}'.format(self.cmd_path, cpu, sample_id, mem_sort_bam, mem_sort_hit_bam)
        self.run_command('qc', cmd)

    def set_output(self):
        sample_id = self.option("sample_id")
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        file0 = '/{}.mem.sort.hit.vcf.tab'.format(sample_id)
        file1 = '/{}.qc'.format(sample_id)
        file2 = '/{}.ref.dp.xls'.format(sample_id)
        os.link(self.work_dir + file0, self.output_dir + file0)
        os.link(self.work_dir + file1, self.output_dir + file1)
        os.link(self.work_dir + file2, self.output_dir + file2)
        self.logger.info("设置文件夹路径成功")

        qc_path = self.output_dir + '/' + sample_id + '.qc'
        tab_path = self.output_dir + '/' + sample_id + '.mem.sort.hit.vcf.tab'
        qc = dict()
        with open(qc_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n').split(':')
                qc[line[0]] = line[1]
        line_count = 0
        db = self.api.api('medical.paternity_test_v2')
        if os.path.getsize(tab_path) > 0:
            if float(qc['dp1']) > 10:
                self.logger.info('{}该样本非异常样本'.format(sample_id))
                file_name = sample_id + ".tab"
                if '-F' in sample_id:                             # 如果该样本为爸爸将tab文件放入查重父本库
                    if not os.path.exists(self.ref_data + "/" + file_name):
                        os.link(self.output_dir + file0, self.ref_data + "/" + file_name)
                        self.logger.info("参考库中没有{},开始移动该样本".format(sample_id))
                    else:
                        self.logger.info("参考库中有{}".format(sample_id))
                if not os.path.exists(self.ref_all + "/" + file_name):
                    os.link(self.output_dir + file0, self.ref_all + "/" + file_name)
                    self.logger.info("在tab_all中备份该样本：{}".format(sample_id))
                else:
                    self.logger.info("tab_all中已有该样本：{}".format(sample_id))
            table_name = 'sg_sample_tab' if float(qc['dp1']) > 10 else 'sg_sample_tab_problem'
            if db.check_existence(table_name, {'sample_id': sample_id}):
                pass
                # self.set_error('可能样本{}重名，请检查！'.format(sample_id))
                # raise Exception('可能样本重名，请检查！')
            # else:
            line_count = db.add_tab(tab_path, table_name)
            if table_name == 'sg_sample_tab':
                db.update_sample_tab(sample_id, self.option('board_batch'), self.option('batch_id'), 'dcpt')
            if db.check_existence('sg_sample_qc', {'sample_id': sample_id}):
                pass
                # self.set_error('可能样本{}重名，请检查！'.format(sample_id))
                # raise Exception('可能样本{}重名，请检查！'.format(sample_id))
            # else:
            db.dcpt_qc(qc_path, sample_id, line_count, self.option('board_batch'))
        else:
            db.dcpt_qc(qc_path, sample_id, line_count, self.option('board_batch'), False)
        db.update_analysis_status(self.option('batch_id'), 'snp', self.option("is_update"))   # modified by hongdong 20171218

    def run_command(self, cmd_name, cmd):
        self.logger.info('开始运行{}'.format(cmd_name))
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行{}完成".format(cmd_name))
        else:
            if command.return_code is None:
                command.rerun()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("重新运行{}完成".format(cmd_name))
                elif command.return_code is None:
                    self.logger.info("重新运行{}返回码仍然为None".format(cmd_name))
                else:
                    self.logger.info('Run Return Code: {}'.format(command.return_code))
                    raise Exception("{}运行出错!".format(cmd_name))
            else:
                self.logger.info('Run Return Code: {}'.format(command.return_code))
                raise Exception("{}运行出错!".format(cmd_name))

    def run_single_command(self, cmd):
        exitcode = subprocess.call(cmd, shell=True)
        if exitcode == 0:
            self.logger.info("运行成功：" + cmd)
        else:
            exitcode = subprocess.call(cmd, shell=True)
            if exitcode == 0:
                self.logger.info("重运行成功：" + cmd)
            else:
                raise Exception("运行失败：" + cmd)

    def run(self):
        super(Bam2tabDcTool, self).run()
        sample_id = self.option("sample_id")
        self.run_bam2tab()
        with open('{}.mem.sort.hit.vcf.tab'.format(sample_id), 'r') as f:
            tab_line_count = len(["" for line in f])
        self.logger.info("第一次得到的tab文件行数为：{}".format(tab_line_count))
        if tab_line_count < 1000:
            self.run_bam2tab(rerun=True)
        self.qc()
        self.set_output()
        self.end()
