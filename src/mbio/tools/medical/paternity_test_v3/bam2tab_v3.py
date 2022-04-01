# -*- coding: utf-8 -*-
# __author__ = 'hongyu.chen'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import subprocess
import datetime
import time
import random


class Bam2tabV3Agent(Agent):
    def __init__(self, parent):
        super(Bam2tabV3Agent, self).__init__(parent)
        options = [
            {"name": "sample_id", "type": "string"},                 # 输入F/M/S的fastq文件的样本名
            {"name": "bam_path", "type": "string"},                  # bam文件的路径
            {"name": "ref_fasta", "type": "string"},  # hg38.chromosomal_assembly/ref.fa
            {"name": "targets_bedfile", "type": "infile", "format": "paternity_test.rda"},
            {"name": "batch_id", "type": "string"},                    # 拆分主表id
            {"name": "board_batch", "type": "string"},                  # 版号170928_TPNB500180_0112_AHMWJKAFXX
            {"name": "cpu_number", "type": "int", "default": 4},
            {"name": "memory", "type": "int", "default": 10},
            {"name": "split_type", "type": "string", "default": "PE"},
            {"name": "is_update", "type": "string", "default": "true"}  # add by hd 20180204
        ]
        self.add_option(options)
        self.step.add_steps("Bam2tabV3")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self.queue = "gada"

    def stepstart(self):
        self.step.Bam2tabV3.start()
        self.step.update()

    def stepfinish(self):
        self.step.Bam2tabV3.finish()
        self.step.update()

    def check_options(self):
        if not self.option("sample_id"):
            raise OptionError("必须输入样本编号")
        if not self.option("bam_path"):
            raise OptionError("必须输入bam文件的所在路径")
        if not self.option("ref_fasta"):
            raise OptionError("必须输入参考基因组序列fasta文件")
        if not self.option('targets_bedfile'):
            raise OptionError('必须提供target_bedfile文件')
        if not self.option("batch_id"):
            raise OptionError("必须输入batch_id")
        if not self.option("board_batch"):
            raise OptionError("必须输入board_batch")
        if self.option("split_type") not in ["PE", "SE"]:
            raise OptionError("split_type 必须为PE或SE")
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
        super(Bam2tabV3Agent, self).end()


class Bam2tabV3Tool(Tool):
    def __init__(self, config):
        super(Bam2tabV3Tool, self).__init__(config)
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
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_v3_fm"
        self.ref_all = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_all"

        self.threshold = 0.005                              # 对单端亲子多突变碱基进行过滤的阈值

    def bam2vcf(self):
        sample_id = self.option("sample_id")
        bam_path = self.option("bam_path")
        ref = self.option("ref_fasta")
        bed = self.option("targets_bedfile").prop["path"]
        if not os.path.exists(self.script_absolute_path):
            os.mkdir(self.script_absolute_path)
        now_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        random_num = str(random.randint(1000, 10000))
        full_path = self.script_absolute_path + '/' + sample_id + '_bam2vcf_dcpt.sh' + now_time + random_num
        path = self.script_path + '/' + sample_id + '_bam2vcf_dcpt.sh' + now_time + random_num
        cmd_pool = list()
        with open(full_path, 'w') as out:
            header_line = '#!/bin/bash'

            cmd1 = 'samtools mpileup -d 200000 -Q 20 -t INFO/AD -Auvf {} -l {} {} -o {}.bcf'.format(ref, bed, bam_path,
                                                                                                    sample_id)
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
        ref = self.option("ref_fasta")
        if not os.path.exists(self.script_absolute_path):
            os.mkdir(self.script_absolute_path)
        now_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        random_num = str(random.randint(1000, 10000))
        full_path = self.script_absolute_path + '/' + sample_id + '_vcf2tab_dcpt.sh' + now_time + random_num
        path = self.script_path + '/' + sample_id + '_vcf2tab_dcpt.sh' + now_time + random_num
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

            cmd5 = "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP\\t%DP4\\t%AD\\n' {}.mem.sort.hit.dedup.vcf" \
                   " > {}.vcf.tab".format(sample_id, sample_id)
            cmd_pool.append(cmd5)

            out.write('{}\n'.format(header_line))
            for i in range(0, len(cmd_pool)):
                out.write('{}\n'.format(cmd_pool[i]))
        self.run_single_command('chmod +x {}'.format(full_path))
        return path, full_path

    def vcf2tab(self):
        """
        1、根据dp4的四个值算出ref_dp和alt_dp
        2、当样品为单端样品时，根据AD和设定的阈值，对多个突变碱基进行过滤
        3、挑选出 dp > 3 的chrY位点
        :return:
        """
        sample_id = self.option("sample_id")
        with open(sample_id + ".vcf.tab", "r") as f, open(sample_id + ".mem.sort.hit.vcf.tab", "w") as out:
            for line in f:
                chrome, pos, ref, alt, dp, dp4, ad = line.rstrip("\n").split("\t")
                alt_list = alt.split(",")
                ad_list = [float(i) for i in ad.split(",")]
                dp4_list = [int(i) for i in dp4.split(",")]
                if len(alt_list) + 1 != len(ad_list) and alt != ".":  # 当alt不为"."时，ad的值应该比alt的值多一个
                    print line
                    raise Exception("{}.tab：AD与ALT不匹配".format(sample_id))
                # 过滤多个突变碱基
                # if len(alt_list) > 1 and self.option("split_type") == "SE":  # 仅在alt有多个且为单端样品时，对多个alt进行过滤
                if False:                                                      # 跳过对单端样品多突变碱基的过滤
                    all_dp = sum(ad_list)
                    new_list = list()
                    for i in range(len(alt_list)):
                        if ad_list[i + 1] / all_dp >= self.threshold:  # 过滤掉深度占比小于阈值的突变
                            new_list.append(alt_list[i])
                    if len(new_list) == 0:  # 如果alt碱基全被过滤，将原先的alt修改为"."
                        alt = "."
                    else:  # 否则，将过滤后的碱基赋给alt
                        alt = ",".join(new_list)
                # 计算ref_dp和alt_dp
                ref_dp = dp4_list[0] + dp4_list[1]
                alt_dp = dp4_list[2] + dp4_list[3]
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_id, chrome, pos, ref, alt, dp, ref_dp,
                                                                    alt_dp))

        # 筛选chrY位点
        with open(sample_id + ".mem.sort.hit.vcf.tab", "r") as f1, open(sample_id + '.chrY.tab', 'w') as out1:
            for line in f1:
                name, chrome, pos, ref, alt, dp, ref_dp, alt_dp = line.rstrip("\n").split("\t")
                if chrome == "chrY" and int(ref_dp) + int(alt_dp) > 3:
                    out1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrome, pos, ref, alt, dp, ref_dp, alt_dp))

    def run_bam2tab(self):
        cmd = self.bam2vcf()
        cmd_name = 'bam2vcf'
        self.run_command(cmd_name, cmd[0])
        os.system('rm {}'.format(cmd[1]))                 # 清理程序生成的临时sh脚本

        cmd1 = self.deal_with_vcf()
        cmd_name1 = 'deal_with_vcf'
        self.run_command(cmd_name1, cmd1[0])
        os.system('rm {}'.format(cmd1[1]))                # 清理程序生成的临时sh脚本

        self.vcf2tab()

    def qc(self):
        sample_id = self.option("sample_id")
        bam_path = self.option("bam_path")
        bam_dir = os.path.dirname(bam_path)
        mem_sort_bam = os.path.join(bam_dir, '{}.mem.sort.bam'.format(sample_id))
        mem_sort_hit_bam = os.path.join(bam_dir, '{}.mem.sort.hit.bam'.format(sample_id))
        cpu = self.option("cpu_number")
        cmd = '{}/dcpt_qc_v3.sh {} {} {} {}'.format(self.cmd_path, cpu, sample_id, mem_sort_bam, mem_sort_hit_bam)
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
        db = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        if os.path.getsize(tab_path) > 0:
            if float(qc['dp1']) > 10:
                self.logger.info('{}该样本非异常样本'.format(sample_id))
                file_name = sample_id + ".tab"
                if '-F' in sample_id or '-M' in sample_id:                # 将父本、母本放入查重库
                    if not os.path.exists(self.ref_data + "/" + file_name):
                        os.link(self.output_dir + file0, self.ref_data + "/" + file_name)
                        self.logger.info("参考库中没有{},开始移动该样本".format(sample_id))
                    else:
                        self.logger.info("参考库中有{}".format(sample_id))
                time.sleep(random.randint(1, 30))  # 随机停避免多个tool同时复制胎儿tab造成的文件冲突2019111 by hd
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
                db.update_sample_tab(sample_id, self.option('board_batch'), self.option('batch_id'), 'wqcf')
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
        super(Bam2tabV3Tool, self).run()
        self.run_bam2tab()
        self.qc()
        self.set_output()
        self.end()
