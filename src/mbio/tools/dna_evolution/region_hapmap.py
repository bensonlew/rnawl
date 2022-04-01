# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# last_modify: 2018.08.29

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import time
from bson.objectid import ObjectId


class RegionHapmapAgent(Agent):
    """
    一个性状的三个csv，运行三次这个tool
    glm_csv_path mlm_csv_path farmcpu_csv_path
    按照一个性状的三个方法并行
    distance ld_r2二选一 ，p_value q_value二选一，互组搭配运行程序
    """
    def __init__(self, parent):
        super(RegionHapmapAgent, self).__init__(parent)
        options = [
            {"name": "recode_vcf_path", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "anno_summary_path", "type": "string"},
            {"name": "csv_path", "type": "string"},
            {"name": "pop_path", "type": "string"},
            {"name": "distance", "type": "int"},    # , "default": 200000
            {"name": "ld_r2", "type": "float"},     # , "default": 0.8
            {"name": "p_value", "type": "float"},   # , "default": 0.05
            {"name": "q_value", "type": "float"},   # , "default": 0.05
            {"name": "task_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "is_wgs_result", "type": "string"},
            {"name": "file_path", "type": "string"}  # 对象存储路径
        ]
        self.add_option(options)
        self.step.add_steps("region_hapmap")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.region_hapmap.start()
        self.step.update()

    def stepfinish(self):
        self.step.region_hapmap.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('recode_vcf_path'):
            raise OptionError('必须输入:recode_vcf_path')
        if not self.option('anno_summary_path'):
            raise OptionError('必须输入:anno_summary_path')
        if not self.option('csv_path'):
            raise OptionError('必须输入:csv_path')
        if not self.option('pop_path'):
            raise OptionError('必须输入:pop_path')
        if not self.option('p_value') and not self.option('q_value'):
            raise OptionError('p:{}/q_value:{}有且只能输入一个'.format(self.option('p_value'), self.option('q_value')))
        if not self.option('distance') and not self.option('ld_r2'):
            raise OptionError('distance:{}/ld_r2:{}有且只能输入一个'.format(self.option('distance'), self.option('ld_r2')))
        if not self.option('task_id'):
            raise OptionError('必须输入task_id')
        if not self.option("main_id"):
            raise OptionError("请设置main_id")
        if not self.option("file_path"):
            raise OptionError("请设置file_path")
        if not self.option("is_wgs_result"):
            raise OptionError("请设置is_wgs_result")

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(RegionHapmapAgent, self).end()


class RegionHapmapTool(Tool):
    """
    /mnt/ilustre/users/sanger-dev/app/bioinfo/dna_evolution/plink/plink
    """
    def __init__(self, config):
        super(RegionHapmapTool, self).__init__(config)
        self._version = "v1.0.1"
        # self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.R_path = "program/R-3.3.3/bin/Rscript"
        self.plink_path = "bioinfo/dna_evolution/plink/plink"
        self.plink_more_path = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/plink/plink"
        self.vcftools_path = self.config.SOFTWARE_DIR + "/bioinfo/dna_evolution/vcftools"
        self.region_snp_path = self.config.PACKAGE_DIR + "/dna_evolution/region.snp.pl"
        self.ld_region_path = self.config.PACKAGE_DIR + "/dna_evolution/ld_region.pl"
        self.gwas_region_path = self.config.PACKAGE_DIR + "/dna_evolution/gwas.region.pl"
        self.region_stat_path = self.config.PACKAGE_DIR + "/dna_evolution/region.stat.pl"
        self.ldheatmap_vcfsplit_path = self.config.PACKAGE_DIR + "/dna_evolution/LDheatmap.vcfsplit.pl"
        self.LDheatmap_vcf2matrix_path = self.config.PACKAGE_DIR + "/dna_evolution/LDheatmap.vcf2matrix.pl"
        self.LDheatmap_path = self.config.PACKAGE_DIR + "/dna_evolution/LDheatmap.R"
        self.anno_summary_path = os.path.join((self.config.SOFTWARE_DIR + "/database/dna_geneome/"),
                                              self.option("anno_summary_path"))
        if self.option("is_wgs_result") == "yes":
            self.anno_summary_path = os.path.join((self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/"),
                                                  self.option("anno_summary_path"))
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.para_fly = 'program/parafly-r2013-01-21/bin/bin/ParaFly'  # 并行投递cmd
        self.threshold, self.adjust = '', ''
        if self.option('p_value'):
            self.threshold = self.option('p_value')
            self.adjust = 0
        else:
            self.threshold = self.option('q_value')
            self.adjust = 1
        name = os.path.basename(self.option('csv_path')).strip().split('.')
        self.output_name = name[-3] + "." + name[-2]
        self.vcfsplit_name = name[-3] + "_" + name[-2]
        self.file_path = []
        self.analysis = name[-2]
        self.trait = name[-3]
        self.name_list = []

    def run_region_snp(self):
        """
        第一步
        perl /mnt/ilustre/users/minghao.zhang/newmdt/ZMH_Pipeline/Gwas_v2.0/bin/bin/region.snp.pl
        -input ./TRAIT3_GLM.csv -output ./select.snp -threshold 0.05 -adjust 1
        trait.GLM.select.snp /GLM/FarmCPU/MLM
        存在路径select_snp下
        """
        if not os.path.exists(self.work_dir + "/select_snp"):
            os.mkdir(self.work_dir + "/select_snp")
        cmd = "{} {}".format(self.perl_path, self.region_snp_path)
        cmd += " -input {}".format(self.option('csv_path'))
        cmd += " -threshold {}".format(self.threshold)
        cmd += " -adjust {}".format(self.adjust)
        cmd += " -output {}".format(self.work_dir + "/select_snp/" + self.output_name + ".select.snp")
        self.logger.info('开始运行region_snp:{}提取'.format(self.output_name))
        self.run_cmd(cmd, "region_snp_{}".format(str(self.output_name).lower()))

    def run_plink(self):
        """
        第二步
        /mnt/ilustre/users/dna/.env/bin/plink
        --file /mnt/ilustre/users/qingmei.cui/newmdt/sanger/7.pop/hapmap-hyplotype/pop
        --r2 --ld-snp-list /select.snp --ld-window-kb 500 --ld-window 99999
        --out /mnt/ilustre/users/qingmei.cui/newmdt/sanger/7.pop/hapmap-hyplotype/plink_ld/select.snp(.ld)
        """
        select_snp_path = self.work_dir + "/select_snp/" + self.output_name + ".select.snp"
        cmd = "{} --file {}".format(self.plink_path, self.option('pop_path'))
        cmd += " --r2 --ld-snp-list {}".format(select_snp_path)
        cmd += " --ld-window-kb 500 --ld-window 99999"
        cmd += " --out {}".format(select_snp_path)
        self.logger.info('开始运行plink:{}提取'.format(self.output_name))
        # self.run_cmd(cmd, "plink_{}".format(str(self.output_name).lower()))
        command = self.add_command("plink_{}".format(str(self.output_name).lower()), cmd, ignore_error=True).run()
        # command = self.run_cmd(cmd, "plink_" + str(self.output_name).lower())
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("plink完成！")
        elif command.return_code == 1:
            if self.check_rerun():
                self.logger.info("返回码为1，重新运行一次")
                command.rerun()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("重新运行一次成功！")
                else:
                    self.set_error("plink出错！", code="34503601")
                    raise Exception("plink出错！")
            else:
                self.set_error("plink出错！", code="34503602")
                raise Exception("plink出错！")

    def run_ld_region(self):
        """
        第三步可选
        perl ld_region.pl -input ../3.plink_ld/select.snp.ld -output select.region
        """
        select_snp_ld_path = self.work_dir + "/select_snp/" + self.output_name + ".select.snp.ld"
        cmd = "{} {}".format(self.perl_path, self.ld_region_path)
        cmd += " -input {}".format(select_snp_ld_path)
        cmd += " -R2 {}".format(self.option('ld_r2'))
        cmd += " -output {}".format(self.work_dir + "/" + self.output_name + ".select.region")
        self.logger.info('开始运行ld_region:{}注释统计:'.format(self.output_name))
        self.run_cmd(cmd, "ld_region_" + str(self.output_name).lower())

    def run_gwas_region(self):
        """
        第三步可选
        perl gwas.region.pl -input MVP.TRAIT2.GLM.csv -output region -distance 200000 -threshold 0.05 -adjust 1
        """
        cmd = "{} {}".format(self.perl_path, self.gwas_region_path)
        cmd += " -input {}".format(self.option('csv_path'))
        cmd += " -distance {}".format(self.option('distance'))
        cmd += " -threshold {}".format(self.threshold)
        cmd += " -adjust {}".format(self.adjust)
        cmd += " -output {}".format(self.work_dir + "/" + self.output_name + ".select.region")
        self.logger.info('开始运行gwas_region:{}注释统计:'.format(self.output_name))
        self.run_cmd(cmd, "gwas_region_{}".format(str(self.output_name).lower()))

    def run_region_stat(self):
        """
        第四步,用region.stat.pl统计第三步的基因组个数，p值，snp/indel个数
        """
        cmd = "{} {}".format(self.perl_path, self.region_stat_path)
        cmd += " -region {}".format(self.work_dir + "/" + self.output_name + ".select.region")
        cmd += " -p {}".format(self.option('csv_path'))
        cmd += " -vcf {}".format(self.option('recode_vcf_path').prop['path'])
        cmd += " -ann {}".format(self.anno_summary_path)
        cmd += " -out {}".format(self.output_dir + "/" + self.output_name + ".select.region.xls")
        self.logger.info('开始运行region_stat:{}统计:'.format(self.output_name))
        self.run_cmd(cmd, "region_stat_{}".format(str(self.output_name.lower())))

    def run_ldheatmap_vcfsplit(self):
        """
        第五步
        perl LDheatmap.vcfsplit.pl -i pop.recode.vcf -r  select.region -o ./
        """
        if not os.path.exists(os.path.join(self.work_dir, self.vcfsplit_name)):
            os.mkdir(os.path.join(self.work_dir, self.vcfsplit_name))
        cmd = "{} {} -i {} -r {} -o {}".format(self.perl_path, self.ldheatmap_vcfsplit_path,
                                               self.option('recode_vcf_path').prop['path'],
                                               (self.work_dir + "/" + self.output_name + ".select.region"),
                                               os.path.join(self.work_dir, self.vcfsplit_name))
        self.logger.info('开始运行ldheatmap_vcfsplit:{}统计:'.format(self.output_name))
        self.run_cmd(cmd, "ldheatmap_vcfsplit_{}".format(str(self.output_name.lower())))

    def run_get_data(self):
        """
        获得生成单倍体图的数据
        """
        cmd_list = []
        vcf_list = os.path.join(os.path.join(self.work_dir, self.vcfsplit_name), "vcf.list")
        with open(vcf_list, 'r')as fr:
            lines = fr.readlines()
            for line in lines:
                temp = line.strip().split('/')
                chr_name = temp[-1].strip().split('.')
                vcf_file_path = os.path.join(os.path.join(self.work_dir, self.vcfsplit_name), (chr_name[0] + '.LD.vcf'))
                with open(vcf_file_path, 'r')as r:
                    lines = r.readlines()
                    count = len(lines)
                    if count > 2:
                        self.name_list.append(chr_name[0])
        for i in self.name_list:
            if not os.path.exists(os.path.join(self.output_dir, i)):
                os.mkdir(os.path.join(self.output_dir, i))
                vcf_path = os.path.join(os.path.join(self.work_dir, self.vcfsplit_name), (i + '.LD.vcf'))
                out_path = os.path.join(os.path.join(self.output_dir, i), (i + '.data'))
                cmd = "{} {} -i {} -o {}".format(os.path.join(self.config.SOFTWARE_DIR, self.perl_path),
                                                 self.LDheatmap_vcf2matrix_path, vcf_path, out_path)
                cmd_list.append(cmd)
        self.run_cmd_more(cmd_list, "get_data")

    def run_heatmap(self):
        """
        生成热图。
        """
        cmd_list = []
        for i in self.name_list:
            input_path = os.path.join(os.path.join(self.output_dir, i), (i + '.data'))
            with open(input_path, "r") as f:
                lines = f.readlines()
                try:
                    line = lines[0]
                    tmp = line.strip().split("\t")
                    if len(tmp) > 1:
                        output_path = os.path.join(os.path.join(self.output_dir, i), i)
                        cmd = "{} {} --input {} --output {}".format(os.path.join(self.config.SOFTWARE_DIR, self.R_path),
                                                                    self.LDheatmap_path, input_path, output_path)
                        cmd_list.append(cmd)
                except:
                    pass
        self.run_cmd_more(cmd_list, "heatmap")

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("{}运行失败".format(cmd_name))
            raise Exception("{}运行失败".format(cmd_name))

    def run_cmd_more(self, cmd_list, cmd_name):
        """
        将多个cmd命令并行执行
        """
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.para_fly, cmd_file, 2, wrong_cmd)
        command = self.add_command("more_" + cmd_name, cmd_more, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}完成！".format(cmd_name))
        else:
            self.set_error("{}出错！".format(cmd_name), code="34503602")
            raise Exception("{}出错！".format(cmd_name))

    def set_db(self):
        """
        单倍体图导表
        """
        time_start = time.time()
        region_id = ObjectId(self.option("main_id"))
        self.logger.info("开始sg_heatmap的导表！")
        api = self.api.api("dna_evolution.hapmap")
        task_id = self.option('task_id')
        api.add_sg_haploid(main_id=region_id, file_path=self.option("file_path"), analysis=self.analysis,
                           trait=self.trait, task_id=task_id, area=self.name_list)
        time_end = time.time()
        self.logger.info("设置sg_heatmap的导表成功！")
        self.logger.info("导表花费{}s".format(time_end-time_start))
        self.end()

    def linkfile(self, orifile_path, targetfile_path):
        if os.path.exists(targetfile_path):
            os.remove(targetfile_path)
        else:
            os.link(orifile_path, targetfile_path)

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))

    def run(self):
        """
        运行
        :return:
        """
        vcf_list = os.path.join(os.path.join(self.work_dir, self.vcfsplit_name), "vcf.list")
        super(RegionHapmapTool, self).run()
        self.run_region_snp()
        self.run_plink()
        if self.option('ld_r2'):
            self.run_ld_region()
        else:
            self.run_gwas_region()
        self.run_region_stat()
        self.run_ldheatmap_vcfsplit()
        self.logger.info("----------------end---------------")
        self.logger.info(vcf_list)
        with open(vcf_list, "r")as fr:
            lines = fr.readlines()
            if lines == []:
                self.logger.info("文件为空!")
                self.logger.info("----------------1---------------")
                self.end()
            else:
                self.logger.info("----------------2---------------")
                self.run_get_data()
                self.run_heatmap()
                self.set_db()
                # self.end()
