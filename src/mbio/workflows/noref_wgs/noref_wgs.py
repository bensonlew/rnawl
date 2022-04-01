# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'

"""无参变异检测基础工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import os
import datetime
import json
import time
import shutil


class NorefWgsWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        version = 1.0.0
        lasted modifed by HONGDONG 20180418
        写在前面：所有的样本名必须是下划线，不能有中划线（A8_10），有不同批次的时候（A8_10-1）
        """
        self._sheet = wsheet_object
        super(NorefWgsWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'in_fastq', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # 输入的fq文件夹
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 分组方案表--我们这边导入
            {"name": 'enzyme_method', 'type': 'string', "default": 'single'},
            {"name": "enzyme_site", "type": "string"},
            {"name": "analysis_method", "type": "string", "default": "ipyrad"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.api_base = self.api.api("noref_wgs.api_base")   # 导表的基础模块
        self.norefwgsbase = self.api.api("noref_wgs.norefwgs_base")   # 导表的基础模块
        self.fastp = self.add_module("noref_wgs.fastp")  # fastp_质控
        self.ipyrad = self.add_module("noref_wgs.ipyrad")  # ipyrad
        self.uniform = self.add_module("noref_wgs.uniform")  # uniform均一化
        self.stacks = self.add_module("noref_wgs.stacks")  # stacks
        self.step.add_steps("qc_stat", "uniform")
        self.enzyme_method = 'GBS' if self.option('enzyme_method').lower() != "single" else 'RAD'
        self.target_dir = ""  # 文件上传到磁盘对应的路径
        self.fastq_list = ''
        self.group_file = ''
        self.result_path = ''
        self.cutseq = ''

    def check_options(self):
        if not self.option("in_fastq").is_set:
            raise OptionError("please input raw data fastq!", code="15500313")
        if self.option("enzyme_method"):
            if self.option('enzyme_method') not in ['single', 'double']:
                raise OptionError("enzyme_method必须是single or double！", code="15500314")
        else:
            raise OptionError("必须要输入enzyme_method参数！", code="15500315")
        if self.option("analysis_method"):
            if self.option("analysis_method") not in ['ipyrad', 'stacks']:
                raise OptionError('分析方法必须是ipyrad or stacks！', code="15500316")
        else:
            raise OptionError('必须输入analysis_method参数！', code="15500317")
        if self.option("analysis_method") == 'ipyrad' and not self.option('enzyme_site'):
            raise OptionError("当分析方法是ipyrad的时候，必须要输入酶切位点参数！", code="15500318")
        return True

    def run_fastp(self):
        """
        计算所有fastq文件的qc相关文件，并进行统计
        :return:
        """
        self.fastp.set_options({
            "sample_list": self.fastq_list
        })
        self.fastp.on("end", self.set_output, "qc_stat")
        self.fastp.on("start", self.set_step, {'start': self.step.qc_stat})
        self.fastp.on("end", self.set_step, {'end': self.step.qc_stat})
        self.fastp.run()

    def run_uniform(self):
        self.uniform.set_options({
            "sample_list": self.fastp.output_dir + "/fastq.list",
            "analysis_method": self.option("analysis_method"),
            "enzyme_method": self.enzyme_method,
            "uniform_length": 120
        })
        self.uniform.on("end", self.set_output, "uniform")
        self.uniform.on("start", self.set_step, {'start': self.step.uniform})
        self.uniform.on("end", self.set_step, {'end': self.step.uniform})
        self.uniform.run()

    def run_ipyrad(self):
        self.ipyrad.set_options({
            "fastq_dir": self.uniform.output_dir,
            "enzyme_method": self.enzyme_method,
            "cutseq": self.cutseq,
        })
        self.ipyrad.on("end", self.set_output, "ipyrad")
        self.ipyrad.on("end", self.set_step, {'end': self.step.ipyrad})
        self.ipyrad.on("start", self.set_step, {'start': self.step.ipyrad})
        self.ipyrad.run()

    def run_stacks(self):
        options = {
            "fastq_list": self.uniform.output_dir + "/fq.list"
        }
        if self.option('group').is_set:
            self.logger.info("分组方案文件：{}".format(self.group_file))
            options.update({"group_list": self.group_file})
        self.stacks.set_options(options)
        self.stacks.on("end", self.set_output, "stacks")
        self.stacks.on("start", self.set_step, {'start': self.step.stacks})
        self.stacks.on("end", self.set_step, {'end': self.step.stacks})
        self.stacks.run()

    def run(self):
        self.get_target_dir()
        self.fastq_list = self.norefwgsbase.fastp_file_list(self.option("in_fastq").prop['path'],
                                                            os.path.join(self.output_dir, "fastp.list"))
        self.norefwgsbase.add_sg_task(self._sheet.member_id, self._sheet.member_type, self._sheet.cmd_id)
        if self.option('group').is_set:
            self.make_group_file()
        self.fastp.on("end", self.run_uniform)
        if self.option('analysis_method').lower() == 'ipyrad':
            self.step.add_steps('ipyrad')
            self.uniform.on("end", self.run_ipyrad)
            # self.ipyrad.on('end', self.end)
        else:
            self.step.add_steps('stacks')
            self.uniform.on("end", self.run_stacks)
            # self.stacks.on('end', self.end)
        self.run_fastp()
        super(NorefWgsWorkflow, self).run()
        # self.start_listener()
        # self.fire("start")
        # self.end()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def move2outputdir(self, olddir, newname):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        :param olddir: 初始路径
        :param newname: 目标路径，可以自定义
        :return:
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="15500303")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        # self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}到{}移动耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        """
        递归移动文件或者文件到指定路径
        :param old_file: 初始路径
        :param new_file: 目的路径
        :return:
        """
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file_ in os.listdir(old_file):
                file_path = os.path.join(old_file, file_)
                new_path = os.path.join(new_file, file_)
                self.move_file(file_path, new_path)

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "qc_stat":
            self.move2outputdir(obj.output_dir, self.output_dir + "/01.fastq_qc")
        if event['data'] == "uniform":
            self.move2outputdir(obj.output_dir, self.output_dir + "/02.uniform")
        if event['data'] == "ipyrad":
            self.move2outputdir(obj.output_dir, self.output_dir + "/02.ipyrad")
            self.end()
        if event['data'] == "stacks":
            self.move2outputdir(obj.output_dir, self.output_dir + "/02.stacks")
            self.end()

    def send_files(self):
        # self.remove_file()
        repaths = [
            [".", "", "基础分析结果文件夹",0,"350001"],
            ["01.fastq_qc", "", "数据质控结果目录",0,"350002"],
            ["01.fastq_qc/json_dir", "", "原始数据统计结果目录",0,"350003"]
        ]
        if self.option('analysis_method').lower() == 'ipyrad':
            repaths.extend([
                ["02.ipyrad", "", "ipyrad分析结果目录",0,"350004"]
            ])
        else:
            repaths.extend([
                ["02.stacks", "", "stacks分析结果目录",0,"350005"],
                ["02.stacks/consensus_stat", "", "consensus分析结果目录",0,"350006"],
                ["02.stacks/genotype", "", "genotype分析结果目录",0,"350007"],
                ["02.stacks/stack_stat", "", "stacks分析结果目录",0,"350008"]
            ])
        regexps = [
            [r"01.fastq_qc/json_dir/.*\.json", "", "各个样本的统计数据结果文件",0,"350009"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # for i in self.get_upload_files():
        #     self.logger.info('upload file:{}'.format(str(i)))

    def run_api(self):
        """
        对所有结果进行导表操作
        :return:
        """
        self.import_specimen_info()
        self.import_tag_info()
        self.import_snp_stat()
        self.set_samples()
        self.import_files_paths()
        self.remove_file()
        self.send_files()

    def import_specimen_info(self):
        self.logger.info("开始进行样本信息导表")
        self.norefwgsbase.add_sg_specimen()
        self.logger.info("样本信息导表成功！")
        if self.option("group").is_set:
            self.logger.info("开始进行样本分组导表")
            self.norefwgsbase.add_sg_specimen_group(self.option("group").prop['path'])
            self.logger.info("样本分组导表成功！")
        self.logger.info("开始进行样本质控信息导表以及碱基含量，碱基错误分布导表")
        self.norefwgsbase.add_sg_specimen_qc(self.output_dir + "/01.fastq_qc/json_dir")
        self.logger.info("进行样本质控信息导表以及碱基含量，碱基错误分布导表成功！")
        self.norefwgsbase.add_sg_software(self._sheet.id)

    def import_tag_info(self):
        """
        酶切片段聚类分析--tags信息统计导表
        以及consensus信息统计导表
        :return:
        """
        enzyme_api = self.api.api("noref_wgs.enzyme_cluster")
        self.logger.info("开始进行聚类结果统计导表")
        main_id = self.norefwgsbase.add_sg_cluster_tag(self._sheet.id, self._sheet.project_sn,
                                                       {"analysis_method": self.option('analysis_method').lower()})
        tag_stat_path = self.result_path + "/ustacks/tags.stat" \
            if self.option('analysis_method') != 'ipyrad' else self.result_path + "/tag_stat.xls"
        enzyme_api.add_sg_cluster_tag_stat(main_id, tag_stat_path, self.option('analysis_method').lower())
        self.logger.info("进行聚类结果统计导表成功--将进行tag深度分布导表")
        dep_path = "{}/tag_depth.xls".format(self.result_path) if \
            self.option("analysis_method") == 'ipyrad' else "{}/consensus_stat/tag_dep".format(self.result_path)
        self.norefwgsbase.add_tag_dep(self._sheet.id, dep_path, main_id, self.option('analysis_method'))
        self.logger.info("tag深度分布导表成功！--将进行consensus信息统计导表！")
        main_id = self.norefwgsbase.add_sg_cluster_consensus(self._sheet.id, self._sheet.project_sn, {"value": "10%"})
        consensus_stat = "{}/consensus_stat.xls".format(self.result_path) if \
            self.option("analysis_method") == 'ipyrad' else "{}/consensus_stat/consensus_stat.txt"\
            .format(self.result_path)
        self.norefwgsbase.add_sg_cluster_consensus_stat(main_id, consensus_stat)
        self.logger.info("consensus信息统计导表成功！--将进行consensus的聚类覆盖度分布图导表！")
        consensus_coverage = "{}/consensus_coverage.xls".format(self.result_path) if \
            self.option("analysis_method") == 'ipyrad' else "{}/consensus_stat/consensus_coverage.txt" \
            .format(self.result_path)
        self.norefwgsbase.add_sg_cluster_consensus_bar(main_id, consensus_coverage)
        self.logger.info("consensus的聚类覆盖度分布图导表成功！")

    def import_snp_stat(self):
        """
        snp位点统计导表
        :return:
        """
        main_id = self.norefwgsbase.add_sg_snp_call({"analysis_method": self.option('analysis_method').lower()})
        snp_call_stat = "{}/snp_stat.xls".format(self.result_path) if \
            self.option("analysis_method") == 'ipyrad' else "{}/stack_stat/snp.stat".format(self.result_path)
        self.norefwgsbase.add_sg_snp_call_stat(main_id, snp_call_stat)
        self.logger.info("snp数据统计导表成功")
        snp_depth = "{}/depth_snp.xls".format(self.result_path) if \
            self.option("analysis_method") == 'ipyrad' else "{}/consensus_stat/depth_snp.xls".format(self.result_path)
        self.norefwgsbase.add_sg_snp_depth(main_id, snp_depth)
        self.logger.info("snp深度分布导表成功")
        snp_density = "{}/tag_snp.xls".format(self.result_path) if \
            self.option("analysis_method") == 'ipyrad' else "{}/consensus_stat/tag_snp.xls".format(self.result_path)
        self.norefwgsbase.add_snp_density(main_id, snp_density)
        self.logger.info("snp密度分布导表成功")

    def import_files_paths(self):
        if self.option("analysis_method") == 'ipyrad':
            vcf_path = self.target_dir + "/02.ipyrad/data.vcf"
            populations_tag = self.target_dir + "/02.ipyrad/populations.tag"
        else:
            vcf_path = self.target_dir + "/02.stacks/genotype/vcf_convert/pop.recode.vcf"
            populations_tag = self.target_dir + "/02.stacks/genotype/tag_generate/populations.tag"
        file_paths = {
            "pop_final_vcf": vcf_path,
            "populations_tag": populations_tag
        }
        self.api_base.update_db_record("sg_task", {"task_id": self._sheet.id}, file_paths)

    def set_samples(self):
        samples = self.norefwgsbase.find_all_sample()
        with open(self.output_dir + "/info.list", "w") as w:
            w.write("PID\t{}\n".format(samples))
            w.write("BID\t{}\n".format(samples))

    def end(self):
        self.run_api()
        self.logger.info("全部导表完成！")
        super(NorefWgsWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        if self.option('enzyme_site'):
            enzyme_site = self.option('enzyme_site').split('_')
            if len(enzyme_site) == 2:
                self.cutseq = ', '.join(enzyme_site)
            else:
                self.cutseq = enzyme_site[0]
        self.target_dir = self._sheet.output.rstrip('/')
        self.result_path = "{}/02.ipyrad/".format(self.output_dir) if \
            self.option("analysis_method") == 'ipyrad' else "{}/02.stacks/".format(self.output_dir)

    def remove_file(self):
        """
        删除不需要上传到磁盘的文件
        :return:
        """
        rm_list = list()
        rm_list.append(self.output_dir + "/01.fastq_qc/fastq")
        rm_list.append(self.output_dir + "/01.fastq_qc/fastq.list")
        rm_list.append(self.output_dir + "/02.uniform")
        rm_list.append(self.output_dir + "/02.stacks/cstacks")
        rm_list.append(self.output_dir + "/02.stacks/sstacks")
        rm_list.append(self.output_dir + "/fastp.list")
        rm_list.append(self.output_dir + "/02.stacks/stack_stat/snp.depth")
        rm_list.append(self.output_dir + "/02.stacks/stack_stat/snp.diff.pdf")
        rm_list.append(self.output_dir + "/02.stacks/stack_stat/snp.diff.png")
        rm_list.append(self.output_dir + "/02.stacks/genotype/genotype")
        if self.option("analysis_method") == "stacks":
            for file_ in os.listdir(self.output_dir + "/02.stacks/ustacks"):
                if file_ != 'tags.stat':
                    rm_list.append(self.output_dir + "/02.stacks/ustacks/" + file_)
        for files in rm_list:
            if os.path.isfile(files):
                os.remove(files)
                self.logger.info("删除文件{}成功！".format(files))
            elif os.path.isdir(files):
                code = os.system("rm -rf {}".format(files))
                if code != 0:
                    self.logger.info("删除文件夹{}失败！".format(files))
                else:
                    self.logger.info("删除文件夹{}成功！".format(files))
            else:
                self.logger.info("文件{}不存在，不用删除！".format(files))

    def make_group_file(self):
        """
        整理我们的分组方案为符合流程运行的格式
        :return:
        """
        self.group_file = self.norefwgsbase.make_group_file(self.option('group').prop['path'], self.work_dir)
