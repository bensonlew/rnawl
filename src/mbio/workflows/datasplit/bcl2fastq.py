# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""数据拆分一次拆分以及一次拆分后的原始数据统计"""

import os
import re
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.api.to_file.datasplit import *


class Bcl2fastqWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Bcl2fastqWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_path", "type": "string"},  # 下机数据路径
            {"name": "sample_sheet", "type": "infile", "format": "datasplit.sample_sheet"},  # csv文件
            {"name": "bases_mask", "type": "string"},  # 测序模式，hiseq4000 为y151,i6nn,y151; miseq为y301,i6,y301
            {"name": "barcode_mismatch", "type": "int", "default": 0},  # --barcode-mismatches，barcode错配数
            {"name": "lib_id", "type": "string"},  # 文件，第一列文库名，第二列对应此次分析的文库id,第三列样本id(逗号分隔)
            {'name': 'status_id', "type": 'string'},  # 分析记录表id
            {"name": "run_type", "type": "string", "default": "auto"},
            {"name": "update_info", "type": "string"},
            # zj 运行类型，auto:运行完文库拆分workflow,运行二次拆分及质控workflow;bcl2fastq:只运行文库拆分workflow
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.step.add_steps('bcl2fastq', 'raw_data_stat', 'sample_split_qc')
        self.bcl2fastq = self.add_tool("datasplit.bcl2fastq_v3")
        self.raw_data_stat = self.add_module("datasplit.raw_data_stat")
        self.lib_id = {}  # {文库名：文库ID}
        self.sample_id = {}  # {文库id：样本ID}
        self.api_datasplit = self.api.api('datasplit_new')

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("data_path"):
            raise OptionError("缺少下机数据路径")
        # base_path = os.path.join(self.option("data_path"), "Data/Intensities/BaseCalls/")
        # if not os.path.exists(base_path):
        #     raise OptionError("下机数据里缺少{}文件夹".format(base_path))
        if not self.option("sample_sheet").is_set:
            raise OptionError("缺少csv表")
        if not self.option("bases_mask"):
            raise OptionError("缺少测序模式")
        if re.match(r"\d+", self.option("bases_mask")):
            raise OptionError("测序模式有问题，请检查")
        if self.option("barcode_mismatch") < 0:
            raise OptionError("参数barcode_mismatch:{}不能小于0，请检查".format(self.option("barcode_mismatch")))
        if self.option("run_type") not in ["auto", "split"]:  # zj
            raise OptionError("运行类型只能是auto/split,请检查是否传错流程")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_bcl2fastq(self):
        self.logger.info("开始运行bcl2fastq")
        self.step.bcl2fastq.start()
        self.bcl2fastq.set_options({
            "data_path": self.option('data_path'),
            "sample_sheet": self.option('sample_sheet'),
            "bases_mask": self.option('bases_mask'),
            "barcode_mismatch": self.option("barcode_mismatch")
        })
        self.bcl2fastq.on('end', self.set_output, 'bcl2fastq')
        self.bcl2fastq.run()
        self.step.bcl2fastq.finish()
        self.step.update()

    def run_raw_data_stat(self):
        self.logger.info("开始进行文库原始序列统计")
        self.step.raw_data_stat.start()
        self.raw_data_stat.set_options({
            "list_file": self.work_dir + '/file_list.txt',
        })
        self.raw_data_stat.on('end', self.set_output, 'raw_data_stat')
        self.raw_data_stat.run()
        self.step.raw_data_stat.finish()
        self.step.update()

    def run_sample_split_qc(self):  # zj
        """
        触发二次拆分及质控的workflow:sample_split_qc
        """
        self.logger.info("开始进行二次拆分及质控")
        self.step.sample_split_qc.start()
        upload_dir = self._sheet.output
        # export_split_qc_params(data=self.option('status_id'), option_name="", dir_path=self.work_dir, bind_obj=None)
        export_split_qc_params_by_bcl2fastq(data=self.option('status_id'), dir_path=self.work_dir, lib_info=self.lib_id,\
                                            sample_info=self.sample_id, bcl2fastq_output=self.output_dir + '/bcl2fastq')
        self.sample_split_qc = self.add_batch(path="datasplit.sample_split_qc", ignore_error=False, batch_type="workflow",\
                               upload_dir=upload_dir, up_api_name=None, import_report_data=False, import_report_after_end=True)
        self.sample_split_qc.set_options({
            "project_params": self.work_dir + "/params.json",
            "status_id": self.option('status_id'),
            "update_info": self.option("update_info")
        })
        self.sample_split_qc.on('end', self.end)
        self.sample_split_qc.run()
        self.step.sample_split_qc.finish()
        self.step.update()

    def run(self):
        self.run_bcl2fastq()
        super(Bcl2fastqWorkflow, self).run()

    def get_path_list(self):
        """
        sample_info: 样本路径信息 {样本/文库：[R1_path,R2_path]}
        :return: file_list: 文件，第一列路径，第二列 样本/文库名, 第三列 l/r
        """
        lib_info = {}
        file_list = os.path.join(self.work_dir, 'file_list.txt')
        with open(self.option("sample_sheet").prop['path'])as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split(",")
                lib_info[tmp[1]] = tmp[2]     # 建立文库名与文件夹名称的对应关系 {文件夹名称：文库名}
        all_dir = os.listdir(self.output_dir + '/bcl2fastq')
        with open(file_list, 'w+')as fw:
            for fq_dir in all_dir:
                fq_list = os.listdir(os.path.join(self.output_dir + '/bcl2fastq', fq_dir))
                for fq in fq_list:
                    if re.search(r'_R1_', fq):
                        fq_path = os.path.join(self.output_dir + '/bcl2fastq/' + fq_dir, fq)
                        fw.write(fq_path + '\t' + fq_dir + '\t' + 'l\n')
                    if re.search(r'_R2_', fq):
                        fq_path = os.path.join(self.output_dir + '/bcl2fastq/' + fq_dir, fq)
                        fw.write(fq_path + '\t' + fq_dir + '\t' + 'r\n')
        self.run_raw_data_stat()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def link_has_dir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件一层目录文件夹移动到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        for file_dir in allfiles:
            dir_path = os.path.join(dirpath, file_dir)
            allfqs = os.listdir(dir_path)
            for fq in allfqs:
                new_path = os.path.join(newdir, file_dir)
                if not os.path.exists(new_path):
                    os.mkdir(new_path)
                oldfile = os.path.join(dir_path, fq)
                newfile = os.path.join(new_path, fq)
                if os.path.exists(newfile):
                    os.remove(newfile)
                os.link(oldfile, newfile)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event["data"] == "bcl2fastq":
            chip_name = ''
            for l in os.listdir(self.bcl2fastq.output_dir + '/Reports/html/'):
                path = os.path.join(self.bcl2fastq.output_dir + '/Reports/html/', l)
                if os.path.isdir(path):
                    chip_name = l
            if os.path.exists(self.output_dir + '/laneBarcode.html'):
                os.remove(self.output_dir + '/laneBarcode.html')
            if os.path.exists(self.output_dir + '/lane.html'):
                os.remove(self.output_dir + '/lane.html')
            os.link(self.bcl2fastq.output_dir + '/Reports/html/' + chip_name + '/all/all/all/laneBarcode.html', self.output_dir + '/laneBarcode.html')
            os.link(self.bcl2fastq.output_dir + '/Reports/html/' + chip_name + '/all/all/all/lane.html', self.output_dir + '/lane.html')
            self.link_has_dir(self.bcl2fastq.output_dir + '/Fastq', self.output_dir + '/bcl2fastq')
            # if self._sheet.output.startswith("s3"):
            #     self.upload_to_s3("output", self._sheet.output)  # 提前上传，避免开始样本拆分和质控的时候还没有文件
            # else:
            #     self.add_upload_dir(self.output_dir)
            self.get_path_list()
            self.add_db_path()
            self.logger.info("设置结果目录成功")
        if event["data"] == "raw_data_stat":
            self.link_has_dir(self.raw_data_stat.output_dir, self.output_dir + '/raw_data_stat')
            self.set_db()
            self.logger.info("设置结果目录成功")
            if self.option('run_type') == 'split':
                self.end()  # zj

    def add_db_path(self):
        """
        先更新文库表中的文库路径(一个文库对应一个样本的需要更新样本表raw_path),方便进行二次拆分及质控
        :return:
        """
        self.logger.info("开始更新一次拆分结果文库路径")
        with open(self.option('lib_id'))as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                self.lib_id[tmp[0]] = tmp[1]
                ids = tmp[2].strip().split(",")
                if len(ids) == 1:
                    self.sample_id[tmp[1]] = ids[0]
        if self._sheet.output.startswith("s3"):
            s3_upload_dir = self._sheet.output + "bcl2fastq/"
        else:
            if self._sheet.client01:
                s3_upload_dir = "/mnt/ilustre/data/" + self._sheet.output.split(":")[1] + "/bcl2fastq/"
            else:
                s3_upload_dir = "/mnt/ilustre/tsanger-data/" + self._sheet.output.split(":")[1] + "/bcl2fastq/"
        self.api_datasplit.update_lib_path(self.lib_id, self.output_dir + '/bcl2fastq', self.sample_id, s3_upload_dir)
        self.logger.info("更新一次拆分结果文库路径成功")
        if self.option('run_type') == 'auto':
            self.run_sample_split_qc()  # 触发二次拆分及质控的workflow

    def set_db(self):
        """
        将一次拆分的结果导表
        :return:
        """
        self.logger.info("开始一次拆分结果导表")
        self.api_datasplit.add_flowcell_summary(self.option('status_id'), self.output_dir + '/lane.html', self.output_dir + '/laneBarcode.html')
        self.api_datasplit.add_stat(self.option('status_id'), 'raw_lib', self.lib_id,
                                    self.output_dir + '/raw_data_stat/fastq_stat/fastq_stat.xls',
                                    self.output_dir + '/raw_data_stat/fastq_dup',
                                    self.output_dir + '/raw_data_stat/fastx')
        self.api_datasplit.add_lib_raw_sample_stat(self.option("status_id"), self.option("lib_id"),
                                                   self.output_dir + "/raw_data_stat")

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(Bcl2fastqWorkflow, self).end()
