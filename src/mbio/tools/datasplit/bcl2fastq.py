# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""bcl2fastq 工具 """
import os
import errno
import re
import xml.etree.ElementTree as ET
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from mbio.packages.datasplit.miseq_split import code2index


class Bcl2fastqAgent(Agent):
    """
    bcl2fastq
    version 2.17
    """
    def __init__(self, parent=None):
        super(Bcl2fastqAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'}  # 样本拆分信息表
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检测
        """
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 5
        self._memory = ''


class Bcl2fastqTool(Tool):
    """
    """
    def __init__(self, config):
        super(Bcl2fastqTool, self).__init__(config)
        self._version = 1.0
        self.bcl2fastq_path = "bioinfo/seq/bcl2fastq2-v2.17.1.14/bin/bcl2fastq"
        self.option('sample_info').get_info()
        self.base_mask = ""
        # 由于挂载点的问题，需要将下机文件的路径前缀由/mnt/ilustre换成/mnt/clustre
        self.option('sample_info').prop["file_path"] = re.sub("^\/mnt\/ilustre", "/mnt/clustre", self.option('sample_info').prop["file_path"])
        my_path = self.option('sample_info').prop["file_path"]
        run_info = os.path.join(my_path, "RunInfo.xml")
        if not os.path.exists(my_path):
            raise Exception("数据路径{}不存在".format(my_path))
        if not os.path.exists(run_info):
            raise Exception("RunInfo.xml文件缺失")
        if not self.option('sample_info').check_parent_repeat():
            self.set_error("父样本中的index重复")
            raise Exception("父样本中的index重复")
        if not self.option('sample_info').check_child_repeat():
            self.set_error("属于同一个父样本的子样本中的index重复")
            raise Exception("属于同一个父样本的子样本中的index重复")

    def create_sample_sheet(self):
        """
        生成SampleSheet.csv, 供程序bcl2fastq使用
        """
        sample_sheet_path = os.path.join(self.work_dir, "soft_input", "SampleSheet.csv")
        flag = 0
        # 检测是单barcode模式还是双barcode模式
        for p_id in self.option('sample_info').prop["parent_ids"]:
            index2 = self.option('sample_info').parent_sample(p_id, "index2")
            if index2 != "" and index2 != "null" and index2 != "NULL" and index2 != "NA" and index2 is not None:
                flag = 1
            else:
                if flag == 1:
                    name = self.option('sample_info').parent_sample(p_id, "library_name")
                    raise Exception("文库：{}缺失index2".format(name))
        if flag == 0:
            self.logger.info("检测到单barcode模式， 开始生成sample_sheet表")
            self.create_single_code_sheet(sample_sheet_path)
        else:
            self.logger.info("检测到双barcode模式， 开始生成sample_sheet表")
            self.create_double_code_sheet(sample_sheet_path)

    def create_single_code_sheet(self, sample_sheet_path):
        """
        单barcode模式
        """
        info_path = os.path.join(self.option('sample_info').prop["file_path"], "RunInfo.xml")
        e = ET.parse(info_path).getroot()
        # 解析"RunInfo.xml"文件，获取机器模式中的index长度
        machineIndexLen = 0
        for ele in e.iter("Read"):
            if ele.attrib["Number"] == "2" and ele.attrib['IsIndexedRead'] == "Y":
                machineIndexLen = int(ele.attrib['NumCycles'])
            if ele.attrib["Number"] == "1" and ele.attrib['IsIndexedRead'] == "N":
                sqMethod1 = int(ele.attrib['NumCycles'])
            if ele.attrib["Number"] == "3" and ele.attrib['IsIndexedRead'] == "N":
                sqMethod2 = int(ele.attrib['NumCycles'])
        if machineIndexLen == 0:
            raise Exception("RunInfo.xml解析出错, 无法解析上机时的index长度，请检查{}是否正确".format(info_path))

        # 遍历所有文库， 获取所有文库中最短的barcode长度
        minIndexLen = 9999
        myCode2Index = dict()  # 为了减少数据库的查询次数
        for p_id in self.option('sample_info').prop["parent_ids"]:
            index_code = self.option('sample_info').parent_sample(p_id, "index")
            index = code2index(index_code)[0]
            myCode2Index[index_code] = index
            if len(index) > machineIndexLen:
                myLibraryName = self.option('sample_info').parent_sample(p_id, "library_name")
                raise Exception("文库{}的index长度大于机器模式中的index长度".format(myLibraryName))
            if len(index) < minIndexLen:
                minIndexLen = len(index)

        # 生成bask_mask
        iContent = str(minIndexLen) + "n" * int(machineIndexLen - minIndexLen)
        self.base_mask = "y{},i{},y{}".format(sqMethod1, iContent, sqMethod2)

        with open(sample_sheet_path, 'w+') as w:
            head = "[Data],,,\nLane,Sample_ID,index,Sample_Project\n"
            w.write(head)
            for p_id in self.option('sample_info').prop["parent_ids"]:
                lane = self.option('sample_info').parent_sample(p_id, "lane")
                name = self.option('sample_info').parent_sample(p_id, "sample_id")
                index = self.option('sample_info').parent_sample(p_id, "index")
                index = myCode2Index[index]
                index = index[0:minIndexLen]
                project = self.option('sample_info').parent_sample(p_id, "library_type")
                if project is None:
                    project = "undefine"
                line = lane + "," + name + "," + index + "," + project + "\n"
                w.write(line)

    def create_double_code_sheet(self, sample_sheet_path):
        """
        双barcode模式
        """
        info_path = os.path.join(self.option('sample_info').prop["file_path"], "RunInfo.xml")
        e = ET.parse(info_path).getroot()
        machineIndexLen1 = 0
        machineIndexLen2 = 0
        for ele in e.iter("Read"):
            if ele.attrib["Number"] == "1" and ele.attrib['IsIndexedRead'] == "N":
                sqMethod1 = int(ele.attrib['NumCycles'])
            if ele.attrib["Number"] == "4" and ele.attrib['IsIndexedRead'] == "N":
                sqMethod2 = int(ele.attrib['NumCycles'])
            if ele.attrib["Number"] == "2" and ele.attrib['IsIndexedRead'] == "Y":
                machineIndexLen1 = int(ele.attrib['NumCycles'])
            if ele.attrib["Number"] == "3" and ele.attrib['IsIndexedRead'] == "Y":
                machineIndexLen2 = int(ele.attrib['NumCycles'])
        if machineIndexLen1 == 0 or machineIndexLen2 == 0:
            raise Exception("RunInfo.xml解析出错, 无法解析上机时的index长度，请检查{}是否正确".format(info_path))

        minIndexLen1 = 9999
        minIndexLen2 = 9999
        myCode2Index = dict()
        for p_id in self.option('sample_info').prop["parent_ids"]:
            index_code1 = self.option('sample_info').parent_sample(p_id, "index")
            index1 = code2index(index_code1)[0]
            myCode2Index[index_code1] = index1
            index_code2 = self.option('sample_info').parent_sample(p_id, "index2")
            index2 = code2index(index_code2)[0]
            myCode2Index[index_code2] = index2
            if len(index1) > machineIndexLen1:
                myLibraryName = self.option('sample_info').parent_sample(p_id, "library_name")
                raise Exception("文库{}的index1长度大于机器模式中的index1长度".format(myLibraryName))
            if len(index2) > machineIndexLen2:
                myLibraryName = self.option('sample_info').parent_sample(p_id, "library_name")
                raise Exception("文库{}的index2长度大于机器模式中的index2长度".format(myLibraryName))
            if len(index1) < minIndexLen1:
                minIndexLen1 = len(index1)
            if len(index2) < minIndexLen2:
                minIndexLen2 = len(index2)

        iContent1 = str(minIndexLen1) + "n" * int(machineIndexLen1 - minIndexLen1)
        iContent2 = str(minIndexLen2) + "n" * int(machineIndexLen2 - minIndexLen2)
        self.base_mask = "y{},i{},i{},y{}".format(sqMethod1, iContent1, iContent2, sqMethod2)

        with open(sample_sheet_path, 'w+') as w:
            head = "[Data],,,,\nLane,Sample_ID,index,index2,Sample_Project\n"
            w.write(head)
            for p_id in self.option('sample_info').prop["parent_ids"]:
                lane = self.option('sample_info').parent_sample(p_id, "lane")
                name = self.option('sample_info').parent_sample(p_id, "sample_id")
                index1 = self.option('sample_info').parent_sample(p_id, "index")
                index2 = self.option('sample_info').parent_sample(p_id, "index2")
                index1 = myCode2Index[index1]
                index2 = myCode2Index[index2]
                index1 = index1[0:minIndexLen1]
                index2 = index1[0:minIndexLen2]
                project = self.option('sample_info').parent_sample(p_id, "library_type")
                if project is None:
                    project = "undefine"
                line = lane + "," + name + "," + index1 + "," + index2 + "," + project + "\n"
                w.write(line)

    def bcl2fastq(self):
        """
        运行bcl2fastq
        """
        basecall = os.path.join(self.option('sample_info').prop['file_path'], "Data/Intensities/BaseCalls")
        output_dir = os.path.join(self.work_dir, "soft_output")
        sample_sheet_path = os.path.join(self.work_dir, "soft_input", "SampleSheet.csv")
        self.logger.debug(basecall)
        bcl2fastqstr = (self.bcl2fastq_path + " --input-dir " + basecall + " --runfolder-dir " +
                        self.option('sample_info').prop['file_path'] + " --output-dir " + output_dir +
                        " --sample-sheet " + sample_sheet_path +
                        " --barcode-mismatches " + str(self.option('sample_info').prop["index_missmatch"])
                        + " --use-bases-mask " + self.base_mask +
                        " -p 8 -d 8"
                        )
        if self.option('sample_info').prop["ignore_missing_bcl"]:
            bcl2fastqstr = bcl2fastqstr + " --ignore-missing-bcl"
        log_path = os.path.join(self.work_dir, "sh.log")
        with open(log_path, 'wb') as w:
            my_str = "bcl2fastq"
            w.write(my_str.center(79, "#"))
            w.write("\n")
            w.write(bcl2fastqstr)
            w.write("\n")
        bcl2fastqcmd = self.add_command("bcl2fastq", bcl2fastqstr)
        self.logger.info("开始运行bcl2fastq")
        self.logger.debug(bcl2fastqstr)
        bcl2fastqcmd.run()
        self.wait(bcl2fastqcmd)
        self.logger.debug(bcl2fastqcmd.return_code)
        if bcl2fastqcmd.return_code == 0:
            self.logger.info("bcl2fastq运行成功")
        else:
            self.logger.info("bcl2fastq运行失败")
            raise OSError("bcl2fastq运行失败")

    def make_ess_dir(self):
        """
        为软件bcl2fastq的运行创建必要的运行目录
        """
        input_dir = os.path.join(self.work_dir, "soft_input")
        output_dir = os.path.join(self.work_dir, "soft_output")
        dir_list = [input_dir, output_dir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def run(self):
        super(Bcl2fastqTool, self).run()
        self.make_ess_dir()
        self.create_sample_sheet()
        self.bcl2fastq()
        self.end()
