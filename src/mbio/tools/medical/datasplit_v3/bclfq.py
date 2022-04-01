# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# created at 20171111

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from bs4 import BeautifulSoup
from biocluster.api.file.lib.transfer import MultiFileTransfer
from collections import defaultdict
import datetime
import re
import os


class BclfqAgent(Agent):
    """
    测序下机数据拆分程序
    lasted modified by hongdong@20180821
    """
    def __init__(self, parent):
        super(BclfqAgent, self).__init__(parent)
        options = [
            {"name": "split_tab", "type": "infile", "format": "medical.tab"},
            {"name": "data_dir", "type": "infile", "format": "medical.run_dir"},
            {"name": "split_type", "type": "string"},  # PE/SE
            {"name": "Barcode_sum", "type": "outfile", "format": "medical.html"},
            {"name": "sample_dir", "type": "outfile", "format": "medical.sample_dir"},
            {"name": "batch_id", "type": "string"},
            {"name": "sanger_type", 'type': "string", "default": "sanger"},   # 临时使用
            {"name": 'indextype', 'type': 'string', 'default': 'single'}
        ]
        self.add_option(options)
        self.queue = "gada"

    def check_options(self):
        """
        参数检查
        """
        if not self.option('split_tab'):
            raise OptionError("缺少必要参数split_tab:样本拆分表")
        if not self.option('data_dir'):
            raise OptionError("缺少必要参数data_dir:下机数据文件夹")
        if self.option('split_type') not in ['PE', 'SE']:
            raise OptionError("参数split_type设置不正确")
        self.logger.info("tab表头：{}".format(self.option('split_tab').rows[0]))
        for x in ['sample_name', 'index', 'department', 'analysis_type']:
            if x not in self.option('split_tab').rows[0]:
                raise OptionError("样本拆分表中缺少表头字段：{}".format(x))

    def set_resource(self):
        """
        设置运行资源
        """
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        super(BclfqAgent, self).end()


class BclfqTool(Tool):
    def __init__(self, config):
        super(BclfqTool, self).__init__(config)
        self.bcl2fq = "bioinfo/seq/bcl2fastq2-v2.17.1.14/bin/bcl2fastq"
        self.sample_sheet = os.path.join(self.work_dir, "sample_sheet.csv")
        self.indextype = self.option("indextype")
        if self.indextype == 'single':
            # 单index 单端与双端
            # self.bases_mask = 'y76,i6n'   # 接头6碱基
            self.bases_mask = 'y76,i8'   # 接头8碱基
            if self.option('split_type') == 'PE':
                # self.bases_mask = 'y76,i6n,y76'
                self.bases_mask = 'y76,i8,y76'
        else:
            # 双index 单端与双端
            # self.bases_mask = 'y76,i6nn,nnnnnnnn'
            self.bases_mask = 'y76,i8,nnnnnnnn'
            if self.option('split_type') == 'PE':
                self.logger.info("come here!")
                # self.bases_mask = 'y76,i6nn,nnnnnnnn,y76'  # 亲子变为双index时，这里要改成y76,i6nn,i6nn,y76
                self.bases_mask = 'y76,i8,nnnnnnnn,y76'
        self.sample_names = []
        self.sample_tab = []
        self.zero_samples = []
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')

    def creat_sample_sheet(self):
        """
        创建拆分程序表--单index
        """
        idx = [self.option('split_tab').rows[0].index(i) for i in ['sample_name', 'sample_name', 'index',
                                                                   'department', "analysis_type"]]
        colums = self.option('split_tab').get_colums(idx)
        with open(self.sample_sheet, "w") as f:
            f.write("[Data],,,\nSample_ID,Sample_Name,index,Sample_Project")
            for row in zip(*colums)[1:]:
                if row[4] in ['pass']:
                    continue
                else:
                    f.write('\nSample_' + ','.join(row[:-1]))

    def creat_sample_sheet_doubleindex(self):
        """
        创建拆分程序表--双index
        """
        idx = [self.option('split_tab').rows[0].index(i) for i in ['sample_name', 'sample_name', 'index', 'index2', 'department']]
        colums = self.option('split_tab').get_colums(idx)
        with open(self.sample_sheet, "w") as f:
            f.write("[Data],,,\nSample_ID,Sample_Name,index,index2,Sample_Project")
            for row in zip(*colums)[1:]:
                f.write('\nSample_' + ','.join(row))

    def run_bclfq(self):
        """
        运行bcl2fastq
        """
        self.logger.info(self.option('data_dir').prop['path'])
        cmd = self.bcl2fq + \
            ' -i ' + self.option('data_dir').prop['path'] + '/Data/Intensities/BaseCalls/' + \
            ' -o ' + self.work_dir + \
            ' --sample-sheet ' + self.sample_sheet + \
            ' --use-bases-mask ' + self.bases_mask + \
            ' --ignore-missing-bcl ' + \
            ' -R ' + self.option('data_dir').prop['path'] + \
            ' -p 10 -r 4 -w 4 -d 2 --barcode-mismatches 0'
        self.logger.info(cmd)
        command = self.add_command("bcl2fq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行完成bcl2fastq:" + cmd)
        else:
            self.set_error("运行bcl2fastq失败！")
            raise Exception("运行bcl2fastq失败！")

    def merge_sample(self):
        """
        合并一个样本不同lane的fq
        """
        id_s = self.option('split_tab').rows[0].index('sample_name')
        self.sample_names = self.option('split_tab').colums[id_s][1:]
        sample_path = os.path.join(self.output_dir, 'Samples/')
        if os.path.exists(sample_path):
            os.system('rm -rf ' + sample_path)
        os.mkdir(sample_path)
        for sample in self.sample_names:
            sample_idx = self.sample_names.index(sample) + 1
            titles = self.option('split_tab').rows[0]
            id_d, id_t = [titles.index(i) for i in ['department', 'analysis_type']]
            sample_split_path = os.path.join(self.work_dir, self.option('split_tab').colums[id_d][sample_idx],
                                             'Sample_'+sample)
            analysis_type = self.option('split_tab').colums[id_t][sample_idx].upper()
            if re.search('DCPT', analysis_type) or re.search('WQCF', analysis_type):
                analysis_type = 'PT'
            sample_merge_dir = sample_path + "/" + analysis_type
            if not os.path.exists(sample_merge_dir):
                os.mkdir(sample_merge_dir)
            sample_merge_path1 = os.path.join(sample_merge_dir, sample + '_R1.fastq.gz')
            sample_merge_path2 = os.path.join(sample_merge_dir, sample + '_R2.fastq.gz')
            os.system('cat {}/{}_*_R1_*.fastq.gz > {}'.format(sample_split_path, sample, sample_merge_path1))
            # 修复bug，当单端的时候，只合并R1就好了，不对R2进行合并 add by hongdong 20180303 line 130
            if self.option('split_type') == 'PE':
                os.system('cat {}/{}_*_R2_*.fastq.gz > {}'.format(sample_split_path, sample, sample_merge_path2))

    def add_zj_group_to_html(self, input_html, output_html):
        """
        读取bcl2fastq生成的html文件，在结果表中加入"Zj group"（杂交组号）列, 将结果表按照lane、杂交组号、样品名排列
        :param input_html:
        :param output_html:
        :return:
        """
        idx = [self.option('split_tab').rows[0].index(i) for i in ['sample_name', 'zj_group']]
        colums = self.option('split_tab').get_colums(idx)
        zj_group = dict(zip(*colums)[1:])

        soup = BeautifulSoup(open(input_html), "html.parser")
        
        new_table = soup.new_tag("table")
        new_table["border"] = "1"
        new_table["id"] = "ReportTable"
        lane_results = defaultdict(list)
        lane_undetermined = dict()
        all_results = list()

        for child in soup.find_all("table")[2].children:
            if child == "\n":
                continue
            sample = child.contents[5].string
            lane = child.contents[1].string
            if sample == "Sample":
                new_tag = soup.new_tag("th")
                new_tag.string = "Zj group"
                child.insert(5, "\n")
                child.insert(5, new_tag)
                new_table.append(child)                 # 保存标题
            elif sample == "Undetermined":
                new_tag = soup.new_tag("td")
                new_tag.string = "unknown"
                new_tag["style"] = "font-style:italic"
                child.insert(5, "\n")
                child.insert(5, new_tag)
                lane_undetermined[lane] = child
            elif sample in zj_group:
                new_tag = soup.new_tag("td")
                new_tag.string = zj_group[sample]
                child.insert(5, "\n")
                child.insert(5, new_tag)
                lane_results[lane].append((zj_group[sample], sample, child))
            else:
                new_tag = soup.new_tag("td")
                new_tag.string = "unknown"
                new_tag["style"] = "font-style:italic"
                child.insert(5, "\n")
                child.insert(5, new_tag)
                lane_results[lane].append(("unknown", sample, child))
        
        for lane in sorted(lane_results.keys()):
            results = sorted(lane_results[lane], key=lambda x: (x[0], x[1]))      # 按照杂交组号、样品名排序
            all_results.extend([i[2] for i in results])
            all_results.append(lane_undetermined[lane])
        for row in all_results:
            new_table.append(row)
        soup.find_all("table")[2].replace_with(new_table)                         # 替换原先的表

        with open(output_html, "w") as out:
            out.write(soup.prettify())

    def check_seqsnum(self):
        """
        检查序列数
        """
        self.option('Barcode_sum').check()
        self.logger.info("2")
        self.sample_tab = self.option('Barcode_sum').tab_list[2]
        # self.logger.info("self.sample_tab{}".format(self.sample_tab))
        self.zero_samples = []
        for smp in self.sample_tab:
            if len(smp) == 12:   # modified by hd 20180103 smp有异常的情况，长度小于8会报错
                if smp[8] == '0':
                    self.zero_samples.append(smp[2])
            elif len(smp) == 7:
                if smp[-1] == '0':
                    self.zero_samples.append(smp[2])
            else:
                self.logger.info("拆分结果表有值为空")
        self.logger.info(self.zero_samples)
        if 5 > len(self.zero_samples) > 0:
            self.set_error("发现序列数为0样本:" + ','.join(self.zero_samples))
            raise Exception("发现序列数为0样本：" + ','.join(self.zero_samples))
        elif len(self.zero_samples) >= 5:
            self.set_error("发现序列数为0样本:{}等！".format(','.join(self.zero_samples[0:5])))
            raise Exception("发现序列数为0样本:{}等！".format(','.join(self.zero_samples[0:5])))

    def set_output(self):
        """
        设置输出文件
        """
        chip_name = ''
        for l in os.listdir(self.work_dir + '/Reports/html/'):
            if re.match(r'^\w+$', l):
                chip_name = l
        # html_name = os.path.basename(str(self.option('data_dir').prop['path']).rstrip()) + datetime.datetime.now().\
        #     strftime("%Y%m%d%H%M%S") + "_laneBarcode.html"
        if os.path.exists(self.output_dir + '/laneBarcode.html'):
            os.remove(self.output_dir + '/laneBarcode.html')
        self.add_zj_group_to_html(self.work_dir + '/Reports/html/' + chip_name + '/all/all/all/laneBarcode.html',
                                  self.output_dir + '/laneBarcode.html')
        # os.link(self.work_dir + '/Reports/html/' + chip_name + '/all/all/all/laneBarcode.html',
        #         self.output_dir + '/laneBarcode.html')
        # tar_path = "{}/data_split/{}/"\
        #     .format(Config().get_project_region_bucket(project_type="pt_v3").rstrip("/"),
        #             os.path.basename(str(self.option('data_dir').prop['path']).rstrip('/')))
        # if Config().RGW_ENABLE:  # 对象存储
        #     transfer = MultiFileTransfer()
        #     transfer.add_upload(self.output_dir + "/laneBarcode.html",
        #                         tar_path, base_path=os.path.dirname(self.output_dir))
        #     transfer.perform()
        #     mongo_path = tar_path + 'output/laneBarcode.html'
        # else:
        #     if self.option('sanger_type') == 'sanger':
        #         target_path = '/mnt/ilustre/data/rerewrweset/MEDfiles/medical.data_split/'
        #     else:
        #         target_path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/medical.data_split/'
        #     if os.path.exists(target_path + '/' + html_name):
        #         os.remove(target_path + '/' + html_name)
        #     os.link(self.output_dir + '/laneBarcode.html', os.path.join(target_path, html_name))
        #     mongo_path = "rerewrweset/MEDfiles/medical.data_split/" + html_name
        # self.api.api('medical.paternity_test_v3.paternity_test_v3').update_datasplit_html(mongo_path,
        #                                                                                   self.option('batch_id'))
        self.option('Barcode_sum').set_path(self.work_dir + '/Reports/html/' + chip_name +
                                            '/all/all/all/laneBarcode.html')  # 使用原始的html而不是add_zj_group_to_html方法修改后的
        # self.option('Barcode_sum').set_path(self.output_dir + '/laneBarcode.html')
        self.option('sample_dir').set_path(self.output_dir + '/Samples/')

    def check_file_number(self):
        basecall = self.option('data_dir').prop['path'] + '/Data/Intensities/BaseCalls/'
        file_numbers = ['330', '324', '322', '320', '294', '256', '177', '172', '171', '168', '167', '166', '146',
                        '122', '84', '170', '338']  # 从以往下机BaseCalls文件夹L001 - L004中统计出来的文件个数
        for i in ["L001", "L002", "L003", "L004"]:
            lane_dir = os.path.join(basecall, i)
            if not os.path.isdir(lane_dir):
                self.set_error("{}:文件夹不存在！".format(lane_dir))
                raise Exception("{}:文件夹不存在！".format(lane_dir))
            count = len(os.listdir(lane_dir))
            if str(count) not in file_numbers:
                self.set_error("{}文件夹下文件数为{}，与预先设定的文件数({})不符"
                               .format(i, str(count), '|'.join(file_numbers)))

    def run(self):
        """
        运行
        modified by hongdong @ 20171212
        """
        super(BclfqTool, self).run()
        pt_api = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        self.check_file_number()
        # 因为当前亲子还是上的是单index，所以依旧使用creat_sample_sheet（）去生成拆分表，
        # 当亲子变成双端后，这里要修改成使用creat_sample_sheet_doubleindex()
        if self.indextype == 'single':
            self.creat_sample_sheet()
        else:
            self.creat_sample_sheet()
            # self.creat_sample_sheet_doubleindex()
        self.run_bclfq()
        self.merge_sample()
        self.set_output()
        self.check_seqsnum()
        pt_api.add_datasplit_info(os.path.basename(str(self.option('data_dir').prop['path']).rstrip()),
                                  self.option('sample_dir').prop['path'])
        pt_api.update_datasplit_status(self.option('batch_id'))
        self.end()
