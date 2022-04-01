# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# created at 20171111

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import datetime
import re
import os


class BclfqAgent(Agent):
    """
    测序下机数据拆分程序
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
            {"name": "sanger_type", 'type': "string", "default": "sanger"}   # 临时使用
        ]
        self.add_option(options)

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
        self.bases_mask = 'y76,i6n'
        if self.option('split_type') == 'PE':
            self.bases_mask = 'y76,i6n,y76'
        self.sample_names = []
        self.sample_tab = []
        self.zero_samples = []
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')

    def creat_sample_sheet(self):
        """
        创建拆分程序表
        """
        # colums = [self.option('split_tab').colums[x] for x in ['sample_name', 'sample_name', 'index', 'department']]
        idx = [self.option('split_tab').rows[0].index(i) for i in ['sample_name', 'sample_name', 'index', 'department']]
        colums = self.option('split_tab').get_colums(idx)
        with open(self.sample_sheet, "w") as f:
            f.write("[Data],,,\nSample_ID,Sample_Name,index,Sample_Project")
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
            if re.search('DCPT', analysis_type):
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
            # r1_list = []
            # r2_list = []
            # file_list = os.listdir(sample_split_path)
            # for p in file_list:
            #     if re.match('{}_(.*)_R1_([0-9].*).fastq.gz'.format(sample), p):
            #         p = os.path.join(sample_split_path, p)
            #         r1_list.append(p)
            #     elif re.match('{}_(.*)_R2_([0-9].*).fastq.gz'.format(sample), p):
            #         p = os.path.join(sample_split_path, p)
            #         r2_list.append(p)
            #     else:
            #         pass
            # r1_list.sort()
            # cat_files1 = ' '.join(r1_list)
            # self.logger.info(cat_files1)
            # os.system('cat ' + cat_files1 + ' > ' + sample_merge_path1)
            # if self.option('split_type') == 'PE':
            #     r2_list.sort()
            #     cat_files2 = ' '.join(r2_list)
            #     os.system('cat ' + cat_files2 + ' > ' + sample_merge_path2)

    def check_seqsnum(self):
        """
        检查序列数
        """
        self.option('Barcode_sum').check()
        self.logger.info("2")
        self.sample_tab = self.option('Barcode_sum').tab_list[2]
        self.logger.info("self.sample_tab{}".format(self.sample_tab))
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
        if self.option('sanger_type') == 'sanger':
            target_path = '/mnt/ilustre/data/rerewrweset/MEDfiles/medical.data_split/'
        else:
            target_path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/medical.data_split/'
        chip_name = ''
        for l in os.listdir(self.work_dir + '/Reports/html/'):
            if re.match(r'^\w+$', l):
                chip_name = l
        html_name = os.path.basename(str(self.option('data_dir').prop['path']).rstrip()) + datetime.datetime.now().\
            strftime("%Y%m%d%H%M%S") + "_laneBarcode.html"
        self.logger.info("-")
        if os.path.exists(self.output_dir + '/laneBarcode.html'):
            os.remove(self.output_dir + '/laneBarcode.html')
        if os.path.exists(target_path + '/' + html_name):
            os.remove(target_path + '/' + html_name)
        self.logger.info("--")
        # os.symlink(self.work_dir + '/Reports/html/' + chip_name + '/all/all/all/laneBarcode.html',
        #            self.output_dir + '/laneBarcode.html')
        os.link(self.work_dir + '/Reports/html/' + chip_name + '/all/all/all/laneBarcode.html',
                self.output_dir + '/laneBarcode.html')
        self.logger.info("---")
        os.link(self.output_dir + '/laneBarcode.html', os.path.join(target_path, html_name))
        self.logger.info("-----")
        self.api.api('medical.paternity_test_v2').update_datasplit_html("medical.data_split/" + html_name,
                                                                        self.option('batch_id'))
        self.option('Barcode_sum').set_path(self.output_dir + '/laneBarcode.html')
        self.option('sample_dir').set_path(self.output_dir + '/Samples/')
        self.logger.info("1")

    def check_file_number(self):
        basecall = self.option('data_dir').prop['path'] + '/Data/Intensities/BaseCalls/'
        file_numbers = ['330', '324', '322', '320', '294', '256', '177', '172', '171', '168', '167', '166', '146',
                        '122', '84']  # 从以往下机BaseCalls文件夹L001 - L004中统计出来的文件个数
        for i in ["L001", "L002", "L003", "L004"]:
            lane_dir = os.path.join(basecall, i)
            if not os.path.exists(lane_dir):
                self.set_error("{}:文件夹不存在！".format(lane_dir))
                raise Exception("{}:文件夹不存在！".format(lane_dir))
            count = len(os.listdir(lane_dir))
            if str(count) not in file_numbers:
                self.set_error("{}文件夹下文件数为{}，与预先设定的文件数({})不符".format(i, str(count), '|'.join(file_numbers)))
                raise Exception("{}文件夹下文件数为{}，与预先设定的文件数({})"
                                "不符".format(i, str(count), '|'.join(file_numbers)))

    def run(self):
        """
        运行
        modified by hongdong @ 20171212
        """
        super(BclfqTool, self).run()
        self.check_file_number()
        self.creat_sample_sheet()
        self.run_bclfq()
        self.merge_sample()
        self.set_output()
        self.check_seqsnum()
        self.api.api('medical.paternity'
                     '_test_v2').add_datasplit_info(os.path.basename(str(self.option('data_dir').prop['path']).rstrip())
                                                    , self.option('sample_dir').prop['path'])
        self.api.api('medical.paternity_test_v2').update_datasplit_status(self.option('batch_id'))
        self.end()
