# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

'''医学检验所-亲子鉴定数据拆分流程'''
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import re
from biocluster.wpm.client import worker_client as WC
import datetime
import json
from mainapp.models.mongo.submit.paternity_test_mongo import PaternityTest as PT
from bson import ObjectId
from biocluster.config import Config

class PtDatasplitWorkflow(Workflow):
    """
    名称：亲子鉴定数据拆分流程
    作用：完成下机数据的拆分，合并以及分组(wq\ws\undetermined)
    author：zhouxuan
    last_modified: 2017.04.26
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PtDatasplitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "message_table", "type": "infile", "format": "paternity_test.tab"},  # 拆分需要的数据表
            # {"name": "data_dir", "type": "infile", "format": "paternity_test.tab"},  # 拆分需要的下机数据压缩文件夹
            {"name": "data_dir", "type": "string"},

            {"name": "family_table", "type": "infile", "format": "paternity_test.tab"},  # 亲子鉴定的家系表
            {"name": "customer_table", "type": "infile", "format": "paternity_test.tab"},  # 产前筛查家系表

            {"name": "pt_data_split_id", "type": "string"},
            {"name": "member_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.data_split = self.add_tool("paternity_test.data_split")
        self.merge_fastq = self.add_tool("paternity_test.merge_fastq")
        self.set_options(self._sheet.options())
        self.tools = []
        self.sample_name_wq = []
        self.sample_name_ws = []
        self.sample_name_un = []
        self.data_dir = ''
        self.wq_dir = ''
        self.ws_dir = ''
        self.un_dir = ''
        self.done_wq = ''
        self.done_ws = ''
        self.done_data_split = ''
        self.ws_single = ''
        self.api_read_tab = self.api.tab_file
        self.check_file = self.api.sg_paternity_test

    def check_options(self):
        '''
        检查参数设置
        '''

        if not self.option("message_table"):
            raise OptionError("缺少拆分需要的数据表")
        else:
            self.message_table = os.path.basename(self.option('message_table').prop['path'])
            self.message_no_time = ('_').join(self.message_table.split('_')[0:-1])
            data_name = self.option("data_dir").split('/')[-2]
            if self.message_no_time != data_name:
                raise OptionError('拆分表的文件名不正确，必须和下机数据文件夹的名称一致')
        if not self.option("data_dir"):
            raise OptionError("缺少拆分需要的下机数据")
        if not self.option("family_table"):
            if not self.option("customer_table"):
                raise OptionError("缺少家系表不进行任何分析")
        return True

    def run_data_split(self):
        self.data_split.set_options({
            "message_table": self.option('message_table'),
            "data_dir": self.option('data_dir'),
            'ws_single': self.ws_single,
        })
        if self.ws_single == 'false':
            if self.pt_sample_name:   # add by hd 20180129 解决掉不是单端，但是没有亲子样本的，导致运行失败的问题
                self.data_split.on('end', self.run_merge_fastq_wq)
            else:
                self.data_split.on('end', self.run_merge_fastq_ws)
        else:
            self.data_split.on('end', self.run_merge_fastq_ws)
        self.data_split.run()

    def db_customer(self):  # 不做任何判断全部都是导表
        db_customer = self.api.pt_customer
        if self.option("family_table").is_set:
            self.logger.info("开始导入pt家系表")
            db_customer.add_pt_customer(main_id=self.option('pt_data_split_id'),
                                        customer_file=self.option('family_table').prop['path'])
            self.logger.info("pt家系表导入完成")

        # self.logger.info('更新胎儿为重送样时相对应的家系表中的受理日期')  # modify 20170706
        # time = os.path.basename(self.option('message_table').prop['path']).split('-')[0]
        # year = time[0:4]
        # mon = time[4:6]
        # day = time[6:]
        # report_time = datetime.datetime(int(year), int(mon), int(day), 0, 0)
        # accept_time = report_time - datetime.timedelta(days=3)  # 拆分表格的日期为上机日期不是分析日期所以要少减一日
        # if len(str(accept_time.month)) == 1:
        #     ti = str(accept_time.year) + '-0' + str(accept_time.month)
        # else:
        #     ti = str(accept_time.year) + '-' + str(accept_time.month)
        # if len(str(accept_time.day)) == 1:
        #     ti = ti + '-0' + str(accept_time.day)
        # else:
        #     ti = ti + '-' + str(accept_time.day)
        # self.logger.info('time:{}'.format(ti))
        # with open(self.option('message_table').prop['path'], 'r') as m:
        #     for line in m:
        #         line = line.strip().split('\t')
        #         # 如果是胎儿重上机不更新信息依旧用之前的信息（重送样或者爸爸妈妈的信息）
        #         if re.match('WQ([0-9]{8,})-(S)(.*)(T)([0-9])', line[3]):
        #             continue
        #         else:
        #             if re.match('WQ([0-9]{8,})-(SC)(.*)', line[3]):  # 胎儿重送样的样本名称为-这个是系统中生成的一般会标记为SC1
        #                 family_id = line[3].split('-')[0]
        #                 self.logger.info('存在重送样的胎儿样本——{}'.format(line[3]))
        #                 db_customer.update_pt_family(family_id, ti)
        # self.logger.info('更新胎儿为重送样时相对应的家系表中的受理日期完成')  # modify 20170706

        self.logger.info("导入样本类型信息")
        db_customer.add_sample_type(self.option('message_table').prop['path'])

        if self.option("customer_table").is_set:
            self.logger.info("开始导入nipt家系表")
            self.api_nipt = self.api.nipt_analysis
            file = self.option('customer_table').prop['path']
            self.api_nipt.nipt_customer(file)
            self.logger.info("nipt家系表导入完成")

    def judge_sample_type(self, file_path):
        n = 0
        self.pt_sample_name = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if re.match('WQ([0-9]{2,})-(M|F|S)(.*)', line[3]):
                    self.pt_sample_name.append(line[3])
                if re.match('(.*)WQ(.+)', line[3]):
                    n += 1
                elif line[2] != "nipt":   # add by hd 20180129 修复产筛单端与双端判断错误
                    n += 1
                else:
                    continue
        if n == 0:
            self.ws_single = 'true'
        else:
            self.ws_single = 'false'

    def judge_sample_name(self):  # 判断样本是否重命名 20170721
        for i in self.pt_sample_name:
            x = self.api_read_tab.tab_exist(i)  # 判断是否重名
            if x:
                self.logger.error('请确认{}样本是否重名'.format(i))
                raise Exception('请确认{}样本是否重名'.format(i))
            else:
                continue

    # def family_search(self):  # 判断样本是否存在于家系表中(包括重上机样本和一般样本)
    #     paternity_sample = []
    #     pt_customer = self.api.pt_customer
    #     tab_file = self.api.tab_file
    #     with open(self.option('message_table').prop['path'], 'r') as m:
    #         for line in m:
    #             line = line.strip().split('\t')
    #             if re.match('WQ([0-9]{8,})-(M|F|S)(.*)', line[3]):
    #                 paternity_sample.append(line[3])
    #     for i in paternity_sample:  # 把所有样本照常入库
    #         tab_file.add_pt_tab(i, self.option('pt_data_split_id'))
    #     self.family_id = tab_file.family_unanalysised()  # tuple
    #     check_sample = []
    #     for i in range(len(self.family_id)):
    #         check_sample.append(self.family_id[i][0])
    #         check_sample.append(self.family_id[i][1])
    #     for i in check_sample:
    #         m = re.match('WQ([0-9]{8,})-(M|F)(.*)', i)  # 家系表里面只有爸爸和妈妈，所以此处只匹配父本和母本的名称
    #         pt_number = 'WQ' + m.group(1)
    #         if 'T' in m.group(3):
    #             n = re.match('(.*)-T(.*)', m.group(3))
    #             name = m.group(2) + n.group(1)
    #         else:
    #             name = m.group(2) + m.group(3)
    #         if 'M' in i:
    #             result = pt_customer.family_search(family_id=pt_number, mom_id_=name)
    #         else:
    #             result = pt_customer.family_search(family_id=pt_number, dad_id_=name)
    #         if result == 'False':
    #             self.logger.info('{}该样本命名有问题，家系表信息中不存在该样本'.format(i))
    #             raise Exception('{}该样本命名有问题，家系表信息中不存在该样本'.format(i))
    #         else:
    #             self.logger.error('{}该样本存在于家系表中'.format(i))

    def get_sample(self):
        self.sample_name_wq = []
        self.sample_name_ws = []
        self.sample_name_un = []
        self.data_dir = self.data_split.output_dir + "/MED"
        sample_name = os.listdir(self.data_dir)
        for j in sample_name:
            p = re.match('Sample_WQ([0-9]{2,})-(M|F|S)(.*)', j)  # 20170703 修改匹配规则
            q = re.match('Sample_WS([0-9]{2,})(.*)', j)
            if p:
                self.sample_name_wq.append(j)
            elif q:
                self.sample_name_ws.append(j)
            else:
                self.sample_name_un.append(j)

    def run_merge_fastq_wq(self):
        self.get_sample()
        n = 0
        self.tools = []
        self.wq_dir = os.path.join(self.output_dir, "wq_dir")
        if not os.path.exists(self.wq_dir):
            os.mkdir(self.wq_dir)
        for i in self.sample_name_wq:
            merge_fastq = self.add_tool("paternity_test.merge_fastq")
            self.logger.info(i)
            merge_fastq.set_options({
                "sample_dir_name": i,
                "data_dir": self.data_dir,
                "result_dir": self.wq_dir,
            })
            self.tools.append(merge_fastq)
            n += 1
        if len(self.tools) > 1:
            if len(self.sample_name_ws) == 0 and len(self.sample_name_un) != 0:
                self.on_rely(self.tools, self.run_merge_fastq_un)
            elif len(self.sample_name_ws) == 0 and len(self.sample_name_un) == 0:
                self.on_rely(self.tools, self.end)
            else:
                self.on_rely(self.tools, self.run_merge_fastq_ws)
        else:
            if len(self.sample_name_ws) == 0 and len(self.sample_name_un) != 0:
                self.tools[0].on('end', self.run_merge_fastq_un)
            elif len(self.sample_name_ws) == 0 and len(self.sample_name_un) == 0:
                self.tools[0].on('end', self.end)
            else:
                self.tools[0].on('end', self.run_merge_fastq_ws)
        for tool in self.tools:
            tool.run()

    def run_merge_fastq_ws(self):
        if self.sample_name_ws == []:
            self.get_sample()
        if self.option("family_table").is_set and self.ws_single != 'true':
            self.run_wq_wf()  # 启动亲子鉴定流程和导表工作
        n = 0
        self.tools = []
        self.ws_dir = os.path.join(self.output_dir, "ws_dir")
        if not os.path.exists(self.ws_dir):
            os.mkdir(self.ws_dir)
        for i in self.sample_name_ws:
            merge_fastq = self.add_tool("paternity_test.merge_fastq")
            self.logger.info(i)
            merge_fastq.set_options({
                "sample_dir_name": i,
                "data_dir": self.data_dir,
                "result_dir": self.ws_dir,
                "ws_single": self.ws_single,
            })
            self.tools.append(merge_fastq)
            n += 1
        if len(self.tools) > 1:
            if len(self.sample_name_un) != 0:
                self.on_rely(self.tools, self.run_merge_fastq_un)
            else:
                self.on_rely(self.tools, self.end)
        else:
            if len(self.sample_name_un) != 0:
                self.tools[0].on('end', self.run_merge_fastq_un)
            else:
                self.tools[0].on('end', self.end)
        for tool in self.tools:
            tool.run()

    def run_wq_wf(self):  # 亲子鉴定流程
        if self.wq_dir == '':
            self.logger.info("wq_dir不存在，不能进行亲子鉴定流程")
        else:
            self.logger.info("给pt_batch传送数据路径")
            mongo_data = [
                ('batch_id', ObjectId(self.option('pt_data_split_id'))),
                ("type", "pt"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("status", "start"),
                ("member_id", self.option('member_id'))
            ]
            main_table_id = PT().insert_main_table('sg_analysis_status', mongo_data)
            update_info = {str(main_table_id): 'sg_analysis_status'}
            update_info = json.dumps(update_info)
            if self.done_data_split == "true":
                value = "False"
            else:
                value = "True"
            data = {
                'stage_id': 0,
                'UPDATE_STATUS_API': self._update_status_api(),
                "id": 'pt_batch' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
                "type": "workflow",
                "name": "paternity_test.pt_dedup",
                "instant": False,
                "IMPORT_REPORT_DATA": True,
                "IMPORT_REPORT_AFTER_END": False,
                "options": {
                    "member_id": self.option('member_id'),
                    "fastq_path": self.wq_dir,
                    "cpu_number": 8,
                    "ref_fasta": Config().SOFTWARE_DIR + "/database/human/hg38.chromosomal_assembly/ref.fa",
                    "targets_bedfile": Config().SOFTWARE_DIR + "/database/human/pt_ref/snp.chr.sort.3.bed",
                    "ref_point": Config().SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda",
                    "err_min": 9,  # 11
                    "batch_id": self.option('pt_data_split_id'),
                    "dedup_num": 10,
                    "update_info": update_info,
                    "direct_get_path": value
                }
            }
            WC().add_task(data)
            self.logger.info("亲子鉴定数据拆分结束，pt_batch流程开始")
        self.done_wq = "true"

    def run_ws_wf(self):
        if self.ws_dir == '':
            self.logger.info("ws_dir不存在，不能进行产前筛查流程")
        else:
            self.logger.info("给nipt的workflow传送数据")
            mongo_data = [
                ('batch_id', ObjectId(self.option('pt_data_split_id'))),
                ("type", "nipt"),
                ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                ("status", "start"),
                ("member_id", self.option('member_id'))
            ]
            main_table_id = PT().insert_main_table('sg_analysis_status', mongo_data)
            update_info = {str(main_table_id): 'sg_analysis_status'}
            update_info = json.dumps(update_info)
            if self.done_data_split == "true":
                nipt_value = "False"
            else:
                nipt_value = "True"
            data = {
                'stage_id': 0,
                'UPDATE_STATUS_API': self._update_status_api(),
                "id": 'nipt_batch' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
                "type": "workflow",
                "name": "nipt.nipt",
                "instant": False,
                "IMPORT_REPORT_DATA": True,
                "IMPORT_REPORT_AFTER_END": False,
                "options": {
                    "fastq_path": self.ws_dir,
                    "batch_id": self.option('pt_data_split_id'),
                    'member_id': self.option('member_id'),
                    "bw": 10,
                    "bs": 1,
                    "ref_group": 1,
                    "update_info": update_info,
                    "single": self.ws_single,
                    'sanger_type': self.option('data_dir').split(":")[0],
                    "direct_get_path": nipt_value,
                }
            }
            WC().add_task(data)
            self.logger.info("亲子鉴定数据产筛部分结束，nipt_batch流程开始")
        self.done_ws = "true"

    def _update_status_api(self):
        name = self.option('data_dir').split(":")[0]
        if name == "sanger" or name == "i-sanger":
            return 'pt.med_report_update'
        else:
            return 'pt.med_report_tupdate'

    def run_merge_fastq_un(self):
        if self.option("customer_table").is_set:
            self.logger.info("启动产前筛查流程")
            self.run_ws_wf()  # 进行nipt分析
        if self.done_wq != "true" and self.option("family_table").is_set:
            if self.ws_single != 'true':
                self.run_wq_wf()
        n = 0
        self.tools = []
        self.un_dir = os.path.join(self.output_dir, "undetermined_dir")
        if not os.path.exists(self.un_dir):
            os.mkdir(self.un_dir)
        for i in self.sample_name_un:
            merge_fastq = self.add_tool("paternity_test.merge_fastq")
            self.logger.info(i)
            merge_fastq.set_options({
                "sample_dir_name": i,
                "data_dir": self.data_dir,
                "result_dir": self.un_dir
            })
            self.tools.append(merge_fastq)
            n += 1
        if len(self.tools) > 1:
            self.on_rely(self.tools, self.end)
        else:
            self.tools[0].on('end', self.end)
        for tool in self.tools:
            tool.run()

    def run(self):
        """
        判断这组数据是不是已经跑过拆分了，如果数据库中已有，说明已经有wq和ws以及undetermined的路径了
        判断路径是否存在，如果存在给self.wq_dir等赋值，如果不存在，直接重跑
        :return:
        """
        self.db_customer()  # 家系表导表，不管是否做过拆分导表都进行一下
        self.judge_sample_type(self.option('message_table').prop['path'])  # 判断ws是否为单端
        db_customer = self.api.pt_customer
        if self.option('family_table').is_set:
            db_customer.family_search(self.pt_sample_name)  # 判断这些样本能组成的家系是否存在家系信息
        self.logger.info("开始检查家系")
        db_customer.get_urgency_sample(self.option('message_table').prop['path'], self.option('pt_data_split_id'))  # 用于获取加急样本，并导表
        dir_list = db_customer.get_wq_dir(self.option('data_dir').split(":")[1] + '-' + self.message_table)
        # 上述记录拆分表是为了再次拆分的时候，上传表格改变就会重新拆分
        self.logger.info(dir_list)
        if len(dir_list) == 3 and (os.path.exists(dir_list[0]) or os.path.exists(dir_list[1])):
            self.wq_dir = dir_list[0]
            self.ws_dir = dir_list[1]
            self.un_dir = dir_list[2]
            self.start_listener()
            self.end()
        else:
            self.logger.info('开始进行拆分-------------------')
            if self.option('family_table').is_set:
                self.judge_sample_name()  # 确认样本是否重命名只能在拆分的时候进行，当第二次调用的时候可能样本已经入库了 所以不能再进行检查
            self.done_data_split = "true"   # 本次workflow是否进行数据拆分，true为进行
            self.run_data_split()
            super(PtDatasplitWorkflow, self).run()

    def set_output(self, event):  # 暂时无用
        obj = event["bind_object"]
        if event['data'] == "data_split":
            self.linkdir(obj.output_dir, self.output_dir + "/data_split")
        if event['data'] == "merge_fastq":
            wq_dir = os.path.join(self.output_dir, "wq_dir")
            self.wq_dir = wq_dir
            ws_dir = os.path.join(self.output_dir, "ws_dir")
            undetermined_dir = os.path.join(self.output_dir, "undetermined_dir")
            if not os.path.exists(wq_dir):
                os.mkdir(wq_dir)
            if not os.path.exists(ws_dir):
                os.mkdir(ws_dir)
            if not os.path.exists(undetermined_dir):
                os.mkdir(undetermined_dir)
            file_name = os.listdir(obj.output_dir)
            m = re.match('WQ([0-9].*)-(.*)', file_name[0])  # wq
            n = re.match('WS(.*)', file_name[0])  # ws
            if m:
                self.linkdir(obj.output_dir, wq_dir)
            else:
                if n:
                    self.linkdir(obj.output_dir, ws_dir)
                else:
                    self.linkdir(obj.output_dir, undetermined_dir)

    def end(self):
        self.logger.info("医学流程数据拆分结束")
        if self.done_data_split == "true":  # true表示这是第一次进行拆分
            self.logger.info("开始导入拆分结果路径")
            db_customer = self.api.pt_customer
            db_customer.add_data_dir(self.option('data_dir').split(":")[1] + '-' + self.message_table,
                                     self.wq_dir, self.ws_dir, self.un_dir)
        if self.done_wq != "true" and self.option('family_table').is_set:
            if self.ws_single != 'true':
                self.run_wq_wf()
        if self.option('customer_table').is_set and self.done_ws != "true":
            self.run_ws_wf()
        super(PtDatasplitWorkflow, self).end()

    def linkdir(self, dirpath, dirname):  # 暂时无用
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
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
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))
