# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# created at 20171110

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import datetime


class DataSplitWorkflow(Workflow):
    """
    医学数据拆分流程（目前包括亲子、产筛）
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DataSplitWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "split_tab", "type": "infile", "format": "medical.tab"},
            {"name": "batch_id", "type": "string"},  # 数据拆分主表id
            {"name": "board_batch", "type": "string"},  # 数据拆分板号
            {"name": "split_type", "type": "string"},  # PE/SE
            {"name": "member_id", "type": "string"},  # 用户ID
            {"name": "update_info", "type": "string"},  # 状态更新
            {"name": "project_types", "type": "string", "default": "all"}  # 项目类型，用于区分拆分完之后
            # 激发那个项目的工作流
        ]
        self.add_option(options)
        self.bclfq = self.add_tool("medical.datasplit.bclfq")
        self.manage = self.add_module("medical.med_manage")
        self.set_options(self._sheet.options())
        self.is_split = None
        self.sample_dir = None
        self.wq_sample = None  # 存储的所有的亲子样本id
        self.logger.info("设置好了参数")

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("split_tab"):
            raise OptionError("缺少拆分需要的数据表")
        return True

    def run_bclfq(self):
        """
        运行bcl2fq tool
        """
        self.bclfq.set_options({
            "split_tab": self.option('split_tab'),
            # "data_dir": "/mnt/clustre/upload/nextseq1/" + self.option('board_batch'),
            "data_dir": self.option('board_batch'),  # 这里是路径，其他地方是板号
            "split_type": self.option('split_type'),
            "batch_id": self.option('batch_id'),
            "sanger_type": "sanger" if self._sheet.client == 'client01' else 'tsanger'  # 这里传递运行环境，测试还是正式，临时用
        })
        self.logger.info("开始第一个tool：bclfq")
        self.bclfq.run()

    def manage_works(self):
        """
        开始各产品流程分析
        """
        self.manage.set_options({
            "member_id": self.option('member_id'),
            "batch_id": self.option('batch_id'),
            "board_batch": self.option('board_batch'),
            "sample_tab": self.option('split_tab'),
            "sample_dir": self.bclfq.option("sample_dir") if not self.is_split else self.sample_dir,
            "project_types": self.option("project_types")
        })
        self.manage.run()

    def set_output(self):
        if os.path.exists(self.output_dir + "/data"):
            os.remove(self.output_dir + "/data")
        self.linkdir(self.bclfq.output_dir,
                     self.output_dir + "/data-{}".format(os.path.basename(str(self.option('board_batch')).rstrip())))
        # html_name = os.path.basename(str(self.option('board_batch')).rstrip()) + datetime.datetime.now().\
        #     strftime("%Y%m%d%H%M%S") + "_laneBarcode.html"
        # os.link(self.bclfq.output_dir + "/laneBarcode.html", self.output_dir + "/" + html_name)
        # repaths = [
        #     [html_name, "", "laneBarcode统计文件"]
        # ]
        # self.add_upload_dir(self.output_dir)
        # self.add_upload_dir(self.output_dir).add_relpath_rules(repaths)
        # self.api.api('medical.paternity_test_v2').update_datasplit_html("medical.data_split/" + html_name,
        #                                                                 self.option('batch_id'))

    def linkdir(self, dirpath, dirname):
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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        """
        开始运行拆分
        modified by hongdong @ 20171212
        运行逻辑：
        首先检测该板子是否已经拆分成功了，如果拆分成功了就会跳过拆分的tool，一般情况都是没有进行拆分的，
        该情况较多运行在测试的过程中,不需要一直等拆分过程
        """
        self.is_split, self.sample_dir, self.wq_sample = self.api.api('medical.paternity_test_v2').find_fastq_info(
            os.path.basename(str(self.option('board_batch')).rstrip()))
        self.logger.info("is_split{}{}".format(self.is_split, self.sample_dir))
        if self.wq_sample:
            self.api.api('medical.paternity_test_v2').add_sg_family(self.wq_sample)
        else:
            self.logger.info("wq_sample列表为空，不进行sg_family的导表操作！")
        if self.is_split:
            self.manage.on('end', self.end)
            self.manage_works()
        else:
            self.bclfq.on('end', self.manage_works)
            self.manage.on('end', self.end)
            self.run_bclfq()
        super(DataSplitWorkflow, self).run()

    def end(self):
        self.set_output()
        super(DataSplitWorkflow, self).end()

