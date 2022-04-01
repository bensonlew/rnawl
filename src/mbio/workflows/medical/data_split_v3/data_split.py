# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
# created at 20171110
# lasted modifeid by hd@20200807

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from biocluster.api.file.lib.transfer import MultiFileTransfer
import os
import re
import datetime


class DataSplitWorkflow(Workflow):
    """
    医学数据拆分流程（目前包括亲子、产筛）
    lasted modified by hongdong@20180821
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
            {"name": "project_types", "type": "string", "default": "all"},  # 项目类型，用于区分拆分完之后
            # 激发那个项目的工作流
            {"name": 'indextype', 'type': 'string', 'default': 'single'},
            {"name": "datatype", "type": "string", 'default': "notwailai"},  # 是否为外来数据notwailai and iswailai
            {"name": "aliyun_path", 'type': "string"},  # 数据在aliyun上面的路径
        ]
        self.add_option(options)
        self.bclfq = self.add_tool("medical.datasplit_v3.bclfq")
        self.manage = self.add_module("medical.paternity_test_v3.med_manage")
        self.set_options(self._sheet.options())
        self.is_split = None
        self.sample_dir = None
        self.wq_sample = None  # 存储的所有的亲子样本id
        self.to_s3 = True     # 是否将拆分表传到对象存储
        self.indextype = self.option("indextype")
        self.ossscript = os.path.join(Config().SOFTWARE_DIR, "bioinfo/medical/ossutil64")
        self.osskey = os.path.join(Config().SOFTWARE_DIR, "bioinfo/medical/osskey")
        self.raw_fastq = ""

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
            "sanger_type": "sanger" if self._sheet.client == 'client01' else 'tsanger',  # 这里传递运行环境，测试还是正式，临时用
            "indextype": self.indextype 
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

    def transfer_html_file(self):
        html_name = "laneBarcode.html"
        project_bucket = Config().get_project_region_bucket(project_type="pt_v3").rstrip("/")
        tar_path = "{}/data_split/{}/" \
            .format(project_bucket, os.path.basename(str(self.option('board_batch')).rstrip('/')))

        if Config().RGW_ENABLE and self.to_s3:  # 对象存储
            transfer = MultiFileTransfer()
            transfer.add_upload(self.bclfq.output_dir + "/laneBarcode.html",
                                tar_path, base_path=os.path.dirname(self.bclfq.output_dir))
            transfer.perform()
            mongo_path = tar_path + 'output/laneBarcode.html'
        else:
            target_path = '/mnt/ilustre/users/sanger-dev/tsanger/med_data/{}/'.format(os.path.basename(
                self.option('board_batch')))
            if not os.path.exists(target_path):
                os.makedirs(target_path)
            if os.path.exists(target_path + '/' + html_name):
                os.remove(target_path + '/' + html_name)
            os.link(self.bclfq.output_dir + '/laneBarcode.html', os.path.join(target_path, html_name))
            mongo_path = "{}/".format(os.path.basename(self.option('board_batch'))) + html_name
        self.api.api('medical.paternity_test_v3.paternity_test_v3').update_datasplit_html(mongo_path,
                                                                                          self.option('batch_id'))

    def run(self):
        """
        开始运行拆分
        modified by hongdong @ 20171212
        modified by hd@20200807  增加外来数据的线上分析功能
        运行逻辑：
        首先检测该板子是否已经拆分成功了，如果拆分成功了就会跳过拆分的tool，一般情况都是没有进行拆分的，
        该情况较多运行在测试的过程中,不需要一直等拆分过程
        """
        pt_api = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        self.is_split, self.sample_dir, self.wq_sample = pt_api.find_fastq_info(
            os.path.basename(str(self.option('board_batch')).rstrip()))
        self.logger.info("is_split{}{}".format(self.is_split, self.sample_dir))
        if self.wq_sample:
            pt_api.add_sg_family(self.wq_sample)
        else:
            self.logger.info("wq_sample列表为空，不进行sg_family的导表操作！")
        if self.option("datatype") == "iswailai":
            self.is_split = True
            self.wailaidata_pro()
            self.manage.on('end', self.end)
            self.manage_works()
        else:
            if self.is_split:
                self.manage.on('end', self.end)
                self.manage_works()
            else:
                self.bclfq.on('end', self.manage_works)
                self.bclfq.on('end', self.transfer_html_file)
                self.manage.on('end', self.end)
                self.run_bclfq()
        super(DataSplitWorkflow, self).run()

    def end(self):
        self.set_output()
        super(DataSplitWorkflow, self).end()

    def wailaidata_pro(self):
        """
        外来数据处理流程，主要功能是下载fastq并完成文件名字检查+重名后，移动到指定文件夹位置
        :return:
        """
        if self.needdownlad():
            self.download_data()
        self.make_samples_name()
        self.check_sample_fastq_num()
        pt_api = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        pt_api.update_datasplit_status(self.option('batch_id'))

    def download_data(self):
        """
        下载外来数据，从aliyun上面下载
        ./ossutil64 cp -r -f  oss://delivery-data/s1500/Project_s1500g01008_6Samples_20200712_1594548853/
        Sample_R20028074-WQHZ0003-WQ2001210-M/ ./test/
        :return:
        """
        self.logger.info("开始下载文件夹：{}".format(self.option("aliyun_path")))
        if not self.option("aliyun_path"):
            self.set_error("需要提供阿里云路径！")
        cmd = "{} cp -r -f {} {} --config-file {} --output-dir {} --checkpoint-dir {} > {}"\
            .format(self.ossscript, self.option("aliyun_path").strip().rstrip("/") + "/",
                    os.path.join(self.work_dir, "ossdata"), self.osskey, os.path.join(self.work_dir, "ossutil_output"),
                    os.path.join(self.work_dir, "checkpoint_dir"),
                    os.path.join(self.work_dir, "download.log"))
        self.logger.info("cmd: {}".format(cmd))
        code = os.system(cmd)
        if code == 0:
            self.raw_fastq = os.path.join(self.work_dir, "ossdata")
            self.logger.info("下载文件夹成功：{}".format(self.option("aliyun_path")))
        else:
            self.set_error("下载文件夹失败：{}；解决方案："
                           "mnt-ilustre-users-yixuezhuanhua-raw-data文件夹中创建板号，并手动下载fastq，"
                           "下载完成后创建download_success文件！".format(self.option("aliyun_path")))

    def needdownlad(self):
        """
        检查是否需要下载阿里云上面的数据。
        逻辑是:
        1. 检查/mnt/ilustre/users/yixuezhuanhua/raw-data中是否有对应板子的信息
        2. 有板子信息，并且文件夹中有download_success文件
        3. 符合上述条件就跳过下载
        :return:
        """
        if os.path.exists(self.option('board_batch')) and os.path.isdir(self.option('board_batch')):
            if os.path.exists(os.path.join(self.option('board_batch'), "download_success")):
                self.raw_fastq = self.option('board_batch')
                self.logger.info("检测到fastq已经下载了，不用再次下载！")
                return False
        return True

    def make_samples_name(self):
        """
        重名了样本名字，并检查样本名字是否与拆分表中一致，这步检测可以避免掉报告组输入的aliyun路径是否正确
        :return:
        """
        opd = self.raw_fastq
        target_sample = os.path.join(self.work_dir, "Samples/PT/")
        if os.path.exists(target_sample):
            os.system('rm -rf {}/*'.format(target_sample))
        else:
            os.makedirs(target_sample)
        for m in os.listdir(opd):
            if re.match('Sample_*', m) and os.path.isdir(os.path.join(opd, m)):
                for n in os.listdir(os.path.join(opd, m)):
                    if not re.match(r".*\.gz\.md5$", n):
                        sample_name = '-'.join(n.split('_')[0].split('-')[2:])
                        if sample_name not in self.wq_sample:
                            self.set_error("aliyun中样本{}不再拆分表中，检查aliyun路径是否正确".format(sample_name))
                        sample_dir = os.path.join(opd, m)
                        result = re.match(r".*(_R1\.fastq\.gz|_R2\.fastq\.gz)$", n)
                        if result:
                            os.link(os.path.join(sample_dir, n),
                                    os.path.join(target_sample, sample_name + result.group(1)))
        self.logger.info("所有样本重命名成功！")
        self.sample_dir = os.path.join(self.work_dir, "Samples/")
        self.logger.info("sample_dir: {}".format(self.sample_dir))

    def check_sample_fastq_num(self):
        """
        检查样本的fastq个数是否正确，目前外送样本默认都是双端数据
        :return:
        """
        fastqs = os.listdir(os.path.join(self.sample_dir, "PT"))
        for sample in self.wq_sample:
            if "".join([sample, "_R1.fastq.gz"]) not in fastqs:
                self.set_error("样本{}缺少：{}文件--请检查从aliyun下载的数据是否缺少--同时检查阿里云中数据是否已经删除!"
                               .format(sample, "".join([sample, "_R1.fastq.gz"])))
            if "".join([sample, "_R2.fastq.gz"]) not in fastqs:
                self.set_error("样本{}缺少：{}文件--请检查从aliyun下载的数据是否缺少--同时检查阿里云中数据是否已经删除!"
                               .format(sample, "".join([sample, "_R2.fastq.gz"])))
