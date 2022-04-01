# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v3.functions import workfuncdeco
import os
import shutil
import datetime
from collections import OrderedDict
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import re
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class DownloadDetailSeqWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DownloadDetailSeqWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'seq_detail_dir', "type": "infile", "format": "ref_rna_v2.common_dir"},
            {'name': 'geneset_extract', "type": "infile", "format": "ref_rna_v2.common"},
            # {'name': 'seq_type', 'type': 'string', 'default': 'PE'},  # ['PE', 'SE']
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'extract_info', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None},
            {'name': 'upload_offline', 'type': 'bool', 'default': False}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.download = self.add_module('ref_rna_v3.download')
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/07 Seq_Download')
        self.download = self.add_tool("ref_rna_v3.extract_detail_seq")

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(DownloadDetailSeqWorkflow, self).send_log(data)

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        # if not self.option('bam_file_list').is_set:
        #     self.set_error("bam文件路径必须提供")
        # if self.option('seq_type') not in ['PE', 'SE']:
        #     self.set_error('序列类型必须为PE/SE')
        if not self.option("extract_info"):
            self.set_error('必须提供提取内容的相关信息')

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    @workfuncdeco
    def run(self):
        # self.get_run_log()
        self.run_extract()
        super(DownloadDetailSeqWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_seq_extract", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    @workfuncdeco
    def run_extract(self):
        self.step.add_steps('download')
        self.extract_info = json.loads(self.option("extract_info"), object_pairs_hook=OrderedDict)
        options = {
            'seq_detail_dir': self.option("seq_detail_dir"),
            'geneset_extract': self.option("geneset_extract"),
            'level': self.option("level"),
            'extract_info': self.option("extract_info")
        }
        self.download.set_options(options)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        self.download.on('end', self.set_output)
        self.download.run()



    def set_output(self):
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        for file in os.listdir(self.download.output_dir):
            shutil.move(os.path.join(self.download.output_dir, file), self.output_dir)
        if self.option("upload_offline"):
            self.transfer()
        else:
            self.end()

    def transfer(self):
        self.logger.info("开始向线下服务器传递结果文件，请耐心等待")
        time_now = datetime.datetime.now()
        target_dir = "{}_{}".format(self.option("task_id"), time_now.strftime('%Y%m%d_%H%M%S'))
        self.offline_results = os.path.join(self.work_dir, target_dir)
        if os.path.isdir(self.offline_results):
            shutil.rmtree(self.offline_results)
        os.mkdir(self.offline_results)
        for result_dir in os.listdir(self.output_dir):
            shutil.move(os.path.join(self.output_dir, result_dir), self.offline_results)
        cmd = "scp -r -i ~/.ssh/id_rsa {} dongmei.fu@192.168.10.46:/mnt/ilustre/centos7users/meng.luo/sanger_data".format(
            self.offline_results)
        try:
            code = os.system(cmd)
            if code == 0:
                self.logger.info("命令{}执行成功！".format(cmd))
            else:
                self.logger.info("命令{}执行失败！".format(cmd))
                self.set_error("向线下服务器传递数据失败")
        except:
            self.set_error("向线下服务器传递数据失败")
        self.end()

    @workfuncdeco
    def end(self):
        # if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
        #     os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        # os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["07 Seq_Download", "", "序列下载结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "序列下载结果目录", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
            ['./extract_gene_seqs', '', '基因序列', 0],
            ['./extract_trans_seqs', '', '转录本序列', 0],
            ['./extract_cds_seqs', '', 'CDS序列', 0],
            ['./extract_pep_seqs', '', '氨基酸序列', 0]
        ])
        super(DownloadDetailSeqWorkflow, self).end()