# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v3.functions import workfuncdeco
import os
import shutil
import datetime
from collections import OrderedDict
import json


class ExtractSeqWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExtractSeqWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'bam_file_list', 'type': 'infile', 'format': 'ref_rna_v3.common'},
            {'name': 'seq_type', 'type': 'string', 'default': 'PE'},  # ['PE', 'SE']
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'extract_info', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string', 'default': None},
            {'name': 'upload_offline', 'type': 'bool', 'default': True}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.download = self.add_module('ref_rna_v3.download')

    @workfuncdeco
    def check_options(self):
        for k, v in self.sheet.options().items():
            self.logger.debug('{} = {}'.format(k, v))
        if not self.option('bam_file_list').is_set:
            self.set_error("bam文件路径必须提供")
        if self.option('seq_type') not in ['PE', 'SE']:
            self.set_error('序列类型必须为PE/SE')
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
        self.run_download()
        super(ExtractSeqWorkflow, self).run()

    @workfuncdeco
    def run_download(self):
        self.step.add_steps('download')
        self.extract_info = json.loads(self.option("extract_info"), object_pairs_hook=OrderedDict)
        samples = set()
        for key in self.extract_info:
            if self.extract_info[key]:
                for value in self.extract_info[key]:
                    samples.add(value)
        s3_file_list = os.path.join(self.work_dir, 's3_file_list')
        with open(self.option('bam_file_list').prop['path'], 'r') as f, open(s3_file_list, "w") as w:
            for line in f:
                sample = line.strip().split("/")[-1].replace('.bam', '')
                if sample in samples:
                    w.write(line)
        options = {
            's3_file_list': s3_file_list
        }
        self.download.set_options(options)
        self.download.on('start', self.set_step, {'start': self.step.download})
        self.download.on('end', self.set_step, {'end': self.step.download})
        rely = list()
        if 'mapped' in self.extract_info and self.extract_info['mapped']:
            self.mapped_extract = self.add_module('ref_rna_v3.extract_seq')
            self.download.on('end', self.run_extract_mapped)
            rely.append(self.mapped_extract)
        if 'unmapped' in self.extract_info and self.extract_info['unmapped']:
            self.unmapped_extract = self.add_module('ref_rna_v3.extract_seq')
            self.download.on('end', self.run_extract_unmapped)
            rely.append(self.unmapped_extract)
        self.logger.info("Rely: {}".format(rely))
        if len(rely) == 1:
            if 'mapped' in self.extract_info and self.extract_info['mapped']:
                self.mapped_extract.on('end', self.set_output)
            else:
                self.unmapped_extract.on('end', self.set_output)
        elif len(rely) >= 2:
            self.on_rely(rely, self.set_output)
        else:
            self.download.on('end', self.set_output)
        self.download.run()

    def run_extract_mapped(self):
        self.step.add_steps('mapped_extract')
        bam_list = os.path.join(self.download.work_dir, 's3_file_list')
        samples = set(self.extract_info['mapped'])
        mapped_list = os.path.join(self.work_dir, 'mapped_list')
        with open(bam_list, 'r') as f, open(mapped_list, "w") as w:
            for line in f:
                sample = line.strip().split("/")[-1].replace('.bam', '')
                if sample in samples:
                    w.write(line)
        self.mapped_extract.set_options({
            'input_file': mapped_list,
            'seq_type': self.option('seq_type'),
            'extract_type': 'mapped',
        })
        self.mapped_extract.on('start', self.set_step, {'start': self.step.mapped_extract})
        self.mapped_extract.on('end', self.set_step, {'end': self.step.mapped_extract})
        self.mapped_extract.run()

    def run_extract_unmapped(self):
        self.step.add_steps('unmapped_extract')
        bam_list = os.path.join(self.download.work_dir, 's3_file_list')
        samples = set(self.extract_info['unmapped'])
        unmapped_list = os.path.join(self.work_dir, 'unmapped_list')
        with open(bam_list, 'r') as f, open(unmapped_list, "w") as w:
            for line in f:
                sample = line.strip().split("/")[-1].replace('.bam', '')
                if sample in samples:
                    w.write(line)
        self.unmapped_extract.set_options({
            'input_file': unmapped_list,
            'seq_type': self.option('seq_type'),
            'extract_type': 'unmapped',
        })
        self.unmapped_extract.on('start', self.set_step, {'start': self.step.unmapped_extract})
        self.unmapped_extract.on('end', self.set_step, {'end': self.step.unmapped_extract})
        self.unmapped_extract.run()

    def set_output(self):
        shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        if 'bam' in self.extract_info and self.extract_info['bam']:
            bam_dir = os.path.join(self.output_dir, "bam")
            os.mkdir(bam_dir)
            for result in os.listdir(self.download.output_dir):
                if result.replace(".bam", "") in self.extract_info['bam']:
                    os.link(os.path.join(self.download.output_dir, result), os.path.join(bam_dir, result))
        if 'mapped' in self.extract_info and self.extract_info['mapped']:
            mapped_dir = os.path.join(self.output_dir, "mapped")
            os.mkdir(mapped_dir)
            for result in os.listdir(self.mapped_extract.output_dir):
                os.link(os.path.join(self.mapped_extract.output_dir, result), os.path.join(mapped_dir, result))
        if 'unmapped' in self.extract_info and self.extract_info['unmapped']:
            unmapped_dir = os.path.join(self.output_dir, "unmapped")
            os.mkdir(unmapped_dir)
            for result in os.listdir(self.unmapped_extract.output_dir):
                os.link(os.path.join(self.unmapped_extract.output_dir, result), os.path.join(unmapped_dir, result))
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
        cmd = "scp -r -i ~/.ssh/id_rsa {} dongmei.fu@10.2.4.236:/mnt/ilustre/users/chun.luo/bam_data_for_prok_rna".format(
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
        super(ExtractSeqWorkflow, self).end()