# -*- coding: utf-8 -*-
# __author__: HONGDONG
# modified: 20190312

import re
import os
import json
from biocluster.agent import Agent
from biocluster.tool import Tool
from bson.objectid import ObjectId


class SeqDealAgent(Agent):
    """
    SeqDeal分析
    """
    def __init__(self, parent):
        super(SeqDealAgent, self).__init__(parent)
        options = [
            {"name": "sample_region_ids", "type": "string", "required": True},
            {"name": "scaffold_ids", "type": "string", "required": True},
            {"name": "seq_path", "type": "infile", 'format': 'bsa.dir', "required": True},  # 序列的路径
            {"name": "update_info", "type": "string", "required": True},
            {"name": "main_id", "type": "string", "required": True},
            {"name": "possname", "type": "string", "required": True},
            {"name": "target_path", 'type': "string"},
            {"name": "project_type", "type": "string", 'default': "dna_wgs_v2"}
        ]
        self.add_option(options)

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SeqDealAgent, self).end()


class SeqDealTool(Tool):
    def __init__(self, config):
        super(SeqDealTool, self).__init__(config)

    def per_make_params(self):
        sample_region_ids = self.option("sample_region_ids").split(',')
        scaffold_ids = self.option('scaffold_ids').split(',')
        possname = json.loads(self.option('possname'))   # "{\"region01\": \"chr3:-20000\"}"
        for i in range(0, len(sample_region_ids)):
            tt = sample_region_ids[i].split(':')[0].split("_")
            name = "_".join(tt[:-1])   # 样本名字 YC_bulk
            region = sample_region_ids[i].split('_')[-1]  # 区域名字 region01
            # name, region = sample_region_ids[i].split('_')
            if region in possname.keys() or region == 'unmapping':
                if region == 'unmapping':
                    denovo_seq_file = os.path.join(self.option('seq_path').prop['path'],
                                                   '_'.join([name, "unmapping.denovo.scafSeq"]))
                else:
                    file_name = possname[region]
                    denovo_seq_file = os.path.join(self.option('seq_path').prop['path'],
                                                   '_'.join([name, file_name + ".denovo.scafSeq"]))
                    if not os.path.exists(denovo_seq_file):
                        if re.match('(.*):(.*)-$', file_name):  # chr1:1-
                            region = file_name.split('-')[0]
                        elif re.match("(.*):-(.*)$", file_name):  # chr1:-2000
                            region = ':'.join(
                                [file_name.split(':')[0], '-'.join(['1', file_name.split(':')[1].split('-')[1]])])
                        elif re.match('(.*):-$', file_name):  # chr1:-
                            region = file_name.split(':-')[0]
                        else:
                            self.set_error("region{}不合法！".format(region))
                        denovo_seq_file = os.path.join(self.option('seq_path').prop['path'],
                                                       '_'.join([name, region + ".denovo.scafSeq"]))
                self.logger.info("denovo111:{}".format(denovo_seq_file))
                self.get_fasta(denovo_seq_file, scaffold_ids[i], sample_region_ids[i])
            else:
                self.set_error('{} not in {}'.format(region, possname.keys()))

    def get_fasta(self, denovo_seq_file, scaffold_id, sample_region):
        """
        根据scaffold_id去获取对应序列
        /mnt/ilustre/users/sanger-dev/workspace/20180510/Assembly_wgs_test_0510094133_9523_7275/output/soap_denovo
        :return:
        """
        self.logger.info("开始从{}，获取序列:{}".format(denovo_seq_file, scaffold_id))
        if not os.path.exists(os.path.join(self.output_dir, 'seq_dir')):
            os.mkdir(os.path.join(self.output_dir, 'seq_dir'))
        seq_path = os.path.join(self.output_dir, "seq_dir/{}_{}.fa".format(sample_region, scaffold_id))
        if not os.path.exists(denovo_seq_file):
            self.set_error("文件%s不存在！" % denovo_seq_file)
        with open(denovo_seq_file, "r") as r, open(seq_path, "w") as w:
            data = r.readlines()
            # self.logger.info("dd:{}".format(len(data)))
            n = 0
            for line in data:
                # self.logger.info("1:{}".format(line))
                m = re.match(r'^>{}'.format(scaffold_id), line)
                if m:
                    # self.logger.info("2:{}".format('^>{}'.format(self.option('scaffold_id'))))
                    w.write(line.strip().split(' ')[0] + '\n')
                    break
                n += 1
            # self.logger.info("3:{}".format(n))
            while True:
                if len(data) == n + 1:
                    self.logger.info("cc:{}".format(n))
                    break
                if not re.match(r'^>.*', data[n + 1]):
                    # self.logger.info("4:{}".format(data[n+1]))
                    w.write(data[n + 1])
                    n += 1
                else:
                    break
        self.logger.info("获取序列:{}成功".format(scaffold_id))

    def compression_result(self):
        os.system('tar -czf seq_dir.tar.gz {}'.format(os.path.join(self.output_dir, 'seq_dir')))

    def set_db(self):
        self.logger.info("开始进行导表")
        api = self.api.api("wgs_v2.api_base")
        if self.option("project_type") != "dna_wgs_v2":
            api.project_type = self.option("project_type")
        api.update_db_record('sg_sequence', {'_id': ObjectId(self.option('main_id'))},
                             {'seq_path': self.option('target_path') + "/seq_dir"})

    def run(self):
        super(SeqDealTool, self).run()
        self.per_make_params()
        if self.option("main_id"):
            self.set_db()
        self.end()
