# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import shutil
import glob
import re
import pickle
import time, datetime
import psutil
import shlex
import json
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class ProteinsetStringPictureWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetStringPictureWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "proteinset_list", "type": "string"},
            {"name": "useblast", "type": "string"},
            {"name": "species", "type": "int", "default": 0},
            {"name": "identity", "type": "float", "default": 98.0},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "infile", "format": "labelfree.common_dir"},
            {"name": "type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + '/labelfree/get_string_picture.py'
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/5_Proteinset/06_StringPic')
        self.inter_dirs = []

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
        super(ProteinsetStringPictureWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        # super(ProteinsetStringPictureWorkflow, self).run()
        self.run_string_picture()
        self.set_db()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        for file in glob.glob(os.path.join(self.output_dir, '*', '*')):
            shutil.copy(file, os.path.join(self.output_dir, os.path.basename(file)))
        try:
            shutil.rmtree(os.path.dirname(file))
        except:
            self.logger.info('数据爬取应该是失败了')
            return
        api_proteinset = self.api.api('dia.proteinset')

        self.logger.info("开始进行string_picture的导表")
        api_proteinset.add_string_picture(self.option('main_table_id'), self.output_dir)
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = api_proteinset.db["sg_proteinset_string_picture"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        for file in os.listdir(self.output_dir):
            if file.endswith('.svg'):
                svg = file
                break
        try:
            graph_dir = os.path.join(self.workflow_output, svg)
        except:
            self.logger.info('没能生成图片')
            self.end()
        conn.update({"_id": record_id}, {"$set": {'graph_dir': graph_dir}}, upsert=True)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析",0],
            ["5_Proteinset/06_StringPic", "", "蛋白网络分析（静态图）", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集STRING数据库爬取结果目录"],
            ["*.annotation.xls", "", "node注释信息"],
            ["*.bitscore.xls", "", "蛋白间关联度得分"],
            ["*.interaction.xls", "", "蛋白间关系文件"],
        ])
        super(ProteinsetStringPictureWorkflow, self).end()

    def run_string_picture(self):
        # 为了保证同一时间只能运行一个爬虫任务
        time_file = self.config.PACKAGE_DIR + '/labelfree/string.time'
        if not os.path.exists(time_file):
            with open(time_file, 'w') as tw:
                pickle.dump(datetime.datetime.now(), tw)
        wait_start = datetime.datetime.now()
        while 1:
            tr = open(time_file)
            last_time = pickle.load(tr)
            now_time = datetime.datetime.now()
            if not last_time:
                if (now_time - wait_start).seconds < 720:
                    time.sleep(300)
                    continue
                else:
                    break
            # if (now_time-last_time).seconds < 720:
            if (now_time - last_time).seconds < 20:
                self.logger.info('需要睡5分钟')
                time.sleep(300)
            else:
                break
        with open(time_file, 'w') as tw:
            # pickle.dump(datetime.datetime.now(), tw)
            pickle.dump('', tw)
        if self.option('species') == -1 or not self.option('species'):
            specie = self.get_most_common_specie(self.option('origin_result').prop['path'], )
        else:
            specie = self.option('species')
        useblast = 'no'
        if self.option('useblast').lower() == 'yes':
            useblast = 'yes'
        list_dir = os.path.join(self.work_dir, 'list_dir')
        if os.path.exists(list_dir):
            shutil.rmtree(list_dir)
        os.makedirs(list_dir)
        os.link(self.option('proteinset_list'), os.path.join(list_dir, os.path.basename(self.option('proteinset_list'))))
        shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)
        cmd = self.python_path + ' ' + self.script
        cmd += ' -specie ' + str(specie)
        cmd += ' -list_path ' + self.work_dir
        cmd += ' -vsstring ' + os.path.join(self.option('origin_result').prop['path'], "blast_xml/string.xml")
        cmd += ' -identity ' + str(self.option('identity'))
        # cmd += ' -max_num ' + str(self.option('max_num'))
        cmd += ' -useblast ' + useblast
        cmd += ' -out ' + self.output_dir
        self.logger.info(cmd)
        # os.system(cmd)
        stdout = open(os.path.join(self.work_dir, 'scrapy.log') ,'w')
        fn = stdout.fileno()
        p = psutil.Popen(shlex.split(cmd), stdout=fn)
        begin = datetime.datetime.now()
        while 1:
            time.sleep(20)
            self.logger.info(p.status())
            if p.status() == 'zombie':
                break
                # pass
            else:
                if (datetime.datetime.now() - begin).seconds > 600:
                    self.logger.info('运行时间太长，已经被杀掉')
                    p.terminal()
                    break
        stdout.close()
        with open(time_file, 'w') as tw:
            pickle.dump(datetime.datetime.now(), tw)

    def get_most_common_specie(self, xml):
        from Bio.Blast import NCBIXML
        from collections import Counter
        species = list()
        with open(xml, 'r') as blast_r:
            records = NCBIXML.parse(blast_r)
            for rec in records:
                for align in rec.alignments:
                    for hsp in align.hsps:
                        hit = align.hit_id
                        des = align.hit_def
                        if u'.' in des and '|' not in des:
                            hit = des
                        ident = float(hsp.identities)
                        hit_len = float(hsp.align_length)
                        if ident / hit_len * 100 >= self.option('identity'):
                            # specie = re.split('.', hit, maxsplit=1)[0]
                            specie = hit.split('.')[0]
                            # self.logger.info(hit)
                            species.append(specie)
        if not species:
            return 9096
        most_s = Counter(species).most_common()[0][0]
        return int(most_s)