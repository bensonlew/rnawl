# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180409
import json
import re
import os
import unittest

from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class DrawCircosWorkflow(Workflow):
    """
    交互分析：circos图
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DrawCircosWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "target_trans", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "target_cis", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "diff_exp", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "rna_type", "type": "infile", "format": "lnc_rna.lnc_common"},
            {"name": "top_ref_num", "type": "int", "default": 10},
            # {"name": "diff_group", "type": "string"},

            {"name": "lncrna_gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
            {"name": "mrna_gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},

            {"name": "ref_fa_fai", "type": "infile", "format": "lnc_rna.lnc_common"},

            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "project_type", "type": "string"},
            {"name": "task_id", "type": "string", "default": None},
            # {"name": "geneset_type", "type": "string", "default": 'T'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.circos_tool = self.add_tool('lnc_rna.geneset.draw_circos')
        self._sheet.output = self._sheet.output.replace('interaction_results/circos',
                                                        'interaction_results/04 GeneSet/11 Circos')
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
        super(DrawCircosWorkflow, self).send_log(data)

    def check_options(self):
        for name in (
                'target_trans', 'target_cis', 'diff_exp', 'rna_type', 'ref_fa_fai', 'lncrna_gtf', 'mrna_gtf', 'ref_fa_fai'):
            if not self.option(name).is_set:
                # get_size
                raise OptionError('必须提供 %s w文件' % name)
            if self.option(name).get_size() <= 1:
                raise OptionError('%s文件为空，请检查' % name)
        return True

    def run_circos(self):
        """
        {"name": "gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
        :return:
        """
        options = {
            'target_trans': self.option('target_trans'),
            'target_cis': self.option('target_cis'),
            'diff_exp': self.option('diff_exp'),
            'rna_type': self.option('rna_type'),
            'ref_fa_fai': self.option('ref_fa_fai'),
            'top_ref_num': self.option('top_ref_num'),
            'geneset_type': 'T',
            "lncrna_gtf": self.option('lncrna_gtf'),
            "mrna_gtf": self.option('mrna_gtf'),
        }
        self.circos_tool.set_options(options)
        self.circos_tool.on("end", self.set_output, "circos_tool")
        self.circos_tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'circos_tool':
            self.linkdir(obj.output_dir, self.output_dir)
        self.set_db()

    def linkdir(self, dirpath, dirname):
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
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("设置circos的结果路径")
        sample_id = []
        for file_ in os.listdir(self.output_dir):
            m = re.match(r"(.*)\.png$", file_)
            if m:
                sample_id.append(m.group(1))
        api = self.api.api("lnc_rna.api_base")
        # if self.option("project_type"):
        #     api._project_type = self.option("project_type")
        target_dir = self._sheet.output.rstrip('/') + '/'
        api.update_db_record("sg_circos", ObjectId(self.option("main_id")),
                             insert_dict={"circos_result_path": target_dir, "sample_list": sample_id})
        self.end()

    def run(self):
        self.get_run_log()
        self.run_circos()
        super(DrawCircosWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_circos", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/11 Circos", "", "染色体分布circos图", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "染色体分布circos图文件", 0],
            ['./circos.pdf', '', '染色体分布circos图pdf', 0],
            ['./circos.png', '', '染色体分布circos图png', 0],
            ['./circos.svg', '', '染色体分布circos图svg', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        # result_dir.add_regexp_rules([
        #     ["", "", ""]
        # ])
        super(DrawCircosWorkflow, self).end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        '''
        This is test for the workflow. Just run this script to do test.
        '''

        def test(self):
            from biocluster.wsheet import Sheet
            import random
            data = {
                'id': 'draw_circos_workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                'type': 'workflow',
                'name': 'lnc_rna.report.draw_circos',
                'options': {
                    # {"name": "target_trans", "type": "infile", "format": "lnc_rna.lnc_common"},
                    'target_trans': '/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                                    'trans_targets.txt',
                    # {"name": "target_cis", "type": "infile", "format": "lnc_rna.lnc_common"},
                    'target_cis': '/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                                  'cis_targets.txt',
                    # {"name": "diff_exp", "type": "infile", "format": "lnc_rna.lnc_common"},
                    'diff_exp': '/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                                'diff_exp.txt',
                    # {"name": "rna_type", "type": "infile", "format": "lnc_rna.lnc_common"},
                    'rna_type': '/mnt/ilustre/users/sanger-dev/workspace/20190401/DrawCircos_lnc_rna_8841_1244/'
                                'seq_id_type.txt',
                    # {"name": "top_ref_num", "type": "int", "default": 10},
                    'top_ref_num': 10,
                    # {"name": "gtf", "type": "infile", "format": "lnc_rna.lnc_gtf"},
                    'gtf': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-cufflinks/'
                           'output/NewTranscripts/ref_and_new.gtf',
                    # {"name": "ref_fa_fai", "type": "infile", "format": "lnc_rna.lnc_common"},
                    'ref_fa_fai': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens'
                                  '/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa.fai',
                    # {"name": "main_id", "type": "string"},
                    'main_id': '5ca1c67b17b2bf12ac7e36bf',
                    # {"name": "update_info", "type": "string"},
                    'update_info': json.dumps({'5ca1c67b17b2bf12ac7e36bf': "sg_circos"}),
                    # {"name": "project_type", "type": "string"}
                }
            }
            wsheet = Sheet(data=data)
            wf = DrawCircosWorkflow(wsheet)
            wf.sheet.id = 'lnc_rna'
            wf.sheet.project_sn = 'lnc_rna'
            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_DATA = False
            wf.run()


    unittest.main()
