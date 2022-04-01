# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import shutil
from bson.objectid import ObjectId
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetEnrichWorkflow(Workflow):
    '''
    geneset_list: geneset.list
    all_list: all.list
    go_list: go.list
    kegg_table: gene_kegg_table.xls
    geneset_kegg: multi_geneset_list.txt
    add_info: add_info.txt
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'geneset_id', 'type': 'string'},
            {'name': 'type', 'type': 'string'},
            {'name': 'anno_type', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'geneset_type', 'type': 'string', 'default': 'G'},
            {'name': 'method', 'type': 'string'},
            {'name': 'geneset_list', 'type': 'string'},
            {'name': 'all_list', 'type': 'string'},
            # anno_type == 'go'
            {'name': 'go_list', "type": 'string'},
            # anno_type == 'kegg'
            {'name': 'kegg_table', 'type': 'string'},
            {'name': 'geneset_kegg', 'type': 'string'},
            {'name': 'add_info', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            # essential key for updating information
            {"name": "update_info", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "kegg":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/06 KEGG_Enrich')
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/05 GO_Enrich')
        self.inter_dirs = []
        if self.option('anno_type') == 'go':
            self.go_enrich = self.add_tool('small_rna.geneset.go_enrich')
        elif self.option('anno_type') == 'kegg':
            self.kegg_rich = self.add_tool('small_rna.geneset.kegg_rich')
            self.kegg_class = self.add_tool('small_rna.geneset.kegg_class')
        self.api_geneset = self.api.api('small_rna.geneset_enrich')

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
        super(GenesetEnrichWorkflow, self).send_log(data)


    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self.sheet.options():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''
        self.get_run_log()
        if self.option('anno_type') == 'go':
            self.go_enrich.on('end', self.set_output)
            self.run_go_enrich()
        elif self.option('anno_type') == 'kegg':
            self.kegg_rich.on('end', self.run_kegg_class)
            self.kegg_class.on('end', self.set_output)
            self.run_kegg_rich()
        super(GenesetEnrichWorkflow, self).run()

    def get_run_log(self):
        if self.option('anno_type') == 'kegg':
            table = "sg_geneset_kegg_enrich"
        elif self.option('anno_type') == 'go':
            table = "sg_geneset_go_enrich"
        get_run_log = GetRunLog("small_rna", table=table, main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_go_enrich(self):
        opts = {
            'diff_list': self.option('geneset_list'),
            'go_list': self.option('go_list'),
            'method': self.option('method'),
        }
        self.go_enrich.set_options(opts)
        self.go_enrich.run()

    def run_kegg_rich(self):
        opts = {
            'diff_list': self.option('geneset_list'),
            'correct': self.option('method'),
            'kegg_table': self.option('kegg_table'),
            'kegg_version': self.option('kegg_version'),
            'add_info': self.option('add_info')
        }
        self.kegg_rich.set_options(opts)
        self.kegg_rich.run()

    def run_kegg_class(self):
        opts = {
            'geneset_kegg': self.option('geneset_kegg'),
            'task_id': self.option('task_id'),
            'kegg_table': self.option('kegg_table'),
            'geneset_id': self.option('geneset_id'),
            'kegg_version': self.option('kegg_version'),
            'background_links': self.option('add_info'),
            'type': self.option('type'),
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()

    def set_output(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        if self.option('anno_type') == 'go':
            os.link(
                os.path.join(self.go_enrich.output_dir, 'go_enrich_geneset.xls'),
                os.path.join(self.output_dir, 'go_enrich_geneset.xls')
            )
        elif self.option('anno_type') == 'kegg':
            os.link(
                os.path.join(self.kegg_rich.output_dir, 'geneset.list.DE.list.check.kegg_enrichment.xls'),
                os.path.join(self.output_dir, 'kegg_enrichment.xls')
            )
            os.mkdir(os.path.join(self.output_dir, 'pathways'))
            for pdf in glob.glob(os.path.join(self.kegg_class.output_dir, 'pathways/*.pdf')):
                shutil.copy(pdf, os.path.join(self.output_dir, 'pathways', os.path.basename(pdf)))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        '''
        export result in output_dir of workflow to api
        '''
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        output_file = glob.glob('{}/*.xls'.format(self.output_dir))[0]
        result_file = os.path.join(self.get_workflow_output_dir(), os.path.basename(output_file))
        # check anno_type and adapt different strategies
        if self.option('anno_type') == 'go':
            self.api_geneset.add_go_enrich_detail(go_enrich_id=self.option('main_id'),
                                                  go_enrich_dir=output_file)
            self.api_geneset.update_db_record('sg_geneset_go_enrich', self.option("main_id"), result_file=result_file)
        elif self.option('anno_type') == 'kegg':

            self.api_geneset.add_kegg_enrich_detail(kegg_enrich_id=self.option('main_id'),
                                                    kegg_enrich_table=output_file,
                                                    geneset_list_path=self.option('geneset_list'),
                                                    all_list_path=self.option('all_list'))
            self.api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_id"), result_file=result_file)
            png_dir = os.path.join(self.kegg_class.output_dir, 'pathways')
            self.api_geneset.add_kegg_enrich_pic(main_table_id=self.option('main_id'),
                                                 level_path=output_file,
                                                 png_dir=png_dir)
            graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
            self.api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option('main_id'), graph_dir=graph_dir)
        self.logger.debug('output_file: {}'.format(output_file))
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if self.option('anno_type') == 'go':
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/05 GO_Enrich", "", "靶基因GO功能富集", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', '靶基因GO富集分析文件', 0],
                ['go_enrich_geneset.xls', '', '靶基因GO富集分析统计表 ', 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        elif self.option('anno_type') == 'kegg':
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/06 KEGG_Enrich", "", "靶基因KEGG富集分析", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', '靶基因KEGG富集分析文件', 0],
                ['pathways', '', 'KEGG富集通路图', 0],
                ['kegg_enrichment.xls', '', '靶基因KEGG富集分析统计表 ', 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
            result_dir.add_regexp_rules([
                [r'pathways/map.*\.html', '', 'KEGG通路html文件',0],
                [r'pathways/map.*\.png', '', 'KEGG通路图片png',0],
            ])
        super(GenesetEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
