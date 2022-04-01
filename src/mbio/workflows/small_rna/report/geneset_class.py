# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from biocluster.workflow import Workflow
import os
from biocluster.core.exceptions import OptionError
import re
from bson.objectid import ObjectId
import shutil
import pandas as pd
import glob
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetClassWorkflow(Workflow):
    '''
    geneset_cog: cog_class_table.xls
    geneset_go: go_class_table.xls
    kegg_table: gene_kegg_table.xls
    kegg_level: gene_kegg_level_table.xls
    geneset_kegg: multi_geneset_list.txt
    add_info: add_info.txt
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetClassWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'geneset_id', 'type': 'string'},
            {'name': 'type', 'type': 'string'},
            {'name': 'anno_type', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'geneset_type', 'type': 'string', 'default': 'G'},
            # anno_type == 'cog'
            {'name': 'geneset_cog', 'type': 'string'},
            # anno_type == 'go'
            {'name': 'geneset_go', 'type': 'string'},
            # anno_type == 'kegg'
            {'name': 'kegg_table', 'type': 'string'},
            {'name': 'kegg_level', 'type': 'string'},
            {'name': 'geneset_kegg', 'type': 'string'},
            {'name': 'add_info', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            # essential key for updating information
            {'name': 'update_info', 'type': 'string'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "cog":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/02 COG_Annotation')
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/03 GO_Annotation')
        if self.option("anno_type") == "kegg":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/04 KEGG_Annotation')
        self.inter_dirs = []
        self.api_geneset = self.api.api('small_rna.geneset_class')

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
        super(GenesetClassWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self.sheet.options():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        if self.option('anno_type') == 'kegg' and len(self.option('geneset_id').split(',')) > 2:
            raise OptionError('number of geneset_id can not be greater than two when anno_type is kegg')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''
        if self.option('anno_type') == 'kegg':
            self.logger.info('start run tool when anno_type is {}'.format(self.option('anno_type')))
            self.tool = self.add_tool('small_rna.geneset.kegg_class')
            self.tool.on('end', self.set_output)
            self.get_run_log()
            self.run_tool()
            super(GenesetClassWorkflow, self).run()
        else:
            self.logger.info('skip run tool when anno_type is {}'.format(self.option('anno_type')))
            self.start_listener()
            self.fire('start')
            self.get_run_log()
            self.set_output()

    def get_run_log(self):
        if self.option('anno_type') == 'kegg':
            table = "sg_geneset_kegg_class"
        elif self.option('anno_type') == 'cog':
            table = "sg_geneset_cog_class"
        elif self.option('anno_type') == 'go':
            table = "sg_geneset_go_class"
        get_run_log = GetRunLog("small_rna", table=table, main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        opts = {
            'geneset_kegg': self.option('geneset_kegg'),
            'kegg_table': self.option('kegg_table'),
            'geneset_id': self.option('geneset_id'),
            'background_links': self.option('add_info'),
            'type': self.option('type'),
            'task_id': self.option('task_id')
        }
        if self.option('anno_type') == 'kegg':
            opts.update({"kegg_version": self.option('kegg_version')})
        self.tool.set_options(opts)
        self.tool.run()

    def set_output(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        if self.option('anno_type') == 'kegg':
            pass
        elif self.option('anno_type') == 'go':
            source = self.option('geneset_go')
        elif self.option('anno_type') == 'cog':
            source = self.option('geneset_cog')
        if self.option('anno_type') != 'kegg':
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        '''
        export result in output_dir of workflow to api
        '''
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        # check anno_type and adapt different strategies
        if self.option('anno_type') == 'kegg':
            output_file = os.path.join(self.tool.output_dir, 'kegg_stat.xls')
            self.api_geneset.add_kegg_regulate_new(main_table_id=self.option('main_id'),
                                                   geneset_id=self.option('geneset_id'),
                                                   kegg_stat_xls=output_file,
                                                   gene_kegg_level_table_xls=self.option('kegg_level'),
                                                   work_dir=self.work_dir,
                                                   geneset_type=self.option('geneset_type'))
            pathway_dir = os.path.join(self.tool.output_dir, 'pathways')
            self.api_geneset.add_kegg_regulate_pic(main_table_id=self.option('main_id'),
                                                   level_path=self.option('kegg_level'),
                                                   png_dir=pathway_dir)
            main_id = ObjectId(self.option('main_id'))
            self.workflow_output_tmp = self._sheet.output
            if re.match(r'tsanger:', self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
            else:
                self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
            graph_dir = os.path.join(self.workflow_output, 'pathways')
            conn = self.api_geneset.db['sg_geneset_kegg_class']
            conn.update({'_id': main_id}, {'$set': {'graph_dir': graph_dir}}, upsert=True)
        elif self.option('anno_type') == 'go':
            output_file = self.option('geneset_go')
            self.api_geneset.add_geneset_go_class(go_class_table=output_file, main_id=self.option('main_id'))
        elif self.option('anno_type') == 'cog':
            output_file = self.option('geneset_cog')
            self.api_geneset.add_geneset_cog_class(cog_class_table=output_file, main_id=self.option('main_id'))
        self.logger.debug('output_file: {}'.format(output_file))
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        if self.option('anno_type') == 'kegg':
            shutil.rmtree(os.path.join(self.tool.output_dir, 'ko'))
            shutil.copyfile(
                os.path.join(self.work_dir, 'kegg_analysis_of_anotate'),
                os.path.join(self.tool.output_dir, 'kegg_detail.xls')
            )
            statis = pd.read_table(os.path.join(self.work_dir, 'kegg_statistic'), sep='\t', header=0)
            statis_new = statis.drop(['kegg_id', 'geneset_type', 'geneset_id'], axis=1)
            statis_new.to_csv(os.path.join(self.work_dir, 'kegg_statistics'), sep='\t', index=False)
            shutil.copyfile(
                os.path.join(self.work_dir, 'kegg_statistics'),
                os.path.join(self.tool.output_dir, 'kegg_statistics.xls')
            )
            shutil.copyfile(
                os.path.join(self.tool.output_dir, 'kegg_statistics.xls'),
                os.path.join(self.output_dir, 'kegg_statistics.xls')
            )
            shutil.copyfile(
                os.path.join(self.tool.output_dir, 'kegg_detail.xls'),
                os.path.join(self.output_dir, 'kegg_detail.xls')
            )
            os.mkdir(os.path.join(self.output_dir, 'pathways'))
            for pdf in glob.glob(os.path.join(self.tool.output_dir, 'pathways/*.pdf')):
                shutil.copy(pdf, os.path.join(self.output_dir, 'pathways', os.path.basename(pdf)))
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/04 KEGG_Annotation", "", "靶基因KEGG功能注释", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', '靶基因KEGG功能注释文件', 0],
                ['pathways', '', '靶基因KEGG通路图', 0],
                ['kegg_statistics.xls', '', '靶基因KEGG分类统计表 ', 0],
                ['kegg_detail.xls', '', 'KEGG注释详情表', 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
            result_dir.add_regexp_rules([
                [r'pathways/map.*\.html', '', 'KEGG通路html文件', 0],
                [r'pathways/map.*\.png', '', 'KEGG通路图片png', 0],
            ])
        elif self.option('anno_type') == 'go':
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录",0],
                ["04 GeneSet/03 GO_Annotation", "", "靶基因GO功能注释", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', '靶基因GO功能注释文件', 0],
                ['go_class_table.xls', '', '靶基因GO分类统计表 ', 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        elif self.option('anno_type') == 'cog':
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
                os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
            os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录",0],
                ["04 GeneSet/02 COG_Annotation", "", "靶基因COG功能注释", 0]
            ]
            result_dir.add_relpath_rules([
                ['.', '', '靶基因COG功能注释文件', 0],
                ['cog_class_table.xls', '', '靶基因COG分类统计表 ', 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        super(GenesetClassWorkflow, self).end()
        # self.logger.debug(self.get_upload_files())
        # self.set_end()
        # self.fire('end')
        # self.end_unfinish_job()
        # self._upload_result()
        # self._import_report_data()
        # self._update('set_end')
        # self.step.finish()
        # self.step.update()
        # self.logger.info('succeed in running workflow')
        # self._save_report_data()
