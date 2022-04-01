# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import os
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json


class TaxonCompareWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """reads 物种注释比较分析, 目前支持pcoa"""
        self._sheet = wsheet_object
        super(TaxonCompareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "main_col", "type": "string"},
            {"name": "table", "type": "infile", "format": "meta_genomic.taxon_dir"},
            {"name": "col", "type": "string"},
            {"name": "name2id", "type": "string"},
            {"name": "level_id", "type": "int"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "analysis_type", "type": "string", "default": 'pcoa'},
            {"name": "distance_method", "type": "string", "default": "bray_curtis"},
            {"name": "diff_check", "type": "string", "default": "none"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
            {"name": "permutation", "type": "int", "default": 999},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.compare = self.add_module("meta.beta_diversity.beta_diversity")
        self.filter_table = self.add_tool("metagenomic.table_select")
        self.ellipse = self.add_tool("graph.ellipse")
        self.skip_ellipse = True

    def run(self):
        self.filter_table.on("end", self.run_compare)
        group_detail = json.loads(self.option("group_detail"))
        if len(group_detail.keys()) > 1:
            for k in group_detail.keys():
                if len(group_detail[k]) > 2:
                    self.skip_ellipse = False
                    break
        if not self.skip_ellipse:
            self.compare.on('end', self.run_ellipse)
            self.ellipse.on('end', self.set_db)
        else:
            self.compare.on('end', self.set_db)

        self.run_filter_table()
        super(TaxonCompareWorkflow, self).run()

    def run_filter_table(self):
        samples = json.loads(self.option("name2id")).keys()
        select_col = [self.option("col")] + samples
        opts = {
            "cols": json.dumps(select_col),
            "table": self.option("table").get_level(self.option("level_id"))
        }
        self.filter_table.set_options(opts)
        self.filter_table.run()

    def run_compare(self):
        otutable = self.filter_table.option("out_table").path
        opts = {
            'analysis': self.option('analysis_type'),
            'otutable': otutable,
        }
        if self.option('analysis_type') in ['pcoa', 'nmds', 'dbrda']:
            opts['dis_method'] = self.option('distance_method')
        if self.option('analysis_type') in ['pca', 'pcoa', 'nmds']:
            opts['group'] = self.option('group')
            if not self.skip_ellipse:
                opts['ellipse'] = "T"
        if self.option("diff_check") != "none":
            opts["diff_test_method"] = self.option("diff_check")
            opts["change_times"] = str(self.option("permutation"))
        self.compare.set_options(opts)
        self.compare.run()

    def run_ellipse(self):
        opts = {}
        if self.option("group").is_set:
            opts['group_table'] = self.option("group")
        pc_map = {'pca': "/Pca/pca_sites.xls", 'pcoa': "/Pcoa/pcoa_sites.xls",
                  'dbrda': '/Dbrda/db_rda_sites.xls', 'nmds': '/Nmds/nmds_sites.xls',
                  'rda_cca': '/Rda'
                  }
        opts['analysis'] = self.option('analysis_type')
        opts['pc_table'] = self.compare.output_dir + pc_map[self.option('analysis_type')]
        self.ellipse.set_options(opts)
        self.ellipse.run()

    def set_db(self):
        for f in os.listdir(self.compare.output_dir):
            path1 = os.path.join(self.compare.output_dir, f)
            if os.path.isdir(path1):
                for f2 in os.listdir(path1):
                    self.link(os.path.join(path1, f2))
            else:
                self.link(path1)
        api_beta = self.api.api('metagenomic.beta_diversity')
        api_beta.add_beta_diversity(self.output_dir, self.option('analysis_type'),
                                    main=None, main_id=str(self.option('main_id')),
                                    main_col=self.option("main_col"),
                                    diff_check=self.option("diff_check"))
        if not self.skip_ellipse:
            self.link(os.path.join(self.ellipse.work_dir, "ellipse_out.xls"))
            api_beta.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls', str(self.option('main_id')))
        self.logger.info('运行self.end')
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "taxon_compare")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": "taxon_compare_pcoa",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        sdir = self.add_upload_dir(self.output_dir)
        repaths = [
            [".", "", "PCoA分析结果目录", 0, ""],
            ["pcoa_eigenvalues.xls", "xls", "矩阵特征值", 0, "120121"],
            ["pcoa_eigenvaluespre.xls", "xls", "特征解释度百分比",0,"120280"],
            ["pcoa_sites.xls", "xls", "样本坐标表", 0, "120123"],
            ["PCoA.pdf", "pdf", "PCoA图"],
            ["PCoA_box.pdf", "pdf", "PCoA箱线图"]
        ]
        regexps = [
            [r'%s.*\.xls$' % self.option('distance_method'), 'xls', '样本距离矩阵文件', 0, "120156"],
        ]
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(TaxonCompareWorkflow, self).end()
