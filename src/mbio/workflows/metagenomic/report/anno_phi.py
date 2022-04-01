# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2018/11/22'

from biocluster.workflow import Workflow
import os
from comfun import ComfunWorkflow
from mbio.packages.metabolome.common import link_file, link_dir
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class AnnoPhiWorkflow(ComfunWorkflow):
    """
    phi注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AnnoPhiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {'name': 'save_pdf', 'type': 'int', 'default': 0}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.sth_tool = self.add_tool('metabolome.annotation.anno_keggc')  # edit tool
        # self.sth_module = self.add_module('metabolome.annotation.anno_keggp')  # edit module
        # self.step.add_steps("keggp", "keggc", "overview")  # edit steps

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run anno workflow")
        self.run_filter(self.set_db)
        # self.edit_steps.on("end", self.set_db)
        # self.run_sth()
        super(AnnoPhiWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        # self.sheet.id = "tsg_31796"
        # self.sheet.project_sn = "188_5b5acb3018914"
        link_dir(self.filter.output_dir, self.output_dir)
        link_file(self.filter.recal_abu_tool.option('gene_anno_table').path, self.output_dir + '/gene_phi_anno.xls')
        api_path = self.api.api("metagenomic.mg_anno_phi")
        # main_id = api_path.add_anno_phi({"":""}, "5a54a80bedcb254cbab9d45e", self.output_dir + "/gene_phi_anno.xls")
        api_path.add_anno_phi_host(self.option('main_table_id'), self.output_dir + '/phi_host_profile.xls')
        api_path.add_anno_phi_pathogen(self.option('main_table_id'), self.output_dir + '/phi_pathogen_profile.xls')
        api_path.add_anno_phi_phenotype(self.option('main_table_id'), self.output_dir + '/phi_phenotype_profile.xls')
        api_path.add_anno_phi_detail(self.option('main_table_id'), self.output_dir + '/gene_phi_anno.xls')
        api_path.update_anno_file(self.option('main_table_id'), os.path.join(self._sheet.output, 'gene_phi_anno.xls'))
        if self.option("save_pdf"):
        # self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        # if self.pdf_status:
            name = get_name(self.option("main_table_id"), "anno_phi")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "anno_phi")
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.option("save_pdf"):
        # if self.pdf_status:
            # if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PHI结果目录",0,"120258"],
            ["gene_phi_anno.xls", "xls", "PHI功能注释结果信息表", 0, "120259"],
            ["phi_host_profile.xls", "xls", "各样品Host丰度表", 0, "120260"],
            ["phi_pathogen_profile.xls", "xls", "各样品Pathogen丰度表", 0, "120261"],
            ["phi_phenotype_profile.xls", "xls", "各样品Phenotype丰度表", 0, "120262"],
            ["phi_protein_profile.xls", "xls", "各样品Protein丰度表", 0, "120263"],
            ["phi_phenotype_bar.pdf", "pdf", "phenotype基因丰度图"]
        ])
        super(AnnoPhiWorkflow, self).end()

    # def run_sth(self):
    #     opts = {
    #         {"name": "anno_list", "type": "string", "default": "keggc,keggp,overveiw"},  # str/int/float/bool...
    #         {"name": "metab_table", "type": "infile", "format": "metabolome.metab_desc"},  # infile
    #         {"name": "overview_ko", "type": "outfile", "format": "sequence.profile_table"}  # outfile
    #     }
    #     self.sth_tool.set_options(opts)
    #     self.sth_tool.run()
