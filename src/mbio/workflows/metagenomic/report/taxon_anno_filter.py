# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
# last_modifiy =  20210304

from biocluster.workflow import Workflow
import json
import os


class TaxonAnnoFilterWorkflow(Workflow):
    """
    宏基因组个性化注释
    """

    def __init__(self, wsheet_object):
        super(TaxonAnnoFilterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_meth", "type": "string", "default": "metaphlan3"},  # metaphlan3 kraken2
            {"name": "anno_file", "type": "infile", "format": "meta_genomic.taxon_dir"},
            {"name": "samples", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "main_col", "type": "string"},
            {"name": "level_filter", "type": "string"},
            {"name": "sample_filter", "type": "string"},
            {"name": "abu_filter", "type": "string"},
            {"name": "name2id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.anno_filter = self.add_tool("metagenomic.taxon_anno_filter")

    def check_options(self):
        if self.option('anno_meth') not in ['metaphlan3', 'kraken2']:
            self.set_error("物种注释功能目前只支持 ['metaphlan3', 'kraken2']")
        if not self.option('anno_file'):
            self.set_error("缺少输入文件")
        return True

    def run_filter(self):
        name2id = json.loads(self.option("name2id"))
        id2name = {v: k for k, v in name2id.items()}
        samples = [id2name[i] for i in self.option("samples").split(',')]
        opts = {
            "anno_file": self.option("anno_file"),
            "samples": json.dumps(samples),
            "level_filter": self.option("level_filter"),
            "sample_filter": self.option("sample_filter"),
            "abu_filter": self.option("abu_filter"),
        }
        self.anno_filter.set_options(opts)
        self.anno_filter.run()

    def run(self):
        self.anno_filter.on("end", self.set_db)
        self.run_filter()
        super(TaxonAnnoFilterWorkflow, self).run()

    def set_db(self):
        """
        保存结果output，导mongo数据库
        """
        api = self.api.api('metagenomic.common_api')
        outfiles = os.listdir(self.anno_filter.output_dir)
        key_map = json.loads(self.option("name2id"))
        key_map.update({
            "Domain": 'd__' , "Kingdom": 'k__', 'Phylum': "p__", 'Class': "c__",
            'Order': "o__", 'Family': "f__", 'Genus': "g__", 'Species': "s__"
        })
        for f in outfiles:
            file_path = os.path.join(self.anno_filter.output_dir, f)
            self.logger.info(file_path)
            api.add_detail(file_path, self.option('main_col') + '_detail',
                           self.option('main_table_id'), 'anno_id',
                           main_table=self.option('main_col'),
                           key_map=key_map,
                           update_main={"anno_file": self.sheet.output})
        self.end()

    def end(self):
        for f in os.listdir(self.anno_filter.output_dir):
            self.link(os.path.join(self.anno_filter.output_dir, f))
        up_dir = self.add_upload_dir(self.output_dir)
        up_dir.add_relpath_rules([
            ['.', '', '结果文件夹', 0, '0000'],
            ['*.txt', 'txt', '结果', 0, '0000'],
        ])
        super(TaxonAnnoFilterWorkflow, self).end()
