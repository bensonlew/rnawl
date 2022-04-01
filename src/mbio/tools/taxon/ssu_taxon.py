# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config


class SsuTaxonAgent(Agent):
    """
    author: shenghe
    last_modify: 2017.05.12
    """

    def __init__(self, parent):
        super(SsuTaxonAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml_dir, align.blast.blast_table_dir"},  # 输入文件
            {"name": "taxon_out", "type": "outfile", "format": "annotation.nr.nr_taxon"},  # 输出结果文件 没有定义内容
            {"name": "ssu_db", 'type': 'string', 'default': 'None'}  # 数据库版本有 silva119/silva123/silva128
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须设置输入文件")
        if self.option('ssu_db') == 'None':
            raise OptionError("必须设置ssu数据库比对版本")
        else:
            if self.option('ssu_db') not in ['silva119', 'silva123', 'silva128']:
                raise OptionError('ssu_db必须为silva128、silva123或者silva119:{}'.format(self.option('ssu_db')))
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
            ])
        super(SsuTaxonAgent, self).end()


class SsuTaxonTool(Tool):
    def __init__(self, config):
        super(SsuTaxonTool, self).__init__(config)
        self._version = "1.0"


    def run(self):
        """
        运行
        :return:
        """
        super(SsuTaxonTool, self).run()
        table_file = self.work_dir + '/temp_blastable.xls'
        if self.option("blastout").format == 'align.blast.blast_xml':
            self.option("blastout").convert2table(table_file)
        else:
            table_file = self.option('blastout').path
        try:
            self.taxon_ssu(table_file, self.option('ssu_db'), self.output_dir + "/ssu_taxon.xls")
            self.end()
        except Exception:
            import traceback
            self.logger.debug(traceback.format_exc())
            self.set_error('注释查询出错！')

    def taxon_ssu(self, table_file, dbversion, outfile):
        query_ids = []
        hit_ids = []
        for i in open(table_file):
            values = i.split('\t')
            query_ids.append(values[5])
            hit_ids.append(values[10])
        query_hit = zip(query_ids, hit_ids)
        # db = self.config.biodb_mongo_client.sanger_biodb
        db = Config().get_mongo_client(mtype="meta", ref=True)[Config().get_mongo_dbname("meta", ref=True)]
        coll_name = dbversion + '_ssu_taxon'
        collection = db[coll_name]
        outw = open(outfile, "w")
        for query, hit in query_hit:
            find = collection.find_one({"_id": hit}, {"taxon": 1})
            outw.write(query + '\t' + find["_id"] + '\t' + find["taxon"] + '\n')
