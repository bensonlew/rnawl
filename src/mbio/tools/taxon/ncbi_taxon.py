# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import threading
import re
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon



class NcbiTaxonAgent(Agent):
    """
    gi2taxon.py v1.0
    author: shenghe
    last_modify: 2016.6.21
    """

    def __init__(self, parent):
        super(NcbiTaxonAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml, align.blast.blast_table"},  # 输入文件
            {"name": "out_type", "type": "int", "default": 0},  # type为0或1,0为未整理的物种信息，1为整理后的物种信息 add by zouxuan
            {"name": "taxon_out", "type": "outfile", "format": "taxon.seq_taxon,annotation.nr.nr_taxon"},  # 输出结果文件
            {"name": "blastdb", 'type': 'string', 'default': 'None'}  # 输入文件的blast比对类型，必须为nr或者nt
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须设置输入文件", code="34200201")
        if self.option('blastdb') == 'None':
            raise OptionError("必须设置输入文件的blast比对类型", code="34200202")
        else:
            if self.option('blastdb') not in ['nr', 'nt', "nt_v20200604",'nt_v20210917']:
                raise OptionError('blastdb必须为nr或者nt:%s', variables=(self.option('blastdb')), code="34200203")
        if self.option('out_type') not in [0, 1]:
            raise OptionError("out_type必须是0或1", code="34200204")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '12G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
        ])
        super(NcbiTaxonAgent, self).end()


class NcbiTaxonTool(Tool):
    def __init__(self, config):
        super(NcbiTaxonTool, self).__init__(config)
        self._version = "1.0"

    class result_thread(threading.Thread):
        def __init__(self, func, *argu, **kwargu):
            super(NcbiTaxonTool.result_thread, self).__init__()
            self.func = func
            self.argu = argu
            self.kwargu = kwargu
            self.result = None

        def run(self):
            self.result = self.func(*self.argu, **self.kwargu)

    def run(self):
        """
        运行
        :return:
        """
        super(NcbiTaxonTool, self).run()
        if self.option("blastout").format == 'align.blast.blast_xml':
            if self.option('blastdb') in ['nt_v20200604','nt_v20210917']:
                from mbio.packages.align.blast.xml2table import xml2table
                table = xml2table(self.option('blastout').path, self.work_dir + '/temp_blastable.xls', hit_id="yes")
            else:
                from mbio.packages.align.blast.xml2table import xml2table
                table = xml2table(self.option('blastout').path, self.work_dir + '/temp_blastable.xls')
        else:
            table = self.option('blastout').path
        if self.option('blastdb') == 'nr':
            db = 'prot'
        else:
            db = 'nucl'
        if self.process_run(table, db):
            # run_gitaxon = self.result_thread(self.process_run, table, db)
            # run_gitaxon.start()
            # run_gitaxon.join()
            # if run_gitaxon.result:
            if self.option('out_type') == 1:
                if self.new_taxon():
                    self.option('taxon_out', self.output_dir + '/query_taxons.xls')
                else:
                    self.set_error('注释结果转换出错！', code="34200201")
            else:
                self.option('taxon_out', self.output_dir + '/query_taxons_detail.xls')
            self.end()
        else:
            self.set_error('注释查询出错！', code="34200201")

    def process_run(self, fp, db):
        """单独的线程运行查询分类注释"""
        if db == 'prot':
            from mbio.packages.taxon.gi2taxon_update import gi_taxon  # import进入了sqlite3的对象，这个对象不可以跨线程使用
            query = self.filter_query(fp, db)
            id_taxons = gi_taxon(set(query.values()), db)
        else:
            if self.option('blastdb') in ['nt']:
                from mbio.packages.taxon.accession2taxon import accession_taxon  # import进入了sqlite3的对象，这个对象不可以跨线程使用
            elif self.option('blastdb') in ['nt_v20210917']:
                from mbio.packages.taxon.accession2taxon_v20210917 import accession_taxon
            else:
                from mbio.packages.taxon.accession2taxon_update import accession_taxon
            query = self.filter_query(fp, db)
            id_taxons = accession_taxon(set(query.values()))  # modified by zouxuan,add accession2taxon
        for i in query:
            query[i] = (query[i], id_taxons[query[i]])
        with open(self.output_dir + '/query_taxons_detail.xls', 'w') as w:
            for item in query.iteritems():
                if item[1][1]:
                    w.write(item[0] + '\t' + item[1][0] + '\t' + item[1][1] + '\n')
        return True

    def filter_query(self, fp, db):
        """生成query对gi的"""
        query = dict()
        openfp = open(fp)
        openfp.readline()
        for i in openfp:
            line_sp = i.split('\t')
            if line_sp[5] in query:
                continue
            else:
                hitname = line_sp[10].split('|')
                if len(hitname) != 1:  #新版本只有accession号
                    if len(hitname) < 3 or hitname[0] != 'gi' or (not re.match(r'^\d+$', hitname[1])):
                        self.set_error('输入文件中不是nr库比对结果,不含有gi信息:%s', variables=(line_sp[10]), code="34200203")
                if db == 'prot':
                    query[line_sp[5]] = hitname[1]
                else:
                    query[line_sp[5]] = hitname[0]  # modified by zouxuan,add get accession
        return query

    def new_taxon(self):
        '''
         注释结果转换 add by zouxuan
        :return:
        '''
        tax = mg_taxon()
        self.logger.info("start new_taxon(detail_to_level)")
        try:
            tax.detail_to_level(self.output_dir + '/query_taxons_detail.xls', self.output_dir)
        except Exception as e:
            self.set_error("new_taxon(detail_to_level) failed", code="34200204")
            raise Exception("new_taxon(detail_to_level) failed {}".format(e))
        return True

"""
class TestFunction(unittest.TestCase):

    This is test for the tool. Just run this script to do test.

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Taxon",
            "type": "tool",
            "name": "taxon.ncbi_taxon",
            "instant": True,
            "options": dict(
                blastout="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/metag_v2/nr/blast_tmp/fasta_10_vs_nr.xml",
                blastdb="nr",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
"""
