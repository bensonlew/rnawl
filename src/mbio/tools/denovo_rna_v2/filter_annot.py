# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import pandas as pd
__author__ = 'gdq'


class FilterAnnotAgent(Agent):
    """
    filter xml result
    """
    def __init__(self, parent):
        super(FilterAnnotAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'xml', 'format': 'denovo_rna_v2.blast_xml'},
            {'type': 'infile', 'name': 'hmm', 'format': 'denovo_rna_v2.common'},
            {'default': 'xml', 'type': 'string', 'name': 'types'},
            {'default': 0.001, 'type': 'float', 'name': 'evalue'},
            {'default': 0, 'type': 'float', 'name': 'identity'},
            {'default': 0, 'type': 'float', 'name': 'similarity'},
            {'type': 'outfile', 'name': 'outxml', 'format': 'denovo_rna_v2.blast_xml'},
            {'default': None, 'type': 'string', 'name': 'exclude_taxon'},
            {'type': 'outfile', 'name': 'outtable', 'format': 'denovo_rna_v2.common'},
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("types") == "xml":
            if not self.option("xml").is_set:
                raise OptionError("必须设置参数query", code = "32003501")
        elif self.option("types") == "hmm":
            if not self.option("hmm").is_set:
                raise OptionError("必须设置参数hmm", code = "32003502")
        else:
            raise OptionError('method %s必须为xml或hmm', variables = (self.option('types')), code = "32003503")

        if not (self.option('evalue') >= 0 and self.option('evalue') <= 0.001) :
            raise OptionError('E-value值设定必须为[0-0.001]之间：%s', variables = (str(self.option('evalue'))), code = "32003504")
        if not (self.option('identity') >= 0 and self.option('identity') <= 100) :
            raise OptionError('Identity值设定必须为[0-100]之间：%s', variables = (str(self.option('identity'))), code = "32003505")
        if not (self.option('similarity') >= 0 and self.option('similarity') <= 100) :
            raise OptionError('similarity值设定必须为[0-100]之间：%s', variables = (str(self.option('similarity'))), code = "32003506")

        return True


    def set_resource(self):
        self._cpu = 1
        mem = 5
        if self.option('xml').is_set:
            file_size = float(os.path.getsize(self.option('xml').prop['path'])) / 1024 / 1024
            mem = int(float(file_size)/1024 * 18) + 2
        mem = max(mem, 32)
        self._memory = "{}G".format(str(mem))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ])
        result_dir.add_regexp_rules([
            [r".+_vs_.+\.xml", "xml", "blast比对输出结果，xml格式"],
            [r".+_vs_.+\.xls", "xls", "hmm比对输出结果，表格(制表符分隔)格式"],
            ])
        # print self.get_upload_files()

        super(FilterAnnotAgent, self).end()


class FilterAnnotTool(Tool):
    """
    filter xml result
    """
    def __init__(self, config):
        super(FilterAnnotTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
 
    def run_filterxml(self):
        basename = os.path.basename(self.option('xml').prop['path']) + ".filter.xml"
        if self.option("exclude_taxon"):
            self.option('xml').filter_blast_xml(basename, self.option('evalue'), self.option('identity'), self.option('similarity'), filter_taxon_ids=self.option("exclude_taxon"))
        else:
            self.option('xml').filter_blast_xml(basename, self.option('evalue'), self.option('identity'), self.option('similarity'))

    def run_filter_hmm(self):
        self.logger.info("hmm 文件 {}".format(os.path.basename(self.option('hmm').prop['path'])))
        basename = os.path.basename(self.option('hmm').prop['path'])
        pfam_table = pd.read_table(self.option('hmm').prop['path'], header=0)
        pfam_select = pfam_table[pfam_table['DomainE-Value'] < self.option('evalue')]
        pfam_select.to_csv(basename + '.filter.xls', sep="\t",index=False)

    def set_output(self):
        if(self.option('types') == "xml"):
            filter_xml = os.path.basename(self.option('xml').prop['path']) + ".filter.xml"
            if os.path.exists(os.path.join(self.output_dir, filter_xml)):
                os.remove(os.path.join(self.output_dir, filter_xml))
            os.link(filter_xml, os.path.join(self.output_dir, filter_xml))
            self.option('outxml', os.path.join(self.output_dir, filter_xml))
        elif(self.option('types') == "hmm"):
            filter_xml = os.path.basename(self.option('hmm').prop['path']) + ".filter.xls"
            if os.path.exists(os.path.join(self.output_dir, filter_xml)):
                os.remove(os.path.join(self.output_dir, filter_xml))
            os.link(filter_xml, os.path.join(self.output_dir, filter_xml))
            self.option('outtable', os.path.join(self.output_dir, filter_xml))

    def run(self):
        super(FilterAnnotTool, self).run()
        if(self.option('types') == "xml"):
            self.run_filterxml()
        elif(self.option('types') == "hmm"):
            self.run_filter_hmm()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/annotation/output/blast_xml/'
        data = {
            "id": "filter_annot" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.filter_annot",
            "instant": False,
            "options": dict(
                xml=test_dir + "Trinity_vs_swissprot.xml",
                types="xml",
                evalue=0.000001,
                identity=90,
                similarity=90
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

        data = {
            "id": "filter_hm" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.filter_annot",
            "instant": False,
            "options": dict(
                hmm=test_dir + "pfam_domain",
                types="hmm",
                evalue=0.0001,
                identity=90,
                similarity=90
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
