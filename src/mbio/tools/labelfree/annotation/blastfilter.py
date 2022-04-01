# coding=utf-8
import os
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
__author__ = 'liubinxu'


class BlastfilterAgent(Agent):
    """
    根据阈值过滤blast结果，用于注释重运行
    version 1.0
    author: 刘彬旭
    last_modify: 2017.12.17
    """
    def __init__(self, parent):
        super(BlastfilterAgent, self).__init__(parent)
        options = [
            {"name": "blast", "type": "infile", "format": "denovo_rna_v2.blast_xml"},
            {"name": "hmm", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "method", "type": "string", "default": "blast"},
            # 比对数据库方法blast或hmm
            {"name": "evalue", "type": "float", "default": 1e-3},
            # evalue过滤阈值
            {"name": "identity", "type": "float", "default": 0},
            # identity过滤阈值
            {"name": "similarity", "type": "float", "default": 0},
            # similarity过滤阈值
            {"name": "outxml", "type": "outfile", "format": "denovo_rna_v2.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "denovo_rna_v2.common"},  # 输出格式为6时输出
        ]
        self.add_option(options)
        # self.step.add_steps('blast_filter')
        # self.on('start', self.step_start)
        # self.on('end', self.step_end)

    def step_start(self):
        self.step.blast_filter.start()
        self.step.update()

    def step_end(self):
        self.step.blast_filter.finish()
        self.step.update()

    def check_options(self):
        if self.option("method") == "blast":
            if not self.option("blast").is_set:
                raise OptionError("必须设置参数query")
        elif self.option("method") == "hmm":
            if not self.option("hmm").is_set:
                raise OptionError("必须设置参数hmm")
        else:
            raise OptionError('method {}必须为blast或hmm'.format(self.option('method')))

        if not (self.option('evalue') >= 0 and self.option('evalue') <= 0.001) :
            raise OptionError('E-value值设定必须为[0-0.001]之间：{}'.format(self.option('evalue')))
        if not (self.option('identity') >= 0 and self.option('identity') <= 100) :
            raise OptionError('Identity值设定必须为[0-100]之间：{}'.format(self.option('identity')))
        if not (self.option('similarity') >= 0 and self.option('similarity') <= 100) :
            raise OptionError('similarity值设定必须为[0-100]之间：{}'.format(self.option('similarity')))

        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '1G'

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
        super(BlastfilterAgent, self).end()


class BlastfilterTool(Tool):
    def __init__(self, config):
        super(BlastfilterTool, self).__init__(config)

    def run_filter_blast(self):
        """
        过滤blast结果

        """
        basename = os.path.basename(self.option('blast').prop['path']) + ".filter.xml"
        self.option('blast').filter_blast_xml(basename, self.option('evalue'), self.option('identity'), self.option('similarity'))

    def run_filter_hmm(self):
        self.logger.info("hmm 文件 {}".format(os.path.basename(self.option('hmm').prop['path'])))
        basename = os.path.basename(self.option('hmm').prop['path'])
        pfam_table = pd.read_table(self.option('hmm').prop['path'], header=0)
        pfam_select = pfam_table[pfam_table['PfamEnd'] < self.option('evalue')]
        pfam_select.to_csv(basename + '.filter.xls', sep="\t",index=False)
    def get_inviro(self):
        cmd = 'env'
        cmd_name = 'transrate'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()

    def set_output(self):
        if(self.option('method') == "blast"):
            filter_xml = os.path.basename(self.option('blast').prop['path']) + ".filter.xml"
            if os.path.exists(filter_xml):
                os.remove(filter_xml)
            os.link(filter_xml, os.path.join(self.output_dir, filter_xml))
            self.option('outxml', os.path.join(self.output_dir, filter_xml))
        elif(self.option('method') == "hmm"):
            filter_xml = os.path.basename(self.option('hmm').prop['path']) + ".filter.xls"
            if os.path.exists(filter_xml):
                os.remove(filter_xml)
            os.link(filter_xml, os.path.join(self.output_dir, filter_xml))
            self.option('outtable', os.path.join(self.output_dir, filter_xml))
    def run(self):
        """
        运行
        :return:
        """
        super(BlastfilterTool, self).run()
        self.get_inviro()
        if(self.option('method') == "blast"):
            self.run_filter_blast()
        elif(self.option('method') == "hmm"):
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
            "id": "Blastfilter" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.annotation.blastfilter",
            "instant": False,
            "options": dict(
                blast=test_dir + "Trinity_vs_swissprot.xml",
                method="blast",
                evalue=0.000001,
                identity=90,
                similarity=90
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

        data = {
            "id": "Blastfilter" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.annotation.blastfilter",
            "instant": False,
            "options": dict(
                hmm=test_dir + "pfam_domain",
                method="hmm",
                evalue=0.000001,
                identity=90,
                similarity=90
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
