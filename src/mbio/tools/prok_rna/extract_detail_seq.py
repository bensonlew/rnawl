# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import json
from collections import OrderedDict
from Bio import SeqIO



class ExtractDetailSeqAgent(Agent):
    """
    This script is used to extracts mapped/unmapped reads from input bam/sam file.
    """
    def __init__(self, parent):
        super(ExtractDetailSeqAgent, self).__init__(parent)
        options = [
            {"name": "seq_detail_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "geneset_extract", "type": "infile", "format": 'ref_rna_v2.common'},
            {"name": "extract_info", "type": "string", "default": None},
        ]
        self.add_option(options)
        self.step.add_steps("extract")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract.finish()
        self.step.update()

    def check_options(self):
        # if not self.option('input_file').is_set:
        #     raise OptionError('SAM/BAM文件必须输入')
        # if self.option('seq_type') not in ["PE", "SE"]:
        #     raise OptionError('测序类型参数输入有误')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(ExtractDetailSeqAgent, self).end()


class ExtractDetailSeqTool(Tool):
    def __init__(self, config):
        super(ExtractDetailSeqTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'seq_extract': os.path.join(self.config.PACKAGE_DIR, 'denovo_rna_v2/seq_extract.py')
        }
        self.extract_info = ""



    def run(self):
        """
        运行
        :return:
        """
        super(ExtractDetailSeqTool, self).run()
        self.run_tool()
        self.set_output()
        self.end()

    def run_tool(self):
        self.extract_info = json.loads(self.option("extract_info"), object_pairs_hook=OrderedDict)
        with open(self.option('geneset_extract').path, "r") as g:
            all_gene = [line.strip() for line in g.readlines()]
        if os.path.exists(os.path.join(self.option('seq_detail_dir').path, 'cds.fa')):
            cds_fa = os.path.join(self.option('seq_detail_dir').path, 'cds.fa')
        elif os.path.exists(os.path.join(self.option('seq_detail_dir').path, 'Sequence_database', 'cds.fa')):
            cds_fa = os.path.join(self.option('seq_detail_dir').path, 'Sequence_database', 'cds.fa')
        else:
            self.set_error("文件cds.fa不存在")
        if os.path.exists(os.path.join(self.option('seq_detail_dir').path, 'cds.faa')):
            cds_faa = os.path.join(self.option('seq_detail_dir').path, 'cds.faa')
        elif os.path.exists(os.path.join(self.option('seq_detail_dir').path, 'Sequence_database', 'cds.faa')):
            cds_faa = os.path.join(self.option('seq_detail_dir').path, 'Sequence_database', 'cds.faa')
        else:
            self.set_error("文件cds.faa不存在")
        cds_fa_dict = dict()
        cds_faa_dict = dict()
        record_fa = SeqIO.parse(cds_fa, 'fasta')
        for i in record_fa:
            cds_fa_dict[i.id] = str(i.seq)
        record_faa = SeqIO.parse(cds_faa, 'fasta')
        for j in record_faa:
            cds_faa_dict[j.id] = str(j.seq)
        if 'cds' in self.extract_info['seq_type']:
            with open(os.path.join(self.output_dir, "gene_seqs.fnn"), "w") as f:
                for gene_id in all_gene:
                    try:
                        f.write(">{}".format(gene_id) + "\n" + cds_fa_dict[gene_id] + "\n")
                    except:
                        pass
        if 'protein' in self.extract_info['seq_type']:
            with open(os.path.join(self.output_dir, "gene_seqs.faa"), "w") as f:
                for gene_id in all_gene:
                    try:
                        f.write(">{}".format(gene_id) + "\n" + cds_faa_dict[gene_id] + "\n")
                    except:
                        pass

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # try:
        #     output_files = glob.glob(self.work_dir + "/*.fq")
        #     for file_path in output_files:
        #         os.link(file_path, os.path.join(self.output_dir, os.path.basename(file_path)))
        # except Exception as e:
        #     self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/medical_transcriptome/test_files'
        data = {
            "id": "ExpPca" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "medical_transcriptome.extract_detail_seq",
            "instant": False,
            "options": dict(
                seq_detail_dir ="/mnt/ilustre/users/sanger-dev/workspace/20201104/DownloadDetailSeq_gbnijnfcuslii6rieegpbdgpt5_9045_87/remote_input/seq_detail_dir/SequenceDetail/",
                geneset_extract = "/mnt/ilustre/users/sanger-dev/workspace/20201104/DownloadDetailSeq_gbnijnfcuslii6rieegpbdgpt5_9045_87/geneset_extract_gene.list",
                level = "G",
                extract_info = json.dumps({"seq_type":["transcript", "cds", "pep"]}).replace('"', '\\"'),
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
