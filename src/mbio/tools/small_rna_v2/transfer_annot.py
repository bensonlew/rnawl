# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
import shutil
import tarfile
import unittest
from biocluster.config import Config
from biocluster.agent import Agent
from biocluster.tool import Tool
from mainapp.models.mongo.denovo_rna_v2 import DenovoRnaV2
from mainapp.models.mongo.ref_rna_v2 import RefRnaV2
from mbio.packages.ref_rna_v2.functions import tryforgood
from biocluster.api.file.lib.transfer import MultiFileTransfer

class TransferAnnotAgent(Agent):
    """
    last_modify: 2017.07.27
    get annotation file from ref_rna_v2 denovo_rna_v2
    """

    def __init__(self, parent):
        super(TransferAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'task_id', 'type': "string", 'default': None},
            {'name': 'task_type', 'type': "string", 'default': 'denovo_rna_v2'},
            {'name': 'outdir', 'type': 'outfile', 'format': 'small_rna.common_dir'},
            {'name': 'genome', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'assembly_gtf', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'assembly_fa', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'assembly_gene2trans', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'assembly_t2g', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'anno_detail', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'anno_class', 'type': 'outfile', 'format': 'small_rna.common_dir'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(TransferAnnotAgent, self).end()


class TransferAnnotTool(Tool):
    def __init__(self, config):
        super(TransferAnnotTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.python = 'miniconda2/bin/python'
        self.cog_summary_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/cog_summary.py')
        self.down_py =  self.config.PACKAGE_DIR + '/rna/download_froms3.py'

    def run(self):
        super(TransferAnnotTool, self).run()
        self.get_annot_dir()
        '''
        if os.path.exists(os.path.basename(self.annot_dir.rstrip("/")) + '/'):
            pass
        else:
        '''
        self.download_file()
        self.cog_summary()
        self.set_output()
        self.end()

    def get_annot_dir(self):
        if self.option("task_type") == "denovo_rna_v2":
            denovo_rna_v2 = DenovoRnaV2(None)

            task_search = {"$regex": str(self.option("task_id"))}
            stat_info = denovo_rna_v2.get_main_info_by_record("sg_annotation_stat", task_id=task_search, type="origin")
            task_info = denovo_rna_v2.get_main_info_by_record("sg_task", task_id=task_search)
            self.logger.info("task_info {}, stat_info {}".format(task_info, stat_info))
            self.annot_dir = stat_info["result_dir"]
            self.assemble_fa = task_info["assemble_fa"]
            self.assemble_t2g = task_info["assemble_t2g"]

        elif self.option("task_type") == "ref_rna_v2":
            ref_rna_v2 = RefRnaV2(None)
            task_search = {"$regex": str(self.option("task_id"))}
            stat_info = ref_rna_v2.get_main_info_by_record("sg_annotation_stat", task_id=task_search, type="origin")
            task_info = ref_rna_v2.get_main_info_by_record("sg_task", task_id=task_search)
            self.logger.info("task_info {}, stat_info {}".format(task_info, stat_info))
            self.annot_dir = stat_info["result_dir"]
            if "assemble_fa" in task_info:
                self.assemble_fa = task_info["assemble_fa"]
                self.assemble_gtf = task_info["as_gtf"]
                self.assemble_t2g = task_info["assemble_t2g"]
            else:
                genome_id = task_info["genome_id"]
                genome_info = ref_rna_v2.get_main_info_by_record("sg_genome_db", genome_id=genome_id)
                self.assemble_fa = self.config.SOFTWARE_DIR + "/database/" + genome_info["dna_fa"]
                self.assemble_gtf = self.config.SOFTWARE_DIR + "/database/" + genome_info["gtf"]
                self.assemble_t2g = self.config.SOFTWARE_DIR + "/database/" + genome_info["g2t2p"]
            self.genome_fa = self.config.SOFTWARE_DIR + "/" + task_info["ref_genome"].split("/app/")[1]
        else:
            self.set_error("不支持此类型项目")
        if not self.annot_dir.endswith("/"):
            self.annot_dir += "/"

    def download_file(self):

        if os.path.exists(os.path.basename(self.annot_dir.rstrip("/")) + '/'):
             shutil.rmtree(os.path.basename(self.annot_dir.rstrip("/")))


        cmd = "{} {} {} {}".format(self.python_path,
                                   self.down_py,
                                   self.annot_dir,
                                   os.path.basename(self.annot_dir.rstrip("/")) + '/')
        command = self.add_command("download", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("download 文件正确")
        else:
            self.set_error("download 文件错误", code = "32001003")

        self.download_from_s3(self.assemble_fa,  os.path.basename(self.assemble_fa))
        self.option("assembly_fa", os.path.join(self.work_dir, os.path.basename(self.assemble_fa)))
        self.download_from_s3(self.assemble_t2g,  os.path.basename(self.assemble_t2g))
        self.option("assembly_t2g", os.path.join(self.work_dir, os.path.basename(self.assemble_t2g)))

        if self.option("task_type") == "ref_rna_v2":
            annot_dir = os.path.join(self.work_dir, os.path.basename(self.annot_dir.rstrip("/")), "allannot_class")
        else:
            annot_dir = os.path.join(self.work_dir, os.path.basename(self.annot_dir.rstrip("/")))
        self.option("anno_class", annot_dir)
        self.option("assembly_gene2trans", annot_dir + "/all_tran2gene.txt")
        self.option("anno_detail", annot_dir +  "/all_annot.xls")

        if self.option("task_type") == "ref_rna_v2":
            self.download_from_s3(self.assemble_gtf,  os.path.basename(self.assemble_gtf))
            self.option("assembly_gtf", os.path.join(self.work_dir, os.path.basename(self.assemble_gtf)))
            self.option("genome", self.genome_fa)
        self.wait()

    def cog_summary(self):
        if self.option("task_type") == "denovo_rna_v2":
            annot_dir = self.work_dir + '/Annotation'
        else:
            annot_dir = self.work_dir + '/Annotation/allannot_class'
        self.txpt_summary_tsv =  annot_dir + '/cog/summary.T.tsv'
        self.gene_summary_tsv =  annot_dir + '/cog/summary.G.tsv'
        self.cog_tsv = annot_dir + '/cog/cog_list_tran.xls'
        self.trans2gene = annot_dir + '/all_tran2gene.txt'
        with open(self.trans2gene, 'r') as fin, open(self.work_dir + "/longest_t2g", 'w') as fo:
            for line in fin:
                cols = line.strip().split("\t")
                if cols[2] == "yes":
                    fo.write("{}\t{}\n".format(cols[0], cols[1]))

        cmd = '{} {}'.format(self.python, self.cog_summary_py)
        cmd += ' --t2g {}'.format(self.work_dir + "/longest_t2g")
        cmd += ' --cog {}'.format(self.cog_tsv)
        cmd += ' --txpt {}'.format(self.txpt_summary_tsv)
        cmd += ' --gene {}'.format(self.gene_summary_tsv)
        cmd_name = 'run_cog_summary'

        command = self.add_command(cmd_name, cmd).run()
        self.wait()

        if command.return_code == 0:
            self.logger.info("cog summary 成功")
        else:
            self.set_error("cog summary错误")


        # self.run_code(cmd_name, cmd)



    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run %s, abord', variables=(n), code="33710702")


    def set_output(self):
        pass
        '''
        self.logger.info(self.option('indir').path)

        for name in os.listdir(self.option('indir').path):
            src = os.path.join(self.option('indir').path, name)
            dst = os.path.join(self.output_dir, name)
            self.logger.info(name)
            self.logger.info(src)
            self.logger.info(dst)
            if os.path.isfile(src):
                self.link2outputdir(src, dst)
            elif os.path.isdir(src):
                self.move2outputdir(src, name)
        self.option('outdir').set_path(self.output_dir)
        '''


    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        self.logger.info("start download {} {}".format(path, to_path))

        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            os.system('ln -s {} {}'.format(path.rstrip("/"), to_path.rstrip("/")))
        else:
            try:
                transfer2 = MultiFileTransfer()
                self.logger.info("start from {} to {}".format(path, to_path))
                transfer2.add_download(path, to_path, base_name=path)
                self.logger.info("files is {}".format(transfer2._files.items()))
                transfer2.perform()
            except:
                return False
        return to_path




class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.

                'task_id': 'tsg_38042',
                'task_type': 'ref_rna_v2'
                'task_id': 'tsg_38052',
                'task_type': 'denovo_rna_v2'
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'transfer_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'small_rna_v2.transfer_annot',
            'instant': False,
            'options': {
                'task_id': 'tsg_38042',
                'task_type': 'ref_rna_v2'

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
