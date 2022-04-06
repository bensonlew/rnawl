# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
import json
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class TfbsPredictAgent(Agent):
    """
    tfbs_predict description
    """
    def __init__(self, parent):
        super(TfbsPredictAgent, self).__init__(parent)
        options = [
            {'name': 'motifs_user', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'tf_selected', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'gtf', 'type': 'string', 'default': None, 'format': None},
            {'name': 'genome', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'gtfParser', 'type': 'string', 'default': '/mnt/ilustre/users/isanger/sg-users/deqing/mbio/packages/transcription_factor/parseGTF.pl', 'format': 'None'},
            {'name': 'bedtools', 'type': 'string', 'default': '/mnt/ilustre/users/isanger/app/bioinfo/seq/bedtools-2.25.0/bin/bedtools', 'format': 'None'},
            {'name': 'genes', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'fimo', 'type': 'string', 'default': '/mnt/ilustre/users/isanger/app/bioinfo/rna/meme_4.12.0/bin/fimo', 'format': 'None'},
            {'name': 'seqdb', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'thresh', 'type': 'float', 'default': '0.0001', 'format': 'None'},
            {'name': 'qv_thresh', 'type': 'int', 'default': '0', 'format': 'None'},
            {'name': 'exp_matrix', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'geneid2tfid', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'corr_cutoff', 'type': 'float', 'default': '0.5', 'format': 'None'},
            {'name': 'corr_pvalue', 'type': 'float', 'default': '0.05', 'format': 'None'},
            {'name': 'corr_use_padjust', 'type': 'int', 'default': '0', 'format': 'None'},
            {'name': 'backward', 'type': 'int', 'default': '450', 'format': 'None'},
            {'name': 'forward', 'type': 'int', 'default': '50', 'format': 'None'},
            {'name': 'motif_species', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'motif_species_type', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'organism', 'type': 'string', 'default': 'None', 'format': 'None'},
            {'name': 'genome_id', 'type': 'string', 'default': 'None', 'format': 'None'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 3
        if len(self.option("tf_selected").split(",")) > 15:
            self._memory = "{}G".format('30')
        else:
            self._memory = "{}G".format('20')


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfbsPredict"]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
        ])
        """
        super(TfbsPredictAgent, self).end()


class TfbsPredictTool(Tool):
    """
    tfbs_predict description
    """
    def __init__(self, config):
        super(TfbsPredictTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'miniconda2/bin/python'
        self.tfbs_predict = self.config.PACKAGE_DIR + '/transcription_factor/tfbs_predict.py'
        self.gtfParser = self.config.PACKAGE_DIR + '/transcription_factor/parseGTF.pl'
        self.fimo = software_dir + '/bioinfo/rna/meme_4.12.0/bin/fimo'
        self.bedtools = software_dir + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_tfbs_predict(self):
        cmd = '{} {} '.format(self.python_path, self.tfbs_predict)
        # prepare input motif file
        motif_files = list()
        geneid_tfid_new = list()
        if self.option("motifs_user").lower() != "none":
            motif_files.append(self.option("motifs_user"))
        if self.option("tf_selected").lower() != "none":
            geneid_tfid_pair_list = self.option("tf_selected").split(',')
            self.logger.info(geneid_tfid_pair_list)
            geneid_tfid = [x.split("|") for x in geneid_tfid_pair_list]
            for g, tf in geneid_tfid:
                motif_db = self.config.SOFTWARE_DIR + "/database/TFDB/{}_tf_motif/".format(self.option("motif_species_type"))
                motif_path = motif_db + '{}.meme'.format(tf)
                if os.path.exists(motif_path):
                    motif_files.append(motif_path)
                    geneid_tfid_new.append([g,tf])
                else:
                    if '.' in tf:
                        tmp_list = tf.split(".")
                        for i in range(0,len(tmp_list)):
                            tmp_x = '.'.join(tmp_list[0:len(tmp_list)-i-1])
                            motif_path = motif_db + '{}.meme'.format(tmp_x)
                            if os.path.exists(motif_path):
                                motif_files.append(motif_path)
                                geneid_tfid_new.append([g,tmp_x])
                                continue

            # replace original geneid2tfid
            with open(self.option("geneid2tfid"), 'w') as fw:
                fw.write('gene_id\ttf_id\n')
                for g, t in geneid_tfid_new:
                    fw.write("{}\t{}\n".format(g, t))
        motif_path = self.work_dir + "/query_motif.meme"
        self.merge_meme(motif_files, motif_path, has_header=False)
        # prepare genome and annot file
        if self.option("genome_id") == "None":
            with open(self.config.SOFTWARE_DIR + '/database/Genome_DB_finish/annot_species.v2.json') as f:
                annot_dict = json.load(f)
            gtf_file = self.config.SOFTWARE_DIR + '/database/Genome_DB_finish/' + annot_dict[self.option('organism')]['gtf']
            if self.option('gtf'):
                all_gtf = self.work_dir + '/all.gtf'
                os.system('cat {} {} > {}'.format(gtf_file, self.option('gtf'), all_gtf))
                gtf_file = all_gtf
            genome = self.config.SOFTWARE_DIR + '/database/Genome_DB_finish/' + annot_dict[self.option('organism')]['dna_fa']
        else:
            db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
            col = db["sg_genome_db"]
            db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
            genome_info = col.find_one({"genome_id": self.option("genome_id")})
            gtf_file = os.path.join(db_path, genome_info["gtf"])
            genome = os.path.join(db_path, genome_info["dna_fa"])
        #
        cmd += '-{} {} '.format("motifs", motif_path)
        if self.option("seqdb") != "none":
            cmd += '-{} {} '.format("seqdb", self.option("seqdb"))
        else:
            cmd += '-{} {} '.format("gtf", gtf_file)
            cmd += '-{} {} '.format("genome", genome)
            cmd += '-{} {} '.format("gtfParser", self.gtfParser)
            cmd += '-{} {} '.format("bedtools", self.bedtools)
            if self.option("genes").lower() not in ['all', 'none', 'refall']:
                cmd += '-{} {} '.format("genes", self.option("genes"))
            cmd += '-{} {} '.format("fimo", self.fimo)
            cmd += '-{} {} '.format("backward", self.option("backward"))
            cmd += '-{} {} '.format("forward", self.option("forward"))
        cmd += '-{} {} '.format("thresh", self.option("thresh"))
        cmd += '-{} {} '.format("qv_thresh", self.option("qv_thresh"))
        if self.option("exp_matrix").lower() != "none":
            cmd += '-{} {} '.format("exp_matrix", self.option("exp_matrix"))
        cmd += '-{} {} '.format("geneid2tfid", self.option("geneid2tfid"))
        cmd += '-{} {} '.format("corr_cutoff", self.option("corr_cutoff"))
        cmd += '-{} {} '.format("corr_pvalue", self.option("corr_pvalue"))
        cmd += '-{} {} '.format("corr_use_padjust", self.option("corr_use_padjust"))
        cmd_name = 'tfbs_predict'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708301")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33708302")

    def set_output(self):
        target_files = glob.glob(self.work_dir + '/tf*.xls')
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(TfbsPredictTool, self).run()
        self.run_tfbs_predict()
        self.set_output()
        self.end()

    @staticmethod
    def merge_meme(file_list, out_path, has_header=True):
        """cat files"""
        if not file_list:
            return
        with open(out_path, 'w') as fw:
            fr = open(file_list[0])
            fw.write(fr.read())
            fr.close()
            for each in file_list[1:]:
                fr = open(each)
                if has_header:
                    _ = fr.readline()
                fw.write(fr.read())
                fr.close()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "TfbsPredict" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna.tfbs_predict",
            "instant": False,
            "options": dict(
                motifs_user="None",
                # seqdb="None",
                thresh="0.0001",
                qv_thresh="0",
                exp_matrix="None",
                tfid2geneid="None",
                corr_cutoff="0.5",
                corr_pvalue="0.05",
                corr_use_padjust="0",
                backward="450",
                forward="50",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
