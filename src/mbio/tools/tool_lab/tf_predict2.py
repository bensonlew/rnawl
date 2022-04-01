# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import csv
# import pandas as pd
__author__ = 'liubinxu'


class TfPredict2Agent(Agent):
    """
    tf_predict  使用JASPAR Profile Inference 预测蛋白中的转录因子
    """
    def __init__(self, parent):
        super(TfPredict2Agent, self).__init__(parent)
        options = [
            {'name': 'taxon', 'type': 'string', 'default': 'all', 'format': 'None'},
            # "fungi", "insects", "nematodes", "plants", "vertebrates"
            {'name': 'output', 'type': 'string', 'default': 'tf_predict.xls', 'format': 'None'},
            {'name': 'seqfile', 'type': 'infile', 'default': 'None', 'format': 'sequence.fasta'},
            {'name': 'thread', 'type': 'int', 'default': 10},
            {'name': 'evalue', 'type': 'float', 'default': 1e-5},
            {'name': 'organism', 'type': 'string', 'default': 'all'},
            {'name': 'rost', 'type': 'int', 'default': 5}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = int(self.option('thread'))
        self._memory = "{}G".format('20')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfPredict2"]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
        ])
        """
        super(TfPredict2Agent, self).end()


class TfPredict2Tool(Tool):
    """
    tf_predict description
    """
    def __init__(self, config):
        super(TfPredict2Tool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.jaspar_profile = software_dir + '/bioinfo/rna/JASPAR-profile-inference-master/infer_profile.py'
        self.jaspar_db =  software_dir + '/bioinfo/rna/JASPAR-profile-inference-master/'
        self.jaspar_file_dir = software_dir + "/database/JASPAR2020"
        self.hmmscan = software_dir + '/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries'
        self.blast = software_dir + '/bioinfo/align/ncbi-blast-2.3.0+/bin'
        self.python3 = 'program/miniconda3/bin/python3'
        self.python3_lib = software_dir + '/program/miniconda3/lib'

        self.set_environ(PATH=self.hmmscan)
        self.set_environ(PATH=self.blast)
        self.set_environ(LD_LIBRARY_PATH=self.python3_lib)
        self.jaspar_dict =dict()
        self.jump_finished = True

    def read_table(self):
        jaspar_files = glob.glob(self.jaspar_file_dir + '/*.xls')
        for each in jaspar_files:
            with open(each, 'r') as f:
                for dic in csv.DictReader(f, delimiter='\t'):
                    # print "dic is", dic
                    self.jaspar_dict[dic["Matrix ID"]] = dic


    def split_fasta(self):
        self.split_fasta_list = list()
        with open(self.option("seqfile").prop['path']) as f:
            seq_id = 0
            split_num = 0
            self.split_fasta_list.append(os.path.basename(f.name) + '_' + str(split_num))
            out_name = os.path.basename(f.name) + '_' + str(split_num)
            file_split = {split_num: open(out_name, 'w')}
            for line in f:
                if line.startswith('>'):
                    split_num = seq_id/3000
                    seq_id = seq_id + 1
                    if split_num not in file_split.keys():
                        self.split_fasta_list.append(os.path.basename(f.name) + '_' + str(split_num))
                        file_split[split_num-1].close()
                        out_name = os.path.basename(f.name) + '_' + str(split_num)
                        file_split[split_num] = open(out_name,'w')
                else:
                    line = line.replace('*', '').replace('.', '')

                file_split[split_num].write(line)
            else:
                file_split[split_num].close()



    def run_tf_predict(self, fa):
        num = fa.split("_")[-1]
        cmd = '{} {} '.format(self.python3, self.jaspar_profile)
        cmd += '{} {} '.format("--fasta-file", fa)
        cmd += '{} {} '.format("--taxon", self.option("taxon"))
        cmd += '{} {} '.format("--threads", self.option("thread"))
        cmd += '{} {} '.format("--rost", self.option("rost"))
        cmd += '{} {} '.format("--output-file", self.option("output") + '_' + num)
        cmd += '--dummy-dir ./tmp '
        cmd += '--files-dir {}files/ '.format(self.jaspar_db)
        cmd += '--models-dir {}models/ '.format(self.jaspar_db)
        print cmd
        cmd_name = "tf_profile_{}".format(str(num))
        print cmd_name
        if self.jump_finished and os.path.exists(cmd_name + '.o'):
            return
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
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        # target_file = self.work_dir + '/tf_predict.xls'
        target_files = glob.glob(self.work_dir + '/tf_predict.xls*')
        target_detail_file = self.work_dir + '/tf_predict_detail.xls'
        is_header =True
        with open(target_detail_file, 'w') as fo:
            for target_file in target_files:
                with open(target_file, 'r') as f:
                    header = f.readline()
                    if is_header == True:
                        fo.write(header.strip("\n") + "\tClass" + "\tFamily" + "\tUniprot ID\thit_link\n")
                        is_header = False
                    else:
                        pass
                    for line in f:
                        cols = line.strip("\n").split("\t")
                        ja_class = self.jaspar_dict[cols[2]].get("Class", "")
                        ja_family = self.jaspar_dict[cols[2]].get("Family", "")
                        ja_uniprot = self.jaspar_dict[cols[2]].get("Uniprot ID", "")
                        link = "http://jaspar.genereg.net/matrix/{}/".format(cols[2])
                        ja_species = self.jaspar_dict[cols[2]].get("Species", "").lower().split(";")
                        if self.option("organism").lower() not in ['all', 'none', 'unknown']:
                            if self.option("organism").lower().replace("_", " ") in ja_species:
                                fo.write(line.strip("\n") + "\t" + ja_class + "\t" + ja_family + "\t" + ja_uniprot + "\t" + link + "\n")
                        else:
                            fo.write(line.strip("\n") + "\t" + ja_class + "\t" + ja_family + "\t" + ja_uniprot + "\t" + link + "\n")


        for each in [target_detail_file]:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(TfPredict2Tool, self).run()
        self.read_table()
        self.split_fasta()
        for fa in self.split_fasta_list:
            self.run_tf_predict(fa)
        # self.run_tf_predict()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "TfPredict2" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.tf_predict2",
            "instant": False,
            "options": dict(
                taxon="plants",
                seqfile="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/soft/JASPAR-profile-inference-master/test.pep.fa"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
