# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
import os
import shutil
import subprocess
from mbio.packages.tool_lab.common import down_seq_files
import HTMLParser

class AniWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(AniWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {'name': 'method', 'type': 'string', 'default': 'ani'},
            {'name': 'ani_method', 'type': 'string', 'default': 'ANIm'}, #ANIm,ANIb,ANIblastall,TETRA
            {'name': 'ani_cluster', 'type': 'string', 'default': 'average'},
            {'name': 'evalue', 'type': 'string', 'default': '1e-3'},
            {'name': 'identity', 'type': 'string', 'default': '30'},
            {'name': 'alignment_length', 'type': 'string', 'default': '70'},
            {'name': 'aai_cluster', 'type': 'string', 'default': 'average'},
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        if self.option('project_data'):
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            for i in self.project_data['specimens']:
                self.samples[i['id']] = i['name']
            assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
            if os.path.exists(assemble_dir):
                shutil.rmtree(assemble_dir)
            (self.assemble_dir, self.analysis_type) = down_seq_files(self.project_data['my_type'],
                                                                     self.project_data['db_version'], assemble_dir,
                                                                     self._sheet.option("project_task_id"),
                                                                     self.samples)

    #def check_options(self):
    #    if not self.option("fasta_dir").is_set:
    #        raise OptionError("请传入序列文件夹！", code="")

    def run(self):
        if self.option("project_data"):
            self.fasta_sort_dir = self.download_file()
        else:
            self.fasta_dir_sort()
        self.run_ani()
        super(AniWorkflow, self).run()

    def run_ani(self):
        """
        Ani/AAi分析
        :return:
        """
        if self.option("method") == "ani":
            self.ani = self.add_module("tool_lab.genome_ani")
            opts = {
                "seq_dir": self.fasta_sort_dir,
                "method": self.option("ani_method"),
                "linkage": self.option("ani_cluster")
            }
            self.ani.set_options(opts)
            self.ani.on('end', self.set_db)
            self.ani.run()
        else:
            self.aai = self.add_module("tool_lab.genome_aai")
            if self.option("project_data"):
                file_ext = os.listdir(self.fasta_sort_dir)[0].split(".")[-1]
            else:
                file_ext = "fasta"
            opts = {
                "seq_dir": self.fasta_sort_dir,
                "evalue": self.option("evalue"),
                "identity": self.option("identity"),
                "aln_len": self.option("alignment_length"),
                "file_ext": file_ext,
                "linkage": self.option("aai_cluster")
            }
            self.aai.set_options(opts)
            self.aai.on('end', self.set_db)
            self.aai.run()


    def set_db(self):
        api_main = self.api.api("tool_lab.ani")
        if self.option("method") == "ani":
            with open(self.ani.output_dir + "/all.ANI_summary.xls") as f, open(self.ani.output_dir + "/ani_summary.xls","w") as t:
                raw_data = f.readlines()
                t.write(raw_data[0])
                for i in raw_data[1:]:
                    t.write(i.split("\t")[0])
                    for x in i.strip().split("\t")[1:]:
                        t.write("\t" + str(float(x) * 100))
                    t.write("\n")
            sumary_file = self.ani.output_dir + "/ani_summary.xls"
            cluster_file = self.ani.output_dir + "/hcluster.tre"
        else:
            sumary_file = self.aai.output_dir + "/aai_summary.xls"
            cluster_file = self.aai.output_dir + "/hcluster.tre"
        api_main.add_ani_detail(main_id=self.option('main_id'), detail_file=sumary_file, cluster_file=cluster_file)
        self.end()

    def fasta_dir_sort(self):
        self.fasta_sort_dir = self.work_dir + "/fasta_dir_sort"
        if os.path.exists(self.fasta_sort_dir):
            shutil.rmtree(self.fasta_sort_dir)
        os.mkdir(self.fasta_sort_dir)
        if os.path.exists(self.option("fasta_dir").prop['path'] + '/' + "list.txt"):
            sample_dict = {}
            with open(self.option("fasta_dir").prop['path'] + '/' + "list.txt") as t:
                for i in t.readlines()[1:]:
                    if i:
                        if len(i.split("\t")) < 2:
                            raise OptionError("list文件应为两列数据！", code="")
                        if i.split("\t") not in sample_dict.keys():
                            sample_dict[i.split("\t")[0]] = [i.split("\t")[1].strip()]
                        else:
                            sample_dict[i.split("\t")[0]].append(i.split("\t")[1].strip())
        else:
            raise OptionError("list文件不存在！", code="")
        for sample in sample_dict.keys():
            if len(sample_dict[sample]) > 1:
                cat_file = self.fasta_sort_dir + '/' + sample.replace(".","_") + '.fasta'
                os.system("cat {} >> {}".format(" ".join(sample_dict[sample]), cat_file))
            else:
                os.link(self.option("fasta_dir").prop['path'] + '/' + sample_dict[sample][0],self.fasta_sort_dir + '/' + sample.replace(".","_") + '.fasta')

    def download_file(self):
        """
        download file from s3
        :return:
        """
        path = ''
        if self.analysis_type in ['complete']:
            assemble_dir2 = os.path.join(self.work_dir, "assemble_dir2")
            if os.path.exists(assemble_dir2):
                shutil.rmtree(assemble_dir2)
            os.mkdir(assemble_dir2)
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir2, self.samples[sample] + ".fna")
                assemble_dir3 = os.path.join(self.assemble_dir, self.samples[sample])
                dir_list = os.listdir(assemble_dir3)
                for file2 in dir_list:
                    os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))
            path = assemble_dir2
        elif self.analysis_type in ['uncomplete']:
            path = self.assemble_dir
        return path

    def end(self):
        if self.option("method") == "ani":
            result_dir = self.add_upload_dir(self.ani.output_dir)
            os.remove(self.ani.output_dir + "/all.ANI_summary.xls")
            result_dir.add_relpath_rules([
                [".", "dir", "ANI分析结果目录", 0],
                ["./ani_summary.xls", "xls", "ANI分析结果表", 0],
                ["./hcluster.tre", "tre", "ANI分析聚类树", 0]
            ])
        else:
            result_dir = self.add_upload_dir(self.aai.output_dir)
            result_dir.add_relpath_rules([
                [".", "dir", "AAI分析结果目录", 0],
                ["./aai_summary.xls", "xls", "AAI分析结果表", 0],
                ["./hcluster.tre", "tre", "AAI分析聚类树", 0]
            ])
        super(AniWorkflow, self).end()