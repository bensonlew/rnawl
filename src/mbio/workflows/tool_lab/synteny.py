# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagbin.common_function import link_dir
import os
import subprocess
import shutil
from mbio.packages.tool_lab.common import down_seq_files
import HTMLParser

class SyntenyWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self.samples = {}
        self._sheet = wsheet_object
        super(SyntenyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "ref", "type": "infile", "format": "sequence.fasta"},
            {'name': 'method', 'type': 'string', 'default': 'nucmer'}, #nucmer/promer
            {'name': 'project_data', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.mummer = self.add_tool('tool_lab.mummer')
        self.plot = self.add_tool('tool_lab.synteny_plot')
        if self.option('project_data'):
            self.project_data = eval(HTMLParser.HTMLParser().unescape(self._sheet.option('project_data')))
            self.logger.info(self._sheet.option('project_data'))
            self.logger.info(self.project_data)
            for i in self.project_data['specimens']:
                self.samples[i] = self.project_data['specimens'][i]
            assemble_dir = os.path.join(self._sheet.work_dir, "assemble_dir")
            if os.path.exists(assemble_dir):
                shutil.rmtree(assemble_dir)
            (self.assemble_dir, self.analysis_type) = down_seq_files(self.project_data['my_type'],
                                                                     self.project_data['db_version'], assemble_dir,
                                                                     self._sheet.option("project_task_id"),
                                                                     self.samples)

    def check_options(self):
        if self.option("project_data"):
            pass
        elif not self.option("query").is_set and not self.option("ref").is_set:
            raise OptionError("请传入查询基因组序列和参考基因组序列！", code="")

    def run(self):
        self.get_seq()
        self.run_mummer()
        super(SyntenyWorkflow, self).run()

    def run_plot(self):
        options = {
            'result': self.mummer.output_dir + "/mummer.output.xls",
            'ref': self.ref_name,
            'samples': self.query_name,
            'seq_dir': self.work_dir + "/fasta_dir/",
        }
        self.plot.set_options(options)
        self.plot.on('end', self.set_db)
        self.plot.run()

    def run_mummer(self):
        options = {
            'mummer': self.option('method'),
            'ref': self.ref_name,
            'samples': self.query_name,
            'seq_dir': self.work_dir + "/fasta_dir/",
        }
        self.mummer.set_options(options)
        self.mummer.on('end', self.set_db)
        self.mummer.run()

    def get_seq(self):
        if self.option("project_data"):
            self.download_file()
            query = self.project_data['sample']["query"]
            self.query_name = query
            ref = self.project_data['sample']["ref"]
            self.ref_name =  ref
            if os.path.exists(self.work_dir + "/fasta_dir"):
                shutil.rmtree(self.work_dir + "/fasta_dir")
            os.mkdir(self.work_dir + "/fasta_dir")
            os.link(os.path.join(self.work_dir, "assemble_dir", query + ".fna"), self.work_dir + "/fasta_dir/" + self.query_name+".fna")
            os.link(os.path.join(self.work_dir, "assemble_dir", ref + ".fna"), self.work_dir + "/fasta_dir/" + self.ref_name+".fna")
        else:
            if os.path.exists(self.work_dir + "/fasta_dir"):
                shutil.rmtree(self.work_dir + "/fasta_dir")
            os.mkdir(self.work_dir + "/fasta_dir")
            self.query_name = self.option("query").prop["path"].split("/")[-1].split(".")[0]
            self.logger.info(self.query_name)
            self.ref_name = self.option("ref").prop["path"].split("/")[-1].split(".")[0]
            if self.query_name == self.ref_name:
                self.query_name = self.query_name + "_query"
                self.ref_name = self.ref_name + "_ref"
            os.link(self.option("query").prop["path"], self.work_dir + "/fasta_dir/" + self.query_name + ".fna")
            os.link(self.option("ref").prop["path"], self.work_dir + "/fasta_dir/" + self.ref_name + ".fna")

    def set_db(self):
        mm_api = self.api.api('tool_lab.synteny')
        if self.option('main_id'):
            main_id = self.option('main_id')
        else:
            main_id = mm_api.add_syteney_main()
        self.logger.info('开始导表！')
        filepath = os.path.join(self.mummer.output_dir, 'mummer.output.xls')
        os.rename(filepath,self.work_dir + "/tmp_mummer.xls")
        with open(self.work_dir + "/tmp_mummer.xls") as f, open(filepath,"w") as t:
            data = f.readlines()
            t.write(data[0])
            for i in data[1:]:
                if i.strip().split("\t")[9] == i.strip().split("\t")[10]:
                    pass
                else:
                    t.write(i)
        os.remove(self.work_dir + "/tmp_mummer.xls")
        mm_api.add_detail(filepath, 'synteny_detail', main_id,  'synteny_id', type="synteny")
        #mm_api.add_detail(filepath, 'synteny_detail', main_id, 'synteny_id', type="synteny_circle")
        #mm_api.add_detail(filepath, 'synteny_detail', main_id, 'synteny_id', type="dotplot")
        mongo_keys = {0: 'chr', 1: 'len'}
        len_dir = {}
        for sp in os.listdir(self.mummer.work_dir):
            if sp.endswith('.seq_len.xls'):
                if sp.split('.seq_len.xls')[0] == self.query_name or sp.split('.seq_len.xls')[0] in self.query_name:
                   len_dir[self.query_name] = sp
                elif sp.split('.seq_len.xls')[0] == self.ref_name or sp.split('.seq_len.xls')[0] in self.ref_name:
                    len_dir[self.ref_name] = sp
                else:
                    self.logger.info('{}len文件没有对应样本！'.format(sp))
        for sampel in [self.query_name,self.ref_name]:
            sp = len_dir[sampel]
            s = sp.split('.seq_len.xls')[0]
            fp = os.path.join(self.mummer.work_dir, sp)
            mm_api.add_detail(fp, 'synteny_detail',
                                main_id, 'synteny_id', mongo_keys=mongo_keys,
                                tag_key=['sp', ], tag_value=[s, ], header=False, type="band_synteny")
        mm_api.updata_status(main_id)
        self.logger.info('导表结束！')
        os.link(filepath, self.output_dir + '/syteney.xls')
        """
        link_dir(self.plot.output_dir,self.output_dir)
        remote_dir = self.sheet.output
        self.logger.info(self.sheet.output)
        if not remote_dir:
            remote_dir = self.output_dir
        mm_api = self.api.api('tool_lab.common')
        mm_api.add_synteny(name="synteny",dual_synteny_path=remote_dir+"/dual_synteny.PNG", circle_path=remote_dir+"/circle.PNG",
                           dot_path=remote_dir+"/dot.PNG",main_id=self.option('main_id'))
        """
        self.end()

    def download_file(self):
        """
        download file from s3
        :return:
        """
        if self.analysis_type in ['complete']:
            assemble_dir = os.path.join(self.work_dir, "assemble_dir")
            for sample in self.samples.keys():
                sample_path = os.path.join(assemble_dir, self.samples[sample] + ".fna")
                assemble_dir3 = os.path.join(assemble_dir, self.samples[sample])
                dir_list = os.listdir(assemble_dir3)
                for file2 in dir_list:
                    os.system("cat {} >> {}".format(os.path.join(assemble_dir3, file2), sample_path))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        repaths = [
            [".", "dir", "共线性分析结果目录", 0],
            ['./dual_synteny.PNG', 'png', '共性性分析线图', 0],
            ['./circle.PNG', 'png', '共性性分析圈图', 0],
            ['./dot.PNG', 'png', '共性性分析点图', 0],
        ]
        result_dir.add_relpath_rules(repaths)
        super(SyntenyWorkflow, self).end()
