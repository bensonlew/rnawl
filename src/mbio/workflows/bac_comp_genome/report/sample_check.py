# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modifies 20191014


import os,shutil
import json
import re
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
import gevent
from biocluster.file import download, exists
from mbio.packages.bac_comp_genome.common_function import get_fasta, get_sample_from_tree
from biocluster.core.exceptions import OptionError, FileError


class SampleCheckWorkflow(Workflow):
    """
    比较基因组文件检查的模块，主要是对gff文件处理，并导表，上传gff文件路径
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SampleCheckWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "type", "type": "string", "default": "up_load"},  # up_load、on_line
            {"name": "task_id", "type": "string"},  #任务id
            {"name": "mapping_file", "type": "string"}, #前端传过来的文件
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.file_path = self._sheet.output
        self.gff =[]
        self.modules =[]


    def check_file(self):
        gff2 = []
        gff = {}
        seq = {}
        sample_info = []
        download(self.option("mapping_file"), self.work_dir + "/mapping_file.txt")
        with open(self.work_dir + "/mapping_file.txt") as f:
            #tmp = f.readline().strip().encode('utf-8')
            #data = json.loads(tmp)
            data = json.load(f, encoding='utf-8')
            for key, value in data.items():
                for file in value:
                    name = str(file['alias'])
                    if name.endswith((".gff")):
                        gff2.append(file)
                        gff[name] = file
                    elif name.endswith((".fasta", ".fa", ".fna")):
                        seq[name] = file
                    elif name in ["list.txt"]:
                        sample_info.append(file)
        self.logger.info(gff)
        self.logger.info(seq)
        if len(sample_info) >1:
            raise OptionError('该文件夹下有多个list.txt文件！')
        elif len(sample_info) == 1:
            download(sample_info[0]['file_path'], self.work_dir + "/" + sample_info[0]['alias'])
        else:
            raise OptionError('该文件夹下没有list.txt文件！')
        with open (self.work_dir + "/" + sample_info[0]['alias'], "r") as f,open (self.work_dir + "/sample_info.xls", "w") as g:
            lines = f.readlines()
            g.write(lines[0].strip()+"\tseq_path\tgff_path\n")
            sampless = []
            for line in lines[1:]:
                lin = line.strip().split("\t")
                sampless.append(lin[0])
                gff_path = ''
                seq_path = ''
                self.logger.info(lin[6])
                if lin[3] in ["Complete"]:
                    if lin[4] == "-":
                        raise OptionError('该{}样品需在list.txt提供location！'.format(lin[0]))
                if lin[6] != "-":
                    if lin[3] in ['Draft']:
                        if lin[6] not in gff.keys():
                            raise OptionError('该文件夹下的缺少{}文件！'.format(lin[6]))
                        else:
                            gff_path = gff[lin[6]]['file_path']
                    elif lin[3] in ['Complete']:
                        if not re.search(";", lin[6]):
                            gff_path = gff[lin[6]]['file_path']
                        else:
                            list = []
                            for file in lin[6].split(";"):
                                seq1 = gff[file]['file_path']
                                list.append(seq1)
                            gff_path = ";".join(list)
                else:
                    gff_path = "-"
                if lin[3] in ['Draft']:
                    seq_path = seq[lin[5]]['file_path']
                elif lin[3] in ['Complete']:
                    if not re.search(";", lin[5]):
                        seq_path = seq[lin[5]]['file_path']
                    else:
                        list1 = []
                        for file in lin[5].split(";"):
                            seq1 = seq[file]['file_path']
                            list1.append(seq1)
                        seq_path = ";".join(list1)
                g.write("{}\t{}\t{}\n".format("\t".join(lin), seq_path, gff_path))
            if len(sampless) != len(set(sampless)):
                raise OptionError('该文件夹下样品名称重复！')
        return gff2


    def check_file2(self):
        gff2 = []
        download(self.option("mapping_file"), self.work_dir + "/mapping_file.txt")
        with open(self.work_dir + "/mapping_file.txt") as f,open (self.work_dir + "/sample_info.xls", "w") as g:
            data = json.load(f)
            g.write("sample\torganism\tstrain\tseq_status\tlocation\tfna\tgff\tfna_path\tgff_path\n")
            sampless =[]
            for line in data:
                gff = ''
                gff_path = ''
                seq = ''
                seq_path = ''
                location = ''
                if line['seq_status'] in ['Draft']:
                    seq = line['seq_alias']
                    seq_path = line['seq']
                    location = "-"
                    if line['giff'] != None:
                        gff = line['giff_alias']
                        gff_path = line['giff']
                        gff2.append({"alias": line['giff_alias'], "file_path": line['giff']})
                    else:
                        gff = "-"
                        gff_path = "-"
                elif line['seq_status'] in ['Complete']:
                    seq1 = []
                    seq_path = []
                    location = []
                    gff = []
                    gff_path = []
                    for file in line['seq_type']:
                        seq1.append(file['seq_alias'])
                        seq_path.append(file['seq'])
                        location.append(file['genome_type'])
                        if file['giff'] != None:
                            gff.append(file['giff_alias'])
                            gff_path.append(file['giff'])
                            gff2.append({"alias": file['giff_alias'], "file_path": file['giff']})
                    seq = ";".join(seq1)
                    seq_path = ";".join(seq_path)
                    location = ";".join(location)
                    if len(gff) >=1:
                        gff = ";".join(gff)
                        gff_path = ";".join(gff_path)
                    else:
                        gff = "-"
                        gff_path = "-"
                g.write("{}\n".format("\t".join([line['sample'],line['organism'],line['strain'],line['seq_status'],location,seq,gff,seq_path,gff_path])))
                sampless.append(line['sample'])
            if len(sampless) != len(set(sampless)):
                raise OptionError('该文件夹下样品名称重复！')
        return gff2

    def get_gff_seq(self):
        gff = []
        with open(self.work_dir + "/sample_info.xls","r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if lin[6] != "-":
                    if re.search(";", lin[6]):
                        for i in zip(lin[5].split(";"),lin[6].split(";"),lin[7].split(";"),lin[8].split(";")):
                            gff.append(i)
                    else:
                        gff.append((lin[5],lin[6],lin[7],lin[8]))
        return gff


    def run_gff(self):
        self.gff = self.get_gff_seq()
        self.logger.info(self.gff)
        for file in self.gff:
            self.logger.info(file)
            gff = self.add_tool("bac_comp_genome.gff_check")
            download(file[2], gff.work_dir + "/" + file[0])
            download(file[3], gff.work_dir + "/" + file[1])
            gff.set_options({
               "gff": gff.work_dir + "/" + file[1],
               "fa": gff.work_dir + "/" + file[0],
            })
            self.modules.append(gff)
        if len(self.modules) >1:
            self.on_rely(self.modules, self.run_gff_stat)
        elif len(self.modules) == 1:
            self.modules[0].on("end", self.run_gff_stat)
        for module in self.modules:
            module.run()

    def run_gff_stat(self):
        if os.path.exists(self.output_dir + "/gff"):
            shutil.rmtree(self.output_dir + "/gff")
        os.mkdir(self.output_dir + "/gff")
        with open(self.output_dir + "/all.gff_stat.xls", "w") as f:
            f.write("name\tcds_num\ttrna_num\trrna_num\gff_status\n")
            if len(self.modules) > 1:
                for module in self.modules:
                    for i in os.listdir(module.output_dir):
                        if re.search(".gff", i):
                            os.link(module.output_dir + "/" + i, self.output_dir + "/gff/" + i)
                        else:
                            with open(module.output_dir + "/" + i, "r") as g:
                                lines =g.readlines()
                                f.write(lines[1])
            elif len(self.modules) == 1:
                for i in os.listdir(self.modules[0].output_dir):
                    if re.search(".gff", i):
                        os.link(self.modules[0].output_dir + "/" + i, self.output_dir + "/gff/" + i)
                    else:
                        with open(self.modules[0].output_dir + "/" + i, "r") as g:
                            lines = g.readlines()
                            f.write(lines[1])
        self.set_db()

    def run(self):
        if self.option("type") in ['up_load']:
            self.gff2 = self.check_file()
            if len(self.gff2) >= 1:
                self.run_gff()
            else:
                gevent.spawn_later(5, self.set_db)
        elif self.option("type") in ['on_line']:
            self.gff2 = self.check_file2()
            if len(self.gff2) >= 1:
                self.run_gff()
            else:
                gevent.spawn_later(5, self.set_db)
        super(SampleCheckWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_path = self.api.api("bac_comp_genome.common_api")
        sample_id = self.option("main_id")
        path = self.file_path + "gff/"
        if len(self.gff2) >= 1:
            self.logger.info("sssssssssssss")
            self.logger.info(self.gff2)
            api_path.add_gff_detail(self.output_dir + "/all.gff_stat.xls", self.work_dir + "/sample_info.xls", sample_id, "sample_detail", path)
        else:
            api_path.add_sample_detail(self.work_dir + "/sample_info.xls", "sample_detail", sample_id, "sample,species,strain,seq_status,g_location,seq_file,gff_file,seq_path,gff_path",main_name="sample_id", insert_d={"cds_no":"-", "rrna_no":"-", "trna_no":"-"})
        self.end()

    def end(self):
        repaths = [
            [".", "", "",0],
        ]
        regexps = [
            [r'.*\.xls$', 'xls', '',0],
            [r'.*\.stat$', 'stat', '',0]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(SampleCheckWorkflow, self).end()
