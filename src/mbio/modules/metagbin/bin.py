# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# last_modify: 2019.01.10

import os,shutil,re
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir

class BinModule(Module):
    """
     宏基因组Binning模块
    """
    def __init__(self, work_id):
        super(BinModule, self).__init__(work_id)
        options = [
            {"name": "minContig", "type": "string","default":"1000"},  #metabat2最小contigs
            {"name": "contig_fa", "type": "infile", "format": "sequence.fasta"},  # metabat2输入文件contigs.fa
            {'name': 'specimen_info', 'type': 'infile', 'format': 'meta_genomic.specimen_info'},
            {"name": "metabat_depth", "type": "infile", "format": "sequence.profile_table"},
            {"name": "maxbin_depth", "type": "infile", "format": "sequence.profile_table"},
            {"name": "bam_dir", "type": "infile", "format": "metagbin.bam_dir"},
            {"name": "sofware_bin", "type": "string", "default": "meatbat"},
            {"name": "bin_dir", "type": "outfile", "format": "sequence.fasta_dir"},
        ]
        self.add_option(options)
        self.metabat = self.add_tool('metagbin.metabat')
        self.maxbin = self.add_tool('metagbin.maxbin')
        self.concoct = self.add_module('metagbin.concoct')
        self.dastools = self.add_tool('metagbin.dastools')
        self.metabat_dasinput = self.add_tool('metagbin.dastool_input')
        self.maxbin_dasinput = self.add_tool('metagbin.dastool_input')
        self.concoct_dasinput = self.add_tool('metagbin.dastool_input')
        self.step.add_steps('metabat', 'maxbin', 'concoct', 'dastools')
        self.tools =[]

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("contig_fa").is_set:
            raise OptionError("必须提供contigs序列文件!")
        if self.option("sofware_bin") == "":
            raise OptionError("必须提供bin的软件！")
        if self.option('sofware_bin') != "":
            if re.search(',',self.option('sofware_bin')):
                sof_list = self.option('sofware_bin').split(",")
                for sof in sof_list:
                    self.sof =sof
                    if self.sof in ["metabat"]:
                        self.tools.append(self.metabat_dasinput)
                    if self.sof in ["maxbin"]:
                        self.tools.append(self.maxbin_dasinput)
                    if self.sof in ["concoct"]:
                        self.tools.append(self.concoct_dasinput)
        return True

    def run(self):
        super(BinModule, self).run()
        if self.option('sofware_bin') != "":
            if re.search(',',self.option('sofware_bin')):
                self.on_rely(self.tools, self.run_das_tools)
                sof_list = self.option('sofware_bin').split(",")
                for sof in sof_list:
                    self.sof =sof
                    if self.sof in ["metabat"]:
                        self.run_metabat()
                    if self.sof in ["maxbin"]:
                        self.run_maxbin()
                    if self.sof in ["concoct"]:
                        self.run_concoct()
            else:
                sof = self.option('sofware_bin')
                if sof in ["metabat"]:
                    self.run_metabat()
                if sof in ["maxbin"]:
                    self.run_maxbin()
                if sof in ["concoct"]:
                    self.run_concoct()

    def run_metabat(self):
        if int(self.option("minContig")) >= 1500:
            mincontig =str(self.option("minContig"))
        else:
            mincontig = str(1500)
        opts = {
            "minContig": mincontig,
            "contig_fa": self.option("contig_fa"),
            "depth_file": self.option("metabat_depth"),
        }
        self.metabat.set_options(opts)
        if re.search(',', self.option('sofware_bin')):
            self.metabat.on("end",self.run_metabat_dasinput)
            self.metabat.run()
        else:
            self.metabat.on('end', self.set_output, "metabat")
            self.metabat.run()

    def run_maxbin(self):
        opts = {
            "minContig":self.option("minContig"),
            "contig_fa": self.option("contig_fa"),
            "depth_file": self.option("maxbin_depth"),
        }
        self.maxbin.set_options(opts)
        if re.search(',', self.option('sofware_bin')):
            self.maxbin.on("end", self.run_maxbin_dasinput)
            self.maxbin.run()
        else:
            self.maxbin.on('end', self.set_output, "maxbin")
            self.maxbin.run()

    def run_concoct(self):
        opts = {
            "minContig":self.option("minContig"),
            "contig_fa": self.option("contig_fa"),
            "bam_dir": self.option("bam_dir"),
            'specimen_info': self.option("specimen_info")
        }
        self.concoct.set_options(opts)
        if re.search(',', self.option('sofware_bin')):
            self.concoct.on("end", self.run_concoct_dasinput)
            self.concoct.run()
        else:
            self.concoct.on('end', self.set_output, "concoct")
            self.concoct.run()

    def run_metabat_dasinput(self):
        opts = {
            "sof": 'metabat',
            'bin_dir':self.metabat.option("metabat_bin"),
        }
        self.metabat_dasinput.set_options(opts)
        self.metabat_dasinput.run()

    def run_maxbin_dasinput(self):
        opts = {
            "sof": 'maxbin',
            'bin_dir': self.maxbin.option("maxbin_bin"),
        }
        self.maxbin_dasinput.set_options(opts)
        self.maxbin_dasinput.run()

    def run_concoct_dasinput(self):
        opts = {
            "sof": 'concoct',
            'bin_dir': self.concoct.option("concoct_bin"),
        }
        self.concoct_dasinput.set_options(opts)
        self.concoct_dasinput.run()

    def run_das_tools(self):
        opts= {
            "sofware_bin":self.option('sofware_bin'),
            "contig_fa":self.option('contig_fa'),
        }
        sof_list = str(self.option('sofware_bin')).split(",")
        if len(sof_list) >=2:
            for sof in sof_list:
                if sof in ["metabat"]:
                    opts['metabat_in']=self.metabat_dasinput.option('out')
                elif sof in ["maxbin"]:
                    opts['maxbin_in']=self.maxbin_dasinput.option('out')
                elif sof in ["concoct"]:
                    opts['concoct_in']=self.concoct_dasinput.option('out')
        self.dastools.set_options(opts)
        self.dastools.on('end', self.set_output, "dastools")
        self.dastools.run()

    def set_output(self, event):
        obj = event['bind_object']
        if not re.search(',', self.option('sofware_bin')):
            if event['data'] == 'metabat':
                if os.path.exists(obj.output_dir + '/metabat_bin'):
                    link_dir(obj.output_dir + '/metabat_bin', self.output_dir + '/bin')
                self.option("bin_dir",self.output_dir + '/bin')
                self.end()
            if event['data'] == 'maxbin':
                if os.path.exists(obj.output_dir + '/maxbin_bin'):
                    link_dir(obj.output_dir + '/maxbin_bin', self.output_dir + '/bin')
                self.option("bin_dir", self.output_dir + '/bin')
                self.end()
            if event['data'] == 'concoct':
                if os.path.exists(obj.output_dir + '/concoct_bin'):
                    link_dir(obj.output_dir + '/concoct_bin', self.output_dir + '/bin')
                self.option("bin_dir", self.output_dir + '/bin')
                self.end()
        else:
            if event['data'] == 'dastools':
                if os.path.exists(obj.output_dir + '/bin'):
                    link_dir(obj.output_dir + '/bin', self.output_dir + '/bin')
                self.option("bin_dir", self.output_dir + '/bin')
                self.end()

    def end(self):
        super(BinModule, self).end()
