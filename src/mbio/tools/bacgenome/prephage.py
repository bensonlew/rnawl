# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.4.3

import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class PrephageAgent(Agent):
    def __init__(self, parent):
        super(PrephageAgent, self).__init__(parent)
        options = [
            {"name": "prot_seq", "type": "infile", "format": "sequence.fasta"},  # 基因蛋白文件
            {"name": "scaf_seq", "type": "infile", "format": "sequence.fasta"},  # scaffold序列文件
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"}, # 基因的gff3文件
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},  # 注释总览表
            {'name': 'sample_name', "type": "string"},  # 样本名
            {"name": "analysis", "type": "string", "default": "uncomplete"},  ###流程分析模式complete，uncomplete
            {"name": "ana_type", "type": "string", "default": "bacgenome"},  ###流程分析模式bacgenome，bac_comparative
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("prot_seq").is_set:
            raise OptionError("必须添加基因的蛋白的序列文件！")
        if not self.option("scaf_seq").is_set:
            raise OptionError("必须添加基因组scaffold的序列文件")
        if not self.option("gene_gff").is_set:
            raise OptionError("必须添加基因组的gff3的文件！")

    def set_resource(self):
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        super(PrephageAgent, self).end()


class PrephageTool(Tool):
    def __init__(self, config):
        super(PrephageTool, self).__init__(config)
        self.prot_seq = self.option("prot_seq").prop['path']
        self.scaf_seq = self.option("scaf_seq").prop['path']
        self.gene_gff = self.option("gene_gff").prop['path']
        self.sample_name = self.option("sample_name")
        self.type = self.option("analysis")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.perl_script1 = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.phage_finder = "/bioinfo/Genomic/Sofware/phage_finder_v2.1/bin/Phage_Finder_v2.1.pl"
        self.lib = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/lib/tRNAscan-SE'
        self.set_environ(PERL5LIB=self.lib)
        self.aragorn = "/bioinfo/Genomic/Sofware/aragorn1.2.38/aragorn"
        self.blastp = "/bioinfo/Genomic/Sofware/ncbi-blast-2.2.28+/bin/blastp"
        self.blastdb = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/phage_finder_v2.1/DB/phage_10_02_07_release.db'
        self.id_phagename = self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/phage_finder_v2.1/DB/id_phagename.xls'
        self.tRNA = "/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/tRNAscan-SE"
        self.hmm_search = "/bioinfo/Genomic/Sofware/phage_finder_v2.1/bin/HMM3_searches.sh"
        self.trna_out = self.work_dir + "/" + 'tRNAscan.out'
        self.ncbi_out = self.work_dir + "/" + 'ncbi.out'
        self.tmrna = self.work_dir + "/" + 'tmRNA_aragorn.out'
        self.phage_info = self.work_dir + "/" + 'phage_finder_info.txt'
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/lib/tRNAscan-SE:"+self.config.SOFTWARE_DIR + "/program/perl-5.24.0/lib"
        self.set_environ(PERL5LIB=self.perl5path)

    def run_hmm_search(self):
        cmd = '{} {}'.format(self.hmm_search, self.prot_seq)
        self.logger.info(cmd)
        command = self.add_command("run_hmm_search", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_hmm_search运行完成")
        else:
            self.set_error("run_hmm_search运行出错!", code="31402701")

    def run_blastp(self):
        cmd = '{} -db {} -outfmt 6 -evalue 0.001 -query {} -out {} -task blastp -num_threads 5 -num_alignments 4'.format(self.blastp, self.blastdb, self.prot_seq, self.ncbi_out)
        self.logger.info(cmd)
        command = self.add_command("run_blastp", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_blastp运行完成")
        else:
            self.set_error("run_blastp运行出错!", code="31402702")

    def run_trna(self):
        cmd = '{} -B -L -o {} {}'.format(self.tRNA, self.trna_out, self.scaf_seq)
        self.logger.info(cmd)
        command = self.add_command("run_trna", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_trna运行完成")
        else:
            self.set_error("run_trna运行出错!", code="31402703")

    def run_aragorn(self):
        cmd = '{} -m -o {} {}'.format(self.aragorn, self.tmrna, self.scaf_seq)
        self.logger.info(cmd)
        command = self.add_command("run_aragorn", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_aragorn运行完成")
        else:
            self.set_error("run_aragorn运行出错!", code="31402704")

    def run_phage_to_info(self):
        if self.option("ana_type") in ['bacgenome']:
            cmd = '{} {}phage_to_info.pl {} {} {} {}'.format(self.perl_path, self.perl_script, self.type, self.gene_gff,
                                                         self.scaf_seq,
                                                         self.phage_info)
        else:
            cmd = '{} {}phage_to_info.pl {} {} {}'.format(self.perl_path, self.perl_script1, self.gene_gff,
                                                             self.scaf_seq,
                                                             self.phage_info)
        self.logger.info(cmd)
        command = self.add_command("run_phage_to_info", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_phage_to_info运行完成")
        else:
            self.set_error("run_phage_to_info运行出错!", code="31402705")

    def run_prephage(self):
        cmd = '{} -b {} -t {} -i {} -r {} -n {} -A {} -S'.format(self.phage_finder,self.work_dir ,self.ncbi_out, self.phage_info,
                                                           self.trna_out, self.tmrna, self.scaf_seq)
        self.logger.info(cmd)
        command = self.add_command("run_prephage", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_prephage运行完成")
        else:
            self.set_error("run_prephage运行出错!", code="31402706")

    def run_phage_stat(self):
        self.anno =self.option("anno").prop['path']
        des =self.work_dir + '/PFPR_tab.txt'
        if os.path.exists(des):
            cmd = '{} {}get-prephage_gene.pl {} {} {} {}'.format(self.perl_path, self.perl_script, self.type,
                                                                 des, self.anno,
                                                                 self.sample_name)
            cmd += ' {} '.format(self.id_phagename) # modify by ysh in 20190409 for add phage detail name
            self.logger.info(cmd)
            command = self.add_command("run_phage_stat", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_phage_stat运行完成")
            else:
                self.set_error("run_phage_stat运行出错!", code="31402707")

    def run_comp_prephage(self):
        des = self.work_dir + '/PFPR_tab.txt'
        if os.path.exists(des):
            cmd = '{} {}prephage_info.pl {} {} {} {}'.format(self.perl_path, self.perl_script1, self.type,
                                                                 des, self.gene_gff,
                                                                 self.sample_name)
            cmd += ' {} '.format(self.id_phagename)
            self.logger.info(cmd)
            command = self.add_command("run_comp_prephage", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run_comp_prephage运行完成!")
            else:
                self.set_error("run_comp_prephage运行出错!")

    def tiqu_fasta(self):
        if os.path.exists(self.work_dir + '/' + self.sample_name + '_prephage.fna'):
            os.remove(self.work_dir + '/' + self.sample_name + '_prephage.fna')
        cmd = '{} {}get_phage_fasta.pl {} {} {} {} {}'.format(self.perl_path, self.perl_script1, self.sample_name, "prephage",  self.scaf_seq, self.work_dir + '/' + self.sample_name + '.summary.xls', self.work_dir + '/' + self.sample_name + '_prephage.fna')
        self.logger.info(cmd)
        command = self.add_command("tiqu_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("tiqu_fasta运行完成!")
        else:
            self.set_error("tiqu_fasta运行出错!")

    def set_output(self):
        if self.option("ana_type") in ['bacgenome']:
            path = os.path.join(self.output_dir, 'prephage')
            if not os.path.exists(path):
                os.mkdir(path)
            files = os.listdir(path)
            for file in files:
                os.remove(path + file)
            num = self.get_size(self.work_dir + '/' + self.sample_name + '.summary.xls')
            self.logger.info(num)
            if num > 1:
                if os.path.exists(path + '/' + self.sample_name + '_prephage_summary.xls'):
                    os.remove(path + '/' + self.sample_name + '_prephage_summary.xls')
                os.link(self.work_dir + '/' + self.sample_name + '.summary.xls',
                        path + '/' + self.sample_name + '_prephage_summary.xls')
                if os.path.exists(path + '/' + self.sample_name + '_prephage_detail.xls'):
                    os.remove(path + '/' + self.sample_name + '_prephage_detail.xls')
                os.link(self.work_dir + '/' + self.sample_name + '.detail.xls',
                        path + '/' + self.sample_name + '_prephage_detail.xls')
                if os.path.exists(path + '/' + self.sample_name + '.stat.xls'):
                    os.remove(path + '/' + self.sample_name + '.stat.xls')
                os.link(self.work_dir + '/' + self.sample_name + '.stat.xls',
                        path + '/' + self.sample_name + '.stat.xls')
            if os.path.exists(self.work_dir + '/PFPR.con'):
                if os.path.exists(path + '/' + self.sample_name + '_prephage.fna'):
                    os.remove(path + '/' + self.sample_name + '_prephage.fna')
                os.link(self.work_dir + '/PFPR.con', path + '/' + self.sample_name + '_prephage.fna')
        else:
            num = self.get_size(self.work_dir + '/' + self.sample_name + '.summary.xls')
            if num > 1:
                if os.path.exists(self.output_dir + '/' + self.sample_name + '_prephage_summary.xls'):
                    os.remove(self.output_dir + '/' + self.sample_name + '_prephage_summary.xls')
                os.link(self.work_dir + '/' + self.sample_name + '.summary.xls',
                        self.output_dir + '/' + self.sample_name + '_prephage_summary.xls')
                if os.path.exists( self.work_dir + '/' + self.sample_name + '_prephage.fna'):
                    if os.path.exists(self.output_dir + '/' + self.sample_name + '_prephage.fna'):
                        os.remove(self.output_dir + '/' + self.sample_name + '_prephage.fna')
                    os.link(self.work_dir  + '/' + self.sample_name + '_prephage.fna', self.output_dir + '/' + self.sample_name + '_prephage.fna')

    def run(self):
        super(PrephageTool, self).run()
        if self.option("ana_type") in ['bacgenome']:
            self.run_hmm_search()
            self.run_blastp()
            self.run_trna()
            self.run_aragorn()
            self.run_phage_to_info()
            self.run_prephage()
            self.run_phage_stat()
            self.set_output()
            self.end()
        else:
            self.run_hmm_search()
            self.run_blastp()
            self.run_trna()
            self.run_aragorn()
            self.run_phage_to_info()
            self.run_prephage()
            self.run_comp_prephage()
            self.tiqu_fasta()
            self.set_output()
            self.end()


    def get_size(self,file):
        with open(file,'r') as f:
            lines = f.readlines()
            num =len(lines)
        return num