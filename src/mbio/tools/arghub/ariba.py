# -*- coding: utf-8 -*- # __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class AribaAgent(Agent):
    def __init__(self, parent):
        super(AribaAgent, self).__init__(parent)
        options = [
            {'name': 'read1', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'read2', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'nthread', 'type': 'int', 'default': 8},
            {'name': 'db_type', 'type': 'string', 'default': 'core'},  # core or plus
            {'name': 'out_prefix', 'type': 'string', 'default': 'ariba_out'},
            {'name': 'output', 'type': 'outfile', 'format': 'sequence.profile_table'}
        ]
        self.add_option(options)

    def check_options(self):
        if not all([self.option('read1').is_set, self.option('read1').is_set]):
            raise OptionError('必须输入双端reads')

    def set_resource(self):
        self._memory = '16G'
        self._cpu = self.option('nthread')


class AribaTool(Tool):
    def __init__(self, config):
        super(AribaTool, self).__init__(config)
        self.ariba_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/arghub/ariba/'  # 待修改
        #self.ariba_path = '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/apps/ariba/'
        self.preref = self.config.SOFTWARE_DIR + '/database/ArgHub/ariba/' + self.option('db_type')
        self.read1 = self.option('read1').prop['path']
        self.read2 = self.option('read2').prop['path']
        self.set_environ(PATH='{0}/bin:{0}'.format(self.ariba_path))

    def run(self):
        super(AribaTool, self).run()
        self.run_ariba()
        self.set_output()
        self.end()

    def run_ariba(self):
        #cmd = 'export PATH={0}/bin:{0}:$PATH && ariba run --threads {1} {2} {3} {4} {5}'
        cmd = '{0}/bin/ariba run --force --verbose --threads {1} {2} {3} {4} {5}'
        cmd = cmd.format(self.ariba_path, self.option('nthread'), self.preref,
                         self.read1, self.read2, self.option('out_prefix'))
        command = self.add_command('run_ariba', cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('ariba runs done')
        else:
            self.set_error('wrong in running ariba')

    def get_db(self, arghub_ids):
        db_table_path = os.path.join(self.config.SOFTWARE_DIR, 'database/ArgHub/tables/arghub_gene.txt')
        db_table = pd.read_csv(db_table_path, header=0, sep='\t')
        db_table = db_table[db_table['arghub_id'].isin(arghub_ids)]
        annot_col = [
            'arghub_id', 'gene_name', 'db_type', 'seq_type', 'anti_type', 'mutations',
            'anti_name', 'resistance_mechanism', 'has_mutation', 'db_source', "anti_class"
        ]
        return db_table[annot_col]

    def get_variant_type(self, tb):
        v_p = (tb['has_mutation'] == 1) & (tb['db_source'] == 'core') & (tb['seq_type'] == "coding")
        vl_p = (tb['has_mutation'] == 1) & (tb['db_source'] != 'core') & (tb['seq_type'] == "coding")
        v_r = (tb['has_mutation'] == 1) & (tb['db_source'] == 'core') & (tb['seq_type'] != "coding")
        vl_r = (tb['has_mutation'] == 1) & (tb['db_source'] != 'core') & (tb['seq_type'] != "coding")
        tb["variant_type"] = "Homologou"
        tb["variant_type"][v_p] = "Varient (Protein)"
        tb["variant_type"][vl_p] = "Varient-like (Protein)"
        tb["variant_type"][v_r] = "Varient (rDNA)"
        tb["variant_type"][vl_r] = "Varient-like (rDNA)"
        return tb

    def set_output(self):
        self.logger.info('set output')
        result_file = self.work_dir + '/' + self.option('out_prefix') + '/' + 'report.tsv'
        with open(result_file, 'r') as r:
            r.readline()
            infos = {}
            infoss = {}
            mutations = {}
            for line in r:
                li = line.strip().split('\t')
                arg_id = ':'.join(li[1].split('_'))
                # reads_num reads_depth reads_cov has_mutation mutaion_type 
                infos[arg_id] = '{}\t{}\t{}\t{}\t{}\t{}'.format(li[5], li[12], float(li[8])/float(li[7]), li[8], li[7], li[11])
                infoss[arg_id] = [li[5], li[12], float(li[8])/float(li[7]), li[8], li[7], li[11]]
                if arg_id not in mutations:
                    mutations[arg_id] = []
                # mutation_type mutation_effect
                mutations[arg_id].append([li[18], li[19]])
        db_table = self.get_db(infos.keys())
        db_table = self.get_variant_type(db_table.rename(columns={"mutations": "variant_info"}))
        col = ["reads_num", "gene_depth", "reads_cov", "ref_len", "ref_c_len", "ctg_len"]
        for index, c in enumerate(col):
            db_table[c] = db_table["arghub_id"].agg(lambda x: infoss[x][index])
        db_table["mutation_info"] = db_table["arghub_id"].agg(lambda x: ';'.join([m[0] for m in mutations[x]]))
        db_table["mutation_effects"] = db_table["arghub_id"].agg(lambda x: ';'.join([m[1] for m in mutations[x]]))
        db_table.to_csv(self.output_dir + '/report.csv', sep='\t', index=False)
        self.option('output', self.output_dir + '/report.csv')

