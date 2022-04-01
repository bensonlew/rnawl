# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class ArgoutAgent(Agent):
    def __init__(self, parent):
        super(ArgoutAgent, self).__init__(parent)
        options = [
            {'name': 'file_list', 'type': 'infile', 'format': 'arghub.file_list'},
            {'name': 'file_type', 'type': 'string', 'default': 'blast'},  # file_list中的结果类型 
            {'name': 'rename', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'output', 'type': 'outfile', 'format': 'sequence.profile_table'}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('file_list').is_set:
            raise OptionError('找不到输入文件')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'


class ArgoutTool(Tool):
    '''
    对arghub分析产生的结果文件结合数据注释文件进行格式化
    '''
    def __init__(self, config):
        super(ArgoutTool, self).__init__(config)
        self.db_table = self.get_db()
        self.out = self.option('file_type') + '_arg_results.txt'

    def run(self):
        self.logger.info('start outformat #####')
        super(ArgoutTool, self).run()
        try:
            self.run_format()
        except Exception as e:
            self.set_error('wrong in out format {}'.format(e))
        if self.option('rename').is_set:
            self.rename_id()
        self.end()

    def get_db(self):
        db_table_path = os.path.join(self.config.SOFTWARE_DIR, 'database/ArgHub/tables/arghub_gene.txt')
        db_table = pd.read_csv(db_table_path, header=0, sep='\t')
        annot_col = [
            'arghub_id', 'gene_name', 'db_type', 'db_source', 'anti_id', 'anti_type', 
            'anti_class', 'anti_name', 'resistance_mechanism', 'mutations', 'has_mutation', 'seq_type'
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

    def run_format(self):
        raw_table = []
        for f in self.option('file_list').file_list:
            raw_table.append(pd.read_csv(f, header=0, sep='\t'))
        raw_table = pd.concat(raw_table)
        raw_table.to_csv(os.path.join(self.output_dir, 'z.txt'), sep='\t', index=False)
        tb = pd.merge(raw_table, self.db_table, on='arghub_id')
        tb.to_csv(os.path.join(self.output_dir, 'z2.txt'), sep='\t', index=False)
        tb = self.get_variant_type(tb.rename(columns={'mutations': "variant_info"}))
        tb.to_csv(os.path.join(self.output_dir, self.out), sep='\t', index=False)
        self.option('output', os.path.join(self.output_dir, self.out))

    def rename_id(self):
        perl = '/program/perl-5.24.0/bin/perl'
        cmd = perl + " -i.bak -lane 'if ($ARGV[0]){{$i{{$F[0]}} = $F[1];" +\
                " print }} else {{ s/^(\S+)/$i{{$1}}/ if $i{{$F[0]}}; print }}' {} {}"
        cmd = cmd.format(self.option('rename').path, self.option('output').path)
        command = self.add_command('rename_id', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('结果基因id重命名成功')
        else:
            self.set_error('基因组id重命名失败')
