# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import subprocess
from biocluster.core.exceptions import OptionError
import types
import pandas as pd




class FunguildAgent(Agent):
    def __init__(self, parent):
        super(FunguildAgent, self).__init__(parent)
        options = [
            {"name": "taxon_table", "type": "infile","format":"meta.otu.otu_table"},
            {"name": "change_format","type": "string", "default":"True"},
            {"name": "out_detail", "type": "outfile","format":"meta.otu.otu_table"},
            {"name": "out_stat", "type": "outfile","format":"meta.otu.otu_table"},
            {"name": "others", "type" : "float", "default": 0 }

        ]
        self.add_option(options)
        self.step.add_steps('funguild')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.funguild.start()
        self.step.update()

    def step_end(self):
        self.step.funguild.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('taxon_table').is_set:
            raise OptionError('必须提供物种层次的表', code="32706901")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(FunguildAgent, self).end()


class FunguildTool(Tool):
    def __init__(self, config):
        super(FunguildTool, self).__init__(config)
        self._version = '1.0'
        self.python_path = "miniconda2/bin/python"
        self.script1_path = self.config.PACKAGE_DIR + "/meta/guilds_v1_1.py"
        self.taxon_table = self.option('taxon_table').path

        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def change_format_fun(self):
        data = pd.read_table(self.taxon_table,sep='\t',header=0)
        data['taxonomy'] = data['OTU ID'].apply(lambda x : ';'.join(str(x).split(';')[:-1]))
        data['OTU ID'] = data['OTU ID'].apply(lambda x : str(x).split(';')[-1].strip())
        data.to_csv('new_taxon_table.xls',sep='\t', index=False)
        self.taxon_table = 'new_taxon_table.xls'
        self.fungi_db_file = self.config.SOFTWARE_DIR + '/database/meta/funguild/funguild_fungi.json'

    def run_command(self):

        self.pre = self.taxon_table.split('.')[0]
        cmd = self.python_path +' %s -otu %s -db fungi -db_file %s ' % (self.script1_path, self.taxon_table,self.fungi_db_file )
        self.logger.info(cmd)

        # subprocess.check_output(self.config.SOFTWARE_DIR+'/'+cmd, shell=True)

        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Funguild成功")
        else:
            self.set_error("Funguild生成失败", code="32706901")

    def run_stat(self):
        self.guild_file = self.pre+ '.guilds.txt'
        data=pd.read_table(self.guild_file,sep='\t',header=0)
        head_list = data.columns.tolist()
        sample_list = head_list[1:head_list.index('taxonomy')]
        sub_data = data[sample_list+['Guild']]
        #sub_data['Guild'] = sub_data['Guild'].apply(lambda x: str(x).split('-')[0] if  str(x) != '-' else x)
        sub_data['Guild'] = sub_data['Guild'].apply(lambda x: 'unknown' if  (str(x) == '-' or str(x) == "nan") else x)
        data2= sub_data.groupby('Guild').sum()

        self.stat_out = 'FUNGuild_guild.txt'
        data2.to_csv(self.stat_out+'_ori',sep='\t')

        sub_data2 = data2[data2.index != 'unknown']
        unknown_data = data2[data2.index == 'unknown']
        sum_part = float(sub_data2.sum().sum())
        sub_data2['percent'] = sub_data2.apply(lambda  x: x.sum()/sum_part, axis=1 )

        sub_data2.reset_index(inplace=True)
        for index in sub_data2.index:
            # if sub_data2['percent'][index] < self.option('others'):
            #     sub_data2['Guild'][index] = 'others'
            if sub_data2.ix[index,'percent'] < self.option('others'):
                 sub_data2.ix[index, 'Guild'] = 'others'


        sub_data_sum = sub_data2.groupby('Guild').sum()
        sub_data_sum.to_csv('percent_no_unknown.xls' , sep='\t', index=True)

        sub_data_sum.drop(['percent'],inplace=True,axis=1)
        new_all = pd.concat([sub_data_sum,unknown_data])
        new_all.to_csv(self.stat_out , sep='\t', index=True)


    def run(self):
        super(FunguildTool, self).run()
        if self.option('change_format') == 'True':
            self.change_format_fun()
        self.run_command()
        self.run_stat()
        self.set_output()

    def set_output(self):
        detail_data = os.path.join(self.output_dir, 'Funguild.txt')
        linksum_data = os.path.join(self.output_dir, self.stat_out)
        for i in detail_data,linksum_data:
            if os.path.exists(i):
               os.remove(i)
        os.link(self.work_dir + '/'+ self.guild_file, detail_data)
        os.link(self.work_dir + '/'+self.stat_out, linksum_data)
        self.option('out_detail',detail_data)
        self.option('out_stat', linksum_data)
        self.end()



