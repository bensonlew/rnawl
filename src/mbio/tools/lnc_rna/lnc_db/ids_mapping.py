#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/15 10:14
@file    : ids_mapping.py
"""

import csv
import datetime
import json
import os
import time
import unittest
from itertools import chain
import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory
from mbio.files.lnc_rna.lnc_fasta import LncFastaFile


class IdsMappingAgent(Agent):
    def __init__(self, parent):
        """参数说明：
            params 格式: json 字符串
            [
                {
                    'db_name': '该fa文件所属数据库名',
                    'fasta': 'fasta path',
                    'ids_mapping': '自己transcript id 和 参考id 间映射, 文件格式: self_id<TAB>ref_id',
                    'lnc_ids': 'lncrna id list file',
                },
                {
                    ...
                }
            ]

        :param parent:
        """
        super(IdsMappingAgent, self).__init__(parent)
        options = [
            {'name': 'params', 'type': 'string', 'required': True},
            {'name': 'ref_lnc_list', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'ref_db_name', 'type': 'string', 'required': True},
            {'name': 'out_fa', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'out_ids_map', 'type': 'outfile', 'format': 'lnc_rna.common'}
        ]
        self.__num = 5
        self.__memory = '10G'
        self.add_option(options)
        self.step.add_steps("ids_mapping")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.ids_mapping.start()
        self.step.update()

    def step_end(self):
        self.step.ids_mapping.finish()
        self.step.update()

    def check_options(self):
        try:
            self.__num = len(json.loads(self.option('params'))) + 2
            self.__memory = '%sG' % ((self.__num - 1) * 4)
        except (json.decoder.JSONDecodeError, ValueError):
            raise OptionError('params 参数必须为 json 字符串对象')
        return True

    def set_resource(self):
        self.logger.debug('cpu: %s, memory: %s' % (self.__num, self.__memory))
        self._cpu = self.__num
        self._memory = self.__memory

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(IdsMappingAgent, self).end()


class IdsMappingTool(Tool):
    def __init__(self, config):
        super(IdsMappingTool, self).__init__(config)
        self.fasta_formatter = '/bioinfo/seq/fastx_toolkit-0.0.13.2/bin/fasta_formatter'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/fastx_toolkit-0.0.13.2/lib')
        self.__filter_func = None
        self.__convert_chrom = lambda key, line,: line.lstrip('chr')
        self.__out_dict = {}

    def params_parser(self):
        """
            fasta, lnc_ids, outdir, ids_mapping, ref_db, db_name
        :return:
        """
        params_map = {'fasta': '-f', 'lnc_ids': '-l', 'ids_mapping': '-i', 'db_name': '-d', 'ref_lnc_list': '-n'}
        params = json.loads(self.option('params'))
        for dic in params:
            one_run_params = []
            self.logger.debug('params: %s' % json.dumps(dic))
            db_name = dic['db_name']
            ref_name = self.option('ref_db_name')
            dir_name = ''.join(i.capitalize() for i in db_name.replace('.', '_').split('_')) + 'Dir'
            outdir = os.path.join(self.work_dir, dir_name)
            for k, v in dic.items():
                if not v:
                    continue
                one_run_params.extend((params_map[k], v))
            self.__out_dict[db_name] = {'outdir': outdir, 'map_file': db_name + '2' + ref_name + '.txt',
                                           'fa_name': os.path.basename(dic['fasta'])}
            one_run_params.extend(('-r', ref_name))
            one_run_params.extend(('-o', outdir))
            one_run_params.extend(('-n', self.option('ref_lnc_list').path))
            yield db_name, ' '.join(one_run_params)

    def run_cmd(self, cmd_name, cmd, is_wait=True, shell=False):
        self.logger.debug(cmd_name + ': ' + cmd + '%s' % type(cmd_name))
        cmd_obj = self.add_command(str(cmd_name), cmd, shell=shell)
        if shell is True:
            cmd_obj.software_dir = ''
            cmd_obj._start_run_time = int(time.time())
        cmd_obj.run()
        if is_wait is True:
            self._check_stat(cmd_obj)
            return
        return cmd_obj

    def _check_stat(self, *cmd_objs):
        self.wait(*cmd_objs)
        for cmd_obj in cmd_objs:
            if cmd_obj.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_obj.name))
            elif cmd_obj.return_code is None:
                self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_obj.name))
                cmd_obj.rerun()
                self.wait()
                if cmd_obj.return_code is 0:
                    self.logger.info("{} Finished successfully".format(cmd_obj.name))
                else:
                    self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_obj.name, cmd_obj.cmd))

    def merge_result(self):
        files = [os.path.join(dic['outdir'], dic['fa_name']) for _, dic in self.__out_dict.items()]
        date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        outfile = os.path.join(self.output_dir, '%s_merged.fa' % date)
        cmd = "/bin/cat {files} | grep -v '^\s*$' > {outfile}".format(files=' '.join(files), outfile=outfile)
        cmd_obj = self.run_cmd('merge_fas', cmd, is_wait=False, shell=True)

        map_files = [os.path.join(dic['outdir'], dic['map_file']) for _, dic in self.__out_dict.items()]
        for f in map_files[:]:
            if not os.path.isfile(f):
                self.logger.debug('\n' + f + ' is not existent\n')
                del map_files[map_files.index(f)]
        dfs = [pd.read_table(f, sep='\t', header=0, index_col='transcript_id') for f in map_files]
        ids_map_file = os.path.join(self.output_dir, '%s_merged_ids_maping.xls' % date)
        if len(dfs) != 0:
            df = pd.concat(dfs, axis=1, join='outer') if len(dfs) > 1 else dfs[0] # 存在bug: axis 错误
            df.fillna(value='', inplace=True)
            df.index.name = 'transcript_id'
            df.to_csv(ids_map_file, sep='\t', header=True, index=True)
        else:
            with open(ids_map_file, 'w') as out_handler:
                pass
            self.logger.debug('no ids_map file out')
        self.wait(cmd_obj)

        ## 转换fasta格式，统一序列长度
        out_fa = os.path.join(self.output_dir, '%s_merged_format.fa' % date)
        cmd = '{} -i {} -w 100 -o {}'.format(self.fasta_formatter, outfile, out_fa)
        command = self.add_command("fasta_formatter", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd))
        else:
            self.set_error("运行%s出错".format(cmd))

        # 设置输出路径
        self.option('out_fa').set_path(out_fa)
        self.option('out_ids_map').set_path(ids_map_file)

    def run_filter(self):
        cmd_objs = []
        tool_path = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/lnc_db/ids_mapping.py')
        python = 'program/Python/bin/python'
        self.logger.debug('fa filter START =====================')
        for db_name, str_cmd_params in self.params_parser():
            cmd = '{python} {tool} {params}'.format(python=python, tool=tool_path, params=str_cmd_params)
            self.logger.debug('fa filter, '+ cmd)
            cmd_obj = self.run_cmd(db_name + '_filter', cmd, is_wait=False)
            cmd_objs.append(cmd_obj)
        self.wait(*cmd_objs)

    def run(self):
        super(IdsMappingTool, self).run()
        self.run_filter()
        self.merge_result()
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            params = [
                {
                    'db_name': 'ncbi',
                    'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/GCF_000001405.38_GRCh38.p12_rna.fna',
                    'ids_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/ncbi_transcript2ensembl',
                    'lnc_ids': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/filter_dir/GCF_000001405.38_GRCh38.p12_genomic.lnc_ids.list',
                    'ref_lnc_list': ''
                },
                # {
                #     'db_name': 'ncbi_t',
                #     'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/GCF_000001405.38_GRCh38.p12_rna.fna',
                #     'ids_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/ncbi_transcript2ensembl',
                #     'lnc_ids': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/filter_dir/GCF_000001405.38_GRCh38.p12_genomic.lnc_ids.list'
                # },
                # {
                #     'db_name': 'ensembl',
                #     'fasta': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/ensembl/Homo_sapiens.GRCh38.ncrna.fa',
                #     'ids_mapping': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/ensembl/Homo_sapiens.GRCh38.ids_mapping.txt',
                #     'lnc_ids': '/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/2019_02_27/data/filter_dir/Homo_sapiens.GRCh38.95.lnc_ids.list'
                # }
            ]
            data = {
                "id": "IdsMappingTool_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lnc_db.ids_mapping",
                "instant": False,
                "options": dict(
                    ref_db_name='ensembl',
                    params=json.dumps(params)
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
