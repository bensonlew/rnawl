# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from mbio.packages.denovo_rna.assemble.trinity_stat import *
import os
import re
import shutil
import unittest


class BuscoAgent(Agent):
    """
    Trinity拼接
    author: 刘彬旭
    last_modify: 2017.11.9
    """

    def __init__(self, parent):
        super(BuscoAgent, self).__init__(parent)

        options = [
            {"name": "mode", "type": "string",
                "default": "tran"},  # 序列类型 geno tran prot
            # 物种类别目前支持 All(表示Eukaryota) Fungi Animal Plant Protist
            {"name": "species", "type": "string", "default": "All"},
            {"name": "out", "type": "string", "default": "busco_result"},  # 输出目录
            {"name": "cpu", "type": "int", "default": 1},  # busco软件所分配的cpu数量
            {"name": "max_memory", "type": "string",
                "default": '10G'},  # busco软件所分配的最大内存，单位为GB
            {"name": "fa", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "g2t", "type": "string", "default": ""},
            {"name": "busco_result", "type": "outfile",
                "format": "denovo_rna_v2.common"},
            {'name': 'summary_result', 'type': 'outfile',
                'format': 'denovo_rna_v2.common'},
            {'type': 'bool', 'name': 'marker', 'default': False},
            {'name': 'odb9', 'type': 'string'}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fa'):
            raise OptionError('必须设置序列文件 :fasta 格式', code="32002301")
        if self.option("g2t") and self.option('marker'):
            self.option('fa').check_g2t(self.option("g2t"))

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = self.option('cpu')
        self._memory = self.option('max_memory')

    def end(self):
        super(BuscoAgent, self).end()


class BuscoTool(Tool):
    def __init__(self, config):
        super(BuscoTool, self).__init__(config)
        self._version = "v1.0.1"
        software_dir = self.config.SOFTWARE_DIR
        self.python = os.path.join(software_dir, 'program/Python/bin')
        self.set_environ(PATH=self.python)
        self.busco_path = os.path.join(
            software_dir, 'bioinfo/denovo_rna_v2/busco-master-e83a6c94101511484799f9770cdfc148559b136d/')
        self.busco_pydb_path = self.config.SOFTWARE_DIR + \
            '/' + self.busco_path + 'build/lib'
        self.busco_db_path = self.config.SOFTWARE_DIR + '/database/orthodb/all_odb9/'
        self.db = {
            'All': 'eukaryota_odb9',
            'Fungi': 'fungi_odb9',
            'Animals': 'metazoa_odb9',
            'Plants': 'embryophyta_odb9',
            'Protists': 'protists_ensembl',
            'Bacteria': 'bacteria_odb9'
        }
        self.busco_config = ''
        if self.config.SOFTWARE_DIR == '/mnt/ilustre/users/sanger-dev/app':
            self.busco_config = self.config.SOFTWARE_DIR + \
                '/' + self.busco_path + 'config.ini.sanger-dev'
        elif self.config.SOFTWARE_DIR == '/mnt/ilustre/users/sanger-test/app':
            self.busco_config = self.config.SOFTWARE_DIR + \
                '/' + self.busco_path + 'config.ini.sanger-test'
        elif self.config.SOFTWARE_DIR == '/mnt/ilustre/users/isanger/app':
            self.busco_config = self.config.SOFTWARE_DIR + \
                '/' + self.busco_path + 'config.ini.isanger'
        elif self.config.SOFTWARE_DIR == '/mnt/lustre/users/sanger/app':
            self.busco_config = self.config.SOFTWARE_DIR + \
                '/' + self.busco_path + 'config.ini.nsanger'
        else:
            raise OptionError('暂不支持该机器，需重新设置busco配置文件', code="32002302")
        self.set_environ(BUSCO_CONFIG_FILE=self.busco_config)
        self.set_environ(PYTHONPATH=self.busco_pydb_path)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': os.path.join(software_dir, 'program/R-3.3.1/bin/Rscript'),
        }

    def run(self):
        """
        运行
        :return:
        """
        super(BuscoTool, self).run()
        self.run_busco()
        self.set_output()
        self.end()

    def trinity_fasta_marker(self, unigene):
        '''
        获取unigene序列，和相关对应关系列表
        '''
        self.option("fa").set_gene2tran(self.option("g2t"))
        self.option("fa").get_unigene(unigene, 'Trinity')
        return unigene

    def get_unigene(self, unigene):
        '''
        获取unigene序列，和相关对应关系列表
        '''
        # self.option("fasta").set_gene2tran(self.option("gene2trans"))
        self.option("fa").get_unigene_by_marker(unigene, self.option("g2t"))
        return unigene

    def run_busco(self):
        """
        运行trinity软件，进行拼接组装
        """
        seq = self.option('fa').prop['path']
        if self.option('g2t') != "":
            if self.option('marker'):
                seq = self.trinity_fasta_marker('unigene.fasta')
            else:
                seq = self.get_unigene("unigene.fasta")

        if os.path.exists(self.work_dir + '/run_' + self.option('out')):
            shutil.rmtree(self.work_dir + '/run_' + self.option('out'))

        # cmd = '{}scripts/run_BUSCO.py -i {} -l {} -o {} --mode {}'.format(
        #     self.busco_path, seq,
        #     self.busco_db_path + self.db[self.option('species')], self.option('out'), self.option('mode') )
        cmd = '{} {}scripts/run_BUSCO.py -i {} -l {} -o {} --mode {}'.format(self.program['python'],
                                                                             self.busco_path, seq,
                                                                             self.busco_db_path + self.option('odb9'), self.option('out'), self.option('mode'))
        self.logger.info('运行busco软件，进行组装结果评估')
        self.logger.info('busco config file is {}'.format(self.busco_config))

        command = self.add_command("busco_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("busco运行完成")
        else:
            self.set_error("busco运行出错", code="32002303")

    def set_output(self):
        busco_result = self.work_dir + '/run_busco_result/short_summary_busco_result.txt'
        summary_result = os.path.join(self.output_dir, 'summary_result')
        out = os.path.join(self.output_dir, os.path.basename(busco_result))
        if os.path.exists(out):
            os.remove(out)
        os.link(busco_result, out)
        with open(out, 'r') as s, open(summary_result, 'w') as r:
            r.write('name' + '\t' + 'short_name' + '\t' +
                    'persent' + '\t' + 'number' + '\n')
            for line in s.readlines():
                if line.startswith('#') or line.startswith('\n'):
                    continue
                else:
                    if line.strip().startswith('C:'):
                        C_persent = re.findall('C:(.*?)%', line.strip())[0]
                        S_persent = re.findall('S:(.*?)%', line.strip())[0]
                        D_persent = re.findall('D:(.*?)%', line.strip())[0]
                        F_persent = re.findall('F:(.*?)%', line.strip())[0]
                        M_persent = re.findall('M:(.*?)%', line.strip())[0]
                        total = re.findall('n:(.*?)$', line.strip())[0]
                        busco_dict = {'C': C_persent, 'S': S_persent,
                                      'D': D_persent, 'F': F_persent, 'M': M_persent, 'T': total}
                    elif line.strip().endswith('searched'):
                        continue
                    else:
                        num, type_busco = line.strip().split('\t')
                        short_name = re.findall('\((.*)\)', type_busco)[0]
                        r.write(type_busco + '\t' + short_name + '\t' +
                                busco_dict[short_name] + '\t' + num + '\n')
            r.write('Total BUSCO groups searched' + '\t' +
                    'T' + '\t' + '100' + '\t' + busco_dict['T'])
        try:
            self.option('busco_result', out)
            self.logger.info("设置busco结果目录成功")
        except Exception as e:
            self.logger.info("设置busco分析结果目录失败{}".format(e))
            self.set_error("设置busco分析结果目录失败%s", variables=(e), code="32002304")
        self.option('summary_result').set_path(summary_result)


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'busco{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.busco',
            'instant': False,
            'options': {
                'fa': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/Busco/assemble_raw.fasta',
                'odb9': 'eukaryota_odb9',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
