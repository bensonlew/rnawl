# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan, qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import unittest
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig

class KeggAnnotAgent(Agent):
    '''
    last_modify: 2019.04.16
    '''
    def __init__(self, parent):
        super(KeggAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast_xml', 'type': 'infile', 'format': 'ref_rna_v2.blast_xml'},
            {'name': 'longest_t2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'known_ko', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            # ['All', 'Animals', 'Plants', 'Protists']
            {'name': 'taxonomy', 'type': 'string', 'default': None},
            # 通路图链接官网颜色，约定参考基因组为黄色（yellow），新序列为绿色（green）, 两者共有为红色（tomato）
            {'name': 'link_bgcolor', 'type': 'string', 'default': None},
            # 通路图静态图颜色，约定参考基因组为黄色（FFFF00），新序列为绿色（00CD00）
            {'name': 'png_bgcolor', 'type': 'string', 'default': None},
            {'name': 'database', 'type': 'string', 'default': 'kegg'},
            {'name': 'txpt_kegg_table', 'type': 'outfile', 'format': 'ref_rna_v2.kegg_table'},
            {'name': 'gene_kegg_table', 'type': 'outfile', 'format': 'ref_rna_v2.kegg_table'},
            {'name': 'txpt_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {"name": "merge_type", "type": "string", "default": "partial"},
            {"name": "kegg_subtax1", "type": "string", "default": None},
            {"name": "kegg_subtax2", "type": "string", "default": None},
            {"name": "kegg_species", "type": "string", "default": None},
        ]
        self.add_option(options)
        self._memory_increase_step = 100
        self.step.add_steps('kegg_annot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.kegg_annot.start()
        self.step.update()

    def step_end(self):
        self.step.kegg_annot.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        infile_size = os.path.getsize(self.option('blast_xml').prop['path'])
        self._memory = '{}G'.format(int(float(infile_size) / 1024 ** 3 * 8 + 100))

    def end(self):
        super(KeggAnnotAgent, self).end()

class KeggAnnotTool(Tool):
    def __init__(self, config):
        super(KeggAnnotTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=os.path.join(self.config.SOFTWARE_DIR, 'library/lib/lib'))
        self.python = 'program/Python/bin/python'
        self.kegg_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/kegg_annotation_v2.py')
        self.txml2gxml_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/txml2gxml.py')
        self.kegg_known_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/kegg_known.py')
        self.rscript = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.3/bin/Rscript')
        self.map4_r = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/map4.r')
        self.get_venn_list_py = os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/annotation/get_venn_list.py')
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.kegg_map_html = self.kegg_files_dict['html'] + '/'
        self.ko_txt = self.kegg_files_dict['species'] + '/{}.ko.txt'.format(self.option('taxonomy'))

        '''
        if self.option('kegg_version') in ["201909", "202003", "202109"]:

            self.ko_txt = os.path.join(
                self.config.SOFTWARE_DIR, 'database/Annotation/other{}/kegg{}/species/{}.ko.txt'.format(self.option('kegg_version')[0:4],
                                                                                                        self.option('kegg_version'),
                                                                                                        self.option('taxonomy'))
            )
            self.kegg_map_html = os.path.join(self.config.SOFTWARE_DIR, 'database/Annotation/other{}/kegg{}/html'.format(self.option('kegg_version')[0:4],
                                                                                                                         self.option('kegg_version'),))
        else:
            self.ko_txt = os.path.join(
                self.config.SOFTWARE_DIR, 'database/KEGG/species/{}.ko.txt'.format(self.option('taxonomy'))
            )
            self.kegg_map_html = os.path.join(self.config.SOFTWARE_DIR, 'database/KEGG/map_html')
        '''
        self.image_magick_convert = os.path.join(self.config.SOFTWARE_DIR, 'program/ImageMagick/bin/convert')
        if self.option('known_ko').is_set:
            self.txpt_known_ko = os.path.join(self.work_dir, 'known_ko.T.tsv')
            self.gene_known_ko = os.path.join(self.work_dir, 'known_ko.G.tsv')
        else:
            self.txpt_known_ko = 'None'
            self.gene_known_ko = 'None'
        self.prefix = os.path.basename(self.option('blast_xml').path)[:-4]
        self.txpt_xml = os.path.join(self.work_dir, '{}.T.xml'.format(self.prefix))
        self.gene_xml = os.path.join(self.work_dir, '{}.G.xml'.format(self.prefix))
        self.txpt_dir = os.path.join(self.output_dir, 'T')
        self.gene_dir = os.path.join(self.output_dir, 'G')
        self.txpt_kegg_table_tsv = os.path.join(self.txpt_dir, 'kegg_table.tsv')
        self.gene_kegg_table_tsv = os.path.join(self.gene_dir, 'kegg_table.tsv')
        self.txpt_pid_txt = os.path.join(self.work_dir, 'pid.T.txt')
        self.gene_pid_txt = os.path.join(self.work_dir, 'pid.G.txt')
        self.txpt_pathways = os.path.join(self.txpt_dir, 'pathways')
        self.gene_pathways = os.path.join(self.gene_dir, 'pathways')
        self.txpt_pathway_table_tsv = os.path.join(self.txpt_dir, 'pathway_table.tsv')
        self.gene_pathway_table_tsv = os.path.join(self.gene_dir, 'pathway_table.tsv')
        self.txpt_kegg_layer_tsv = os.path.join(self.txpt_dir, 'kegg_layer.tsv')
        self.gene_kegg_layer_tsv = os.path.join(self.gene_dir, 'kegg_layer.tsv')
        self.txpt_list = os.path.join(self.output_dir, '{}.T.list'.format(self.option('database')))
        self.gene_list = os.path.join(self.output_dir, '{}.G.list'.format(self.option('database')))

    def get_limit_ko(self):
        self.client = Config().get_mongo_client(mtype='ref_rna', ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname('ref_rna', ref=True)]
        self.kegg_species = self.mongodb.kegg_species_map

        result = None
        try:
            if self.option("kegg_species"):
                result = self.kegg_species.find_one({'abr': self.option("kegg_species"), "level": "species"})
        except:
            self.logger.info("未找到该分类列表 {}".format(self.option("kegg_species")))

        if not result:
            try:
                if self.option("kegg_subtax2"):
                    result = self.kegg_species.find_one({'abr': self.option("kegg_subtax2"), "level": "classIII"})
            except:
                self.logger.info("未找到该分类列表 {}".format(self.option("kegg_subtax2")))

        if not result:
            try:
                if self.option("kegg_subtax1"):
                    result = self.kegg_species.find_one({'abr': self.option("kegg_subtax1"), "level": "classII"})
            except:
                self.logger.info("未找到该分类列表 {}".format(self.option("kegg_subtax1")))

        if result:
            ko_list = ["ko" + mapid for mapid in result['map_list']]
            with open(self.work_dir + "/ko_limit.txt", 'w') as fo:
                fo.write("\n".join(ko_list))
            self.ko_txt = self.work_dir + "/ko_limit.txt"
            return self.work_dir + "/ko_limit.txt"
        else:
            self.logger.info("未找物种分类列表 ，或数据库连接失败")
            return None



    def run(self):
        super(KeggAnnotTool, self).run()
        if self.option("kegg_subtax1") or self.option("kegg_subtax2"):
            self.get_limit_ko()
        if self.option('known_ko').is_set:
            self.run_kegg_known()
        self.run_txml2gxml()
        self.run_kegg_annotation_t()
        self.run_kegg_annotation_g()
        self.run_get_venn_list_t()
        self.run_get_venn_list_g()
        self.set_output()
        self.end()

    def run_kegg_known(self):
        shutil.copy(self.option('known_ko').path, self.txpt_known_ko)
        cmd = '{} {}'.format(self.python, self.kegg_known_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --ko {}'.format(self.txpt_known_ko)
        cmd += ' --output {}'.format(self.gene_known_ko)
        cmd_name = 'run_kegg_known'
        self.run_code(cmd_name, cmd, block=False)

    def run_txml2gxml(self):
        shutil.copy(self.option('blast_xml').path, self.txpt_xml)
        cmd = '{} {}'.format(self.python, self.txml2gxml_py)
        cmd += ' --t2g {}'.format(self.option('longest_t2g').path)
        cmd += ' --xml {}'.format(self.txpt_xml)
        cmd += ' --output {}'.format(self.gene_xml)
        cmd_name = 'run_txml2gxml'
        self.run_code(cmd_name, cmd)

    def run_kegg_annotation_t(self):
        if os.path.isdir(self.txpt_dir):
            shutil.rmtree(self.txpt_dir)
        os.mkdir(self.txpt_dir)
        if self.option("merge_type") == "refonly":
            self.txpt_xml = "None"
            self.gene_xml = "None"
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python, # 0
            self.kegg_annotation_py, # 1
            self.rscript, # 2
            self.map4_r, # 3
            self.txpt_xml, # 4
            None, # 5
            self.txpt_kegg_table_tsv, # 6
            self.txpt_pid_txt, # 7
            self.txpt_pathways, # 8
            self.txpt_pathway_table_tsv, # 9
            self.txpt_kegg_layer_tsv, # 10
            self.ko_txt, # 11
            self.option('link_bgcolor'), # 12
            self.option('png_bgcolor'), # 13
            self.image_magick_convert, # 14
            self.kegg_map_html, # 15
            self.txpt_known_ko, # 16
            self.option('kegg_version'), #17
        )
        cmd_name = 'run_kegg_annotation_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_kegg_annotation_g(self):
        if os.path.isdir(self.gene_dir):
            shutil.rmtree(self.gene_dir)
        os.mkdir(self.gene_dir)
        if self.option("merge_type") == "refonly":
            self.txpt_xml = "None"
            self.gene_xml = "None"
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            self.python, # 0
            self.kegg_annotation_py, # 1
            self.rscript, # 2
            self.map4_r, # 3
            self.gene_xml, # 4
            None, # 5
            self.gene_kegg_table_tsv, # 6
            self.gene_pid_txt, # 7
            self.gene_pathways, # 8
            self.gene_pathway_table_tsv, # 9
            self.gene_kegg_layer_tsv, # 10
            self.ko_txt, # 11
            self.option('link_bgcolor'), # 12
            self.option('png_bgcolor'), # 13
            self.image_magick_convert, # 14
            self.kegg_map_html, # 15
            self.gene_known_ko, # 16
            self.option('kegg_version'), #17
        )
        cmd_name = 'run_kegg_annotation_g'
        self.run_code(cmd_name, cmd)

    def run_get_venn_list_t(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.txpt_kegg_table_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.txpt_list)
        cmd_name = 'run_get_venn_list_t'
        self.run_code(cmd_name, cmd, block=False)

    def run_get_venn_list_g(self):
        cmd = '{} {}'.format(self.python, self.get_venn_list_py)
        cmd += ' --input {}'.format(self.gene_kegg_table_tsv)
        cmd += ' --database {}'.format(self.option('database'))
        cmd += ' --output {}'.format(self.gene_list)
        cmd_name = 'run_get_venn_list_g'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run %s, abord', variables=(n), code="33711402")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('txpt_kegg_table').set_path(self.txpt_kegg_table_tsv)
        self.option('gene_kegg_table').set_path(self.gene_kegg_table_tsv)
        self.option('txpt_list').set_path(self.txpt_list)
        self.option('gene_list').set_path(self.gene_list)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_ref(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'kegg_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.annotation.kegg_annot',
            'instant': False,
            'options': {
                'blast_xml': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotFilter/output/kegg/blast.filter.xml',
                'longest_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20190413/LncRna_tsg_33844/AnnotClass/AnnotFile/result/longest.t2g.tsv',
                'taxonomy': 'Animals',
                'link_bgcolor': 'yellow',
                'png_bgcolor': 'FFFF00',
                'known_ko': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Test_ref/Ensemble_release_89/KEGG/Oreochromis_niloticus.Orenil1.0.pathway',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_ref')])
    unittest.TextTestRunner(verbosity=2).run(suite)
