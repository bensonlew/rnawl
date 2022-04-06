# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import shutil
import unittest
from mbio.packages.rna.annot_config import AnnotConfig

class KeggAnnotAgent(Agent):
    '''
    last_modify: 2019.02.22
    '''
    def __init__(self, parent):
        super(KeggAnnotAgent, self).__init__(parent)
        options = [
            {'name': 'blast_xml', 'type': 'infile', 'format': 'lnc_rna.blast_xml'},
            {'name': 'known_ko', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'taxonomy', 'type': 'string', 'default': ''},
            # 通路图链接官网颜色，约定参考基因组为黄色（yellow），新序列为绿色(green), 两者共有为红色（tomato）
            {'name': 'link_bgcolor', 'type': 'string', 'default': ''},
            # 通路图静态图颜色，约定为黄色（#FFFF00）和绿色（#00CD00）
            {'name': 'png_bgcolor', 'type': 'string', 'default': ''},
            {"name": "kegg_version", "type": "string", "default": "202003"},
            {'name': 'kegg_table', 'type': 'outfile', 'format': 'lnc_rna.kegg_table'},
            {"name": "kegg_subtax1", "type": "string", "default": None},
            {"name": "kegg_subtax2", "type": "string", "default": None},
            {"name": "kegg_species", "type": "string", "default": None},
        ]
        self.add_option(options)
        self._memory_increase_step = 20
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
        if self.option('blast_xml').is_set:
            self.logger.debug('{} = {}'.format('blast_xml', self.option('blast_xml').prop['path']))
            self.infile_size = os.path.getsize(self.option('blast_xml').prop['path'])
        else:
            raise OptionError('input BLAST XML must be provided')
        self.logger.debug('{} = {}'.format('taxonomy', self.option('taxonomy')))
        if self.option('taxonomy') == '':
            raise OptionError('taxonomy must be in Animals, Plants, Fungi, Protists, Archaea, Bacteria or None')
        self.logger.debug('{} = {}'.format('link_bgcolor', self.option('link_bgcolor')))
        self.logger.debug('{} = {}'.format('png_bgcolor', self.option('png_bgcolor')))
        if self.option('known_ko').is_set:
            self.logger.debug('{} = {}'.format('known_ko', self.option('known_ko').prop['path']))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 8 + 40))

    def end(self):
        super(KeggAnnotAgent, self).end()

class KeggAnnotTool(Tool):
    def __init__(self, config):
        super(KeggAnnotTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=os.path.join(self.config.SOFTWARE_DIR, 'library/lib/lib'))
        self.python = 'miniconda2/bin/python'
        self.kegg_annotation_py = os.path.join(self.config.PACKAGE_DIR, 'rna/annotation/kegg_annotation_v2.py')
        self.rscript = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.3/bin/Rscript')
        self.map4_r = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/map4.r')
        self.kegg_table = os.path.join(self.work_dir, 'kegg_table.xls')
        self.pid_txt = os.path.join(self.work_dir, 'pid.txt')
        self.pathways = os.path.join(self.work_dir, 'pathways')
        self.pathway_table_xls = os.path.join(self.work_dir, 'pathway_table.xls')
        self.kegg_layer_xls = os.path.join(self.work_dir, 'kegg_layer.xls')
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.kegg_map_html = self.kegg_files_dict['html'] + '/'
        self.ko_txt = self.kegg_files_dict['species'] + '/{}.ko.txt'.format(self.option('taxonomy'))
        self.image_magick_convert = os.path.join(self.config.SOFTWARE_DIR, 'program/ImageMagick/bin/convert')
        '''
        self.ko_txt = os.path.join(
            self.config.SOFTWARE_DIR, 'database/KEGG/species/{}.ko.txt'.format(self.option('taxonomy'))
        )

        self.kegg_map_html = os.path.join(self.config.SOFTWARE_DIR, 'database/KEGG/map_html')
        '''
        if self.option('known_ko').is_set:
            self.known_ko = self.option('known_ko').prop['path']
        else:
            self.known_ko = 'None'

    def run(self):
        super(KeggAnnotTool, self).run()
        self.run_kegg_annotation()
        self.set_output()
        self.end()

    def run_kegg_annotation(self):
        cmd = '{} {} {} {} {} {} {} {} {} {} {} {} {} \{} {} {} {} {}'.format(
            self.python, # 0
            self.kegg_annotation_py, # 1
            self.rscript, # 2
            self.map4_r, # 3
            self.option('blast_xml').prop['path'], # 4
            None, # 5
            self.kegg_table, # 6
            self.pid_txt, # 7
            self.pathways, # 8
            self.pathway_table_xls, # 9
            self.kegg_layer_xls, # 10
            self.ko_txt, # 11
            self.option('link_bgcolor'), # 12
            self.option('png_bgcolor'), # 13
            self.image_magick_convert, # 14
            self.kegg_map_html, # 15
            self.known_ko, # 16
            self.option('kegg_version'), #17
        )
        cmd_name = 'run_kegg_annotation'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        kegg_table = os.path.join(self.output_dir, os.path.basename(self.kegg_table))
        if os.path.exists(kegg_table):
            os.remove(kegg_table)
        os.link(self.kegg_table, kegg_table)
        self.option('kegg_table').set_path(kegg_table)
        self.logger.info('succeed in linking {} to {}'.format(self.kegg_table, kegg_table))
        for source in [self.pid_txt, self.pathway_table_xls, self.pathway_table_xls, self.kegg_layer_xls]:
            link_name = os.path.join(self.output_dir, os.path.basename(source))
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        pathways = os.path.join(self.output_dir, os.path.basename(self.pathways))
        if os.path.isdir(pathways):
            shutil.rmtree(pathways)
        shutil.copytree(self.pathways, pathways)
        self.logger.info('succeed in copying {} to {}'.format(self.pathways, pathways))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'kegg_annot_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.annotation.kegg_annot',
            'instant': False,
            'options': {
                'blast_xml': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/annot_filter/output/kegg/blast.filter.xml',
                'taxonomy': 'Animals',
                'link_bgcolor': 'yellow',
                'png_bgcolor': '#FFFF00',
                'known_ko': '/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/KEGG/Homo_sapiens.GRCh38.pathway',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
