# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
# last_modify:2019.06.24

from biocluster.module import Module
import os
import shutil
import glob
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest
import datetime
import time
import pickle
import psutil, shlex


class DiffStringPicturesModule(Module):
    def __init__(self, work_id):
        super(DiffStringPicturesModule, self).__init__(work_id)
        options = [
            {'name': 'diff_path', 'type': 'infile', 'format': 'itraq_and_tmt.common_dir'},
            {'name': 'string_xml', 'type': 'infile', 'format': 'itraq_and_tmt.common'},
            {'name': 'gene_list', 'type': 'string'},
            {'name': 'identity', 'type': 'float', 'default': 98},
            {'name': 'max_num', 'type': 'int', 'default': 300},
            {'name': 'species', 'type': 'int', 'default': 0},
            {'name': 'useblast', 'type': 'string', 'default': 'no'},
        ]
        self.add_option(options)
        self.config = Config()
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + '/itraq_and_tmt/get_string_picture.py'
        self.tools = list()

    def get_most_common_specie(self, xml):
        from Bio.Blast import NCBIXML
        from collections import Counter
        species = list()
        with open(xml, 'r') as blast_r:
            records = NCBIXML.parse(blast_r)
            for rec in records:
                for align in rec.alignments:
                    for hsp in align.hsps:
                        hit = align.hit_id
                        des = align.hit_def
                        if u'.' in des and '|' not in des:
                            hit = des
                        ident = float(hsp.identities)
                        hit_len = float(hsp.align_length)
                        if ident / hit_len * 100 >= self.option('identity'):
                            # specie = re.split('.', hit, maxsplit=1)[0]
                            specie = hit.split('.')[0]
                            # self.logger.info(hit)
                            species.append(specie)
        if not species:
            return None
        # 有的id在string数据库里面没有，只能再加一部判断
        # most_s = Counter(species).most_common()[0][0]
        with open(self.config.SOFTWARE_DIR + '/database/Annotation/all/String/string11.5/ppi_species.v11.5.txt') as ps:
            hit_species = [line.strip().split('\t')[1] for line in ps if line.strip()]
        for spe, count in Counter(species).most_common():
            if spe in hit_species:
                self.logger.info('选择第%s个物种%s 作为优势物种进行ppi分析'%(str(count), spe))
                return int(spe)
        return None

    def run_string_pictures(self):
        if self.option('species') == -1 or not self.option('species'):
            specie = self.get_most_common_specie(self.option('string_xml').prop['path'])
        else:
            specie = self.option('species')
        if not specie:
            return
        self.option('species', specie)
        useblast = 'no'
        if self.option('useblast').lower() == 'yes':
            useblast = 'yes'
        # diff_files = glob.glob(os.path.join(self.option('diff_path').prop['path'], '*_vs_*'))
        # for diff in diff_files:
        #     cmp = os.path.basename(diff).split('_diff.xls')[0]
        #     acc_de = os.path.join(self.work_dir, '%s_all_protein.list'%cmp)
        #     acc_up = os.path.join(self.work_dir, '%s_up_protein.list'%cmp)
        #     acc_down = os.path.join(self.work_dir, '%s_down_protein.list'%cmp)
        #     diff_df = pd.read_csv(diff, sep='\t', index_col=0)
        #     uplist = diff_df[((diff_df['regulate'] == 'up') & (diff_df['significant'] == 'yes'))].index.tolist()
        #     downlist = diff_df[((diff_df['regulate'] == 'down') & (diff_df['significant'] == 'yes'))].index.tolist()
        #     with open(acc_de, 'w') as dew, open(acc_up, 'w') as uw, open(acc_down, 'w') as dw:
        #         dew.write('\n'.join(uplist + downlist))
        #         uw.write('\n'.join(uplist))
        #         dw.write('\n'.join(downlist))
        #     options = dict(
        #         fake="fake",
        #     )
        #     de_tool = self.add_tool("itraq_and_tmt.proteinset.fakefake")
        #     de_tool.set_options(options)
        #     self.tools.append(de_tool)

        gene_list = self.option('gene_list').split(';')
        with open(os.path.join(self.work_dir, 'protein.list'), 'w') as p:
            for i in gene_list:
                p.write(i + '\n')
        if glob.glob(os.path.join(self.output_dir, '*', '*.svg')):
            return
        files = glob.glob(os.path.join(self.work_dir, '*protein.list'))
        cmd = self.python_path + ' ' + self.script
        cmd += ' -specie ' + str(specie)
        cmd += ' -list_path ' + self.work_dir
        if useblast == 'yes':
            cmd += ' -vsstring ' + self.option('string_xml').prop['path']
        cmd += ' -identity ' + str(self.option('identity'))
        cmd += ' -max_num ' + str(self.option('max_num'))
        cmd += ' -useblast ' + useblast
        cmd += ' -out ' + self.output_dir
        self.logger.info(cmd)
        # os.system(cmd)
        stdout = open(os.path.join(self.work_dir, 'scrapy.log'), 'w')
        fn = stdout.fileno()
        p = psutil.Popen(shlex.split(cmd), stdout=fn)
        begin = datetime.datetime.now()
        while 1:
            time.sleep(20)
            self.logger.info(p.status())
            if p.status() == 'zombie':
                break
                # pass
            else:
                use_t = (datetime.datetime.now() - begin).seconds
                if use_t > 60*10*len(files) or use_t > 5 * 60 * 60:
                    self.logger.info('运行时间太长，已经被杀掉')
                    p.terminal()
                    break
        stdout.close()
        # if (datetime.datetime.now() - begin).seconds < 60:
        if not glob.glob(self.output_dir + '/*/*'):
            self.logger.info('爬取失败')
            shutil.rmtree(self.output_dir)
            os.mkdir(self.output_dir)

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        return True

    def set_output(self):
        self.end()

    def run(self):
        # 为了保证同一时间只能运行一个爬虫任务
        time_file = self.config.PACKAGE_DIR + '/itraq_and_tmt/string.time'
        if not os.path.exists(time_file):
            with open(time_file, 'w') as tw:
                pickle.dump(datetime.datetime.now(), tw)
        wait_start = datetime.datetime.now()
        while 1:
            tr = open(time_file)
            last_time = pickle.load(tr)
            now_time = datetime.datetime.now()
            if not last_time:
                if (now_time - wait_start).seconds < 720:
                    time.sleep(300)
                    continue
                else:
                    break

            # if (now_time-last_time).seconds < 720:
            if (now_time-last_time).seconds < 20:
                self.logger.info('需要睡5分钟')
                time.sleep(300)
            else:
                break
        with open(time_file, 'w') as tw:
            # pickle.dump(datetime.datetime.now(), tw)
            pickle.dump('', tw)

        self.run_string_pictures()
        with open(time_file, 'w') as tw:
            pickle.dump(datetime.datetime.now(), tw)
        super(DiffStringPicturesModule, self).run()
        if not self.tools:
            super(DiffStringPicturesModule, self).end()
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(DiffStringPicturesModule, self).end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "diff_string_pictures_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.diff_string_pictures",
            "instant": False,
            "options": dict(
                # diff_path="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/Diff/output",
                gene_list='trpA',
                species='9606',
                # string_xml="/mnt/ilustre/users/sanger-dev/workspace/20190624/Labelfree_tsg_34554/ProteinAnnotation/output/blast_xml/string.xml",
            )
        }

        data['id'] += '_fyt'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
