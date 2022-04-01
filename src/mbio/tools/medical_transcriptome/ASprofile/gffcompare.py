# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
import re
import regex
from mbio.packages.whole_transcriptome.utils import runcmd

class GffcompareAgent(Agent):
    """
    有参转录组gffcompare比较
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2017.02.10
    """
    def __init__(self, parent):
        super(GffcompareAgent, self).__init__(parent)
        options = [
            # {"name": "merged_gtf", "type": "string", "default": ""},  # 拼接合并之后的转录本文件
            {'name': 'merged_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},  # 参考基因的注释文件
            {"name": "gff_gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {'name': 'sample', 'type': 'string'},
            {'name': 'seq_path', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'gene_structure.gtf'},
            {'name': 'tmap', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'old_gtf', 'type': 'outfile', 'format': 'gene_structure.gtf'},
            {'name': 'new_gtf', 'type': 'outfile', 'format': 'gene_structure.gtf'},
            {'name': 'all_gtf', 'type': 'outfile', 'format': 'gene_structure.gtf'}
        ]
        self.add_option(options)
        self.step.add_steps("gffcompare")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gffcompare.start()
        self.step.update()

    def stepfinish(self):
        self.step.gffcompare.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('merged_gtf'):
            raise OptionError('必须输入所有样本拼接合并后gtf文件', code = "33703601")
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf', code = "33703602")
        if not self.option('seq_path'):
            raise OptionError('必须输入参考基因组ref.gasta', code = '33703603')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["gffcmp.annotated.gtf", "gtf", "gffcompare输出的相关文件"]
        ])
        super(GffcompareAgent, self).end()


class GffcompareTool(Tool):
    def __init__(self, config):
        super(GffcompareTool, self).__init__(config)
        # self._version = "v0.9.8.linux_x86_64"
        # self.gffcompare_path = 'bioinfo/rna/gffcompare-0.9.8.Linux_x86_64/'
        # tmp = os.path.join(self.config.SOFTWARE_DIR, self.gffcompare_path)
        # tmp1 = tmp + ":$PATH"
        # self.logger.debug(tmp1)
        # self.set_environ(PATH=tmp1)



        self.file = {
            'merged_gtf': os.path.join(self.work_dir, 'merged.gtf'),
            'tmap': os.path.join(self.work_dir, 'gffcmp.merged.gtf.tmap'),
            'ref_gtf': os.path.join(self.output_dir, 'ref.gtf'),
            'new_gtf': os.path.join(self.output_dir, 'new.gtf'),
            'all_gtf': os.path.join(self.output_dir, 'all.gtf'),
            'ref_fasta': os.path.join(self.output_dir, 'ref.fasta'),
            'new_fasta': os.path.join(self.output_dir, 'new.fasta'),
            'all_fasta': os.path.join(self.output_dir, 'all.fasta'),
        }
        self.program = {
            'gffcompare': 'bioinfo/rna/gffcompare-0.9.8.Linux_x86_64/gffcompare',
            'cuffcompare': 'bioinfo/rna/cufflinks-2.2.1/cuffcompare',
            'python': 'program/Python/bin/python',
            'gffread': 'bioinfo/rna/cufflinks-2.2.1/gffread',
        }
        self.script = {
            'feature_class': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/feature_class.py'),
            'step_code': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/assembly/step_code.py')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(GffcompareTool, self).run()
        self.run_gffcompare()
        self.run_feature_class()
        self.run_gffread_new()
        self.run_gffread_ref()
        self.run_gffread_all()
        # self.run_filter_gtf()
        self.set_output()
        self.end()

    def run_gffcompare(self):
        if os.path.isfile(self.file['merged_gtf']):
            os.remove(self.file['merged_gtf'])
        os.link(self.option('merged_gtf').path, self.file['merged_gtf'])
        cmd = '{} '.format(self.program['gffcompare'])
        cmd += ' -r {}'.format(self.option('ref_gtf').path)
        cmd += ' -s {}'.format(self.option('seq_path').path)
        cmd += ' -V'
        cmd += ' {}'.format(self.file['merged_gtf'])
        runcmd(self, 'run_gffcompare', cmd)

    def run_feature_class(self):
        cmd = '{} {}'.format(self.program['python'], self.script['feature_class'])
        cmd += ' -r {}'.format(self.option('ref_gtf').path)
        cmd += ' -m {}'.format(self.file['merged_gtf'])
        cmd += ' -t {}'.format(self.file['tmap'])
        cmd += ' -o {}'.format(self.output_dir)
        runcmd(self, 'run_feature_class', cmd)


    def run_gffread_ref(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.file['ref_gtf'])
        cmd += ' -g {}'.format(self.option('seq_path').path)
        cmd += ' -w {}'.format(self.file['ref_fasta'])
        runcmd(self, 'run_gffread_ref', cmd)

    def run_gffread_new(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.file['new_gtf'])
        cmd += ' -g {}'.format(self.option('seq_path').path)
        cmd += ' -w {}'.format(self.file['new_fasta'])
        runcmd(self, 'run_gffread_new', cmd)

    def run_gffread_all(self):
        cmd = '{} {} -T'.format(self.program['gffread'], self.file['all_gtf'])
        cmd += ' -g {}'.format(self.option('seq_path').path)
        cmd += ' -w {}'.format(self.file['all_fasta'])
        runcmd(self, 'run_gffread_all', cmd)



    # def run_filter_gtf(self):
    #     """
    #     过滤gffcompare出的异常转录本，新基因存在已知转录本
    #     """
    #     known_tran2gene = {}
    #     known_gtf = self.option('ref_gtf').prop['path']
    #     merge_raw_gtf = self.work_dir + "/gffcmp.annotated.gtf"
    #     merge_gtf = self.work_dir + "/gffcmp.annotated.filter.gtf"
    #
    #     with open(known_gtf, "r") as file:
    #         for line in file:
    #             line = line.strip()
    #             content_m = regex.match(
    #                 r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
    #                 line)
    #             if content_m:
    #                 if 'transcript_id' in content_m.captures(6):
    #                     tran_id = content_m.captures(7)[0]
    #                     gene_id = content_m.captures(10)[0]
    #                 else:
    #                     tran_id = content_m.captures(10)[0]
    #                     gene_id = content_m.captures(7)[0]
    #                 known_tran2gene[tran_id] = gene_id
    #
    #     with open(merge_raw_gtf, "r") as merge_raw_file, open(merge_gtf, "w") as merge_file:
    #         for line in merge_raw_file:
    #             line = line.strip()
    #             content_m = regex.match(
    #                 r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
    #                 line)
    #             if content_m:
    #                 if 'transcript_id' in content_m.captures(6):
    #                     tran_id = content_m.captures(7)[0]
    #                     gene_id = content_m.captures(10)[0]
    #                 else:
    #                     tran_id = content_m.captures(10)[0]
    #                     gene_id = content_m.captures(7)[0]
    #                 if known_tran2gene.has_key(tran_id):
    #                     if known_tran2gene[tran_id] == gene_id:
    #                         merge_file.write(line + "\n")
    #                     else:
    #                         self.logger.info('gffcompare gene_id异常: {}为已知转录本，不属于基因{}'.format(tran_id,gene_id))
    #                 else:
    #                     merge_file.write(line + "\n")
    #             else:
    #                 merge_file.write(line + "\n")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # if os.path.exists(self.output_dir + "/new.gtf"):
        #     os.remove(self.output_dir + "/new.gtf")
        # if os.path.exists(self.output_dir + "/ref.gtf"):
        #     os.remove(self.output_dir + "/ref.gtf")
        self.option('old_gtf').set_path(self.output_dir + "/old.gtf")
        self.option('new_gtf').set_path(self.output_dir + "/new.gtf")
        self.option('all_gtf').set_path(self.output_dir + "/all.gtf")
        self.logger.info("设置拼接比较结果目录成功")

