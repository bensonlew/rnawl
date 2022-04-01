# -*- coding: utf-8 -*-

# __author__ = 'linfang.jin'
# time: 2017/1/16 15:47
import importlib
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.ref_rna.trans_step import *
import re
from mbio.files.sequence.file_sample import FileSampleFile


class AsprofileModule(Module):
    '''
    '''

    def __init__(self, work_id):
        super(AsprofileModule, self).__init__(work_id)

        options = [
            {"name": "sample_transcripts_gtf_dir", "type": "infile", "format": "ref_rna.gene_structure.gtf_dir"},
            {"name": "sample_transcripts_tmap_dir", "type": "infile", "format": "ref_rna.gene_structure.tmap_dir"},
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna.assembly.gtf"},
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)
        self.ref_link = ""
        self.as_tools = []

    def check_options(self):
        if self.option('ref_fa') is None:
            raise OptionError("必须设置参考基因组文件")
        if self.option('sample_transcripts_gtf_dir') is None:
            raise OptionError("必须设置样本转录本文件夹")
        if self.option('ref_gtf') is None:
            raise OptionError('必须设置参考基因组注释文件')
        if self.option('sample_transcripts_tmap_dir') is None:
            raise OptionError('必须设置样本转录本tmap文件夹')
        return True


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def link_ref(self):
        ref_fa = self.option('ref_fa').path
        self.ref_link = self.work_dir + "/" + os.path.basename(ref_fa)
        self.logger.info(self.ref_link)
        if os.path.exists(self.ref_link):
            os.remove(self.ref_link)
        os.symlink(ref_fa, self.ref_link)

    def hdrs_run(self):
        """发送信息给前端"""
        self.link_ref()
        self.step.add_steps('asprofile_hdrs')
        self.step.asprofile_hdrs.start()
        self.step.update()
        self.logger.info("模块开始启动hdrs程序")
        self.asprofile_hdrs = self.add_tool("ref_rna.gene_structure.asprofile_hdrs")
        self.asprofile_hdrs.set_options({
            "ref_fa": self.ref_link
        })

        """绑定下一个将要运行的步骤"""
        self.asprofile_hdrs.on("end", self.finish_update, 'asprofile_hdrs')
        self.asprofile_hdrs.on('end', self.multi_as_run)
        self.asprofile_hdrs.run()

    def get_gtf_tmap_dic(self,gtf_dir,tmap_dir):
        d = {}
        sample_gtfs = os.listdir(gtf_dir)
        sample_tmaps = os.listdir(tmap_dir)
        for gtf in sample_gtfs:
            sample_name = re.match(r'^(\S+)_out.gtf$', gtf.strip()).group(1)
            gtf_path = os.path.join(self.option("sample_transcripts_gtf_dir").path, gtf)
            for tmap in sample_tmaps:
                if re.search(r'.*{}.*'.format(sample_name),tmap):
                    tmap_path = os.path.join(self.option("sample_transcripts_tmap_dir").path, tmap)
                    d[gtf_path] = tmap_path

        return d

    def multi_as_run(self):
        gtf_tmap_dic = self.get_gtf_tmap_dic(self.option('sample_transcripts_gtf_dir').path,self.option('sample_transcripts_tmap_dir').path)
        for gtf_path in gtf_tmap_dic.keys():
            tmap_path = gtf_tmap_dic[gtf_path]
            sample_name = os.path.basename(gtf_path).split('_')[0]
            asprofile_as = self.add_tool('ref_rna.gene_structure.asprofile_as')
            self.step.add_steps('asprofile_as_{}'.format(sample_name))
            asprofile_as.set_options({
                'hdrs': os.path.join(self.work_dir, "AsprofileHdrs/output/ref.fa.hdrs"),
                "ref_gtf": self.option("ref_gtf").path,
                'sample_gtf': gtf_path,
                'sample_tmap': tmap_path,
            })
            step = getattr(self.step, 'asprofile_as_{}'.format(sample_name))
            step.start()
            """绑定下一个将要运行的步骤"""
            asprofile_as.on('end', self.finish_update, 'asprofile_as_{}'.format(sample_name))
            asprofile_as.on('end', self.set_output, 'asprofile_as_{}'.format(sample_name))
            self.as_tools.append(asprofile_as)
        for tool in self.as_tools:
            tool.run()

    def set_output(self):
        self.logger.info('set output')
        for as_tool in self.as_tools:
            outfiles = os.listdir(as_tool.output_dir)
            for f in outfiles:
                if re.match(r'^(\S+).as.*', f.strip()):
                    f_path = os.path.join(as_tool.output_dir, f)
                    target = os.path.join(self.output_dir, f)
                    os.link(f_path, target)
        self.logger.info("set output done")
        self.end()

    def run(self):
        self.hdrs_run()
        super(AsprofileModule, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.', '', '可变剪接事件结果输出目录']
            ]
        )
        result_dir.add_regexp_rules([
            [".as.nr$", '', '样本非冗余可变剪接事件详情表'],
            [".as.summary$", '', '样本非冗余可变剪接事件总结'],
            ['.as$','','样本可变剪接事件详情原始表']

        ])
        super(AsprofileModule,self).end()
