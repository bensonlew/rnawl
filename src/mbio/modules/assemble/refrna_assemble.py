# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.ref_rna.trans_step import *
import re
from mbio.files.sequence.file_sample import FileSampleFile


class RefrnaAssembleModule(Module):
    """
    拼接以及新转录本预测
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.09
    """
    def __init__(self, work_id):
        super(RefrnaAssembleModule, self).__init__(work_id)
        options = [
            {"name": "sample_bam_dir", "type": "infile", "format": "align.bwa.bam_dir"},  # 所有样本的bam文件夹
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            {"name": "assembly_GTF_list.txt", "type": "infile", "format": "assembly.merge_txt"},
            # 所有样本比对之后的bam文件路径列表
            {"name": "cpu", "type": "int", "default": 10},  # 软件所分配的cpu数量
            {"name": "fr_stranded", "type": "string", "default": "fr-unstranded"},  # 是否链特异性
            {"name": "strand_direct", "type": "string", "default": "none"},  # 链特异性时选择正负链
            {"name": "assemble_method", "type": "string", "default": "cufflinks"},  # 选择拼接软件
            {"name": "sample_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的gtf文件
            {"name": "merged_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的合并文件
            {"name": "cuff_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # compare后的gtf文件
            {"name": "new_transcripts_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 新转录本注释文件
            {"name": "new_transcripts_fa", "type": "outfile", "format": "sequence.fasta"},  # 新转录本序列文件
            {"name": "new_gene_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 新基因注释文件
            {"name": "merged_fa", "type": "outfile", "format": "sequence.fasta"},  # 新转录本序列文件
            {"name": "add_code_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 增加class_code后的注释文件
            {"name": "change_id_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 转换ID后的注释文件
            {"name": "change_id_fa", "type": "outfile", "format": "sequence.fasta"}  # 转换ID后的序列文件文件
        ]
        self.add_option(options)
        self.tools = []
        self.step.add_steps('stringtie', 'cufflinks', 'stringtie_merge', 'cuffmerge', 'cuffcompare', 'gffcompare',
                            'new_transcripts', 'refassemble_stat')
        self.sum_tools = []
        self.gtfs = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option('sample_bam_dir'):
            raise OptionError('必须输入样本文件夹，文件夹里的文件为bam格式')
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa')
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf')
        if self.option("fr_stranded") != "fr-unstranded" and self.option("strand_direct") == None:
            raise OptionError("当链特异性时必须选择正负链")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def stringtie_run(self):
        n = 0
        samples = os.listdir(self.option('sample_bam_dir').prop['path'])
        for f in samples:
            f = os.path.join(self.option('sample_bam_dir').prop['path'], f)
            stringtie = self.add_tool('assemble.stringtie')
            self.step.add_steps('stringtie_{}'.format(n))
            stringtie.set_options({
                "sample_bam": f,
                "ref_fa": self.option('ref_fa'),  # 此处不传prop['path']
                "ref_gtf": self.option('ref_gtf'),
                "fr_stranded": self.option("fr_stranded"),
                "strand_direct": self.option("strand_direct"),
            })
            step = getattr(self.step, 'stringtie_{}'.format(n))
            step.start()
            stringtie.on('end', self.finish_update, 'stringtie_{}'.format(n))
            self.tools.append(stringtie)
            self.sum_tools.append(stringtie)
            n += 1
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.stringtie_merge_run)
            self.step.stringtie_merge.start()
            self.step.update()
        for tool in self.tools:
            tool.run()

    def stringtie_merge_run(self):
        self.get_list()
        stringtie_merge = self.add_tool("assemble.stringtie_merge")
        stringtie_merge.set_options({
            "assembly_GTF_list.txt": gtffile_path,
            "ref_fa": self.option('ref_fa'),
            "ref_gtf": self.option('ref_gtf'),
        })
        stringtie_merge.on('end', self.gffcompare_run)
        stringtie_merge.run()
        self.sum_tools.append(stringtie_merge)
        self.step.stringtie_merge.finish()
        self.step.gffcompare.start()
        self.step.update()

    def cufflinks_run(self):
        n = 0
        samples = os.listdir(self.option('sample_bam_dir').prop['path'])
        for f in samples:
            f = os.path.join(self.option('sample_bam_dir').prop['path'], f)
            cufflinks = self.add_tool('assemble.cufflinks')
            self.step.add_steps('cufflinks_{}'.format(n))
            cufflinks.set_options({
                "sample_bam": f,
                "ref_fa": self.option('ref_fa'),
                "ref_gtf": self.option('ref_gtf'),
                "fr_stranded": self.option("fr_stranded"),
            })
            step = getattr(self.step, 'cufflinks_{}'.format(n))
            step.start()
            cufflinks.on('end', self.finish_update, 'cufflinks_{}'.format(n))
            self.tools.append(cufflinks)
            self.sum_tools.append(cufflinks)
            n += 1
        if len(self.tools) == 1:
            self.tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.tools, self.cuffmerge_run)
            self.step.cuffmerge.start()
            self.step.update()
        for tool in self.tools:
            tool.run()

    def cuffmerge_run(self):
        self.get_list()
        cuffmerge = self.add_tool("assemble.cuffmerge")
        cuffmerge.set_options({
            "assembly_GTF_list.txt": gtffile_path,
            "ref_fa": self.option('ref_fa'),
            "ref_gtf": self.option('ref_gtf'),
        })
        cuffmerge.on('end', self.gffcompare_run)
        cuffmerge.run()
        self.sum_tools.append(cuffmerge)
        self.step.cuffmerge.finish()
        self.step.gffcompare.start()
        self.step.update()

    def gffcompare_run(self):
        merged_gtf = ""
        if self.option("assemble_method") == "cufflinks":
            merged_gtf = os.path.join(self.work_dir, "Cuffmerge/output/merged.gtf")
        elif self.option("assemble_method") == "stringtie":
            merged_gtf = os.path.join(self.work_dir, "StringtieMerge/output/merged.gtf")
        gffcompare = self.add_tool("assemble.gffcompare")
        gffcompare.set_options({
             "merged_gtf": merged_gtf,
             "ref_gtf": self.option('ref_gtf'),
         })
        gffcompare.on('end', self.new_transcripts_run)
        gffcompare.run()
        self.sum_tools.append(gffcompare)
        self.step.gffcompare.finish()
        self.step.new_transcripts.start()
        self.step.update()

    def new_transcripts_run(self):
        self.logger.info(self.option('ref_gtf').prop['path'])  # 此处阻塞
        tmap = ""
        merged_gtf = ""
        if self.option("assemble_method") == "cufflinks":
            tmap = os.path.join(self.work_dir, "Cuffmerge/output/cuffcmp.merged.gtf.tmap")
            merged_gtf = os.path.join(self.work_dir, "Cuffmerge/output/merged.gtf")
        elif self.option("assemble_method") == "stringtie":
            tmap = os.path.join(self.work_dir, "StringtieMerge/output/cuffcmp.merged.gtf.tmap")
            merged_gtf = os.path.join(self.work_dir, "StringtieMerge/output/merged.gtf")
        new_transcripts = self.add_tool("assemble.new_transcripts")
        new_transcripts.set_options({
            "tmap": tmap,
            "merged_gtf": merged_gtf,
            "ref_fa": self.option('ref_fa'),
            "ref_gtf": self.option('ref_gtf'),
        })
        new_transcripts.on('end', self.refassemble_stat_run)
        new_transcripts.run()
        self.sum_tools.append(new_transcripts)
        self.step.new_transcripts.finish()
        self.step.refassemble_stat.start()
        self.step.update()

    def refassemble_stat_run(self):
        for tool in self.sum_tools:
            self.linkdir(tool.output_dir, 'assembly_newtranscripts')
        self.refassemble_stat = self.add_tool("assemble.refassemble_stat")
        self.refassemble_stat.set_options({
            "all_files_dir": self.work_dir + '/assembly_newtranscripts',
            "assemble_method": self.option("assemble_method"),
        })
        self.refassemble_stat.on('end', self.set_output)
        self.refassemble_stat.run()
        self.sum_tools.append(self.refassemble_stat)
        self.step.refassemble_stat.finish()
        self.step.update()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                if os.path.exists(newfiles[i]):
                    os.remove(newfiles[i])
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                if os.path.exists(newdir):
                    os.remove(newdir)
                os.link(oldfiles[i], newdir)

    def set_output(self):
        if self.option("assemble_method") == "cufflinks":
            gtf_file = 'Cufflinks'
            merge = 'Cuffmerge'
        else:
            gtf_file = 'Stringtie'
            merge = 'StringtieMerge'
        compare = 'Gffcompare'
        new_transcripts = 'NewTranscripts'
        statistics = 'Statistics'
        gtf_dir = os.path.join(self.output_dir, gtf_file)
        if os.path.exists(gtf_dir):
            shutil.rmtree(gtf_dir)
        os.mkdir(gtf_dir)
        merge_dir = os.path.join(self.output_dir, merge)
        if os.path.exists(merge_dir):
            shutil.rmtree(merge_dir)
        os.mkdir(merge_dir)
        compare_dir = os.path.join(self.output_dir, compare)
        if os.path.exists(compare_dir):
            shutil.rmtree(compare_dir)
        os.mkdir(compare_dir)
        new_transcripts_dir = os.path.join(self.output_dir, new_transcripts)
        if os.path.exists(new_transcripts_dir):
            shutil.rmtree(new_transcripts_dir)
        os.mkdir(new_transcripts_dir)
        statistics_dir = os.path.join(self.output_dir, statistics)
        if os.path.exists(statistics_dir):
            shutil.rmtree(statistics_dir)
        os.mkdir(statistics_dir)
        old_dir = self.work_dir + '/assembly_newtranscripts/'
        for files in os.listdir(old_dir):
            if files.endswith("_out.gtf") or files.endswith("_out.fa"):
                if os.path.exists(gtf_dir + "/" + files):
                    os.remove(gtf_dir + "/" + files)
                os.link(old_dir + files, gtf_dir + "/" + files)
            elif files.endswith("merged.gtf") or files.endswith("merged.fa"):
                if os.path.exists(merge_dir + "/" + files):
                    os.remove(merge_dir + "/" + files)
                os.link(old_dir + files, merge_dir + "/" + files)
            elif files.startswith("cuffcmp."):
                if os.path.exists(compare_dir + "/" + files):
                    os.remove(compare_dir + "/" + files)
                os.link(old_dir + files, compare_dir + "/" + files)
            elif files.startswith("new_transcripts.") or files.startswith("new_genes.") or files.startswith("old_trans.gtf") or files.startswith("old_genes.gtf"):
                if os.path.exists(new_transcripts_dir + "/" + files):
                    os.remove(new_transcripts_dir + "/" + files)
                os.link(old_dir + files, new_transcripts_dir + "/" + files)
        txt_files = os.listdir(self.work_dir + '/RefassembleStat/output')
        for files in txt_files:
            move_files = os.path.join(self.work_dir + '/RefassembleStat/output', files)
            # self.logger.info(move_files)
            # self.logger.info(statistics_dir + "/" + files)
            if os.path.exists(statistics_dir + "/" + files):
                os.remove(statistics_dir + "/" + files)
            os.link(move_files, statistics_dir + "/" + files)
        # self.option("merged_gtf").set_path(merge_dir + '/merged.gtf')
        # self.option("merged_fa").set_path(merge_dir + "/merged.fa")
        self.option("new_transcripts_gtf").set_path(new_transcripts_dir + "/new_transcripts.gtf")
        # self.option("new_transcripts_fa").set_path(new_transcripts_dir + "/new_transcripts.fa")
        self.option("new_gene_gtf").set_path(new_transcripts_dir + "/new_genes.gtf")
        self.option("change_id_gtf").set_path(merge_dir + "/change_id_merged.gtf")
        # self.option("add_code_gtf").set_path(merge_dir + "/add_code_merged.gtf")
        self.option("cuff_gtf").set_path(compare_dir + '/cuffcmp.annotated.gtf')
        # self.option("change_id_fa").set_path(merge_dir + "/change_id_merged.fa")
        self.end()

    def run(self):
        super(RefrnaAssembleModule, self).run()
        if self.option("assemble_method") == "cufflinks":
            self.cufflinks_run()
        elif self.option("assemble_method") == "stringtie":
            self.stringtie_run()

    def get_list(self):
        gtffile_path = os.path.join(self.work_dir, "assembly_gtf.txt")
        global gtffile_path
        with open(gtffile_path, "w+") as w:
            for gtf in self.tools:
                for f in os.listdir(gtf.output_dir):
                    m = re.match(".+\.gtf", f)
                    if m:
                        file_path = os.path.join(gtf.output_dir, f)
                        w.write(file_path + "\n")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("assemble_method") == "stringtie":
            result_dir.add_relpath_rules([
                [r".", "", "结果输出目录"],
                ["Stringtie", "", "拼接后的各样本文件夹"],
                ["StringtieMerge", "", "拼接组装合并之后结果文件夹"],
                ["StringtieMerge/merged.gtf", "gtf", "样本合并之后的注释文件"],
                ["StringtieMerge/merged.fa", "fasta", "样本合并之后的序列文件"],
                ["Gffcompare", "", "进行新转录本预测后的结果文件夹"],
                ["Gffcompare/cuffcmp.annotated.gtf", "", "进行新转录本预测后的结果文件"],
                ["Gffcompare/cuffcmp.merged.gtf.refmap", "", "进行新转录本预测后的结果文件"],
                ["Gffcompare/cuffcmp.merged.gtf.tmap", "", "进行新转录本预测后的结果文件"],
                ["NewTranscripts", "", "新转录本结果文件夹"],
                ["NewTranscripts/new_transcripts.gtf", "gtf", "新转录本注释文件"],
                ["NewTranscripts/new_transcripts.fa", "fa", "新转录本序列文件"],
                ["NewTranscripts/new_genes.gtf", "gtf", "新基因注释文件"],
                ["NewTranscripts/new_genes.fa", "fa", "新基因序列文件"],
                ["NewTranscripts/old_trans.gtf", "gtf", "已知转录本注释文件"],
                ["NewTranscripts/old_genes.gtf", "gtf", "已知基因注释文件"],
                ["Statistics", "", "统计信息的结果文件夹"],
                ["Statistics/code_num.txt", "txt", "class_code统计文件"],
            ])
            result_dir.add_regexp_rules([
                [r"Stringtie/.*_out\.gtf$", "gtf", "样本拼接之后的注释文件"],
                [r"Stringtie/.*_out\.fa$", "fasta", "样本拼接之后序列文件"],
                [r"Statistics/trans_count_stat_.*\.txt$", "txt", "新转录本步长统计文件"],
                [r"Statistics/old_.*\.txt$", "txt", "统计结果文件"],
                [r"Statistics/new_.*\.txt$", "txt", "统计结果文件"],

            ])
        if self.option("assemble_method") == "cufflinks":
            result_dir.add_relpath_rules([
                [r".", "", "结果输出目录"],
                ["Cufflinks", "", "拼接后的各样本文件夹"],
                ["Cuffmerge", "", "拼接组装合并之后结果文件夹"],
                ["Cuffmerge/merged.gtf", "gtf", "样本合并之后的注释文件"],
                ["Cuffmerge/merged.fa", "fasta", "样本合并之后的序列文件"],
                ["Gffcompare", "", "进行新转录本预测后的结果文件夹"],
                ["Gffcompare/cuffcmp.annotated.gtf", "", "进行新转录本预测后的结果文件"],
                ["Gffcompare/cuffcmp.merged.gtf.refmap", "", "进行新转录本预测后的结果文件"],
                ["Gffcompare/cuffcmp.merged.gtf.tmap", "", "进行新转录本预测后的结果文件"],
                ["NewTranscripts", "", "新转录本结果文件夹"],
                ["NewTranscripts/new_transcripts.gtf", "gtf", "新转录本注释文件"],
                ["NewTranscripts/new_transcripts.fa", "fa", "新转录本序列文件"],
                ["NewTranscripts/new_genes.gtf", "gtf", "新基因注释文件"],
                ["NewTranscripts/new_genes.fa", "fa", "新基因序列文件"],
                ["NewTranscripts/old_trans.gtf", "gtf", "已知转录本注释文件"],
                ["NewTranscripts/old_genes.gtf", "gtf", "已知基因注释文件"],
                ["Statistics", "", "统计信息的结果文件夹"],
                ["Statistics/code_num.txt", "txt", "class_code统计文件"],
            ])
            result_dir.add_regexp_rules([
                [r"Cufflinks/.*_out\.gtf$", "gtf", "样本拼接之后的注释文件"],
                [r"Cufflinks/.*_out\.fa$", "fasta", "样本拼接之后序列文件"],
                [r"Statistics/trans_count_stat_.*\.txt$", "txt", "新转录本步长统计文件"],
                [r"Statistics/old_.*\.txt$", "txt", "统计结果文件"],
                [r"Statistics/new_.*\.txt$", "txt", "统计结果文件"],

            ])
        super(RefrnaAssembleModule, self).end()