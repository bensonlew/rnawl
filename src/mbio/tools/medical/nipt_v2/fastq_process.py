## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
# last_modify:20170424

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class FastqProcessAgent(Agent):
    """
    产筛的生信shell部分功能
    处理fastq，得到bed2文件

    参数run_whole为false时，可将shell流程拆分为各个命令跑，可查看中间文件，但效率有所降低
    参数single控制，下机数据是否为单端数据，true时为单端，是单端数据，运行单端流程
    lasted modified by hongdong 20171218
    """
    def __init__(self, parent):
        super(FastqProcessAgent, self).__init__(parent)
        options = [
            {"name": "sample_id", "type": "string"},
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "run_whole", "type": "string", "default": "true"},
            {"name": "single", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.step.add_steps("fastq2bed")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fastq2bed.start()
        self.step.update()

    def stepfinish(self):
        self.step.fastq2bed.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('sample_id'):
            raise OptionError("必须输入样本名")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 12
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [".bed.2", "bed.2", "信息表"],
        ])
        super(FastqProcessAgent, self).end()


class FastqProcessTool(Tool):
    """
    比对生成bed文件
    """
    def __init__(self, config):
        super(FastqProcessTool, self).__init__(config)
        self._version = '1.0.1'
        # self.R_path = 'program/R-3.3.1/bin/'
        self.bed_dir = self.config.SOFTWARE_DIR + "/database/human/hg38_nipt/bed_file/"
        self.script_path = 'bioinfo/medical/scripts/'
        self.java_path = self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin/java'
        self.picard_path = self.config.SOFTWARE_DIR + '/bioinfo/medical/picard-tools-2.2.4/picard.jar'
        self.sam_path = '/bioinfo/align/samtools-1.3.1/samtools'

        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/FastQc')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/bwa-0.7.15')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/bioawk')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/seqtk-master')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.4')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/seq/samblaster-0.1.24')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/Python/bin')  # cutadapt
        
        self.ref1 = self.config.SOFTWARE_DIR + '/database/human/hg38.chromosomal_assembly/ref.fa'
        self.ref = self.config.SOFTWARE_DIR + '/database/human/hg38_nipt/nchr.fa'
        self.bed_ref = self.config.SOFTWARE_DIR + '/database/human/hg38_nipt/nchr.20k.gmn.bed'

    def run_tf(self):
        pre_cmd = '{}nipt_fastq_pre.sh {} {}'.format(self.script_path, self.option("fastq_path").prop['path'],
                                                     self.option('sample_id'))
        self.logger.info(pre_cmd)
        cmd = self.add_command("pre_cmd", pre_cmd).run()
        self.wait(cmd)

        if cmd.return_code == 0:
            self.logger.info("处理接头成功")
        else:
            raise Exception("处理接头出错")

        seq_merge = '{}nipt_merge_align.sh {} {} '.\
            format(self.script_path, self.option("sample_id"), self.ref1)
        self.logger.info(seq_merge)
        cmd = self.add_command("seq_merge", seq_merge).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("seqtk mergepe成功")
        else:
            raise Exception("seqtk mergepe出错")

        cut_adapt = '/bioinfo/medical/cutadapt-1.10-py27_0/bin/cutadapt --format fastq --zero-cap -q 1 --trim-n ' \
                    '--minimum-length 30 --times 7 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {}.cut.fastq  ' \
                    '{}.cutN.fastq'.format(self.option("sample_id"), self.option("sample_id"))
        self.logger.info(cut_adapt)
        cmd = self.add_command("cut_adapt", cut_adapt).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("cutadapt去接头成功")
        else:
            raise Exception("cutadapt去接头出错")

        cut_50 = '{}nipt_cut50.sh {} {}'.format(self.script_path, self.option('sample_id'), self.ref)
        self.logger.info(cut_50)
        cmd = self.add_command("cut_50", cut_50).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("cut_50截取成功")
        else:
            raise Exception("cut_50截取出错")

        sam_cutbam = "{} view -h -@ 10 {}.cut.bam -o {}.temp.cut.sam".format(self.sam_path, self.option('sample_id'),
                                                                             self.option('sample_id'))
        self.logger.info(sam_cutbam)
        cmd = self.add_command("sam_cutbam", sam_cutbam).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam_cutbam成功")
        else:
            raise Exception("sam_cutbam出错")
        file_ = self.option('sample_id') + '.temp.cut.sam'
        with open(file_, 'r') as f:
            lines = f.readlines()
        with open(file_, 'w+') as f_w:
            for line in lines:
                if 'XT:A:U' in line or re.search('^@', line):
                    f_w.write(line)
                else:
                    continue
        sam_cut_uniq = '{} view -@ 10 -bS {}.temp.cut.sam -o {}.cut.uniq.bam'\
            .format(self.sam_path, self.option('sample_id'), self.option('sample_id'))
        self.logger.info(sam_cut_uniq)
        cmd = self.add_command("sam_cut_uniq", sam_cut_uniq).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam_cut_uniq成功")
        else:
            raise Exception("sam_cut_uniq出错") 

        sam_sort = '{} sort -@ 10 {}.cut.uniq.bam -o {}.cut.uniq.sort.bam'.format(self.sam_path,
                                                                                  self.option('sample_id'),
                                                                                  self.option('sample_id'))
        self.logger.info(sam_sort)
        cmd = self.add_command("sam_sort", sam_sort).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam排序成功")
        else:
            raise Exception("sam排序出错")

        picard_cmd = '/program/sun_jdk1.8.0/bin/java -Xmx10g -Djava.io.tmpdir={} -jar {} MarkDuplicates ' \
                     'VALIDATION_STRINGENCY=LENIENT INPUT={}.cut.uniq.sort.bam OUTPUT={}.cut.uniq.sort.md.bam' \
                     ' METRICS_FILE={}.cut.uniq.sort.md.metrics'\
            .format(self.work_dir, self.picard_path, self.option('sample_id'), self.option('sample_id'),
                    self.option('sample_id'))
        self.logger.info(picard_cmd)
        cmd = self.add_command("picard_cmd", picard_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("picard成功")
        else:
            raise Exception("picard出错")

        sam_valid = '{} view -F 1024 -@ 10 -bS {}.cut.uniq.sort.md.bam -o {}.valid.bam'\
            .format(self.sam_path, self.option('sample_id'), self.option('sample_id'))
        self.logger.info(sam_valid)
        cmd = self.add_command("sam_valid", sam_valid).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam_valid成功")
        else:
            raise Exception("sam_valid出错")

        sam_valid_index = '{} index {}.valid.bam'.format(self.sam_path, self.option('sample_id'))
        self.logger.info(sam_valid_index)
        cmd = self.add_command("sam_valid_index", sam_valid_index).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam_valid_index成功")
        else:
            raise Exception("sam_valid_index出错")

        sam_map = '{} view -bF 4 -@ 10 {}.valid.bam -o {}.map.valid.bam'\
            .format(self.sam_path, self.option('sample_id'), self.option('sample_id'))
        self.logger.info(sam_map)
        cmd = self.add_command("sam_map", sam_map).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam_map成功")
        else:
            raise Exception("sam_map出错")

        sam_map_index = '{} index {}.map.valid.bam'.format(self.sam_path, self.option('sample_id'))
        self.logger.info(sam_map_index)
        cmd = self.add_command("re_sam_map_index", sam_map_index).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("sam_map_index成功")
        else:
            raise Exception("sam_map_index出错") 

        bed_qc = '{}nipt_bed.sh {} {} {} {}'\
            .format(self.script_path, self.option('sample_id'), self.bed_ref, self.work_dir,
                    self.option("fastq_path").prop['path'])
        self.logger.info(bed_qc)
        cmd = self.add_command("bed_qc", bed_qc).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("bed文件生成成功")
        else:
            raise Exception("bed文件生成出错")

    def run_continually(self):
        cmd = '{}nipt_script.sh {} {} {} {} {} {} {} {}'\
            .format(self.script_path, self.option('sample_id'), self.work_dir,
                    self.java_path, self.picard_path, self.ref, self.ref1, self.bed_ref,
                    self.option("fastq_path").prop['path'])
        self.logger.info(cmd)
        cmd = self.add_command("total_cmd", cmd).run()
        self.wait(cmd)

        if cmd.return_code == 0:
            self.logger.info("shell部分运行成功")
        else:
            raise Exception("shell部分运行出错")

    def run_single(self):
        cmd = '{}nipt_script_single.sh {} {} {} {} {} {} {} {}' \
            .format(self.script_path, self.option('sample_id'), self.work_dir,
                    self.java_path, self.picard_path, self.ref, self.ref1, self.bed_ref,
                    self.option("fastq_path").prop['path'])
        self.logger.info(cmd)
        cmd = self.add_command("single_cmd", cmd).run()
        self.wait(cmd)

        if cmd.return_code == 0:
            self.logger.info("single-shell部分运行成功")
        else:
            raise Exception("single-shell部分运行出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir)
        for f in results:
            if re.search(r'.*map.valid.bam$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            elif re.search(r'.*qc$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            elif re.search(r'.*map.valid.sam$', f):
                self.logger.info("设置新的GC值")
                gc = self.get_map_gc(self.work_dir + '/' + f)
                with open(self.output_dir + '/new_gc.txt', 'w+') as w:
                    w.write('GC\t' + gc)
                self.logger.info("设置新的GC值成功！")
            elif re.search(r'.*bed.2$', f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                self.logger.info("开始移动bed文件到database库中")
                if os.path.exists(self.bed_dir + f):
                    os.remove(self.bed_dir + f)
                os.link(self.work_dir + '/' + f, self.bed_dir + f)
                self.logger.info("移动bed文件到database库中成功")
        self.logger.info('设置文件夹路径成功')
        if self.api.api('medical.nipt_analysis_v2').check_exist_bed(self.option('sample_id')):
            self.logger.info("样本{}已经在库中存在了，不再进行导表！".format(self.option('sample_id')))
        else:
            self.logger.info(self.output_dir + '/' + self.option('sample_id') + '.bed.2')
            self.api.api('medical.nipt_analysis_v2').add_bed_file(self.output_dir + '/' + self.option('sample_id') +
                                                                  '.bed.2')

    def get_map_gc(self, sam):
        """
        根据Sam文件计算对应的GC值
        :param sam:
        :return:
        """
        list1 = []
        with open(sam, 'rb') as f:
            for line in f:
                l = line.strip().split('\t')
                seq = l[9]
                list1.append(seq)
        a = ''.join(list1)
        length = len(a)
        n = 1
        g = c = 0
        while n < length:
            b = a[n - 1:n]
            if b == 'G':
                g += 1
            elif b == 'C':
                c += 1
            else:
                pass
            n += 1
        gc = g + c
        percent_gc = '{:.3f}'.format(float(gc) / length)
        self.logger.info(percent_gc, gc, length)
        return str(percent_gc)

    def run(self):
        super(FastqProcessTool, self).run()
        if self.option('run_whole') == 'true' and self.option('single') == 'false':
            self.run_continually()
        elif self.option('run_whole') == 'true' and self.option('single') == 'true':
            self.run_single()
        else:
            self.run_tf()
        self.set_output()
        self.end()
