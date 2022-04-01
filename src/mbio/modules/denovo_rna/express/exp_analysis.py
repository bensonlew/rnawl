# -*- coding: utf-8 -*-
# __author__ = 'qiuping'

from biocluster.module import Module
import os
import re
from biocluster.core.exceptions import OptionError


class ExpAnalysisModule(Module):
    def __init__(self, work_id):
        super(ExpAnalysisModule, self).__init__(work_id)
        self.step.add_steps('rsem', 'merge_rsem', 'diff_exp', "gene_corr", "tran_corr", "tran_pca", "gene_pca")
        options = [
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "rsem_fa", "type": "infile", "format": "sequence.fasta"},  # trinit.fasta文件
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq_dir"},  # PE测序，包含所有样本的左端fq文件的文件夹
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq_dir"},  # PE测序，包含所有样本的左端fq文件的文件夹
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq_dir"},  # SE测序，包含所有样本的fq文件的文件夹
            {"name": "gname", "type": "string"},  # 分组方案名称
            {"name": "diff_rate", "type": "float", "default": 0.01},  # 期望的差异基因比率
            {"name": "only_bowtie_build", "type": "bool", "default": False},  # 为true时该tool只建索引
            {"name": "exp_way", "type": "string", "default": "fpkm"},
            {"name": "dispersion", "type": "float", "default": 0.1},  # edger离散值
            {"name": "min_rowsum_counts", "type": "int", "default": 2},  # 离散值估计检验的最小计数值
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "denovo_rna.express.control_table"},  # 对照组文件，格式同分组文件
            {"name": "diff_ci", "type": "float", "default": 0.01},  # 显著性水平
            {"name": "diff_count", "type": "outfile", "format": "denovo_rna.express.express_matrix"},  # 差异基因计数表
            {"name": "diff_fpkm", "type": "outfile", "format": "denovo_rna.express.express_matrix"},  # 差异基因表达量表
            {"name": "diff_list_dir", "type": "outfile", "format": "denovo_rna.express.gene_list_dir"},
            {"name": "all_list", "type": "outfile", "format": "denovo_rna.express.gene_list"},  # 全部基因名称文件
            {"name": "bam_dir", "type": "outfile", "format": "align.bwa.bam_dir"},  # bowtie2的bam格式的比对文件
            {"name": "tran_count", "type": "outfile", "format": "denovo_rna.express.express_matrix"},
            {"name": "gene_count", "type": "outfile", "format": "denovo_rna.express.express_matrix"},
            {"name": "tran_fpkm", "type": "outfile", "format": "denovo_rna.express.express_matrix"},
            {"name": "gene_fpkm", "type": "outfile", "format": "denovo_rna.express.express_matrix"},
            {"name": "diff_list", "type": "outfile", "format": "denovo_rna.express.gene_list"},  # 差异基因名称文件
        ]
        self.add_option(options)
        self.bowtie_build = self.add_tool("denovo_rna.express.rsem")
        self.merge_rsem = self.add_tool("denovo_rna.express.merge_rsem")
        self.diff_exp = self.add_tool("denovo_rna.express.diff_exp")
        self.rsem_tools = []
        self.diff_gene = False
        self.bam_path = self.work_dir + '/bowtie2_bam_dir/'
        # add by qindanhua 161201
        self.sample_num = 0
        self.tool_lists = [self.diff_exp]
        self.gene_corr = self.add_tool("denovo_rna.mapping.correlation")
        self.tran_corr = self.add_tool("denovo_rna.mapping.correlation")
        self.gene_pca = self.add_tool("meta.beta_diversity.pca")
        self.tran_pca = self.add_tool("meta.beta_diversity.pca")
        self.corr_tool_list = [self.gene_corr, self.tran_corr, self.gene_pca, self.tran_pca]
        # self.tool_lists += self.corr_tool_list

    def check_options(self):
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE')
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内')
        if not self.option('only_bowtie_build'):
            if self.option("fq_type") == "PE":
                if not self.option("fq_r").is_set and not self.option("fq_l").is_set:
                    raise OptionError("PE测序时需设置左端序列和右端序列输入文件")
                else:
                    self.sample_num = len(os.listdir(self.option("fq_r").prop['path']))
            if self.option("fq_type") == "SE":
                if not self.option("fq_s").is_set:
                    raise OptionError("SE测序时需设置序列输入文件")
                else:
                    self.sample_num = len(os.listdir(self.option("fq_s").prop['path']))
        if self.option("exp_way") not in ['fpkm', 'tpm']:
            raise OptionError("所设表达量的代表指标不在范围内，请检查")
        if not self.option('control_file').is_set:
            raise OptionError("必须设置输入文件：上下调对照组参考文件")
        if self.option("diff_ci") >= 1 or self.option("diff_ci") <= 0:
            raise OptionError("显著性水平不在(0,1)范围内")
        if self.option("diff_rate") > 1 or self.option("diff_rate") <= 0:
            raise OptionError("期望的差异基因比率不在(0，1]范围内")
        if self.option('group_table').is_set and not self.option('gname'):
            raise OptionError("有分组文件时必须传入分组方案名字")
        if self.option('group_table').is_set and self.option('gname') not in self.option('group_table').prop['group_scheme']:
            raise OptionError("传入分组方案名字不在分组文件内")
        if not isinstance(self.option('only_bowtie_build'), bool):
            raise OptionError('only_bowtie_build只能为bool')
        return True

    # add by qindanhua 161201
    def correlation_run(self):
        self.gene_corr.set_options({
            'fpkm': self.merge_rsem.option('gene_fpkm').prop['path']
        })
        self.tran_corr.set_options({
            'fpkm': self.merge_rsem.option('tran_fpkm').prop['path']
        })
        self.gene_corr.on('end', self.set_step, {'end': self.step.gene_corr, 'start': self.step.gene_corr})
        self.tran_corr.on('end', self.set_step, {'end': self.step.tran_corr, 'start': self.step.tran_corr})
        self.gene_corr.on('end', self.set_output, "gene_correlation")
        self.tran_corr.on('end', self.set_output, "tran_correlation")
        if self.sample_num > 2:
            self.tran_pca.set_options({
                'otutable': self.merge_rsem.option('tran_fpkm').prop['path']
            })
            self.gene_pca.set_options({
                'otutable': self.merge_rsem.option('gene_fpkm').prop['path']
            })
            self.tran_pca.on('end', self.set_step, {'end': self.step.tran_pca, 'start': self.step.tran_pca})
            self.gene_pca.on('end', self.set_step, {'end': self.step.gene_pca, 'start': self.step.gene_pca})
            self.tran_pca.on('end', self.set_output, "tran_correlation")
            self.gene_pca.on('end', self.set_output, "gene_correlation")
        if self.sample_num > 2:
            self.gene_pca.run()
            self.tran_pca.run()
            self.gene_corr.run()
            self.tran_corr.run()
        if self.sample_num == 2:
            self.gene_corr.run()
            self.tran_corr.run()

    def run_bowtie_build(self):
        tool_opt = {
            'fq_type': self.option('fq_type'),
            'rsem_fa': self.option('rsem_fa'),
            'only_bowtie_build': True
        }
        self.bowtie_build.set_options(tool_opt)
        self.bowtie_build.run()

    def rsem_run(self):
        self.step.rsem.start()
        self.step.update()
        if not os.path.exists(self.bam_path):
            os.mkdir(self.bam_path)
        tool_opt = {
            'fq_type': self.option('fq_type'),
            'only_bowtie_build': self.option('only_bowtie_build'),
            'rsem_fa': self.bowtie_build.option('fa_build')
        }
        if self.option('fq_type') == 'SE':
            s_files = os.listdir(self.option('fq_s').prop['path'])
            for f in s_files:
                if re.search(r'fastq$', f):
                    sample = f.split('_sickle_s.fastq')[0]
                    tool_opt.update({'fq_s': self.option('fq_s').prop['path'] + '/' + f})
                    self.rsem = self.add_tool('denovo_rna.express.rsem')
                    print tool_opt
                    # print self.bowtie_build.option('fa_build').prop['path']
                    self.rsem.set_options(tool_opt)
                    if self.sample_num != 1:
                        self.rsem.run()
                    self.rsem_tools.append(self.rsem)
        else:
            r_files = os.listdir(self.option('fq_r').prop['path'])
            l_files = os.listdir(self.option('fq_l').prop['path'])
            for f in r_files:
                if re.search(r'fastq$', f):
                    tool_opt['fq_r'] = self.option('fq_r').prop['path'] + '/' + f
                    sample = f.split('_sickle_r.fastq')[0]
                    for f1 in l_files:
                        if sample in f1:
                            tool_opt['fq_l'] = self.option('fq_l').prop['path'] + '/' + f1
                            self.rsem = self.add_tool('denovo_rna.express.rsem')
                            self.rsem.set_options(tool_opt)
                            if self.sample_num != 1:
                                self.rsem.run()
                            self.rsem_tools.append(self.rsem)
        print self.rsem_tools
        if self.sample_num == 1:
            self.rsem.on("end", self.set_output, 'rsem')
            self.rsem.on('end', self.end)
            self.rsem.run()
        else:
            self.on_rely(self.rsem_tools, self.set_output, 'rsem')
        # self.on_rely(self.rsem_tools, self.merge_rsem_run)
        # self.on_rely(self.rsem_tools, self.set_step, {'end': self.step.rsem, 'start': self.step.merge_rsem})

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def merge_rsem_run(self):
        self.merge_rsem.set_options({
            'rsem_files': self.output_dir + '/rsem',
            'exp_way': self.option('exp_way')
        })
        self.merge_rsem.on('end', self.set_output, 'merge_rsem')
        self.merge_rsem.on('end', self.set_step, {'end': self.step.merge_rsem, 'start': self.step.diff_exp})
        self.merge_rsem.run()

    def diff_exp_run(self):
        tool_opt = {
            'count': self.merge_rsem.option('gene_count'),
            'fpkm': self.merge_rsem.option('gene_fpkm'),
            'dispersion': self.option('dispersion'),
            'min_rowsum_counts': self.option('min_rowsum_counts'),
            'control_file': self.option('control_file'),
            'diff_ci': self.option('diff_ci'),
            'diff_rate': self.option('diff_rate')
        }
        if self.option('group_table').is_set:
            tool_opt['edger_group'] = self.option('group_table')
            tool_opt['gname'] = self.option('gname')
        self.diff_exp.set_options(tool_opt)
        self.diff_exp.on('end', self.set_output, 'diff_exp')
        self.diff_exp.on('end', self.set_step, {'end': self.step.diff_exp})
        self.diff_exp.run()

    def linkdir(self, dirpath, dirname, output_dir):
        """
        link一个文件夹下的所有文件到本module的指定目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'rsem':
            for tool in self.rsem_tools:
                self.linkdir(tool.output_dir, 'rsem', self.output_dir)
                files = os.listdir(tool.work_dir)
                for f in files:
                    if re.search(r'bowtie2\.bam$', f):
                        if os.path.exists(self.bam_path + f):
                            os.remove(self.bam_path + f)
                            os.link(os.path.join(tool.work_dir, f), self.bam_path + f)
                        else:
                            os.link(os.path.join(tool.work_dir, f), self.bam_path + f)
            self.option('bam_dir', self.bam_path)
            if self.sample_num != 1:
                self.set_step({'data': {'end': self.step.rsem, 'start': self.step.merge_rsem}})
                self.merge_rsem.on('end', self.diff_exp_run)
                self.merge_rsem_run()
        elif event['data'] == 'merge_rsem':
            self.linkdir(obj.output_dir, 'rsem', self.output_dir)
            self.option('gene_count', self.merge_rsem.option('gene_count'))
            self.option('gene_fpkm', self.merge_rsem.option('gene_fpkm'))
            self.option('gene_fpkm').get_list(self.work_dir + '/all_list')
            self.option('all_list', self.work_dir + '/all_list')
            self.option('tran_count', self.merge_rsem.option('tran_count'))
            self.option('tran_fpkm', self.merge_rsem.option('tran_fpkm'))
        elif event['data'] == 'diff_exp':
            self.linkdir(obj.output_dir, 'diff_exp', self.output_dir)
            self.diff_gene = obj.diff_gene
            if self.diff_gene:
                self.option('diff_count', obj.option('diff_count'))
                self.option('diff_fpkm', obj.option('diff_fpkm'))
                self.option('diff_list', obj.option('diff_list'))
                self.option('diff_list_dir', obj.option('diff_list_dir'))
            else:
                self.logger.info('此输入文件没有检测到差异基因')
        # add by qindanhua 161205
        elif event['data'] == 'tran_correlation':
            self.linkdir(obj.output_dir, "tran_correlation", self.output_dir)
        elif event['data'] == 'gene_correlation':
            self.linkdir(obj.output_dir, "gene_correlation", self.output_dir)
        else:
            pass

    def run(self):
        super(ExpAnalysisModule, self).run()
        self.bowtie_build.on('end', self.rsem_run)
        self.run_bowtie_build()
        # self.merge_rsem.on('end', self.correlation_run)
        # modify by qindanhua 改变end依赖对象 20170213
        # if self.sample_num < 2:
        #     self.rsem.on('end', self.end)
        if self.sample_num == 2:
            self.tool_lists += [self.gene_corr, self.tran_corr]
            self.merge_rsem.on('end', self.correlation_run)
            self.on_rely(self.tool_lists, self.end)
        else:
            self.tool_lists += self.corr_tool_list
            self.merge_rsem.on('end', self.correlation_run)
            self.on_rely(self.tool_lists, self.end)

    def end(self):
        repaths = [
            [".", "", "表达量分析模块结果输出目录"],
            ["./rsem", "", "rsem分析结果输出目录"],
            ["./diff_exp", "", "edger分析结果输出目录"],
            ["diff_exp/diff_fpkm", "xls", "差异基因表达量表"],
            ["diff_exp/diff_count", "xls", "差异基因计数表"]
        ]
        regexps = [
            [r"rsem/results$", "xls", "rsem结果"],
            [r"rsem/matrix$", "xls", "表达量矩阵"],
            [r"diff_exp/.*_edgr_stat\.xls$", "xls", "edger统计结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(ExpAnalysisModule, self).end()
