# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180414

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os


class AnnovarModule(Module):
    """
    用于对snp indel vcf文件进行注释统计
    """
    def __init__(self, work_id):
        super(AnnovarModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "snp_anno_primary_vcf", "type": "string"},  # snp.anno.primary.vcf
            {"name": "indel_anno_primary_vcf", "type": "string"},  # indel.anno.primary.vcf
            {"name": "snp_anno_genes", "type": 'string'},    # snp.anno.genes.txt
            {"name": "indel_anno_genes", "type": 'string'},  # indel.anno.genes.txt
            {"name": "anno_summary", "type": "string"}  # 该文件每个基因组版本对应一个该文件，用于生成pop summary文件
        ]
        self.add_option(options)
        self.gatk_combine_variants = self.add_tool("wgs.gatk_combine_variants")
        self.anno_count = self.add_tool("wgs.anno_count")
        self.eggnog_anno = self.add_tool("wgs.anno_analysis")
        self.kegg_anno = self.add_tool("wgs.anno_analysis")
        self.go_anno = self.add_tool("wgs.anno_analysis")

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="24500101")
        if not self.option("snp_anno_primary_vcf"):
            raise OptionError("缺少snp_anno_primary_vcf参数", code="24500102")
        if not self.option("indel_anno_primary_vcf"):
            raise OptionError("缺少indel_anno_primary_vcf参数", code="24500103")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        ref_file_name = os.path.basename(self.option("ref_fasta").prop['path']).split(".")[0]
        if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
                or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
            raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24500104")
        return True

    def gatk_combine_variants_run(self):
        self.gatk_combine_variants.set_options({
            "ref_fasta": self.option("ref_fasta").prop['path'],
            "snp_anno_primary_vcf": self.option("snp_anno_primary_vcf"),
            "indel_anno_primary_vcf": self.option("indel_anno_primary_vcf")
        })
        self.gatk_combine_variants.on('end', self.set_output, 'combine_variants')
        self.gatk_combine_variants.on('end', self.anno_count_run)
        self.gatk_combine_variants.run()

    def anno_count_run(self):
        self.anno_count.set_options({
            "snp_anno_genes": self.option("snp_anno_genes"),
            "indel_anno_genes": self.option("indel_anno_genes"),
            "anno_summary": self.option("anno_summary")
        })
        self.anno_count.on("end", self.set_output, 'anno_count')
        self.anno_count.on("end", self.run_go_anno)
        self.anno_count.on("end", self.run_kegg_anno)
        self.anno_count.on("end", self.run_eggnog_anno)
        self.anno_count.run()

    def run_go_anno(self):
        """
        对定位后的区域进行基因注释分析 go anno
        :return:
        """
        self.go_anno.set_options({
            "anno_stat": self.anno_count.output_dir + "/pop.go.stat",
            "anno_type": "go"
        })
        self.go_anno.on("end", self.set_output, "go_anno")
        self.go_anno.run()

    def run_kegg_anno(self):
        """
        对定位后的区域进行基因注释分析 kegg anno
        :return:
        """
        self.kegg_anno.set_options({
            "anno_stat": self.anno_count.output_dir + "/pop.kegg.stat",
            "anno_type": "kegg"
        })
        self.kegg_anno.on("end", self.set_output, "kegg_anno")
        self.kegg_anno.run()

    def run_eggnog_anno(self):
        """
        对定位后的区域进行基因注释分析 eggnog anno
        :return:
        """
        self.eggnog_anno.set_options({
            "anno_stat": self.anno_count.output_dir + "/pop.eggnog.stat",
            "anno_type": "eggnog"
        })
        self.eggnog_anno.on("end", self.set_output, "eggnog_anno")
        self.eggnog_anno.run()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
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
                # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'combine_variants':
            self.linkdir(obj.output_dir, 'combine_variants')
        elif event['data'] == 'anno_count':
            self.linkdir(obj.output_dir, 'anno_count')
        elif event['data'] == 'go_anno':
            self.linkdir(obj.output_dir, 'go_anno')
        elif event['data'] == 'eggnog_anno':
            self.linkdir(obj.output_dir, 'eggnog_anno')
        elif event['data'] == 'kegg_anno':
            self.linkdir(obj.output_dir, 'kegg_anno')
        else:
            pass

    def run(self):
        super(AnnovarModule, self).run()
        self.on_rely([self.eggnog_anno, self.kegg_anno, self.go_anno], self.end)
        self.gatk_combine_variants_run()

    def end(self):
        super(AnnovarModule, self).end()
