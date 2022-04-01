# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180416

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import unittest


class VcfFilterModule(Module):
    """
    用于对pop.variant.vcf进行简并碱基过滤，然后分别对snp与indel进行筛选，过滤与eff分析,最后就是统计分析
    """
    def __init__(self, work_id):
        super(VcfFilterModule, self).__init__(work_id)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "pop_var_vcf", "type": "string"},  # pop_variant_vcf
            {"name": "snpEff_config", "type": "string"}  # # snpEff.config
        ]
        self.add_option(options)
        self.filter_bases = self.add_tool("wgs.filter_bases")
        self.gatk_indel = self.add_tool("ref_rna_v2.gatk_indel")
        self.gatk_snp = self.add_tool("ref_rna_v2.gatk_snp")
        self.snp_eff = self.add_tool("wgs.snpeff")
        self.indel_eff = self.add_tool("wgs.snpeff")
        self.variant_stat1 = self.add_tool("wgs.variant_stat")
        self.variant_stat = self.add_tool("wgs.variant_stat")
        self.variant_qual = self.add_tool("wgs.variant_qual")
        self.variant_qual1 = self.add_tool("wgs.variant_qual")
        self.vcf_stat1 = self.add_tool("wgs.vcf_stat")
        self.vcf_stat = self.add_tool("wgs.vcf_stat")

    def check_options(self):
        if not self.option("ref_fasta").is_set:
            raise OptionError("缺少ref_fasta参数", code="23701704")
        if not self.option("pop_var_vcf"):
            raise OptionError("缺少pop_var_vcf参数", code="23701705")
        if not self.option("snpEff_config"):
            raise OptionError("缺少snpEff_config参数", code="23701706")
        ref_file = os.path.dirname(self.option("ref_fasta").prop['path'])
        #ref_file_name = os.path.split(self.option("ref_fasta").prop['path'])[0]
        #if not os.path.isfile(os.path.join(ref_file, "{}.dict".format(ref_file_name))) \
         #       or not os.path.isfile(os.path.join(ref_file, "{}.fa.fai".format(ref_file_name))):
          #  raise OptionError("参考组配置文件缺少dcit与fa.fai文件", code="24501704")
        return True

    def filter_bases_run(self):
        self.filter_bases.set_options({
            "pop_var_vcf": self.option("pop_var_vcf")
        })
        self.filter_bases.on('end', self.set_output, 'vcf_filter')
        self.filter_bases.run()

    def gatk_indel_run(self):
        self.gatk_indel.set_options({
            "ref_fasta": self.option("ref_fasta").prop['path'],
            "pop_var_vcf": self.filter_bases.output_dir + "/pop.variant.vcf"
        })
        # self.gatk_indel.on("end", self.set_output, 'vcf_filter_')
        self.gatk_indel.on("end", self.indel_eff_run)
        self.gatk_indel.run()

    def gatk_snp_run(self):
        self.gatk_snp.set_options({
            "ref_fasta": self.option("ref_fasta").prop['path'],
            "pop_var_vcf": self.filter_bases.output_dir + "/pop.variant.vcf"
        })
        # self.gatk_snp.on("end", self.set_output, 'vcf_filter_')
        self.gatk_snp.on("end", self.snp_eff_run)
        self.gatk_snp.run()

    def snp_eff_run(self):
        self.snp_eff.set_options({
            "ref_name": "ref",
            "snpEff_config": self.option("snpEff_config"),
            "filter_recode_vcf": self.gatk_snp.output_dir + "/pop.snp.filter.recode.vcf",
            "types": "snp"
        })
        self.snp_eff.on("end", self.set_output, "eff")
        self.snp_eff.on("end", self.variant_stat_snp_run)
        self.snp_eff.on("end", self.variant_qual_snp_run)
        self.snp_eff.on("end", self.vcf_stat_snp_run)
        self.snp_eff.run()

    def indel_eff_run(self):
        self.indel_eff.set_options({
            "ref_name": "ref",
            "snpEff_config": self.option("snpEff_config"),
            "filter_recode_vcf": self.gatk_indel.output_dir + "/pop.indel.filter.recode.vcf",
            "types": "indel"
        })
        self.indel_eff.on("end", self.set_output, "eff")
        self.indel_eff.on("end", self.variant_stat_indel_run)
        self.indel_eff.on("end", self.variant_qual_indel_run)
        self.indel_eff.on("end", self.vcf_stat_indel_run)
        self.indel_eff.run()

    def variant_stat_snp_run(self):
        self.variant_stat1.set_options({
            "anno_primary_vcf": os.path.join(self.snp_eff.output_dir, "snp.anno.primary.vcf"),
            "types": "snp"
        })
        self.variant_stat1.on('end', self.set_output, 'variant_stat')
        self.variant_stat1.run()

    def variant_stat_indel_run(self):
        self.variant_stat.set_options({
            "anno_primary_vcf": os.path.join(self.indel_eff.output_dir, "indel.anno.primary.vcf"),
            "types": "indel"
        })
        self.variant_stat.on('end', self.set_output, 'variant_stat')
        self.variant_stat.run()

    def variant_qual_snp_run(self):
        self.variant_qual.set_options({
            "anno_primary_vcf": os.path.join(self.snp_eff.output_dir, "snp.anno.primary.vcf"),
            "types": "snp"
        })
        self.variant_qual.on('end', self.set_output, 'variant_stat')
        self.variant_qual.run()

    def variant_qual_indel_run(self):
        self.variant_qual1.set_options({
            "anno_primary_vcf": os.path.join(self.indel_eff.output_dir, "indel.anno.primary.vcf"),
            "types": "indel"
        })
        self.variant_qual1.on('end', self.set_output, 'variant_stat')
        self.variant_qual1.run()

    def vcf_stat_snp_run(self):
        self.vcf_stat1.set_options({
            "anno_primary_vcf": os.path.join(self.snp_eff.output_dir, "snp.anno.primary.vcf"),
            "types": "snp"
        })
        self.vcf_stat1.on('end', self.set_output, 'anno_stat')
        self.vcf_stat1.run()

    def vcf_stat_indel_run(self):
        self.vcf_stat.set_options({
            "anno_primary_vcf": os.path.join(self.indel_eff.output_dir, "indel.anno.primary.vcf"),
            "types": "indel"
        })
        self.vcf_stat.on('end', self.set_output, 'anno_stat')
        self.vcf_stat.run()

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
        if event['data'] == 'filter_bases':
            self.linkdir(obj.output_dir, 'filter_bases')
        elif event['data'] == 'vcf_filter_':
            pass
            # self.linkdir(obj.output_dir, 'gatk_indel')
        # elif event['data'] == 'gatk_snp':
        #     self.linkdir(obj.output_dir, 'gatk_snp')
        elif event['data'] == "vcf_filter":
            self.linkdir(obj.output_dir, 'vcf_filter')
        elif event['data'] == "variant_stat":
            self.linkdir(obj.output_dir, 'variant_stat')
        elif event['data'] == "anno_stat":
            self.linkdir(obj.output_dir, 'anno_stat')
        elif event['data'] == "eff":
            self.linkdir(obj.output_dir, 'eff')
        else:
            pass

    def run(self):
        """
        module super放到上面，workflow super放到最下面，第一种方法用on_rely绑定但是不行，暂时不知道原因，权哥讲将绑定放到初始化
        :return:
        """
        super(VcfFilterModule, self).run()
        # self.on_rely([self.gatk_indel, self.gatk_snp], self.snp_eff_run)
        # self.on_rely([self.gatk_indel, self.gatk_snp], self.indel_eff_run)
        # # self.on_rely([self.snp_eff, self.indel_eff], self.end)
        # self.on_rely([self.snp_eff, self.indel_eff], self.variant_stat_snp_run)  # self.variant_stat1
        # self.on_rely([self.snp_eff, self.indel_eff], self.variant_stat_indel_run)  # self.variant_stat
        # self.on_rely([self.snp_eff, self.indel_eff], self.variant_qual_snp_run)  # self.variant_qual
        # self.on_rely([self.snp_eff, self.indel_eff], self.variant_qual_indel_run)  # self.variant_qual1
        # self.on_rely([self.snp_eff, self.indel_eff], self.vcf_stat_snp_run)  # self.vcf_stat1
        # self.on_rely([self.snp_eff, self.indel_eff], self.vcf_stat_indel_run)  # self.vcf_stat
        # self.on_rely([self.variant_stat1, self.variant_stat, self.variant_qual, self.variant_qual1, self.vcf_stat1,
        #               self.vcf_stat], self.end)
        # self.on_rely([self.variant_stat1, self.variant_stat], self.end)
        self.filter_bases.on("end", self.gatk_indel_run)
        self.filter_bases.on("end", self.gatk_snp_run)
        self.on_rely([self.variant_stat1, self.variant_stat, self.variant_qual, self.variant_qual1, self.vcf_stat1,
                      self.vcf_stat], self.end)
        self.filter_bases_run()

    def end(self):
        super(VcfFilterModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp_vcf_filter" + str(random.randint(1, 10000))+"xxx",
            "type": "module",
            "name": "ref_rna_v2.vcf_filter",
            "instant": False,
            "options": dict(
                #ref_dict=test_dir + "/" + "Mus_musculus.GRCm38.dna_rm.toplevel.clean.dict",
                #bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/bamlist_new",
                #call_type="sentieon",
                #ref_fasta=test_dir+"/"+"Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                # scm="complete",
                # scd="correlation",
                # corr_method='pearson',
                #output=None,
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna/output_vcf",
                #ref_genome="customer_mode",
                #combine_vcf=False,
                #des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/biomart/Mus_musculus.GRCm38.biomart_gene.txt",
                #des_type="type1"
                #des="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/biomart/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0.biomart",
                #des_type="type3",
                #ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/gtf/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0.gtf",
                #ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Capsicum_annuum/NCBI/dna/GCF_000710875.1_Pepper_Zunla_1_Ref_v1.0_genomic.fna",
                #input_file="/mnt/ilustre/users/sanger-dev/workspace/20190621/Snp_tsg_34423_8856_3754/CallSnpIndel/GvcfTypingV2/output/pop.variant.vcf"
                pop_var_vcf="/mnt/ilustre/users/sanger-dev/workspace/20190620/Snp_tsg_33912_3063_9782/SamRna/BcftoolVcf/output/variant.vcf",
                #ref_name="ref",
                #pop_var_vcf="/mnt/ilustre/users/sanger-dev/workspace/20190618/Snp_tsg_33912_4050_3785/SnpRna/Annovar/Combine_Variants.vcf",
                #pop_var_vcf="/mnt/ilustre/users/sanger-dev/workspace/20190626/Snp_tsg_33912_5731_9012/CallSnpIndel/output/vcf_call/pop.variant.vcf",
                ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                snpEff_config="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/snpEff.config",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
