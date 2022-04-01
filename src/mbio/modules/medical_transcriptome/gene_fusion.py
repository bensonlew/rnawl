# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200609

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import gevent.subprocess as subprocess
from collections import OrderedDict
from biocluster.config import Config


class GeneFusionModule(Module):
    """
    该Module用于基因融合分析，默认使用方法
    """
    def __init__(self, work_id):
        super(GeneFusionModule, self).__init__(work_id)
        options = [
            {"name": "ref_genome_custom", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "bam_list", "type": "string"},
            {"name": "min_junction_reads", "type": "int", "default": 1},
            # 最小junction_reads数  minimum fusion support = ( # junction_reads + # spanning_frags ) Default: 2   最小融合支持：junction+spanning_frags
            {"name": "min_sum_frags", "type": "int", "default": 2},
            # (minimum of junction reads required if breakpoint  lacks involvement of only reference junctions) 当breakpoint  不在reference junction时，junction的最小数值要求
            {"name": "min_novel_junction_support", "type": "int", "default": 3},
            # minimum number of rna-seq fragments required as fusion evidence if there are no junction reads (default: 5) 如果没有junction_reads只有spanning_frags需要多少的支持才进行考虑
            {"name": "min_spanning_frags_only", "type": "int", "default": 5},
            # minimum FFPM (fusion fragments per million rna-seq frags)  (default: 0.1)
            {"name": "min_FFPM", "type": "float", "default": 0.1},
            {"name": "ref_gtf", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "min_FFPM", "type": "float", "default": 0.1},
            {"name": "task_id", "type": "string"},
            {"name": "circos", "type": "string", "default": None},
            {"name": "id_modify", "type": "string", "default": None},
            {"name": "species", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.sample_info = OrderedDict()
        self.star_fusion_modules=[]
        self.sample_stat = OrderedDict()
        self.ctat_genome_lib_build_dir = ""
        if self.option("species") == "Homo_sapiens":
            self.dfam_path = Config().SOFTWARE_DIR + "/database/DFAM_3.1/human/homo_sapiens_dfam.hmm"
        elif self.option("species") == "Mus_musculus":
            self.dfam_path = Config().SOFTWARE_DIR + "/database/DFAM_3.1/mouse/mus_musculus_dfam.hmm"
        else:
            self.dfam_path =  Config().SOFTWARE_DIR + "/database/DFAM_3.1/common/Dfam.hmm"
        self.pfam_path = Config().SOFTWARE_DIR + "/database/Annotation/other2019/pfam32/Pfam-A.hmm"
        self.make_lib = None



    def check_options(self):
        if not self.option("ref_genome_custom").is_set:
            raise OptionError("缺少ref_genome_custom参数")
        if not self.option("bam_list"):
            raise OptionError("缺少bam_list参数")
        if not self.option("ref_gtf"):
            raise OptionError("缺少ref_gtf参数")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_lib(self):
        genome_path = self.option("ref_genome_custom").prop["path"]
        if os.path.exists(os.path.join(os.path.dirname(genome_path),"star_27")):
            self.logger.info("{}已建立对应的star2.7库并构建融合基因数据库,直接进行分析".format(self.option("species")))
            self.ctat_genome_lib_build_dir = os.path.join(os.path.dirname(genome_path),"star_27","ctat_genome_lib_build_dir")
            self.logger.info("数据库位置是{}".format(self.ctat_genome_lib_build_dir))
            self.star_fusion_run()
        else:
            self.ctat_genome_lib_build_dir = os.path.join(os.path.dirname(genome_path), "star_27",
                                                          "ctat_genome_lib_build_dir")
            self.make_lib=self.add_tool("medical_transcriptome.gene_fusion.make_lib")
            self.make_lib.set_options({
                "ref_fasta":self.option("ref_genome_custom").prop["path"],
                "gtf_path" :self.option("ref_gtf").prop["path"],
                "dfam_path" : self.dfam_path ,
                "pfam_path": self.pfam_path
            })
            self.make_lib.on('end', self.star_fusion_run)
            self.make_lib.run()


    def star_fusion_run(self):
        sample_info = self.get_sample_bam()
        n = 0
        for key in sample_info.keys():
            star_fusion = self.add_module("medical_transcriptome.star_fusion")
            opts={
                "ref_genome_custom": self.option("ref_genome_custom").prop["path"],
                "ref_gtf":self.option("ref_gtf").prop["path"],
                "in_bam": sample_info[key],
                "min_junction_reads":self.option("min_junction_reads"),
                "min_sum_frags": self.option("min_sum_frags"),
                "min_novel_junction_support": self.option("min_novel_junction_support"),
                "min_spanning_frags_only": self.option("min_spanning_frags_only"),
                "min_FFPM": self.option("min_FFPM"),
                "task_id" :self.option("task_id"),
                "sample_name":key
            }
            if self.option("circos"):
                opts.update({"circos":self.option("circos")})
            if self.option("id_modify"):
                opts.update({"id_modify": self.option("id_modify")})
            star_fusion.set_options(opts)
            self.star_fusion_modules.append(star_fusion)
            n += 1
        for j in range(len(self.star_fusion_modules)):
            self.star_fusion_modules[j].on('end', self.set_output, 'star_fusion')
        if self.star_fusion_modules:
            if len(self.star_fusion_modules) > 1:
                self.on_rely(self.star_fusion_modules, self.run_fusion_result_stat)
            elif len(self.star_fusion_modules) == 1:
                self.star_fusion_modules[0].on('end', self.run_fusion_result_stat)
        else:
            self.set_error("star_fusion_modules列表为空！")
        for tool in self.star_fusion_modules:
            gevent.sleep(1)
            tool.run()


    def run_fusion_result_stat(self):
        sample_results = os.listdir(os.path.join(self.output_dir,"star_fusion"))
        for i in self.sample_info:
            sample_name = i
            sample_stat_file = os.path.join(self.output_dir,"star_fusion",i,"star-fusion.fusion_predictions.abridged.tsv")
            sample_fusion_num = len(open(sample_stat_file, 'r').readlines()) #这个写法其实极其不科学，虽然很简单，主要这里文件一般比较小所以可以这么写
            self.sample_stat[sample_name] = sample_fusion_num -1
        with open(os.path.join(self.output_dir,"fusion_stat.txt"),"w") as w:
            w.write("Sample"+"\t"+"Number"+"\n")
            for i in self.sample_stat:
                w.write(i + "\t" + str(self.sample_stat[i]) + "\n")
        self.end()


    def get_sample_bam(self):
        """
        生成样本与bam文件对应的字典文件,后面以bam为单位输入star_fusion中进行分析
        :return:
        """
        sample_info = OrderedDict()
        with open(self.option("bam_list"), "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                if len(temp) == 3:
                    if temp[1].upper() not in ['WGS', 'RAD', 'GBS', 'WES']:
                        self.logger.warn("{}建库类型不合法!".format(temp[1]))
                    sample_info[temp[0]] = temp[2]
                elif len(temp) == 2:
                    sample_info[temp[0]] = temp[1]
                else:
                    sample_info[line.split("/")[-1].split(".bam")[0]] =line.strip()
                    #self.set_error("%sbam_list格式不正确！", variables=(self.option("bam_list")), code="24500418")
        self.logger.info("样本与bam文件对应关系：{}".format(sample_info))
        self.sample_info =  sample_info
        return sample_info



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
        if event['data'] == 'star_fusion':
            self.linkdir(obj.output_dir, 'star_fusion')
        else:
            pass

    def run(self):
        super(GeneFusionModule, self).run()
        self.logger.info("开始运行gene_fusion")
        self.logger.info("首先获取样本信息")
        self.check_lib()
        # self.star_fusion_run()

    def end(self):
        super(GeneFusionModule, self).end()



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
            "id": "gene_fusion" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "medical_transcriptome.gene_fusion",
            "instant": False,
            "options": dict(
                bam_list="/mnt/ilustre/users/sanger-dev/workspace/20200810/Refrna_tsg_38312/RnaseqMapping/output/bamlist",
                ref_genome_custom="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38.p6_ensembl_100_/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa",
                ref_gtf ="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38.p6_ensembl_100_/gtf/Mus_musculus.GRCm38.100.gtf",
                circos ="yes",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()