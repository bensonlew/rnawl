# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 2018.0521
# modified 2018.0524

from mbio.packages.dna_evolution.send_email import SendEmail
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError
from biocluster.core.function import CJsonEncoder, filter_error_info
from biocluster.module import Module
from biocluster.config import Config
import os
import re
import json
import urllib
import urllib2
import random
import hashlib
import time


class GenomeConfigV2Module(Module):
    """
    变异检测-遗传进化流程中，基因组配置的模块
    # /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/modules/wgs/genome_config.py
    # /mnt/ilustre/users/sanger-dev/biocluster/scripts/wgs/genome_config_module_run.py
    lasted modified by hongdong@20180916
    增加大基因组中每条颜色体超过11亿个碱基，也就是一条染色体达到了1G，picard使用不了，这里使用samtools dict去进行处理
    lasted modified by hongdong@20200323
    """
    def __init__(self, work_id):
        super(GenomeConfigV2Module, self).__init__(work_id)
        options = [
            {"name": "reffa", "type": "infile", "format": "sequence.fasta"},       # 更名后ref.fa文件 infile核查
            {"name": "refgff", "type": "string"},      # 更名后ref.gff文件
            {"name": "info", "type": "string"},             # info.log
            {"name": "chromosomelist", "type": "string"},     # 无染色体的时候 可无此文件
            {"name": "update", "type": "string", "default": 0},
            {"name": "need_rename", "type": "string", "default": "true"},  # 是否需要改名字，一般情况是要的，
            # 但是为了兼容以前的数据，这里ref与gff是已经改过名字的
            {"name": "ref_changelog", "type": 'string'},
            {'name': 'anno_summary', 'type': 'string'},
            {'name': 'mail_info', 'type': 'string'},
            {'name': 'target_dir', 'type': 'string'},  # 参考基因组配置需要检查的文件上传到的路径
            {"name": "large_genome", "type": "string", "default": "false"},
            {"name": "client", "type": "string", "default": "client01"}
        ]
        self.add_option(options)
        self.species_path = ''
        self.species = ''
        self.edition = ''
        self.release_time = ''
        self.strain_characteristic = ''
        self.n50 = ''
        self.link = ''
        self.ref_path = ''
        self.abspath = ''
        self.assembly = ''
        self.web_dir = ''
        self.wgsjson = {}
        self.time_edition = ''
        self.genome_size = ''
        self.wgsjson = ''
        self.source = ''
        self.file_path = ''
        self.dict_path = ''
        self.gtf_path = ''
        self.gff_path = ''
        self.change = 0
        self.is_end_by_module = True
        self.gffread_to_fagtf = self.add_tool("wgs.gffread_to_fagtf")
        self.normalize_fa = self.add_tool("wgs.normalize_fa")
        self.bwa_config = self.add_tool("wgs.bwa_config")
        self.samtools_faidx = self.add_tool("wgs.samtools_faidx")
        self.genome_grename = self.add_tool("wgs.genome_grename")
        self.picard_dict = self.add_tool("wgs.picard_dict")
        self.circos_chrlist = self.add_tool("wgs.circos_chrlist")
        # self.snpeff_gffcheck = self.add_tool("wgs.snpeff_gffcheck")
        # self.snpeff_index = self.add_tool("wgs.snpeff_index")
        self.snpeff_index_v2 = self.add_tool("wgs.snpeff_index_v2")
        self.ssr_ref_primer = self.add_module("wgs.ssr_ref_primer")
        self.makeblastdb = self.add_tool("wgs.makeblastdb")
        self.getgenefasta = self.add_tool("wgs.getgenefasta")
        self.gene_anno = self.add_module("wgs.gene_anno_v2")
        self.circos_ref2bit = self.add_tool("wgs.circos_ref2bit")
        self.geneanno_sort = self.add_tool("wgs.geneanno_sort")
        self.gff_check = self.add_tool("wgs.gff_check")
        self.need_anno = True
        self.task_id = ''
        self._response_code, self._response = '', ''
        self._key, self._client = '', ''
        self._binds_id,  self._interface_id, self._env_name = '', '', ''

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置reffa")
        if not self.option("refgff"):
            raise OptionError("请设置refgff")
        if not self.option("info"):
            raise OptionError("请设置info文件")
        if self.option("need_rename") == 'false' and not self.option('ref_changelog'):
            raise OptionError("不需要改名的时候，一定要输入ref_changelog参数")
        if self.option("client") not in ["client01", "client03"]:
            raise OptionError("client:%s只能是client01/client03" % self.option("client"))
        self._client = self.option("client")

    def makedir_info(self):
        """
        info.log:
        cuiArabidopsis_thaliana Columbia0(TAIR10)   https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.3/  NCBI
        GCF_000001735.3 2015.04.16  116M(genome size)   31,248,787(n50) Chromosome(assembly)  web_dir
        self.file_path = Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/cui_test_dir"   # 测试目录 拷贝json文件
        链接输入文件到inputfile；
        self.release_time||self.edition 唯一ID
        """
        sym1 = "_"
        sym2 = "."
        self.file_path = Config().SOFTWARE_DIR + "/database/dna_geneome/"
        with open(self.option("info"), 'r') as f:
            data = f.readlines()
            for line in data:
                temp = line.strip().split('\t')
                self.logger.info(temp)
                if len(temp) >= 10:    # 核查列数
                    self.strain_characteristic = temp[1].lower()
                    self.link = temp[2]
                    self.genome_size = temp[6].upper()
                    self.n50 = temp[7]
                    self.assembly = temp[8].title()
                    self.web_dir = temp[9]
                    self.species = sym1.join(re.split(' |_|-', temp[0])).capitalize()
                    self.source = sym1.join(re.split(' |_|-', temp[3])).upper()
                    self.edition = temp[4].upper()
                    self.release_time = sym2.join(re.split('/|-', temp[5]))
                    self.species_path = os.path.join(self.file_path, self.species, self.source, self.edition,
                                                     self.release_time)
                    self.logger.info("##### 物种绝对路径" + self.species_path)
                    # if not os.path.exists(self.species_path):
                    #     os.makedirs(self.species_path)
                    #     self.logger.info("##### " + self.species_path + " 该物种路径新建成功")
                else:
                    self.logger.info("@@@@@ 请check info.log文件列数--10列，用\\t分割! ! !--当前log文件列数少于10列，"
                                     "正常结束！")
                    self.set_error("请检查info.log文件列数--10列，用\\t分割! --当前log文件列数少于10列流程结束！")
                    self.is_end_by_module = False
                    self.end()

    def openjson(self):
        """
        链接输入文件到目录下
        当输入文件就是更名后的文件时，不需要链接原始输入文件
        """
        self.abspath = os.path.join(self.species, self.source, self.edition, self.release_time)
        self.time_edition = self.release_time + "||" + self.edition
        self.wgsjson = self.get_json(self.file_path + "/wgs_genome.json")
        self.logger.info("##### 物种是" + self.species)
        self.logger.info("##### 物种版本是" + self.time_edition)
        self.logger.info("##### 已有物种名如下：{}".format(self.wgsjson.keys()))
        self.logger.info(self.option("update"))
        if self.species not in self.wgsjson.keys() or\
                (self.species in self.wgsjson.keys() and self.time_edition not in self.wgsjson[self.species].keys()):
            self.change = 1
        elif int(self.option("update")) == 1:
            self.change = 1
        else:
            self.logger.info(self.species + "\t" + self.time_edition + "该物种id无需导入!")
            self.logger.info("无需导入，正常")
            self.is_end_by_module = False
            self.end()
        # b = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        # if not os.path.exists(self.output_dir + "/database/dna_wgs_geneome/z.wgs_json/"):
        #     os.mkdir(Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/z.wgs_json/")
        # self.json_write(Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/z.wgs_json/wgs_genome.json_" + b)
        # # 将新物种的json数据备份到该路径下的json文件里

    def get_json(self, json_path):
        """
        用于解析出json文件，然后可以从中获取参考基因组数据
        :return:
        """
        f = open(json_path, "r")
        json_dict = json.loads(f.read())
        f.close()
        return json_dict

    def run_gff_check(self):
        options = {'gff_file': self.option("refgff")}
        self.gff_check.set_options(options)
        self.gff_check.run()

    def run_genome_grename(self):
        """
        将下载的fa gff进行更名；
        result: ref.fa ref.gff ref.changelog
        """
        opti = {"fasta": self.option("reffa").prop['path'], "gff": self.gff_check.output_dir + '/ref.gff'}
        if self.option("chromosomelist"):
            opti["chromosomelist"] = self.option("chromosomelist")
        self.genome_grename.set_options(opti)
        self.genome_grename.on("end", self.run_normalize_fa)
        self.genome_grename.run()

    def run_normalize_fa(self):
        """
        均一化处理基因组
        """
        # self.fa_link()
        self.normalize_fa.set_options({
            "ref_fa": self.genome_grename.output_dir + "/ref.fa" if self.option('need_rename') == 'true' else
            self.option("reffa").prop['path'],
            'large_genome': self.option('large_genome')
        })
        self.normalize_fa.on("end", self.run_gffread_to_fagtf)
        # self.normalize_fa.on("end", self.hardlink, "normalize_fa")
        self.normalize_fa.run()
        self.logger.info("\n##### normalize_fa运行\n")

    def base_file_link(self):
        """
        将基础数据，如输入的info.log，ref.changelog 等link到module的output中
        :return:
        """
        self.ref_path = os.path.join(self.output_dir, "ref.fa")
        self.os_link(os.path.join(self.normalize_fa.output_dir, "ref.fa"), self.ref_path)
        if self.option('need_rename') == 'true':
            self.os_link(os.path.join(self.genome_grename.output_dir, 'ref.gff'), self.output_dir + "/ref.gff")
            self.os_link(os.path.join(self.genome_grename.output_dir, 'ref.changelog'),
                         self.output_dir + "/ref.changelog")
        else:
            self.os_link(self.option('ref_changelog'), self.output_dir + "/ref.changelog")
            self.os_link(self.gff_check.output_dir + '/ref.gff', self.output_dir + "/ref.gff")
        self.os_link(self.option("info"), self.output_dir + "/info.log")
        if not self.need_anno:
            self.os_link(self.option('anno_summary'), self.output_dir + "/anno.summary")
        self.gff_path = os.path.join(self.output_dir, "ref.gff")

    def os_link(self, origin, target):
        """
        link 文件
        :param origin:
        :param target:
        :return:
        """
        self.checklink(target)
        os.link(origin, target)
        self.logger.info("link文件{}到{}成功！".format(origin, target))

    def run_gffread_to_fagtf(self):
        """
        用gffread处理基因组成cds.fa pro.fa mrna.fa
        """
        self.base_file_link()
        self.gffread_to_fagtf.set_options({
            "fa": os.path.join(self.normalize_fa.output_dir, "ref.fa"),
            "gff_path": self.gff_path
        })
        self.gffread_to_fagtf.on("end", self.hardlink, "gffread_to_fagtf")
        self.gffread_to_fagtf.on("end", self.run_samtools_faidx)
        self.gffread_to_fagtf.on("end", self.uploadfile)
        self.logger.info("111:{}".format(self.need_anno))
        if self.need_anno:   # 需要注释
            self.gffread_to_fagtf.on("end", self.run_gene_anno)
        self.gffread_to_fagtf.run()
        self.logger.info("\n##### gffread_to_fagtf运行\n")

    def run_samtools_faidx(self):
        # self.gtf_path = os.path.join(self.output_dir, "ref.gtf")
        # if self.checklink(self.gtf_path, False):
        #     os.link(os.path.join(self.gffread_to_fagtf.work_dir, "ref.gtf"), self.gtf_path)
        # 链接gtf结束
        self.samtools_faidx.set_options({
            "pop_fa": self.ref_path
        })
        self.samtools_faidx.on("end", self.run_bwa_config)          # bwa配置
        self.samtools_faidx.on("end", self.run_picard_dict)         # picard配置
        # self.samtools_faidx.on("end", self.run_snpeff_gffcheck)     # snpeff配置
        self.samtools_faidx.on("end", self.run_ssr_ref_primer)      # ssr module配置
        self.samtools_faidx.on("end", self.run_makeblastdb)         # makeblastdb
        self.samtools_faidx.on("end", self.run_getgenefasta)        # 生成ref.gene.fa 改成用ref.new.mRNA.fa去注释了
        self.samtools_faidx.on("end", self.run_ref_2bit)            # ref_2bit给基因组浏览器用
        self.samtools_faidx.on("end", self.run_snpeff_index_v2)
        self.samtools_faidx.run()
        self.logger.info("\n##### samtools_faidx运行\n")

    def run_bwa_config(self):       # 1
        """
        self.species_path + "/ref.fa"确定存在！
        """
        self.bwa_config.set_options({
            "ref_fa": self.ref_path,
        })
        self.bwa_config.run()
        self.logger.info("\n##### bwa_config运行\n")

    def run_picard_dict(self):  # 2.1
        """
        check: ref.dict删掉重生成
        """
        self.dict_path = os.path.join(self.output_dir, "ref.dict")
        # self.checklink(self.dict_path)
        self.picard_dict.set_options({
            "fa": self.ref_path,
            'large_genome': self.option('large_genome')
        })
        self.picard_dict.on("end", self.run_circos_chrlist)
        self.picard_dict.run()
        self.logger.info("\n##### picard_dict运行\n")

    def run_circos_chrlist(self):   # 2.2
        """
        生成ref.chrlist
        """
        self.circos_chrlist.set_options({
            "dict": self.dict_path
        })
        self.circos_chrlist.on("end", self.hardlink, "ref.chrlist")
        self.circos_chrlist.run()
        self.logger.info("\n##### circos_chrlist运行\n")

    # def run_snpeff_gffcheck(self):  # 3.1
    #     """
    #     参考基因组snpeff check ref.gff生成genes.gff
    #     """
    #     self.snpeff_gffcheck.set_options({
    #         "refgff": self.gff_path,
    #     })
    #     self.snpeff_gffcheck.on("end", self.run_snpeff_index_v2)
    #     self.snpeff_gffcheck.on("end", self.hardlink, 'genes.gff')
    #     self.snpeff_gffcheck.run()
    #     self.logger.info("\n##### snpeff_gffcheck运行\n")

    def run_snpeff_index_v2(self):     # 3.2
        """
        生成结果在ref.fa所在目录，不需要链接
        """
        self.snpeff_index_v2.set_options({
            "reffa": self.ref_path,
            # "genesgff": os.path.join(self.snpeff_gffcheck.output_dir, "genes.gff")
            "genesgtf": os.path.join(self.gffread_to_fagtf.work_dir, "ref.gtf")
        })
        self.snpeff_index_v2.run()
        self.logger.info("\n##### snpeff_index_v2运行\n")

    def run_ssr_ref_primer(self):   # 4 ssr的module
        """
        结果：ref.ssr.stat  ssr.ref.result.xls 需要硬链接为ssr.stat ssr.ref.result.xls
        """
        self.ssr_ref_primer.set_options({
            "reffa": self.ref_path
        })
        self.ssr_ref_primer.on("end", self.hardlink, "ref.ssr.stat")
        self.ssr_ref_primer.run()
        self.logger.info("\n##### ssr_ref_primer运行\n")

    def run_makeblastdb(self):   # 5
        # dbpath = os.path.join(self.output_dir, "makedbblast")
        # if not os.path.exists(dbpath):
        #     os.mkdir(dbpath)
        options = {
            "pop_fa": self.ref_path,
            "db_name": 'ref'
        }
        self.makeblastdb.set_options(options)
        self.makeblastdb.on("end", self.hardlink, "makeblastdb")
        self.makeblastdb.run()
        self.logger.info("\n##### makeblastdb运行\n")

    def run_ref_2bit(self):
        """
        生成ref2bit
        """
        options = {
            "reffa": self.ref_path,
        }
        self.circos_ref2bit.set_options(options)
        self.circos_ref2bit.on("end", self.hardlink, "ref_2bit")
        self.circos_ref2bit.run()
        self.logger.info("\n##### ref_2bit运行\n")

    def run_geneanno_sort(self):
        """
        将产品线的anno.summary按照chrid排序后存进来
        # /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/wgs/geneanno_sort.py
        """
        self.geneanno_sort.set_options({
            "anno_summary": self.gene_anno.output_dir + "/anno.summary"
        })
        self.geneanno_sort.on("end", self.hardlink, "geneanno_sort")
        self.geneanno_sort.run()
        self.logger.info("\n##### geneanno_sort运行\n")

    def run_getgenefasta(self):
        self.getgenefasta.set_options({
            "reffa": self.ref_path,
            "refgff": self.gff_path,
        })
        # self.getgenefasta.on("end", self.run_gene_anno)
        self.getgenefasta.run()
        self.logger.info("\n##### getgenefasta运行\n")

    def run_gene_anno(self):
        """
        输入文件：ref.gene.fa;先分20份然后被注释
        """
        self.gene_anno.set_options({
            "fasta": os.path.join(self.gffread_to_fagtf.output_dir, "ref.new.mRNA.fa")
        })
        self.gene_anno.on("end", self.hardlink, "gene_anno")
        self.gene_anno.on("end", self.run_geneanno_sort)
        self.gene_anno.run()
        self.logger.info("\n##### gene_anno运行\n")

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

    def hardlink(self, event):
        obj = event['bind_object']
        if event['data'] == 'ref.chrlist':
            self.linkdir(obj.output_dir, self.output_dir)
        elif event['data'] == 'normalize_fa':
            pass
        elif event['data'] == 'gffread_to_fagtf':
            self.linkdir(obj.output_dir, self.output_dir)
        elif event['data'] == 'genes.gff':
            self.linkdir(obj.output_dir, self.output_dir + "/ref")
        elif event['data'] == 'ref.ssr.stat':
            self.checklink(os.path.join(self.output_dir, "ssr.stat"))
            self.checklink(os.path.join(self.output_dir, "ssr.ref.result.xls"))
            os.link(os.path.join(self.ssr_ref_primer.output_dir, "ref.ssr.stat"),
                    os.path.join(self.output_dir, "ssr.stat"))
            os.link(os.path.join(self.ssr_ref_primer.output_dir, "ssr.ref.result.xls"),
                    os.path.join(self.output_dir, "ssr.ref.result.xls"))
        elif event['data'] == 'makeblastdb':
            self.linkdir(obj.output_dir, self.output_dir + "/makedbblast")
            self.logger.info("\nmkblastdb链接成功\n")
        elif event['data'] == 'ref_2bit':
            self.linkdir(obj.output_dir, self.output_dir)
            self.logger.info("ref2bit链接成功\n")
        elif event['data'] == 'gene_anno':
            self.linkdir(obj.output_dir + "/gene_anno/", self.output_dir)
            self.logger.info("anno.summary链接成功\n")
        elif event['data'] == 'geneanno_sort':
            self.linkdir(obj.output_dir, self.output_dir)
            self.logger.info("排序的anno.summary 链接成功\n")
        else:
            pass

    def run_update_json(self):
        self.logger.info("******self.change is " + str(self.change))
        self.wgsjson = self.get_json(self.file_path + "/wgs_genome.json")
        if self.change:
            self.json_write()
        self.end()

    def json_write(self, path=None):
        """
        """
        path_ = self.file_path + "/wgs_genome.json" if not path else path
        self.wgsjson = self.get_json(path_)
        outfile = self.output_dir + "/wgs_genome.json"
        data = {
            'organism_name': self.species,
            'strain_characteristic': self.strain_characteristic,
            'edition': self.edition,
            'release_date': self.release_time,
            'genome_size': self.genome_size,
            'n50': self.n50,
            'assembly': self.assembly,
            'link': self.link,
            'ref': self.abspath + "/ref.fa",
            'gff': self.abspath + "/ref.gff",
            'anno': self.abspath + "/anno.summary",
            'change_log': self.abspath + "/ref.changelog",
            'info_log': self.abspath + "/info.log",
            'ref_fa_amb': self.abspath + "/ref.fa.amb",
            'ref_fa_ann': self.abspath + "/ref.fa.ann",
            'ref_fa_bwt': self.abspath + "/ref.fa.bwt",
            'ref_fa_fai': self.abspath + "/ref.fa.fai",
            'ref_fa_pac': self.abspath + "/ref.fa.pac",
            'ref_fa_sa': self.abspath + "/ref.fa.sa",
            'ref_dict': self.abspath + "/ref.dict",
            'ref_chrlist': self.abspath + "/ref.chrlist",
            'total_chrlist': self.abspath + "/total.chrlist",
            'snpeff_path': self.abspath + "/snpEff.config",
            'ssr_path': self.abspath + "/",
            'web_dir': self.web_dir
        }
        if self.species not in self.wgsjson.keys():
            self.wgsjson[self.species] = {self.time_edition: data}
            self.logger.info("物种存入")
        else:
            self.wgsjson[self.species][self.time_edition] = data
            self.logger.info("物种另一个版本存入")
        with open(outfile, "w") as dump_f:
            json.dump(self.wgsjson, dump_f, sort_keys=True, indent=4)
            self.logger.info("##### json文件成功输出,改写json文件")

    def checklink(self, newfile, rm=True):
        if os.path.exists(newfile):
            if rm:
                os.remove(newfile)
            return True
        else:
            return False

    def run(self):
        super(GenomeConfigV2Module, self).run()
        if self.option("anno_summary"):
            self.check_summary()
            self.need_anno = False
        self.makedir_info()
        self.openjson()
        # self.check_gff()
        self.check_chr_name()
        tasks = [self.bwa_config, self.picard_dict, self.ssr_ref_primer, self.makeblastdb, self.getgenefasta,
                 self.circos_ref2bit, self.circos_chrlist, self.snpeff_index_v2]
        if self.need_anno:
            tasks.append(self.geneanno_sort)
        self.on_rely(tasks, self.run_update_json)
        if self.option("need_rename") == 'false':
            self.gff_check.on("end", self.run_normalize_fa)
        else:
            self.gff_check.on("end", self.run_genome_grename)
        self.run_gff_check()

    def check_summary(self):
        """
        检查anno.summary
        :return:
        """
        self.logger.info('anno summary check is start!')
        with open(self.option('anno_summary'), 'r') as r:
            for line in r:
                if re.match('#', line):
                    pass
                else:
                    temp = line.strip().split('\t')
                    if temp[0].count('|') != 4:
                        self.set_error("anno.summary文件格式不规范，"
                                       "第一列信息必须是gene_name|GeneID|Genbank|transcript_id|protein ")
                    break
        self.logger.info('anno summary check is end!')

    def check_chr_name(self):
        """
        检查chr，sca命名是否ok
        :return:
        """
        if self.option("chromosomelist"):
            with open(self.option("chromosomelist"), 'r') as r:
                for line in r:
                    temp = line.strip().split('\t')
                    if re.match('Chr.*|Sca.*', temp[1]) or re.match('chr0.*|sca0.*', temp[1]):
                        self.set_error("染色体{}命名不规范".format(temp[1]))

    def check_gff(self):
        """
        检查gff的最后一列注释信息中，是否含有冒号和空格，如果有，用下划线替换。暂时只是检查是否有问题，如果有问题，反馈给产品线
        :return:''
        """
        self.logger.info('refgff check is start!')
        with open(self.option("refgff"), 'r') as r:
            for line in r:
                if re.match('#.*', line):
                    pass
                else:
                    temp = line.strip().split('\t')
                    for m in temp[-1].split(';'):
                        if re.match(r'ID=(.*) (.*)', m) or re.match(r'ID=(.*):(.*)', m):
                            self.set_error('gff文件中注释信息ID中有空格或者冒号，请修改后重新运行！')
                        elif re.match(r'Parent=(.*) (.*)', m) or re.match(r'Parent=(.*):(.*)', m):
                            self.set_error('gff文件中注释信息Parent中有空格或者冒号，请修改后重新运行！')
                        elif re.match(r'Name=(.*) (.*)', m) or re.match(r'Name=(.*):(.*)', m):
                            self.logger.info(m)
                            self.set_error('gff文件中注释信息Name中有空格或者冒号，请修改后重新运行！')
        self.logger.info('refgff check is end!')

    def end(self):
        if self.is_end_by_module:
            self.set_dict_config_path()
        super(GenomeConfigV2Module, self).end()

    def set_dict_config_path(self):
        """
        设置ref.dict与snpEff.config文件路径
        :return:
        """
        if not self.checklink(self.output_dir + "/ref.dict", False) or not self.checklink(
                        self.output_dir + "/snpEff.config", False):
            self.set_error("file{} and {} is not exists".format(self.output_dir + "/ref.dict",
                                                                self.output_dir + "/snpEff.config"))
        else:
            self.logger.info("开始设置ref.dict与snpEff.config文件")
            os.rename(self.output_dir + "/ref.dict", self.output_dir + "/ref.dict_bak")
            os.remove(self.output_dir + "/snpEff.config")
            with open(self.output_dir + "/ref.dict_bak", "r") as r, open(self.output_dir + "/ref.dict", "w") as w:
                datas = r.readlines()
                for line in datas:
                    temp = line.strip().split("\t")
                    if re.match("@HD.*", temp[0]):
                        w.write(line)
                    else:
                        path_ = "UR:file:{}".format(os.path.join(self.species_path, "ref.fa"))
                        w.write("{}\t{}\t{}\t{}\t{}\n".format(temp[0], temp[1], temp[2], temp[3], path_))
            with open(self.output_dir + "/snpEff.config", 'w') as w1:
                w1.write("data.dir ={}/\n".format(self.species_path))
                w1.write("ref.genome : ref\n")
            os.remove(self.output_dir + "/ref.dict_bak")
            self.logger.info("设置ref.dict与snpEff.config文件成功！")

    def sendmail(self):
        """
        发送邮件mail_info = {'project_id': 4567, 'email_id': "hongdong.xuan@majorbio.com"}
        :return:
        """
        self.task_id = 'none'
        project_id = 'none'
        email_id = 'qinwen.xue@majorbio.com'
        species = '--'
        version = '--'
        info = json.loads(self.option('mail_info'))
        if 'task_id' in info.keys():
            self.task_id = info['task_id']
        if "project_id" in info.keys():
            project_id = info['project_id']
        if 'email_id' in info.keys():
            email_id = info['email_id']
        if 'species' in info.keys():
            species = info['species']
        if 'version' in info.keys():
            version = info['version']
        a = SendEmail("RNA_bioinfor@majorbio.com", "smtp.qiye.aliyun.com", "RNA_star2017", "RNA_bioinfor@majorbio.com", email_id,
                      "任务id:{}, 物种: {}, 版本: {}的基因组配置信息以完成-请及时检查！"
                      .format(self.task_id, species, version), 465)
        a.send_msg("{}".format("http://www.{}.com/task/project_tasks/project_id/{}.html"
                               .format(self.task_id.split('_')[0], project_id)))
        # a.attachment(self.mapping_stat.output_dir + "/result.stat/Total.mapped.detail.xls")
        a.send_email()
        self.logger.info("邮件发送完成")

    def uploadfile(self):
        """
        上传文件ref.gff文件+ref.new.mRNA.fa文件
        :return:
        """
        if self.option('target_dir'):
            if self.option("need_rename") == 'true':
                gff = self.genome_grename.output_dir + "/ref.gff"
            else:
                gff = self.output_dir + '/ref.gff'
            mrna = self.gffread_to_fagtf.output_dir + '/ref.new.mRNA.fa'
            target = os.path.dirname(self.option('target_dir').rstrip('/')) + '/genome_check/'
            transfer = MultiFileTransfer()
            if os.path.exists(gff):
                transfer.add_upload(gff, target, base_path=os.path.dirname(gff))
            if os.path.exists(mrna):
                transfer.add_upload(mrna, target, base_path=os.path.dirname(mrna))
            transfer.perform()
            self.logger.info("文件上传成功！")
            if self.option('mail_info'):
                self.sendmail()
            if self.task_id != 'none':
                self.send()
        else:
            pass

    def post_data(self):
        data = dict()
        target = os.path.dirname(self.option('target_dir').rstrip('/')) + '/genome_check/'
        region, bucket = target.split('://')
        content = {
            "task": {
                "task_id": self.task_id
            },
            'files': [{"code": "", "path": "ref.gff", "format": "", "description": "", "size": ''},
                      {"code": "", "path": "ref.new.mRNA.fa", "format": "", "description": "", "size": ''}],
            'dirs': [{"code": "", "description": "", "format": "", "region": region, "path": bucket, "size": ""}],
            'base_path': bucket,
            'region': region
        }
        data['sync_task_log'] = json.dumps(content, cls=CJsonEncoder)
        yield urllib.urlencode(data)

    def send(self):
        """
        {'sync_task_log': '{"files": [{"size": "", "path": "ref.gff", "code": "", "description": "", "format": ""},
         {"size": "", "path": "ref.new.mRNA.fa", "code": "", "description": "", "format": ""}],
         "dirs": [{"code": "", "description": "", "format": "", "region": "s3",
         "path": "common/files/m_188/188_5cb3e5db67c6a/tsg_34362/genome_check/", "size": ""}],
         "region": "s3", "task": {"task_id": "tsg_34362"},
         "base_path": "common/files/m_188/188_5cb3e5db67c6a/tsg_34362/genome_check/"}'}
        :return:
        """
        _url = ''
        # if re.match('tsg_.*', self.task_id):
        #     self._key = 'hM4uZcGs9d'
        #     self._client = "client03"
        #     _url = "http://api.tsg.com/task/add_task_log"
        # elif re.match('tsanger_.*', self.task_id):
        #     self._key = 'hM4uZcGs9d'
        #     self._client = "client03"
        #     _url = "http://api.tsanger.com/task/add_task_log"
        # elif re.match('sanger_.*', self.task_id) or re.match('i-sanger_.*', self.task_id) or re.match('majorbio_.*',
        #                                                                                               self.task_id):
        #     self._key = '1ZYw71APsQ'
        #     self._client = "client01"
        #     _url = "http://api.sanger.com/task/add_task_log"
        # else:
        #     self.set_error('任务id，不合法')
        if self._client == "client03":
            # self._key = 'hM4uZcGs9d'
            # _url = "http://api.tsg.com/task/add_task_log"
            _url = "http://apicenter.nsg.com/index/in"
            self._key = "458b97de4c0bb5bf416c8cea208309ed"
            self._binds_id = "5f4324f09b7900009300805c"
            self._interface_id = 57
            self._env_name = "online"
        else:
            # self._key = '1ZYw71APsQ'
            # _url = "http://api.sanger.com/task/add_task_log"
            _url = "http://apicenter.lab.majorbio.com/index/in"
            self._key = "e0bb04dd111155b2e6bc6db26d0e1fef"
            self._binds_id = "5f4324f09b7900009300805c"
            self._interface_id = 57
            self._env_name = "online"
        # self.logger.info(self.post_data())
        for p_data in self.post_data():
            self.logger.info(p_data)
            http_handler = urllib2.HTTPHandler(debuglevel=1)
            https_handler = urllib2.HTTPSHandler(debuglevel=1)
            opener = urllib2.build_opener(http_handler, https_handler)
            urllib2.install_opener(opener)
            post_data = "%s&%s" % (self.get_sig(), p_data)
            request = urllib2.Request(_url, post_data)
            response = urllib2.urlopen(request, timeout=60)
            self._response_code = response.getcode()
            self._response = response.read()
            response.close()
            self.logger.info("Return page:%s" % self._response)
        return self._response

    def get_sig(self):
        nonce = str(random.randint(1000, 10000))
        timestamp = str(int(time.time()))
        x_list = [self._key, timestamp, nonce]
        x_list.sort()
        sha1 = hashlib.sha1()
        map(sha1.update, x_list)
        sig = sha1.hexdigest()
        signature = {
            "client": self._client,
            "nonce": nonce,
            "timestamp": timestamp,
            "signature": sig,
            "binds_id": self._binds_id,
            "interface_id": self._interface_id,
            "env_name": self._env_name
        }
        return urllib.urlencode(signature)
