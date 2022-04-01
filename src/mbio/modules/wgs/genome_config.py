# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0521
# modified 2018.0524

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import os, re
import json
import time, datetime


class GenomeConfigModule(Module):
    """
    样本基因组SSR分析及引物设计
    # /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/modules/wgs/genome_config.py
    # /mnt/ilustre/users/sanger-dev/biocluster/scripts/wgs/genome_config_module_run.py
    """
    def __init__(self, work_id):
        super(GenomeConfigModule, self).__init__(work_id)
        options = [
            {"name": "reffa", "type": "string"},       # 更名后ref.fa文件 infile核查
            {"name": "refgff", "type": "string"},      # 更名后ref.gff文件
            {"name": "info", "type": "string"},             # info.log
            # {"name": "chromosomelist", "type": "string"},     # 无染色体的时候 可无此文件
            {"name": "anno_summary", "type": "string"},     # 需要写一个核查文件
            {"name": "ref_changelog", "type": "string"},
            {"name": "update", "type": "string", "default": 0},
        ]
        self.add_option(options)
        self.species_path = ''
        self.species = ''
        self.edition = ''
        self.release_time = ''
        self.wgsjson = ''
        self.source = ''
        self.bwa_config = self.add_tool("wgs.bwa_config")
        self.samtools_faidx = self.add_tool("wgs.samtools_faidx")
        # self.genome_grename = self.add_tool("wgs.genome_grename")
        self.picard_dict = self.add_tool("wgs.picard_dict")
        self.circos_chrlist = self.add_tool("wgs.circos_chrlist")
        self.snpeff_gffcheck = self.add_tool("wgs.snpeff_gffcheck")
        self.snpeff_index = self.add_tool("wgs.snpeff_index")
        self.ssr_ref_primer = self.add_module("wgs.ssr_ref_primer")
        self.makeblastdb = self.add_tool("wgs.makeblastdb")
        # self.getgenefasta = self.add_tool("wgs.getgenefasta")
        # self.gene_anno = self.add_module("wgs.gene_anno")
        self.circos_ref2bit = self.add_tool("wgs.circos_ref2bit")
        self.geneanno_sort = self.add_tool("wgs.geneanno_sort")

    def check_options(self):
        if not self.option("reffa"):
            raise OptionError("请设置reffa")
        if not self.option("refgff"):
            raise OptionError("请设置refgff")
        if not self.option("info"):
            raise OptionError("请设置info文件")
        if not self.option("anno_summary"):
            raise OptionError("请设置anno_summary")
        if not self.option("ref_changelog"):
            raise OptionError("请设置ref_changelog")

    def makedir_info(self):     # mkdir
        """
        info.log:
        cuiArabidopsis_thaliana Columbia0(TAIR10)   https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.3/  NCBI
        GCF_000001735.3 2015.04.16  116M(genome size)   31,248,787(n50) Chromosome(assembly)  web_dir
        self.file_path = Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/cui_test_dir"   # 测试目录 拷贝json文件
        链接输入文件到inputfile；
        """
        sym1 = "_"
        sym2 = "."
        temp = []
        self.file_path = Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/"
        with open(self.option("info"), 'r') as file:
            data = file.readlines()
            for line in data:
                temp = line.strip().split('\t')
                print(temp)
                print(data)
                print(len(temp))
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
                    self.species_path = os.path.join(self.file_path, self.species, self.source, self.edition, self.release_time)
                    print("##### 物种绝对路径" + self.species_path)
                    if not os.path.exists(self.species_path):
                        os.makedirs(self.species_path)
                        print("##### " + self.species_path + " 该物种路径新建成功")
                else:
                    self.logger.info("@@@@@ 请check info.log文件列数--10列，用\\t分割! ! !")
                    print("log文件列数少于10列，正常结束")
                    self.end()

    def openjson(self):
        """
        链接输入文件到目录下
        不需要链接原始输入文件
        """
        # input_file = self.species_path + "/inputfile"
        # if not os.path.exists(input_file):
        #     os.makedirs(input_file)
        #     print("##### 原输入文件{}链接路径新建成功".format(input_file))
        # for ii in ["reffa", "refgff", "info", "chromosomelist"]:
        #     jj = os.path.basename(self.option(ii))
        #     self.checklink(os.path.join(input_file, jj))
        #     os.link(self.option(ii), os.path.join(input_file, jj))
        self.abspath = os.path.join(self.species, self.source, self.edition, self.release_time)
        self.time_edition = self.release_time + "||" + self.edition
        self.wgsjson = self.get_json(self.file_path + "/wgs_genome.json")
        # b = time.strftime("%Y%m%d",time.localtime(time.time()))
        b = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        os.system("cp {} {}".format(self.file_path + "/wgs_genome.json",
                                    Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/z.wgs_json/wgs_genome.json_" + b))
        print("##### 物种是" + self.species)
        print("##### 物种版本是" + self.time_edition)
        print("##### 已有物种名如下")
        print(self.wgsjson.keys())
        self.change = 0
        print(self.option("update"))
        print("#######\n")
        if self.species not in self.wgsjson.keys() or (self.species in self.wgsjson.keys() and self.time_edition not in self.wgsjson[self.species].keys()):
            self.change = 1
        elif int(self.option("update")) == 1:
            self.change = 1
        else:
            print(self.species + "\t" + self.time_edition + "该物种id无需导入!")
            print("无需导入，正常结束")
            self.end()

    def get_json(self, json_path):
        """
        用于解析出json文件，然后可以从中获取参考基因组数据
        :return:
        """
        f = open(json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

    # def run_genome_grename(self):
    #     """
    #     将下载的fa gff进行更名；
    #     result: ref.fa ref.gff ref.changelog
    #     """
    #     opti = {"fasta": self.option("reffa"), "gff": self.option("refgff")}
    #     if self.option("chromosomelist"):
    #         opti["chromosomelist"] = self.option("chromosomelist")
    #     self.genome_grename.set_options(opti)
    #     self.genome_grename.on("end", self.run_stat_annotation)
    #     self.genome_grename.run()

    def fa_link(self):
        """
        功能：链接和更名
        "ref.fa", "ref.gff", "ref.changelog"：即genome_grename输出文件被链接到基因组目录。
        "info.log"输入文件也被链接到基因组目录。
        """
        # option_list = ["reffa", "refgff", "info", "anno_summary", "ref_changelog"]
        option_list = ["reffa", "refgff", "info", "ref_changelog"]
        name_change = ['ref.fa', 'ref.gff', 'info.log', 'ref.changelog']
        for i in range(len(option_list)):
            print "1", self.option(option_list[i])
            print "2", os.path.join(self.species_path, name_change[i])
            oldpathname = os.path.join(self.option(option_list[i]))
            newpathname = os.path.join(self.species_path, name_change[i])
            print(oldpathname)
            self.checklink(newpathname)
            os.link(oldpathname, newpathname)
            print("##### {}链接到物种目录下成功".format(newpathname))

    def run_samtools_faidx(self):
        self.fa_link()
        self.ref_path = os.path.join(self.species_path, "ref.fa")
        self.gff_path = os.path.join(self.species_path, "ref.gff")
        self.samtools_faidx.set_options({
            "pop_fa": self.ref_path
        })
        self.samtools_faidx.on("end", self.run_bwa_config)          # bwa配置
        self.samtools_faidx.on("end", self.run_picard_dict)         # picard配置
        self.samtools_faidx.on("end", self.run_snpeff_gffcheck)     # snpeff配置
        self.samtools_faidx.on("end", self.run_ssr_ref_primer)      # ssr module配置
        self.samtools_faidx.on("end", self.run_makeblastdb)         # makeblastdb
        # self.samtools_faidx.on("end", self.run_getgenefasta)        # 生成ref.gene.fa
        self.samtools_faidx.on("end", self.run_ref_2bit)            # ref_2bit给基因组浏览器用
        self.samtools_faidx.on("end", self.run_geneanno_sort)
        self.samtools_faidx.run()
        print("\n##### samtools_faidx运行结束\n")

    def run_bwa_config(self):       # 1
        """
        self.species_path + "/ref.fa"确定存在！
        """
        self.bwa_config.set_options({
            "ref_fa": self.ref_path,
        })
        self.bwa_config.run()
        print("\n##### bwa_config运行结束\n")

    def run_picard_dict(self):  # 2.1
        """
        check: ref.dict删掉重生成
        """
        self.dict_path = os.path.join(self.species_path, "ref.dict")
        self.checklink(self.dict_path)
        self.picard_dict.set_options({
            "fa": self.ref_path,
        })
        self.picard_dict.on("end", self.run_circos_chrlist)
        self.picard_dict.run()
        print("\n##### picard_dict运行结束\n")

    def run_circos_chrlist(self):   # 2.2
        """
        生成ref.chrlist
        """
        self.circos_chrlist.set_options({
            "dict": self.dict_path,
        })
        self.circos_chrlist.on("end", self.hardlink, "ref.chrlist")
        self.circos_chrlist.run()
        print("\n##### circos_chrlist运行结束\n")

    def run_snpeff_gffcheck(self):  # 3.1
        """
        参考基因组snpeff check ref.gff生成genes.gff
        """
        self.snpeff_gffcheck.set_options({
            "refgff": self.gff_path,
        })
        self.snpeff_gffcheck.on("end", self.run_snpeff_index)
        self.snpeff_gffcheck.on("end", self.hardlink, 'genes.gff')
        self.snpeff_gffcheck.run()
        print("\n##### snpeff_gffcheck运行结束\n")

    def run_snpeff_index(self):     # 3.2
        """
        生成结果在ref.fa所在目录，不需要链接
        """
        self.snpeff_index.set_options({
            "reffa": self.ref_path,
            "genesgff": os.path.join(self.snpeff_gffcheck.output_dir, "genes.gff")
        })
        self.snpeff_index.run()
        print("\n##### snpeff_index运行结束\n")

    def run_ssr_ref_primer(self):   # 4 ssr的module
        """
        结果：ref.ssr.stat  ssr.ref.result.xls 需要硬链接为ssr.stat ssr.ref.result.xls
        """
        self.ssr_ref_primer.set_options({
            "reffa": self.ref_path
        })
        self.ssr_ref_primer.on("end", self.hardlink, "ref.ssr.stat")
        self.ssr_ref_primer.run()
        print("\n##### ssr_ref_primer运行结束\n")

    def run_makeblastdb(self):   # 5
        dbpath = os.path.join(self.species_path, "makedbblast")
        if not os.path.exists(dbpath):
            os.mkdir(dbpath)
        options = {
            "pop_fa": self.ref_path
        }
        self.makeblastdb.set_options(options)
        self.makeblastdb.on("end", self.hardlink, "makeblastdb")
        self.makeblastdb.run()
        print("\n##### makeblastdb运行结束\n")

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
        print("\n##### ref_2bit运行结束\n")

    def run_geneanno_sort(self):
        """
        将产品线的anno.summary按照chrid排序后存进来
        # /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/wgs/geneanno_sort.py
        """
        self.geneanno_sort.set_options({
            "anno_summary": self.option("anno_summary"),
        })
        self.geneanno_sort.on("end", self.hardlink, "geneanno_sort")
        self.geneanno_sort.run()
        print("\n##### geneanno_sort运行结束\n")

    # def run_getgenefasta(self):
    #     self.getgenefasta.set_options({
    #         "reffa": self.ref_path,
    #         "refgff": self.gff_path,
    #     })
    #     self.getgenefasta.on("end", self.run_gene_anno)
    #     self.getgenefasta.run()
    #     print("\n##### getgenefasta运行结束\n")

    # def run_gene_anno(self):
    #     """
    #     输入文件：ref.gene.fa;先分20份然后被注释
    #     """
    #     self.gene_anno.set_options({
    #         "fasta": os.path.join(self.getgenefasta.output_dir, "ref.gene.fa")
    #     })
    #     self.gene_anno.on("end", self.hardlink, "gene_anno")
    #     self.gene_anno.run()
    #     print("\n##### gene_anno运行结束\n")

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
            self.linkdir(obj.output_dir, self.species_path)
        elif event['data'] == 'genes.gff':
            self.linkdir(obj.output_dir, self.species_path + "/ref")
        elif event['data'] == 'ref.ssr.stat':
            self.checklink(os.path.join(self.species_path, "ssr.stat"))
            self.checklink(os.path.join(self.species_path, "ssr.ref.result.xls"))
            os.link(os.path.join(self.ssr_ref_primer.output_dir, "ref.ssr.stat"), os.path.join(self.species_path, "ssr.stat"))
            os.link(os.path.join(self.ssr_ref_primer.output_dir, "ssr.ref.result.xls"), os.path.join(self.species_path, "ssr.ref.result.xls"))
            print("SSR结果链接成功\n")
        elif event['data'] == 'makeblastdb':
            self.linkdir(obj.output_dir, self.species_path + "/makedbblast")
            # self.checklink(os.path.join(self.species_path, "makedbblast", "ref.nhr"))
            # self.checklink(os.path.join(self.species_path, "makedbblast", "ref.nin"))
            # self.checklink(os.path.join(self.species_path, "makedbblast", "ref.nsq"))
            # if not os.path.exists(os.path.join(self.species_path, "ref")):
            #     os.mkdir(os.path.join(self.species_path, "ref"))
            # os.link(os.path.join(self.makeblastdb.output_dir, "dbname.nhr"), os.path.join(self.species_path, "makedbblast", "ref.nhr"))
            # os.link(os.path.join(self.makeblastdb.output_dir, "dbname.nin"), os.path.join(self.species_path, "makedbblast", "ref.nin"))
            # os.link(os.path.join(self.makeblastdb.output_dir, "dbname.nsq"), os.path.join(self.species_path, "makedbblast", "ref.nsq"))
            print("\nmkblastdb链接成功\n")
        elif event['data'] == 'ref_2bit':
            self.linkdir(obj.output_dir, self.species_path)
            print("ref2bit链接成功\n")
        # elif event['data'] == 'gene_anno':
        #     self.linkdir(obj.output_dir + "/gene_anno/", self.species_path)
        #     print("anno.summary链接成功\n")
        elif event['data'] == 'geneanno_sort':
            self.linkdir(obj.output_dir + "/", self.species_path)
            print("排序的anno.summary 链接成功\n")
        else:
            pass

    def run_update_json(self):
        print("******self.change is " + str(self.change))
        self.wgsjson = self.get_json(self.file_path + "/wgs_genome.json")
        while self.change:
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
                print("物种存入")
            else:
                self.wgsjson[self.species][self.time_edition] = data
                print("物种另一个版本存入")
            with open(self.file_path + "/wgs_genome.json", "w") as dump_f:
                json.dump(self.wgsjson, dump_f, sort_keys=True, indent=4)
                print("##### json文件成功输出,改写json文件")
                print("######################################")
                print("\n" + self.species_path + "\n")
                print("######################################")
                break
        self.end()

    def checklink(self, newfile):
        if os.path.exists(newfile):
            os.remove(newfile)
        else:
            pass
        return True

    def run(self):
        super(GenomeConfigModule, self).run()
        self.makedir_info()
        self.openjson()
        # self.run_genome_grename()
        # self.run_stat_annotation()
        # self.run_getgenefasta()
        self.on_rely([self.bwa_config, self.circos_chrlist, self.snpeff_index, self.ssr_ref_primer,
                      self.makeblastdb, self.geneanno_sort, self.circos_ref2bit], self.run_update_json)
        self.run_samtools_faidx()

    def end(self):
        super(GenomeConfigModule, self).end()
