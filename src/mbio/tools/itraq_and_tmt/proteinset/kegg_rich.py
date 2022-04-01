# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from collections import defaultdict
from bson import ObjectId
import os
import re
import time
import json
from mbio.packages.rna.kegg_html import KeggHtml
from mbio.packages.itraq_and_tmt.kegg_regulate import KeggRegulate
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.config import Config
import subprocess
import tarfile
import shutil

from biocluster.core.exceptions import OptionError
import os
from mbio.packages.rna.annot_config import AnnotConfig


class KeggRichAgent(Agent):
    """
    Kegg富集分析
    version v1.0.1
    author: qiuping
    last_modify: 2016.11.23
    """
    def __init__(self, parent):
        super(KeggRichAgent, self).__init__(parent)
        options = [
            {"name": "kegg_table", "type": "infile", "format": "labelfree.kegg_table"},
            # 只含有基因的kegg table结果文件
            {"name": "diff_list", "type": "infile", "format": "denovo_rna_v2.gene_list"},
            {"name": "proteinset_kegg", "type": "string", "default": None},  # 基因集与基因的对应文件，行数等于基因集的数目
            {"name": "proteinset_id", "type": "string"},
            {"name": "type", "type": "string", "default": "origin"},  # 取最新的注释表还是原来的注释表
            # gene名字文件
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "add_info", "type": "string", "default": None},  # 底图文件地址
            {"name": "task_id", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': '201909'},
            {"name": "png_dir", "type": "string", "default": ""}  # 如果是工作流里的则需要传入这个参数
        ]
        self.add_option(options)
        self.step.add_steps("kegg_rich")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_rich.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_rich.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('kegg_table').is_set:
            raise OptionError('必须设置kegg的pathway输入文件', code = "32503601")
        if self.option('correct').lower() not in ['by', 'bh', 'bonferroni', 'holm']:
            raise OptionError('多重检验校正的方法不在提供的范围内', code = "32503602")
        if not self.option("diff_list").is_set:
            raise OptionError("必须设置输入文件diff_list", code = "32503603")

        # if not self.option("proteinset_kegg"):
        #     raise OptionError("必须设置输入文件proteinset_kegg", code = "32503603")
        return True



    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r"kegg_enrichment.xls$", "xls", "kegg富集分析结果"]
        ])
        super(KeggRichAgent, self).end()


class KeggRichTool(Tool):
    def __init__(self, config):
        super(KeggRichTool, self).__init__(config)
        self._version = "v1.0.1"
        # self.kobas = '/bioinfo/annotation/kobas-2.1.1/src/kobas/scripts/'
        # self.kobas_path = self.config.SOFTWARE_DIR + '/bioinfo/annotation/kobas-2.1.1/src/'
        # self.set_environ(PYTHONPATH=self.kobas_path)
        # self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        # self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        # self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        # self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.python = '/program/Python/bin/'
        # self.all_list = self.option('all_list').prop['gene_list']
        # self.diff_list = self.option('diff_list').prop['gene_list']
        self.script_path = self.config.PACKAGE_DIR + "/itraq_and_tmt/kegg_enrichment.py"

        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.k2e = self.kegg_files_dict['K2enzyme.tab']
        self.brite = self.kegg_files_dict['br08901.txt']

        self.db_path_old = self.config.SOFTWARE_DIR + "/database/KEGG/xml/"
        self.db_path = self.kegg_files_dict["html"]
        self.geneset_ko = list()
        self.proteinset_list = list()
        self.map_dict = {}
        self.species_abr = 'map'
        self.gene2set = dict()
        self.genome_path = self.kegg_files_dict["genome2.xls"]

        self.map_dict = dict()  # 获得背景色连接
        self.front_dict = dict()

    def run(self):
        """
        运行
        :return:
        """
        super(KeggRichTool, self).run()
        # if self.option("diff_stat").is_set:
        #     self.run_kegg_rich()
        # else:
        #     self.run_web_kegg()
        self.run_kegg_rich()
        self.run_identify()
        if self.option("add_info"):
            self.run_add_info()
        if self.option("proteinset_kegg"):
            self.run_html()
        self.end()


    def run_html(self):
        self.get_kegg_png()
        self.get_dicts()
        pathways = self.output_dir + '/pathways'
        if not os.path.exists(pathways):
            os.mkdir(pathways)
        self.generate_new_pics()


    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            if to_path.endswith('/'):
                os.system('ln -s {} {}'.format(path, to_path.rstrip('/')))
        else:
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(path, to_path)
                transfer.perform()
            except:
                self.set_error('file can not find %s', variables=(path), code="32503503")
        if not os.path.exists(to_path):
            self.set_error("文件下载失败")
        elif not os.listdir(to_path):
            self.set_error("下载的文件夹为空")
        else:
            return to_path

    def get_kegg_png(self):
        if self.option('task_id') == 'itraq_and_tmt_default':
            self.logger.info('直接从工作流结果复制kegg图片')
            if os.path.isdir(self.option('png_dir')):
                anno_path = self.work_dir + "/png"
                if os.path.exists(anno_path):
                    shutil.rmtree(anno_path)
                cmd = 'cp -r %s %s' %(self.option('png_dir').rstrip('/'), anno_path)
                os.system(cmd)
            import pandas as pd
            kegg_df = pd.read_csv(self.option('kegg_table').prop['path'], sep='\t').fillna('')
            paths = sum([p.split(';') for p in kegg_df['Paths'] if p], [])
            specie = paths[0]
            self.species_abr = filter(str.isalpha, specie)
            if self.species_abr.lower() == 'ko':
                self.species_abr = 'map'
            return
        proteinset_id = self.option("proteinset_id")
        self.mongo_db = Config().get_mongo_client(mtype="itraq_tmt")[Config().get_mongo_dbname("itraq_tmt")]
        kegg_collection = self.mongo_db["sg_annotation_stat"]
        proteinset_collection = self.mongo_db["sg_proteinset"]
        if "," in proteinset_id:
            proteinset_id = proteinset_id.split(",")[0]
        result = proteinset_collection.find_one({"main_id": ObjectId(proteinset_id)})
        if not result:
            self.set_error('proteinset with main_id: %s was not found', variables = (proteinset_id), code = "32503502")
        task_id = result["task_id"]
        kegg_annot_output = kegg_collection.find_one({"task_id": self.option('task_id'), "type": self.option('type'), 'status': "end"})['result_dir']
        annot_params = kegg_collection.find_one({"task_id": self.option('task_id'), "type": self.option('type'), 'status': "end"})['params']
        params_dict = json.loads(annot_params)
        if params_dict.has_key('kegg_species'):
            self.get_genome_abr(params_dict['kegg_species'])

        kegg_created_ts = kegg_collection.find_one({"task_id": self.option('task_id'), "type": self.option('type'), 'status': "end"})['created_ts']
        if kegg_created_ts < "2021-07-28 18:10:49":
            kegg_pathways = os.path.join(kegg_annot_output, 'kegg', 'pathways')
            if not kegg_pathways.endswith("/"):
                kegg_pathways += "/"
            anno_path = self.work_dir + "/raw_png/"
            self.download_s3_file(kegg_pathways, anno_path)
            shutil.copytree("raw_png", "png")
            if os.path.exists(anno_path + "pathways.tar.gz"):
                tar = tarfile.open(anno_path + "pathways.tar.gz")
                tar.extractall(path=self.work_dir + "/png")
            elif os.path.exists(anno_path):
                pass
            else:
                self.set_error("无法获得kegg注释结果图片信息", code="32503504")
        else:
            # 旧的结果中，pathways是个目录，缺点是文件多传输慢。如今改为打包形式(pathways.tar)，需要下载包并解包。
            kegg_pathways = os.path.join(kegg_annot_output, 'kegg', 'pathways.tar')
            try:
                transfer = MultiFileTransfer()
                transfer.add_download(kegg_pathways, self.work_dir)
                transfer.perform()
            except:
                pass
            if not os.path.exists(os.path.join(self.work_dir, 'pathways.tar')):
                kegg_pathways = os.path.join(kegg_annot_output, '3_anno_KEGG', 'pathways.tar')
                try:
                    transfer = MultiFileTransfer()
                    transfer.add_download(kegg_pathways, self.work_dir)
                    transfer.perform()
                except:
                    self.set_error('无法获得kegg注释结果图片信息,file can not find %s', variables=(kegg_pathways), code="32503503")
            if os.path.exists(os.path.join(self.work_dir, 'pathways.tar')):
                pathways_tar = tarfile.open(os.path.join(self.work_dir, 'pathways.tar'))
                pathways_tar.extractall(path=self.work_dir)
                os.rename(os.path.join(self.work_dir, 'pathways'),os.path.join(self.work_dir, 'png'))

    def get_genome_abr(self, kegg_species):
        with open (self.genome_path, 'r' ) as f:
            lines = f.readlines()
            for line in lines:
                genome_id = re.sub("gn:", "", line.split("\t")[0])
                genome_abr = line.split("\t")[1].split(',')[0]
                genome = line.split("\t")[1].split(';')[-1].strip()
                if kegg_species == genome_abr or kegg_species == genome:
                    self.species_abr = genome_abr
                    break

    def generate_new_pics(self):
        path_ko=self.path_ko
        kos_path=self.output_dir + "/ko"
        out_dir=self.output_dir
        path_list = path_ko.keys()
        cmd1_list = []
        cmd2_list = []
        for path in path_list:
            # ko_path = kos_path + "/" + path
            if 1:
                ko_old = path.replace(self.species_abr, "ko")
                ko = path.replace(self.species_abr, "map")
                if ko.startswith("ko"):
                    ko.replace("ko", "map")
                ko_xml = self.db_path + "/" + ko + ".kgml"
                if self.option("kegg_version") <= "2020":
                    if os.path.exists(self.db_path_old + ko_old + ".xml"):
                        ko_xml = self.db_path_old + ko_old + ".xml"
                    else:
                        ko_xml = self.db_path + ko + ".kgml"

                map_html = KeggHtml(version=self.option('kegg_version'))
                if os.path.exists(self.work_dir + "/png/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/png/" + path + ".html.mark", self.proteinset_ko)
                elif os.path.exists(self.work_dir + "/kegg/" + path + ".html.mark"):
                    map_html.run_gene_set_mark(self.work_dir + "/png/" + path + ".html", out_dir + "/pathways/" + path + ".html.mark", self.gene2set, self.work_dir + "/kegg/" + path + ".html.mark", self.proteinset_ko)
                else:
                    self.logger.info("{}注释结果中没有html 或 mark文件".format(path))

            else:
                pass

        self.logger.info("kegg图片生成完毕")


    def get_dicts(self):
        gene_kegg_gene = dict()
        if self.species_abr == 'map':
            ko_genes, path_ko = self.option('kegg_table').get_pathway_koid(kegg_version=self.option('kegg_version'))
        else:
            ko_genes, path_ko, gene_kegg_gene = self.option('kegg_table').get_pathway_koid2(kegg_version=self.option('kegg_version'))
        proteinset_ko = defaultdict(set)
        regulate_gene = {}
        gene2set = dict()
        proteinset2_ko = list()
        with open(self.option("proteinset_kegg"), "r") as f:
            for line in f:
                line = line.strip().split("\t")
                regulate_gene[line[0]] = line[1].split(",")
                proteinset_ko[line[0]] = []
                self.proteinset_list.append(line[0])
                for key in ko_genes.keys():
                    for gene in regulate_gene[line[0]]:
                        if gene in ko_genes[key]:
                            proteinset_ko[line[0]].append(key)
                proteinset2_ko.append(set(proteinset_ko[line[0]]))
                for gene in regulate_gene[line[0]]:
                    if gene2set.has_key(gene):
                        gene2set[gene].append(line[0])
                    else:
                        gene2set[gene] = [line[0]]

                    if self.species_abr == 'map':
                        pass
                    else:
                        if gene_kegg_gene.has_key(gene):
                            for kegg_genes in gene_kegg_gene[gene]:
                                for kegg_gene in kegg_genes.split(";"):
                                    proteinset_ko[line[0]].append(kegg_gene)
                        else:
                            pass
        for i in gene2set.keys():
            self.gene2set[i] = ",".join(gene2set[i])
        self.logger.info("path is {}".format(";".join(path_ko.keys())))
        self.logger.info("gene2set is {}".format(self.gene2set))

        self.category = proteinset_ko  # 基因集与ko的对应关系
        self.ko_genes= ko_genes  # ko与基因的对应关系
        self.path_ko = path_ko  # path小ko与基因大Ko的对应关系
        self.proteinset_gene = regulate_gene  # 基因集与基因的对应关系
        self.proteinset_ko = proteinset2_ko
        self.gene_kegg_gene = gene_kegg_gene #基因与kegg基因的对应关系

    def get_geneset_type(self, diff_file):
        """
        查看基因集是已知基因集还是已知基因+新基因 刘彬旭
        """
        gene_set_f = open(diff_file, 'r')
        for line in gene_set_f.readlines():
            if line.startswith('MSTR') or line.startswith('TCON') or line.startswith('XLOC'):
                gene_set_f.close()
                return 'all'
        gene_set_f.close()
        return 'known'

    def run_kegg_rich(self):
        """
        运行kobas软件，进行kegg富集分析
        """
        try:
            # self.option('kegg_table').get_kegg_list(self.work_dir, self.all_list, self.diff_list)
            # self.logger.info("kegg富集第一步运行完成")
            self.logger.info("准备gene path:gene_Knumber G2K文件")
            self.option("kegg_table").get_gene2K(self.work_dir)
            self.logger.info("准备gene path:konumber G2K文件")
            self.option("kegg_table").get_gene2path(self.work_dir)
            self.logger.info("准备差异文件")
            # diff_gene, regulate_dict = self.option("diff_list").get_table_info()
            self.option("diff_list").get_stat_file(self.work_dir, self.work_dir + "/gene2K.info")
            self.logger.info("统计背景数量")
            length = os.popen("less {}|wc -l".format(self.work_dir + "/gene2K.info"))
            line_number = int(length.read().strip("\n"))
            self.bgn = line_number - 1  # 去掉头文件
            # self.run_identify()
        except Exception as e:
            self.set_error("kegg富集第一步运行出错:%s", variables = (e), code = "32503604")
            self.logger.info("kegg富集第一步运行出错:{}".format(e))
        # self.run_identify()

    def choose_known_background(self, annot_file):
        """
        如果基因集中仅含有已知基因，修改背景注释为已知基因注释 刘彬旭
        """
        annot_file_new = annot_file + ".known"
        with open(annot_file_new, "w") as newfw, open(annot_file, "r") as anf:
            newfw.write(anf.readline())
            for line in anf.readlines():
                if line.startswith('MSTR') or line.startswith('TCON') or line.startswith('XLOC'):
                    pass
                else:
                    newfw.write(line)
            return annot_file_new

    def check_list(self, ko_file):
        """
        去除diff_list中没有注释信息的数据
        new_file_name为在work_dir中生成的新diff_list文件的绝对路径
        :return:
        """
        file1 = ko_file
        file2 = self.work_dir + "/gene2K.info"
        f1 = open(file1, "r")
        f2 = open(file2, "r")
        lst_2 = [x.split('\t')[0] for x in f2.readlines() if x.split('\t')[1]!="\n"]
        f2.close()
        new_file_name = ko_file + ".check"
        with open(new_file_name, "w") as check_file:
            for item in f1.readlines():
                gene = item.split('\t')[0]
                if gene in lst_2:
                    check_file.write(item)
        f1.close()
        return new_file_name

    def run_identify(self):
        deg_path = self.work_dir + "/" + os.path.basename(self.option('diff_list').prop["path"]) + ".DE.list"
        # kofile = os.path.splitext(os.path.basename(self.option('diff_list').prop['path']))[0]
        kofile = os.path.basename(self.option('diff_list').prop['path']) + ".DE.list.check"
        deg_path = self.check_list(deg_path)

        g2p_path = self.work_dir + "/gene2path.info"
        g2k_path = self.work_dir + "/gene2K.info"
        bgn = self.bgn
        k2e = self.k2e
        brite = self.brite

        # 当基因集仅包含已知基因时只使用已知基因作为背景
        # if self.get_geneset_type(self.option('diff_list').prop['path']) == 'known':
        #     g2p_path = self.choose_known_background(g2p_path)
        #     g2k_path = self.choose_known_background(g2k_path)
        #     length = os.popen("less {}|wc -l".format(self.work_dir + "/gene2K.info.known"))
        #     line_number = int(length.read().strip("\n"))
        #     bgn = line_number - 1
        # correct method:
        correct = self.option("correct")
        correct_method = 3
        if correct.lower() == "bh":
            correct_method = 3
        elif correct.lower() == 'bonferroni':
            correct_method = 1
        elif correct.lower() == 'holm':
            correct_method = 2
        elif correct.lower() == 'by':
            correct_method = 4
        else:
            print('correct method not exist, BH will be used')

        # end of correct method
        cmd_2 = self.python + 'python {} ' \
                              '-deg {} -g2p {} -g2k {} -bgn {} -k2e {} -brite {} --FDR -dn 20 ' \
                              '-correct {}'.format(self.script_path, deg_path, g2p_path, g2k_path,
                                                   bgn, k2e, brite, correct_method)
        self.logger.info('开始运行kegg富集第二步：进行kegg富集分析')
        command_2 = self.add_command("cmd_2", cmd_2).run()
        self.wait(command_2)
        if command_2.return_code == 0:
            self.logger.info("kegg富集分析运行完成")
            self.output_name = kofile + '.kegg_enrichment.xls'
            self.set_output(kofile + '.kegg_enrichment.xls')
            # self.end()
        else:
            self.set_error("kegg富集分析运行出错!", code = "32503605")

    def run_add_info(self):
        background_links = self.option("add_info")
        if not os.path.exists(background_links):
            self.logger.info("不存在背景信息文件，程序退出")
            return
        ko2spe = self.option("kegg_table").get_ko2spegene()
        self.option("kegg_table").get_gene2K(self.work_dir)
        with open(background_links) as r:
            r.readline()
            for line in r:
                tmp = line.strip().split("\t")
                pathway = tmp[0]
                links = tmp[1].split("?")[1].split("/")
                links.pop(0)  # 去除掉map
                dict = {}
                for link in links:
                    ko = link.split("%09")[0]  # K06101
                    color = link.split("%09")[1]  # tomato
                    dict[ko] = color
                self.map_dict[pathway] = dict  # self.map_dict["map05340"] = {"K10887": "yellow"}
        with open(self.output_name, "r") as file, open(self.output_name + "_new", "w") as fw:
            line = file.readline()
            fw.write(line)
            for line in file:
                tmp = line.strip().split("\t")
                links = tmp[9]
                pathway = tmp[3]
                links_ko = links.split("?")[1].split("/")
                links_ko.pop(0)  # 去除掉map
                lnk = links.split("?")[0] + "?map=" + pathway + "&multi_query="
                # self.front_dict[pathway] = []
                # for link in links_ko:
                #     ko = link.split("%09")[0]
                #     self.front_dict[pathway].append(ko)
                #     if pathway in self.map_dict:
                ko_tmp = [x.split("%09")[0] for x in links_ko]  # ko gene_list
                if pathway.startswith("map") or pathway.startswith("ko"):
                    pass
                else:
                    ko_tmp_spe = []
                    for ko_a in ko_tmp:
                        if ko_a in ko2spe:
                            ko_tmp_spe += ko2spe[ko_a]
                    ko_tmp = list(set(ko_tmp_spe))
                if pathway in self.map_dict:
                    for ko in self.map_dict[pathway].keys():  # 对背景色中的所有项进行循环
                        self.logger.info(ko)
                        if ko == "":
                            continue
                        if ko in ko_tmp:  # 当背景色在富集得出的表中时
                            lnk += ko + "+{},{}%0d%0a".format(self.map_dict[pathway][ko], "blue")  # 有背景色，前景色为蓝
                        else:
                            lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                    # tmp[9] = lnk
                    # str_ = "\t".join(tmp) + "\n"
                    # fw.write(str_)
                else:
                    for ko in ko_tmp:
                        if ko == "":
                            continue
                        lnk += ko + "+{},{}%0d%0a".format("white", "blue")  # 只有背景色
                tmp[9] = lnk
                str_ = "\t".join(tmp) + "\n"
                fw.write(str_)
        if os.path.exists(self.output_name + "_bak"):
            os.remove(self.output_name + "_bak")

        shutil.copyfile(self.output_name, self.output_name + "_bak")
        os.remove(self.output_name)
        shutil.copyfile(self.output_name + "_new", self.output_name)
        self.set_output(self.output_name)

    def set_output(self, linkfile):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        try:
            os.link(linkfile, self.output_dir + '/{}'.format(linkfile))
            self.logger.info("设置kegg富集分析结果目录成功")
        except Exception as e:
            self.logger.info("设置kegg富集分析结果目录失败{}".format(e))
            self.set_error("设置kegg富集分析结果目录失败%s", variables = (e), code = "32503606")
