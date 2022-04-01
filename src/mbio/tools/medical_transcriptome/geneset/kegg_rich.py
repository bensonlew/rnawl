# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
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
            {"name": "kegg_table", "type": "infile", "format": "medical_transcriptome.kegg_table"},
            # 只含有基因的kegg table结果文件
            {"name": "diff_list", "type": "infile", "format": "denovo_rna_v2.gene_list"},
            {"name": "kegg_table2", "type": "string"},
            # gene名字文件
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "add_info", "type": "string", "default": None},  # 底图文件地址
            {'name': 'kegg_version', 'type': 'string', 'default': '202007'}
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
            raise OptionError('必须设置kegg的pathway输入文件', code = "33706601")
        if self.option('correct').lower() not in ['by', 'bh', 'bonferroni', 'holm']:
            raise OptionError('多重检验校正的方法不在提供的范围内', code = "33706602")
        if not self.option("diff_list").is_set:
            raise OptionError("必须设置输入文件diff_list", code = "33706603")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
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
        self.script_path = self.config.PACKAGE_DIR + "/ref_rna_v2/kegg_enrichment.py"
        self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))
        self.k2e = self.kegg_files_dict['K2enzyme.tab']
        self.brite = self.kegg_files_dict['br08901.txt']
        self.k2e = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/K2enzyme.tab"
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
        self.end()

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
            self.set_error("kegg富集第一步运行出错:%s", variables = (e), code = "33706604")

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
        lst_2 = [x.split('\t')[0] for x in f2.readlines()]
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
            self.set_error("kegg富集分析运行出错!", code = "33706605")


    def run_add_info(self):
        background_links = self.option("add_info")
        if not os.path.exists(background_links):
            self.logger.info("不存在背景信息文件，程序退出")
            return
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
                if pathway in self.map_dict:
                    for ko in self.map_dict[pathway].keys():  # 对背景色中的所有项进行循环
                        self.logger.info(ko)
                        if ko == "":
                            continue
                        if ko in ko_tmp:  # 当背景色在富集得出的表中时
                            lnk += ko + "+{},{}%0d%0a".format(self.map_dict[pathway][ko], "red")  # 有背景色，前景色为蓝
                        else:
                            lnk += ko + "+{}%0d".format(self.map_dict[pathway][ko])  # 只有背景色
                    # tmp[9] = lnk
                    # str_ = "\t".join(tmp) + "\n"
                    # fw.write(str_)
                else:
                    for ko in ko_tmp:
                        if ko == "":
                            continue
                        lnk += ko + "+{},{}%0d%0a".format("white", "red")  # 只有背景色
                tmp[9] = lnk
                str_ = "\t".join(tmp) + "\n"
                fw.write(str_)
        os.link(self.output_name, self.output_name + "_bak")
        os.remove(self.output_name)
        os.link(self.output_name + "_new", self.output_name)
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
            self.set_error("设置kegg富集分析结果目录失败%s", variables = (e), code = "33706606")
