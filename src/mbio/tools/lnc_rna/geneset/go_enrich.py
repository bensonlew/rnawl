# -*- coding: utf-8 -*-
# __author__ = "shenghe"
import os
import traceback
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.lnc_rna.geneset.go_graph import draw_GO
import subprocess
import unittest
import pandas as pd
from mbio.packages.rna.annot_config import AnnotConfig


class GoEnrichAgent(Agent):
    """
    version v1.0
    author: hesheng
    last_modify: 2016.08.10
    """
    def __init__(self, parent):
        super(GoEnrichAgent, self).__init__(parent)
        options = [
            dict(name="diff_list", type="infile", format="lnc_rna.gene_list"),
            dict(name="go_list", type="infile", format="lnc_rna.go_list"),
            dict(name="pval", type="string", default="0.05"),
            dict(name="go_version", type="string", default="2018"),
            dict(name="method", type="string", default="bonferroni"),
            ]
        self.add_option(options)
        self.step.add_steps("goenrich")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.goenrich.start()
        self.step.update()

    def stepfinish(self):
        self.step.goenrich.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("diff_list").is_set:
            raise OptionError("缺少输入文件:差异基因名称文件" , code = "33706101")
        # if not self.option("all_list").is_set: # edited by shijin
        #     raise OptionError("缺少输入文件:全部基因名称文件")
        if not self.option("go_list").is_set:
            raise OptionError("缺少输入文件:差异基因对应的go_id", code = "33706102")
        method = self.option('method').lower()
        if method == 'bonferroni':
            self.option('method', 'sm_bonferroni')
        elif method == 'bh':
            self.option('method', 'fdr_bh')
        elif method == 'by':
            self.option('method', 'fdr_by')
        elif method == 'holm':
            self.option('method', 'sm_holm')

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
            ])
        result_dir.add_regexp_rules([
            [r"go_enrich_.*", "xls", "go富集结果文件"],
            [r"go_lineage.*", "png", "go富集有向无环图"],
            [r"go_lineage.*", "pdf", "go富集有向无环图"],
            [r"go_lineage.*", "svg", "go富集有向无环图"],
            ])
        super(GoEnrichAgent, self).end()


class GoEnrichTool(Tool):
    """
    """
    def __init__(self, config):
        super(GoEnrichTool, self).__init__(config)
        self.goatools_path = '/bioinfo/annotation/goatools-0.6.5-shenghe'
        self.go_enrich_path = self.goatools_path + '/scripts/find_enrichment.py'
        if self.option("go_version") == "2018":
            self.obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go-basic']
        else:
            self.obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']
        # self.obo = self.config.SOFTWARE_DIR + '/database/GO/go-basic.obo'
        self.set_environ(PYTHONPATH=self.config.SOFTWARE_DIR + self.goatools_path)
        self.python_path = 'miniconda2/bin/python'
        self.out_enrich_fp = self.output_dir + '/go_enrich_' + os.path.splitext(os.path.basename(self.option('diff_list').path))[0] + '.xls'
        self.out_go_graph = self.output_dir + '/go_lineage'
        self.image_magick_path = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/"
        self.out_adjust_graph = self.output_dir + '/adjust_lineage'
        self.class_code_dict = {}
        self.set_environ(FONTCONFIG_PATH=self.config.SOFTWARE_DIR + '/library/fontconfig-2.13.1/etc/fonts')

    def check_list(self):
        """
        去除diff_list中没有注释信息的数据
        new_file_name为在work_dir中生成的新diff_list文件的绝对路径
        :return:
        """
        file1 = self.option("diff_list").prop["path"]
        file2 = self.work_dir + "/all.list"
        f1 = open(file1, "r")
        f2 = open(file2, "r")
        lst_1 = f1.readlines()
        lst_2 = f2.readlines()
        f1.close()
        f2.close()
        new_file_name = self.work_dir + "/" + os.path.basename(file1)
        diff_list = []
        with open(new_file_name, "w") as fw:
            for item in lst_1:
                if item in lst_2:
                    fw.write(item)
                    diff_list.append(item)
        if len(diff_list) == 0:
            self.set_error('该基因集没有注释，请重新选择')
        return new_file_name

    # def get_geneset_type(self, diff_file):
    #     """
    #     查看基因集是已知基因集还是已知基因+新基因 刘彬旭
    #     """
    #     gene_set_f = open(diff_file, 'r')
    #     for line in gene_set_f.readlines():
    #         if line.startswith('MSTR') or line.startswith('TCON') or line.startswith('XLOC'):
    #             gene_set_f.close()
    #             return 'all'
    #     gene_set_f.close()
    #     return 'known'

    # def choose_known_background(self):
    #     """
    #     如果基因集中仅含有已知基因，修改背景注释为已知基因注释 刘彬旭
    #     """
    #     go_list_file = self.option("go_list").prop["path"]
    #     new_go_list = self.work_dir + "/known_go.list"
    #     with open(new_go_list, "w") as ngfw, open(go_list_file, "r") as go_list_f:
    #         for line in go_list_f.readlines():
    #             if line.startswith('MSTR') or line.startswith('TCON') or line.startswith('XLOC'):
    #                 pass
    #             else:
    #                 ngfw.write(line)
    #         return new_go_list

    def run_enrich(self):
        back_ground = self.option('go_list').prop["path"]
        # if self.get_geneset_type(self.option("diff_list").prop["path"]) == 'known':
        #     back_ground = self.choose_known_background()
        cmd0 = "less {}| cut -f1 > {}/all.list".format(back_ground, self.work_dir)
        os.system(cmd0)
        new_file_name = self.check_list()  # edited by shijin 除去背景中不存在的基因
        cmd = self.python_path + ' ' + self.config.SOFTWARE_DIR + self.go_enrich_path + ' '
        cmd = cmd + new_file_name + ' ' + self.work_dir + "/all.list" + ' ' + back_ground
        cmd = cmd + ' --pval ' + self.option('pval') + ' --indent' + ' --method ' + self.option('method') + ' --outfile ' + self.out_enrich_fp
        cmd = cmd + ' --obo ' + self.obo
        command = self.add_command('go_enrich', cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.run_draw_go_graph()
        else:
            self.set_error('goatools计算错误', code = "33706103")

    def run_draw_go_graph(self):
        try:
            self.logger.info('run_draw_go_graph')
            go_pvalue,go_padjust = self.get_go_pvalue_dict()
            self.logger.info('rrrrrrrrrrrrrrrrun_draw_go_graph')
            self.logger.info(go_pvalue)
            self.logger.info('run_draw_go_graphhhhhhhhhhhhhhhhh')
            if go_pvalue:
                # cmd = self.image_magick_path + "convert {} {}".format(self.out_go_graph + ".png", self.out_go_graph + ".pdf")
                draw_GO(go_pvalue, out=self.out_go_graph, obo=self.obo)
                # subprocess.check_output(cmd, shell=True)
            if go_padjust:
                # cmd = self.image_magick_path + "convert {} {}".format(self.out_adjust_graph + ".png", self.out_adjust_graph + ".pdf")
                draw_GO(go_padjust, out=self.out_adjust_graph, obo=self.obo)
                # subprocess.check_output(cmd, shell=True)
            self.end()
        except Exception:
            self.set_error('绘图发生错误:\n%s', variables = (traceback.format_exc()), code = "33706104")

    def get_go_pvalue_dict(self):
        # go2pvalue = {}
        # go2padjust = {}
        # with open(self.out_enrich_fp) as f:
        #     f.readline()
        #     for line in f:
        #         line_sp = line.split('\t')
        #         if line_sp[2] == 'p':
        #             continue
        #         p_bonferroni = float(line_sp[6])
        #         padjust = float(line_sp[9])
        #         go2padjust[line_sp[0]] = padjust
        #         go2pvalue[line_sp[0]] = p_bonferroni
        # tar = sorted(go2pvalue.items(), key=lambda e:e[1], reverse=True)
        # tar_adjust = sorted(go2padjust.items(), key=lambda e:e[1], reverse=True)
        # new_go2padjust = dict(tar_adjust[-10:])
        # new_go2pvalue = dict(tar[-10:])
        # self.logger.info(new_go2pvalue)
        # return new_go2pvalue, new_go2padjust
        go2pvalue = {}
        go2padjust = {}
        # with open(self.out_enrich_fp) as f:
        #     f.readline()
        #     num = 0
        #     for line in f:
        #         num += 1
        #         if num > 10:
        #             break
        #         else:
        #             line_sp = line.split('\t')
        #             p_bonferroni = float(line_sp[6])
        #             padjust = float(line_sp[9])
        #             go2padjust[line_sp[0]] = padjust
        #             go2pvalue[line_sp[0]] = p_bonferroni


        # tar = sorted(go2pvalue.items(), key=lambda e:e[1], reverse=True)
        # tar_adjust = sorted(go2padjust.items(), key=lambda e:e[1], reverse=True)
        # new_go2padjust = dict(tar_adjust[-10:])
        # new_go2pvalue = dict(tar[-10:])
        df = pd.read_table(self.out_enrich_fp, sep='\t', index_col=0)
        df = df[df['enrichment'] != 'p']
        df_10 = df.sort_values(by=['p_uncorrected', 'p_fdr_bh'], ascending=True)[0:10]
        for i in df_10.index:
            go2pvalue[i] = df_10.loc[i]['p_uncorrected']
            go2padjust[i] = df_10.loc[i]['p_fdr_bh']
        new_go2padjust = go2padjust
        new_go2pvalue = go2pvalue
        self.logger.info(new_go2pvalue)
        return new_go2pvalue, new_go2padjust

    def run(self):
        super(GoEnrichTool, self).run()
        self.run_enrich()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'
            data = {
                "id": "GoEnrich" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "denovo_rna_v2.go_enrich",
                "instant": False,
                "options": dict(
                    diff_list=test_dir + "/go_enrich/" + "B_vs_C.deseq2.DE.list",
                    # all_list=test_dir + "/go_enrich/" + "all_gene_list",
                    go_list=test_dir + '/go_enrich/' + "gene_gos.list",
                    pval="0.05",
                    method='bonferroni',
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
