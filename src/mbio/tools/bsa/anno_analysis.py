# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class AnnoAnalysisAgent(Agent):
    """
    bsa分析注释分析
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2018.02.26
    """
    def __init__(self, parent):
        super(AnnoAnalysisAgent, self).__init__(parent)
        options = [
            {"name": "anno_stat", "type": "infile", "format": "bsa.vcf"},  # region.threshold.gene.*.stat文件
            {"name": "anno_type", "type": "string"},  # kegg时需要跑函数，得到pathway的pdf和png文件
        ]
        self.add_option(options)
        self.step.add_steps("anno_analysis")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.anno_analysis.start()
        self.step.update()

    def stepfinish(self):
        self.step.anno_analysis.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('anno_stat'):
            raise OptionError('必须输入各注释统计表',code = "31500101")
        if not self.option('anno_type'):
            raise OptionError('必须输入注释类型',code = "31500102")
        if self.option('anno_type') not in ["kegg", "go", "eggnog"]:
            raise OptionError('注释类型只能是"kegg,go,eggnog"中的一种',code = "31500103")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "4G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AnnoAnalysisAgent, self).end()


class AnnoAnalysisTool(Tool):
    def __init__(self, config):
        super(AnnoAnalysisTool, self).__init__(config)
        self._version = "v1.0.1"
        self.R_path = 'program/R-3.3.1/bin/Rscript '
        self.stat_path = self.config.PACKAGE_DIR + '/bsa/eff-enrich.R'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib')
        self.image_magick = self.config.SOFTWARE_DIR + "/program/ImageMagick/bin/convert"
        self.ko_path = self.config.SOFTWARE_DIR + '/database/KEGG/kegg_2017-05-01/kegg/pathway/ko/'
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoAnalysisTool, self).run()
        with open(self.option('anno_stat').prop['path'])as fr:   # add by wzy 20180329  上游结果为空，工作流修改复杂，在此处修改
            lines = fr.readlines()
            if len(lines) == 0:
                self.direct_output(self.option('anno_type'))  # 判断输入文件是否为空，为空则直接输出空文件，不再运行tool
                self.end()
            else:
                self.run_anno_analysis()
                if self.option('anno_type') == 'kegg':
                    self.get_picture()
                self.set_output()
                self.end()

    def run_anno_analysis(self):
        """
        运行eff-enrich.R,计算滑窗
        """
        base_name = os.path.basename(self.option('anno_stat').prop['path']).split('.stat')[0]
        if self.option('anno_type') == 'eggnog':
            cmd = self.R_path + self.stat_path + ' --input {} --output {} --eggnog'.format(
                self.option('anno_stat').prop['path'], self.work_dir + '/' + base_name + '.final.stat')
        else:
            cmd = self.R_path + self.stat_path + ' --input {} --output {} --top 1'.format(
                self.option('anno_stat').prop['path'], self.work_dir + '/' + base_name + '.final.stat')
        self.logger.info('运行eff-enrich.R,注释统计')
        command = self.add_command("anno_analysis_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("anno_analysis运行完成")
        else:
            self.set_error("anno_analysis运行出错!", code="31500101")

    def get_picture(self):
        """
        根据ko号找到对应的所有ko的pdf和png文件-- 添加如果存在pathway_dir就将该文件夹删除 add by hongdong 20180301
        :return:
        """
        # if os.path.isdir(self.output_dir + '/pathway_dir'):
        #     os.system('rm -r %s' % os.path.join(self.output_dir, '/pathway_dir'))
        cmd_list = []
        os.system("mkdir {}".format(self.output_dir + '/pathway_dir'))
        all_ko = os.listdir(self.ko_path)
        with open(self.option('anno_stat').prop['path'])as fr:
            for line in fr:
                tmp = line.strip().split('\t')
                ko_png = tmp[0] + '.png'
                if ko_png in all_ko:
                    ko_path = os.path.join(self.ko_path, ko_png)
                    os.link(ko_path, self.output_dir + '/pathway_dir/' + ko_png)
                    pdf_path = os.path.join(self.output_dir + '/pathway_dir', tmp[0] + '.pdf')
                    cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + ko_path + ' ' + pdf_path
                    cmd_list.append(cmd)
                else:
                    self.set_error("ko %s png文件在本地数据库 %s 中没有，请核实!",variables = (tmp[0], self.ko_path), code="31500108")
        with open(self.work_dir + "/cmd.list", "w") as fw:
            for i in range(len(cmd_list)):
                fw.write(cmd_list[i] + "\n")
        cmd = self.parafly + " -c {} -CPU 10".format(self.work_dir + "/cmd.list")
        cmd_obj = self.add_command("cmd_list", cmd).run()
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info("cmd list执行成功")
        else:
            self.set_error("cmd list运行出错!", code="31500103")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面  --这里加上output中文件已经存在要删除的判断 modified by HONGDONG 20180301
        :return:
        """
        # for root, dirs, files in os.walk(self.output_dir):
        #     for names in files:
        #         os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        try:
            base_name = os.path.basename(self.option('anno_stat').prop['path']).split('.stat')[0]
            os.link(self.work_dir + '/' + base_name + '.final.stat.detail', self.output_dir + '/' + base_name + '.final.stat.detail')
            self.logger.info("设置anno_analysis分析结果目录成功")

        except Exception as e:
            self.logger.info("设置anno_analysis分析结果目录失败{}".format(e))
            self.set_error("设置anno_analysis分析结果目录失败%s ",variables=(e), code="31500104")
 
    def direct_output(self, anno_type):   # add by wzy 20180329
        """
        判断输入文件是否为空，为空则直接输出空文件，不再运行tool
        """
        self.logger.info("输入文件为空，直接输出空结果文件")
        try:
            base_name = os.path.basename(self.option('anno_stat').prop['path']).split('.stat')[0]
            os.link(self.option('anno_stat').prop['path'], self.output_dir + '/' + base_name + '.final.stat.detail')
            if anno_type == 'kegg':
                os.system("mkdir {}".format(self.output_dir + '/pathway_dir'))
            self.logger.info("设置anno_analysis分析结果目录成功")

        except Exception as e:
            self.logger.info("设置anno_analysis分析结果目录失败{}".format(e))
            self.set_error("设置anno_analysis分析结果目录失败%s", variables=(e), code="31500105")