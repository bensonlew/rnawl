# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last modify by shaohua.yuan 20180403
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re


class MgKeggStatAgent(Agent):
    """
    宏基因kegg注释结果统计tool 调用脚本 meta_kegg_stat.py
    author: zhouxuan
    last_modify: 2018.0402
    last_modify_by:haidong.gu
    """

    def __init__(self, parent):
        super(MgKeggStatAgent, self).__init__(parent)
        options = [
            {"name": "kegg_result_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "reads_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "kegg_profile_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},
            {"name": "align_table_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "xml_file_dir", "type": "infile", "format": "align.blast.blast_xml_dir"}  # 画图用
        ]
        self.add_option(options)
        self._memory_increase_step = 70  # 每次重运行增加内存70G by guhaidong @ 20190220

    def check_options(self):
        if not self.option("kegg_result_dir").is_set:
            raise OptionError("必须设置输入文件夹", code="31203001")
        if not self.option("reads_profile").is_set:
            raise OptionError("必须设置基因丰度表", code="31203002")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '30G'  # 内存5G增加到10G  by GHD @20180427
        # tmp_mem = 10 + 20 * self._rerun_time  # 每次因拼接失败而重运行的内存增加20G by GHD @ 20180402
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('mg_kegg_stat use memory : ' + self._memory)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(MgKeggStatAgent, self).end()


class MgKeggStatTool(Tool):
    def __init__(self, config):
        super(MgKeggStatTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "miniconda2/bin/python"
        # self.python_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/meta_kegg_stat.py'
        self.python_script2 = self.config.PACKAGE_DIR + '/annotation/mg_annotation/kegg_pathway_img_v94.py'
        self.perl = '/program/perl-5.24.0/bin/perl'
        self.gene_profile = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/kegg_gene_profile.pl'
        self.ko_profile = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/kegg_ko_profile.pl'
        self.module_profile = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/kegg_module_profile.pl'
        self.enzyme_profile = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/kegg_enzyme_profile.pl'
        self.pathway_profile = self.config.PACKAGE_DIR + '/metagenomic/anno_profile/kegg_pathway_profile.pl'

        self.sh_path = 'bioinfo/align/scripts/cat.sh'
        self.cat_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/cat_seq.sh' ##add by qingchen.zhang@20190529
        self.anno_result = ''
        self.enzyme_list = ''
        self.module_list = ''
        self.pathway_list = ''
        self.html_path = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007/html/"

    def run(self):
        """
        运行
        :return:
        """
        super(MgKeggStatTool, self).run()
        if self.option("align_table_dir").is_set:
            self.merge_align_table()
        self.merge_table()
        self.run_kegg_stat()
        self.run_kegg_img()
        self.set_output()
        self.end()

    def merge_align_table(self):
        """
        合并比对结果文件
        :return:
        """
        profile_file = os.listdir(self.option('align_table_dir').prop['path'])
        self.align_table = os.path.join(self.output_dir, "kegg_align_table.xls")
        if os.path.exists(self.align_table):
            os.remove(self.align_table)
        cmd = '{}'.format(self.cat_path)
        for i in profile_file:
            file_path = os.path.join(self.option('align_table_dir').prop['path'], i)
            #if kegg_align > 1:
            with open(file_path, 'r') as f:
                line = f.readline()
                if re.search(r'^Score', line):
                    os.system("sed -i '/^Score\t/d ' " + file_path)
            cmd += ' ' + file_path
            #cmd = '{} {} {}'.format(self.sh_path, file_path, self.align_table)
        cmd += ' ' + self.align_table ##add by qingchen.zhang@20190529
        self.logger.info("start cat align")
        command_name = "cat_align"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cat align done")
        else:
            self.set_error("cat align error", code="31203001")
            self.set_error("cat align error", code="31203006")


        xml_file_dir = os.listdir(self.option('xml_file_dir').prop['path'])
        self.logger.info(xml_file_dir)
        self.xml_file = os.path.join(self.output_dir, "kegg_merge.xml")
        if os.path.exists(self.xml_file):
            os.remove(self.xml_file)
        cmd_xml = '{}'.format(self.cat_path)
        for i in xml_file_dir:
            file_xml_path = os.path.join(self.option('xml_file_dir').prop['path'], i)
            #cmd2 = '{} {} {}'.format(self.sh_path, file_path, self.xml_file)
            cmd_xml += ' ' + file_xml_path
        cmd_xml += ' ' + self.xml_file ##add by qingchen.zhang@20190529
        command_name = "cat_xml"
        command2 = self.add_command(command_name, cmd_xml).run()
        self.wait(command2)
        if command2.return_code != 0:
            self.set_error("cat_xml error", code="31203002")
        else:
            self.logger.info("{} done".format(command_name))
        os.system("sed -i '/^Score\t/d ' "+ self.align_table)
        os.system("sed -i '1iScore\tE-Value\tHSP-Len\tIdentity-%\tSimilarity-%\tQuery-Name\tQ-Len\tQ-Begin\tQ-End\tQ-Frame\tHit-Name\tHit-Len\tHsp-Begin\tHsp-End\tHsp-Frame\tHit-Description'" + " %s" %(self.align_table))

    def merge_table(self):
        """
        合并注释结果文件
        :return:
        """
        profile_file = os.listdir(self.option('kegg_result_dir').prop['path'])
        self.anno_result = os.path.join(self.work_dir, "tmp_kegg_anno.xls")
        self.enzyme_list = os.path.join(self.work_dir, "kegg_enzyme_list.xls")
        self.module_list = os.path.join(self.work_dir, "kegg_module_list.xls")
        self.pathway_list = os.path.join(self.work_dir, "kegg_pathway_list.xls")
        if os.path.exists(self.anno_result):
            os.remove(self.anno_result)
        if os.path.exists(self.enzyme_list):
            os.remove(self.enzyme_list)
        if os.path.exists(self.module_list):
            os.remove(self.module_list)
        if os.path.exists(self.pathway_list):
            os.remove(self.pathway_list)
        cmd1 = '{}'.format(self.cat_path)
        cmd2 = '{}'.format(self.cat_path)
        cmd3 = '{}'.format(self.cat_path)
        cmd4 = '{}'.format(self.cat_path)
        suffix = ["anno_result", "enzyme_list", "module_list", "pathway_list"]
        #merge_result = [self.anno_result, self.enzyme_list, self.module_list, self.pathway_list]
        for j in range(0, 4):
            if suffix[j] == "anno_result":
                for i in profile_file:
                    if suffix[j] in i:
                        file_path = os.path.join(self.option('kegg_result_dir').prop['path'], i)
                        #if kegg_number > 1:
                        with open(file_path, 'r') as f:
                            line = f.readline()
                            if re.search(r'^#', line):
                                os.system("sed -i '/^#/d ' " + file_path)
                        cmd1 += ' ' + file_path
                cmd1 += ' ' + self.anno_result
                self.logger.info("start cat_anno")
                command_name = "cat_anno"
                command1 = self.add_command(command_name, cmd1).run()
                self.wait(command1)
                if command1.return_code == 0:
                    self.logger.info("cat_anno done")
                else:
                    self.set_error("cat_anno error", code="31203003")
                    self.set_error("cat_anno error", code="31203007")
            elif suffix[j] == "enzyme_list":
                for i in profile_file:
                    if suffix[j] in i:
                        file_path = os.path.join(self.option('kegg_result_dir').prop['path'], i)
                        cmd2 += ' ' + file_path
                cmd2 += ' ' + self.enzyme_list
                self.logger.info("start cat_enzyme")
                command_name = "cat_enzyme"
                command2 = self.add_command(command_name, cmd2).run()
                self.wait(command2)
                if command2.return_code == 0:
                    self.logger.info("cat_enzyme done")
                else:
                    self.set_error("cat_enzyme error", code="31203003")
                    self.set_error("cat_enzyme error", code="31203007")
            elif suffix[j] == "module_list":
                for i in profile_file:
                    if suffix[j] in i:
                        file_path = os.path.join(self.option('kegg_result_dir').prop['path'], i)
                        cmd3 += ' ' + file_path
                cmd3 += ' ' + self.module_list
                self.logger.info("start cat_module")
                command_name = "cat_module"
                command3 = self.add_command(command_name, cmd3).run()
                self.wait(command3)
                if command3.return_code == 0:
                    self.logger.info("cat_module done")
                else:
                    self.set_error("cat_module error", code="31203003")
                    self.set_error("cat_module error", code="31203007")
            elif suffix[j] == "pathway_list":
                for i in profile_file:
                    if suffix[j] in i:
                        file_path = os.path.join(self.option('kegg_result_dir').prop['path'], i)
                        cmd4 += ' ' + file_path
                cmd4 += ' ' + self.pathway_list
                self.logger.info("start cat_pathway")
                command_name = "cat_pathway"
                command4 = self.add_command(command_name, cmd4).run()
                self.wait(command4)
                if command4.return_code == 0:
                    self.logger.info("cat_pathway done")
                else:
                    self.set_error("cat_pathway error", code="31203003")
                    self.set_error("cat_pathway error", code="31203007")
        os.system("sed -i '/^#/d ' "+ self.anno_result)
        os.system("sed -i '1i#Query\tGene\tKO\tDefinition\tPathway\tEnzyme\tModule\tHyperlink\tIdentity(%)\tAlign_len' " + self.anno_result)

    def run_kegg_stat(self):
        kegg_anno = self.anno_result
        enzyme_list = self.enzyme_list
        module_list = self.module_list
        pathway_list = self.pathway_list
        ko_p = self.perl + " {} {} {} {}".format(self.ko_profile, kegg_anno,
                                                 self.option("reads_profile").path,
                                                 self.output_dir)
        gene_p = self.perl + " {} {} {} {}".format(self.gene_profile, kegg_anno,
                                                   self.option("reads_profile").path,
                                                   self.output_dir)
        enzyme_p = self.perl + " {} {} {} {} {}".format(self.enzyme_profile, kegg_anno,
                                                        self.option("reads_profile").path,
                                                        self.output_dir, enzyme_list)
        module_p = self.perl + " {} {} {} {} {}".format(self.module_profile, kegg_anno,
                                                        self.option("reads_profile").path,
                                                        self.output_dir, module_list)
        pathway_p = self.perl + " {} {} {} {} {}".format(self.pathway_profile, kegg_anno,
                                                         self.option("reads_profile").path,
                                                         self.output_dir, pathway_list)
        all_cmds = []
        all_cmds.append(self.add_command("ko_p", ko_p, ignore_error=True).run())
        all_cmds.append(self.add_command("gene_p", gene_p, ignore_error=True).run())
        all_cmds.append(self.add_command("enzyme_p", enzyme_p, ignore_error=True).run())
        all_cmds.append(self.add_command("module_p", module_p, ignore_error=True).run())
        all_cmds.append(self.add_command("pathway_p", pathway_p, ignore_error=True).run())
        self.wait(*all_cmds)
        errs = []
        for cmd in all_cmds:
            if cmd.return_code != 0:
                errs.append(cmd.name)
        if errs:
            self.set_error("err in {}".format(errs))
        else:
            self.logger.info("kegg_stat succeed")

    def _run_kegg_stat(self):
        """
        统计各功能水平丰度
        :return:
        """
        kegg_anno = self.anno_result
        enzyme_list = self.enzyme_list
        module_list = self.module_list
        pathway_list = self.pathway_list
        cmd = self.python_path + ' {} -k {} -e {} -p {}  -m {} -r {} -o {} '. \
            format(self.python_script, kegg_anno, enzyme_list, pathway_list, module_list,
                   self.option('reads_profile').prop['path'], self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('kegg_stat', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg_stat succeed")
        elif command.return_code in [1, -9]:  # 内存超出是返回值为1,增加-9 by ghd @ 20190129
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("kegg_stat failed", code="31203004")
            self.set_error("kegg_stat failed", code="31203008")

    def run_kegg_img(self):
        """
        kegg注释得到pathway通路图
        :return:
        """
        self.logger.info("start output kegg img")
        # xml_file = self.xml_file
        pathway_file = self.output_dir + "/kegg_pathway_eachmap.xls"
        # cmd2 = "{} {} -i {} -o {} -p {} -ko {} -KO {} -png_file {} -html {}".format(self.python_path, self.python_script2, xml_file, self.output_dir, pathway_file, "#Pathway", "KO_list", "True", self.html_path)
        cmd2 = "{} {} -o {} -p {} -ko {} -KO {} -png_file {} -html {}".format(self.python_path, self.python_script2, self.output_dir, pathway_file, "#Pathway", "KO_list", "True", self.html_path)
        self.logger.info(cmd2)
        command = self.add_command("output_kegg_pathway_img", cmd2).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("output_kegg_pathway_img succeed")
            os.system("mv {}/pathway_img {}".format(self.output_dir, self.work_dir))
            os.system("tar zcf {}/pathway_img.tar.gz pathway_img".format(self.output_dir))
        else:
            self.set_error("output kegg pathway img failed", code="31203005")

    def set_output(self):
        self.logger.info("set_output")
        for f in os.listdir(self.output_dir):  # 删除sed的中间文件
            if f.startswith('sed'):
                fp = os.path.join(self.output_dir, f)
                os.system("rm -f " + fp)
        try:
            self.option("kegg_profile_dir", self.output_dir)
        except Exception as e:
            self.set_error("SET_OUTFILE FAILED %s", variables=(e), code="31203009")
        self.logger.info("OUTPUT RIGHT")
