# -*- coding: utf-8 -*-
# __author__ = 'linna'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class PcaAgent(Agent):
    def __init__(self, parent):
        super(PcaAgent, self).__init__(parent)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "specimen_name", "type": "string", "default": "column"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "scale", "type": "string", "default": "True"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('tooltable').prop['sample_num'] < 3:
            raise OptionError('列数少于3，不可进行分析', code="32702904")
        # if self.option('group_table').is_set:
        #     if self.option('specimen_name') == 'column':
        #         for i in self.option('group_table').prop['sample_name']:
        #             if i not in self.option('tooltable').prop['col_sample']:
        #                 raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702909")
        #     else:
        #         for i in self.option('group_table').prop['sample_name']:
        #             if i not in self.option('tooltable').prop['row_sample']:
        #                 raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702910")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 5
        self._memory = '10G'

    def end(self):
        if self.option('group_table').is_set:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
                ["./pca_importance.xls", "xls", "主成分解释度表"],
                ["./pca_rotation.xls", "xls", "主成分贡献度表"],
                ["./pca_rotation_all.xls", "xls", "主成分贡献度表全部"],
                ["./pca_sites.xls", "xls", "样本坐标表"],
                ["./group.xls", "xls", "分组表"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
                ["./pca_importance.xls", "xls", "主成分解释度表"],
                ["./pca_rotation.xls", "xls", "主成分贡献度表"],
                ["./pca_rotation_all.xls", "xls", "主成分贡献度表全部"],
                ["./pca_sites.xls", "xls", "样本坐标表"],
            ])
        super(PcaAgent, self).end()


class PcaTool(Tool):
    def __init__(self, config):
        super(PcaTool, self).__init__(config)
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.cmd_path = self.config.PACKAGE_DIR + '/statistical/ordination.pl'
        self.script_path = "bioinfo/meta/scripts/beta_diver.sh"
        self.R_path = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.1/bin/R')
        self.R_path1 = '/program/R-3.3.1/bin/'
        self.ellipse_path = self.config.PACKAGE_DIR + '/graph/scripts/Ellipse.R'

    def run(self):
        """
        运行
        """
        super(PcaTool, self).run()
        self.run_ordination()
        self.linkfile()
        # self.set_db()//tool里面不导表，在工作流中导表
        self.end()

    def formattable(self, tablepath):
        """
        转置表格
        """
        imp_file = tablepath
        imp_link = self.work_dir + '/input.txt'
        if os.path.exists(imp_link):
            os.remove(imp_link)
        os.link(imp_file, imp_link)
        this_table = imp_link
        if self.option('specimen_name') != 'column':
            newtable = this_table + '.T'
            self.t_table(this_table, newtable)
            return newtable
        else:
            return this_table

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def run_ordination(self, cmd1='cmd', cmd2='pca'):
        """
        运行ordination.pl计算pca
        """
        if self.option('group_table').is_set:
            group_input = self.option('group_table').path
        else:
            if self.option('specimen_name') == 'column':
                in_table = self.option('tooltable').path + '.T'
                self.t_table(self.option('tooltable').path, in_table)
            else:
                in_table = self.option('tooltable').path
            dis_table = pd.read_table(in_table)
            # 这里没有分组时，调用的脚本不会画椭圆，所以创建一个分组文件，每个样本的分组名都为All
            group_table = pd.DataFrame()
            group_table["Sample_ID"] = dis_table[list(dis_table.columns)[0]]
            group_table["group"] = "All"
            group_table.to_csv(self.work_dir + "/sample.xls", index=0, sep="\t")
            group_input = self.work_dir + "/sample.xls"
        if self.option('scale') == "True":
            # 是否标准化参数传给pl脚本值不对导致不能标准化，这里添加处理  modified by wuqin on 20211027
            scale = "T"
        else:
            scale = "F"
        cmd = '{} {} -type pca -community {} -outdir {} -scale {} -group {}'.format(
            self.perl_path, self.cmd_path, self.formattable(self.option('tooltable').path), self.work_dir,
            scale, group_input)
        self.logger.info('运行ordination.pl程序')
        self.logger.info(cmd)
        cmd1 = self.add_command(cmd1, cmd).run()  # 命令1运行
        self.wait(cmd1)
        if cmd1.return_code == 0:
            self.logger.info('生成 cmd.r 文件成功')
        else:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件', code="32702902")
        cmd_ = self.script_path + ' %s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_)
        cmd2 = self.add_command(cmd2, cmd_).run()  # 命令2 运行
        self.wait(cmd2)
        if cmd2.return_code == 0:
            self.logger.info('pca计算成功')
        else:
            self.logger.info('pca计算失败')
            self.set_error('R程序计算pca失败', code="32702903")
        self.logger.info('运行ordination.pl程序计算pca完成')

    def linkfile(self):
        if self.option('group_table').is_set:
            group_file = self.option('group_table').path
            group_link = self.output_dir + '/group.xls'
            if os.path.exists(group_link):
                os.remove(group_link)
            os.link(group_file, group_link)
            basename = os.path.basename(self.work_dir + '/input.txt')
            basename1 = basename.replace('.', '_')
            imp_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_importance.xls'
            if os.path.exists(imp_file):
                pass
            else:
                imp_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_importance.xls'
            imp_link = self.output_dir + '/pca_importance.xls'
            if os.path.exists(imp_link):
                os.remove(imp_link)
            os.link(imp_file, imp_link)
            rot_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_rotation.xls'
            if os.path.exists(rot_file):
                pass
            else:
                rot_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_rotation.xls'
            rot_link = self.output_dir + '/pca_rotation.xls'
            if os.path.exists(rot_link):
                os.remove(rot_link)
            os.link(rot_file, rot_link)
            all_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_rotation_all.xls'
            if os.path.exists(all_file):
                pass
            else:
                all_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_rotation_all.xls'
            all_link = self.output_dir + '/pca_rotation_all.xls'
            if os.path.exists(all_link):
                os.remove(all_link)
            os.link(all_file, all_link)
            site_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_sites.xls'
            if os.path.exists(site_file):
                pass
            else:
                site_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_sites.xls'
            site_link = self.output_dir + '/pca_sites.xls'
            if os.path.exists(site_link):
                os.remove(site_link)
            os.link(site_file, site_link)
            ellipse_file = self.work_dir + "/pca" + "/" + "ellipse.xls"
            ellipse_link = self.output_dir + '/ellipse.xls'
            if os.path.exists(ellipse_link):
                os.remove(ellipse_link)
            os.link(ellipse_file, ellipse_link)
        else:
            basename = os.path.basename(self.work_dir + '/input.txt')
            basename1 = basename.replace('.', '_')
            imp_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_importance.xls'
            if os.path.exists(imp_file):
                pass
            else:
                imp_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_importance.xls'
            imp_link = self.output_dir + '/pca_importance.xls'
            if os.path.exists(imp_link):
                os.remove(imp_link)
            os.link(imp_file, imp_link)
            rot_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_rotation.xls'
            if os.path.exists(rot_file):
                pass
            else:
                rot_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_rotation.xls'
            rot_link = self.output_dir + '/pca_rotation.xls'
            if os.path.exists(rot_link):
                os.remove(rot_link)
            os.link(rot_file, rot_link)
            all_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_rotation_all.xls'
            if os.path.exists(all_file):
                pass
            else:
                all_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_rotation_all.xls'
            all_link = self.output_dir + '/pca_rotation_all.xls'
            if os.path.exists(all_link):
                os.remove(all_link)
            os.link(all_file, all_link)
            site_file = self.work_dir + '/pca' + '/' + basename1 + '_' + 'pca_sites.xls'
            if os.path.exists(site_file):
                pass
            else:
                site_file = self.work_dir + '/pca' + '/' + basename1 + '_T_' + 'pca_sites.xls'
            site_link = self.output_dir + '/pca_sites.xls'
            if os.path.exists(site_link):
                os.remove(site_link)
            os.link(site_file, site_link)
            ellipse_file = self.work_dir + "/pca" + "/" + "ellipse.xls"
            ellipse_link = self.output_dir + '/ellipse.xls'
            if os.path.exists(ellipse_link):
                os.remove(ellipse_link)
            os.link(ellipse_file, ellipse_link)

    def set_db(self):
        self.logger.info("开始导表")
        api_pca = self.api.api('wgs_v2.pca')
        api_pca.add_pca(self.option('main_id'), self.output_dir)
        self.logger.info("导表结束")