# -*- coding: utf-8 -*-
import os
# import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class EllipseAgent(Agent):
    """
    转录组画置信椭圆
    version 1.0
    author: fuwenyao
    需要R软件
    """
    def __init__(self, parent):
        super(EllipseAgent, self).__init__(parent)
        options = [
            {"name": "pc_table", "type": "string"},
            # {"name": "group_table", "type": "infile", "format": " meta.otu.group_table, toolapps.group_table"},  # 输入的group表格
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},            # 输入的group表格
            {"name": "level", "type": "string", "default": "0.95"},  # 置信度
            {"name": "out", "type": "string", "default": "ellipse_out.xls"},
            {"name": "group_id", "type": "string", "default":''},
            {"name": "analysis", "type": "string", "default":"default"},
            {"name": "meta", "type": "string", "default": ""},  # by houshuang 20190924 用于判断是否为多样性分析
            {"name": "draw_all_sample", "type":"string", "default":"f"}  #画所有样本的椭圆
        ]
        self.add_option(options)


    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("pc_table"):
            raise OptionError("参数pc_table不能为空", code="32301301")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '5G'


class EllipseTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(EllipseTool, self).__init__(config)
        self.R_path = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/'
        self.ellipse_path = self.config.PACKAGE_DIR + '/graph/scripts/Ellipse.R'

    def create_ellipse_table(self):
        """
        调用脚本
        """
        self.logger.info('运行ellipse_cmd')

        infile = self.option('pc_table')
        if self.option('analysis') == 'rda_cca':
            infile = infile + '/rda_sites.xls'
            if not os.path.exists(infile):
                infile = self.option('pc_table') + '/cca_sites.xls'
        if self.option('analysis') == 'nmds':
            new_file = self.work_dir + '/nmds_new.xls'

            with open(infile) as f:
                lines = f.readlines()
                head = lines[0].replace('MDS','NMDS')
                with open(new_file, 'w') as fw:
                    fw.write(head)
                    fw.write(''.join(lines[1:]))
            infile = new_file

        # by houshuang 20190921 >>>
        if not os.path.exists(infile):
            self.set_error("文件不存在:%s", variables=(infile), code="32301302")
        # <<<

        self.outfile = self.work_dir + '/' + self.option('out')
        group_table = ''
        if self.option('group_table') and self.option('group_id') not in ['all', 'ALL', 'All']:
            group_table = self.option('group_table').prop['path']
            with open(self.option('group_table').prop['path']) as f:
                head =f.readline()
                group_detail = {}
                for line in f:
                    line = line.strip()
                    spline = line.split('\t')
                    if spline[1] not in group_detail.keys():
                        group_detail[spline[1]] = [spline[0]]
                    else:
                        group_detail[spline[1]].append(spline[0])
                small_g = {}
                small_list = []
                for k in group_detail.keys():
                    if len(group_detail[k]) < 3:
                        small_g[k] = (group_detail[k])
                        small_list.extend(group_detail[k])
                if small_g != {}:
                    new_group = self.work_dir + '/big_group.xls'
                    with open(new_group,'w') as fw:
                        fw.write(head)
                        for k in group_detail.keys():
                            if k not in small_g.keys():
                                fw.write('\n'.join([i +'\t'+ k for i in group_detail[k]])+'\n')
                    group_table = new_group
                    new_infile = self.work_dir + '/big_group_sites.xls'
                    w_site =open(new_infile, 'w')

                    with open(infile,'r') as f:
                        head = f.readline()
                        w_site.write(head)
                        for line in f:
                            sp = line.split('\t',1)
                            if sp[0] not in small_list:
                                w_site.write(line)
                    infile = new_infile
                    w_site.close()

        cmd = self.R_path + 'Rscript {} -f {} -l {} -o {} '.format(self.ellipse_path, infile, self.option('level'), self.outfile)
        if self.option('meta') != "":
            cmd = self.R_path + 'Rscript {} -f {} -l {} -o {} -m {}'.format(self.ellipse_path, infile, self.option('level'), self.outfile, self.option('meta'))
        if group_table != '':
            lines = open(group_table).readlines()
            open(group_table, 'w').writelines('{}\n'.format(line.strip()) for line in lines)
            cmd += ' -g {}'.format(group_table)
        try:## add by qingchen.zhang ## 目的是为了如果计算不出来置信椭圆，就不进行导表 @ 20201202
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行ellipse_cmd完成")
        except:
            self.logger.info("运行ellipse_cmd出错")

        # command = self.add_command("ellipse_cmd", cmd)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #     self.logger.info("运行ellipse_cmd完成")
        # else:
        #     self.set_error("运行ellipse_cmd出错!", code="32301301")


        if self.option("draw_all_sample")=='t' and  self.option('group_table').is_set:
            cmd = self.R_path + 'Rscript {} -f {} -l {} -o {} '.format(self.ellipse_path, infile, self.option('level'), 'all_sample_ellipse.xls')
            try: ## add by qingchen.zhang ## 目的是为了如果计算不出来置信椭圆，就不进行导表@ 20201202
                self.logger.info(cmd)
                subprocess.check_output(cmd, shell=True)
                self.logger.info("运行ellipse_cmd2完成")
            except:
                self.logger.info("运行ellipse_cmd2出错")


            # command2 = self.add_command("ellipse_cmd2", cmd)
            # command2.run()
            # self.wait(command2)
            # if command2.return_code == 0:
            #     self.logger.info("运行ellipse_cmd2完成")
            # else:
            #     self.set_error("运行ellipse_cmd2出错!")


    def set_output(self):
        """
        将结果文件链接至output
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        if os.path.exists(self.outfile):
            os.link(self.outfile, self.output_dir + '/' + self.option('out'))
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(EllipseTool, self).run()
        self.logger.info("开始运行！")
        self.create_ellipse_table()
        self.logger.info("set out put")
        self.set_output()
        self.end()
