# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuewang'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from biocluster.core.exceptions import OptionError


class AnosimBoxAgent(Agent):
    """
    version 1.0
    author zhaoyuewang
    last_modified:2017.01.18
    """

    def __init__(self, parent):
        super(AnosimBoxAgent, self).__init__(parent)
        options = [
            {"name": "dis_matrix", "type": "infile",
                "format": "meta.beta_diversity.distance_matrix"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "grouplab", "type": "string", "default": ""},
            {"name": "permutations", "type": "int", "default": 999},
        ]
        self.add_option(options)
        self.step.add_steps('distancebox')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.distancebox.start()
        self.step.update()

    def step_end(self):
        self.step.distancebox.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if 10000 >= self.option('permutations') >= 10:
            pass
        else:
            raise OptionError('随机置换次数:%s不在正常范围内[10, 10000]', variables=self.option('permutations'), code="32701601")
        if not self.option('dis_matrix').is_set:
            raise OptionError('必须提供距离矩阵文件', code="32701602")
        else:
            self.option('dis_matrix').get_info()
            samplelist = self.option('dis_matrix').prop['samp_list']
        if not self.option('group').is_set:
            raise OptionError('必须提供分组信息文件', code="32701603")
        else:
            self.option('group').get_info()
            if self.option('grouplab'):
                if self.option('grouplab') not in self.option('group').prop['group_scheme']:
                    raise OptionError('选定的分组方案名:%s在分组文件中不存在', variables=self.option('grouplab'), code="32701604")
            else:
                pass
            if len(samplelist) < len(self.option('group').prop['sample']):
                raise OptionError('分组文件中样本数量:%s多于距离矩阵中的样本数量:%s', variables=(len(self.option('group').prop['sample']),
                                  len(samplelist)), code="32701605")
            for sample in self.option('group').prop['sample']:
                if sample not in samplelist:
                    raise OptionError('分组文件的样本(%s)在距离矩阵的样本中不存在', variables=(sample), code="32701606")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "距离统计和统计检验分析结果目录"],
            ["./fliers_file.xls", "xls", "异常值"],
            ["./box_data.xls", "xls", "箱线图数据"]
        ])
        super(AnosimBoxAgent, self).end()


class AnosimBoxTool(Tool):
    def __init__(self, config):
        super(AnosimBoxTool, self).__init__(config)
        self._version = 'R-3.3.1'
        self.cmd_path = 'program/R-3.3.1/bin/Rscript'
        self.package_path = self.config.SOFTWARE_DIR +'/bioinfo/meta/scripts/anosim_box.r'
        self.grouplab = self.option('grouplab') if self.option('grouplab') else self.option('group').prop['group_scheme'][0]

    def run(self):
        """
        运行
        """
        super(AnosimBoxTool, self).run()
        self.run_box()



    def run_box(self):
        """
        运行anosim_box.r
        """
        cmd = self.cmd_path+' '+self.package_path + ' %s %s %s %s' % \
                                                    (self.option('dis_matrix').path, self.option('group').path,
                                                     self.option('permutations'), self.work_dir)
        self.logger.info('运行anosim_box.r程序')
        box_command = self.add_command('box', cmd)
        box_command.run()
        self.wait(box_command)
        if box_command.return_code == 0:
            self.logger.info('运行anosim_box.r完成')
            self.get_fliers()
            self.logger.info('整理结果文件结束')
            self.end()
        elif box_command.return_code is None:
            self.logger.info('返回吗是none，重新运行一次！')
            re_cmd = self.add_command('re_box', cmd).rerun()
            if re_cmd.return_code == 0:
                self.logger.info("重新运行一次成功！")
                self.get_fliers()
                self.logger.info('整理结果文件结束')
                self.end()
            else:
                self.set_error("重新运行一次出错", code="32701607")
        else:
            self.set_error('运行anosim_box.r出错', code="32701608")

    def get_fliers(self):
        """
        处理fliers.xls和group.xls文件，获得奇值点fliers
        整理数据box_data.xls
        """
        fliers_file = self.work_dir + '/fliers.xls'
        group_file = self.work_dir + '/group.xls'
        fliers = self.work_dir + '/fliers_file.xls'
        new_fliers = self.output_dir + '/fliers_file.xls'
        g = open(group_file)
        f = open(fliers_file)
        n = open(fliers, "w")
        # with open(group_file) as g, open(fliers_file) as f, open(fliers, "w") as n:
        group = []
        n.write("group\tfliers\n")
        for line in g:
            group_member = line.strip("\n").strip("\"")
            group.append(group_member)
        for line in f:
            if re.search("\"", line):
                pass
            else:
                number = line.strip().split("\t")
                new_line = group[int(number[0])] + "\t" + number[1] + "\n"
                n.write(new_line)       # fliers_file存放所有的fliers 格式：group /t XXX
        n.close()
        g.close()
        f.close()

        n = open(fliers)
        dic = dict()
        next(n)
        for line in n:
            key = line.strip().split("\t")[0]
            dic[key] =[]
        n.close()

        n = open(fliers)
        for lines in n:
            member = lines.strip().split("\t")
            if len(member) == 2:
                if dic. has_key(member[0]):
                    dic[member[0]].append(member[1])
                else:
                    pass
            else:
                self.logger.info("work_dir中的fliers_file文件格式有误")   #处理fliers，格式 group /t "XX,YYY"

        data = self.work_dir + "/box_data.xls"
        box_data = self.output_dir + '/box_data.xls'
        with open(data) as d, open(box_data, "a+") as b:
            b.write("names\tmin\tq1\tfliers\tmedian\tq3\tmax\n")
            next(d)
            for line in d:
                line = line.replace("\"", "")
                box_line = line.strip().split("\t")
                name = box_line[0]
                if dic. has_key(name):
                    new_line = "\t".join(box_line[0:3]) + "\t" + ",".join(dic[name]) + "\t" + "\t".join(box_line[3:]) + "\n"
                    b.write(new_line)
                else:
                    b.write("\t".join(box_line[0:3]) + "\t" + "\t" + "\t".join(box_line[3:]) + "\n")

        n.close()

