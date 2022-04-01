# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import subprocess
from biocluster.core.exceptions import OptionError
import types
import glob
from biocluster.config import Config



class MaaslinAgent(Agent):


    def __init__(self, parent):
        super(MaaslinAgent, self).__init__(parent)
        options = [
            {"name": "taxon_table", "type": "infile","format":"meta.otu.otu_table"},
            {"name": "env_table", "type": "infile","format":"meta.otu.otu_table"}
        ]
        self.add_option(options)
        self.step.add_steps('maaslin')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.maaslin.start()
        self.step.update()

    def step_end(self):
        self.step.maaslin.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('taxon_table').is_set:
            raise OptionError('必须提供物种层次的表', code="32707101")

        with open(self.option('taxon_table').path) as f:
            lines = f.readlines()
            if len(lines) < 3:
                raise OptionError('物种类型数目少，无法做分析', code="32707102")


        if not self.option('env_table').is_set:
            raise OptionError('必须提供环境因子表', code="32707103")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        # ])
        super(MaaslinAgent, self).end()


class MaaslinTool(Tool):
    def __init__(self, config):
        super(MaaslinTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        #self.script1_path = self.config.PACKAGE_DIR + "/statistical/regression_analysis.pl"    
        self.script1_path = self.config.PACKAGE_DIR + "/meta/scripts/maaslin_analysis.pl"
        #self.R_path = 'program/R-3.3.1/bin/Rscript'
        conf = Config()
        if conf.wpm_user == 'sanger-dev':
            self.R_path = 'program/R-3.5.1/bin/Rscript'
            LD_LIBRARY_PATH= '{0}/library/gdal-2.0.0/lib:{0}/library/proj/lib:{0}/library/lib:{0}/gcc/5.1.0/lib64:' \
                         '{0}/library/geos-3.4.1/lib:{0}/library/udunits-2.0.4/lib'.format(self.config.SOFTWARE_DIR)
        else:
            self.R_path = 'program/R-3.3.1/bin/Rscript'
            LD_LIBRARY_PATH= '{0}/library/gdal-2.0.0/lib:{0}/library/proj/lib:{0}/library/lib:' \
                         '{0}/library/geos-3.4.1/lib:{0}/library/udunits-2.0.4/lib'.format(self.config.SOFTWARE_DIR)

        self.set_environ(LD_LIBRARY_PATH=LD_LIBRARY_PATH)


    def run_command(self, table):
        cmd = self.perl_path +' %s -data %s -f %s -m maaslin.pcl -o %s ' % (
        self.script1_path, table, self.option('env_table').prop['path'], self.work_dir + '/Maaslin')
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Maaslin的r文件生成成功")
        else:
            self.set_error("Maaslin的r文件生成失败", code="32707101")
            self.set_error("Maaslin的r文件生成失败", code="32707102")

        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/Maaslin.cmd.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_r', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算Maaslin成功")
        else:
            self.set_error("R程序计算Maaslin失败", code="32707103")
            self.set_error("R程序计算Maaslin失败", code="32707104")

    def run(self):
        super(MaaslinTool, self).run()
        input_table, replace_dict = self.replace_table(self.option("taxon_table").prop['path'])
        self.run_command(input_table)
        self.set_output(replace_dict)

    def replace_table(self, table):
        """
        功能：对功能表替换功能名称（因为功能或者物种水平可能含有特殊字符）
        modify: 20200313
        :return:
        """
        self.logger.info("功能表进行物种名称替换！")
        table_path = table
        table_name = os.path.basename(table_path)
        out_table = os.path.join(self.work_dir,table_name)
        spe_num = 1
        replace_dict = {}
        with open(table_path, 'r') as f, open(out_table, 'w') as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                line = line.strip().split("\t")
                species_name = line[0]
                replace_name = "spe{}".format(spe_num)
                if replace_name not in replace_dict.keys():
                    replace_dict[replace_name] = species_name
                spe_num += 1
                w.write(replace_name + "\t" + "\t".join(line[1:]) + "\n")
        return (out_table, replace_dict)

    def replace_back(self, dict, table, outtable):
        """
        功能：根据dict将table的名称替换回来
        :param dict: 前面存的物种与替换名称的对应关系
        :param table: 输入表
        :param outtable: 输出表
        :return:
        """
        self.logger.info("开始对{}进行替换".format(table))
        with open(table, 'r') as f, open(outtable, 'w') as w:
            lines = f.readlines()
            if 'Maaslin.xls' in table:
                new_head_list = []
                head_list = lines[0].strip().split("\t")
                for name in head_list:
                    if name in dict.keys():
                        replace_name = dict[name]
                        new_head_list.append(replace_name)
                    else:
                        new_head_list.append(name)
                w.write("\t".join(new_head_list) + '\n')
                for line in lines[1:]:
                    w.write(line)
            else:
                w.write(lines[0])
                if 'Maaslin.txt' in table:
                    for line in lines[1:]:
                        line = line.strip().split("\t")
                        species_name = line[1]
                        if species_name in dict.keys():
                            replace_name = dict[species_name]
                            w.write(line[0] + '\t' +replace_name + "\t" + "\t".join(line[2:]) + "\n")
                else:
                    for line in lines[1:]:
                        line = line.strip().split("\t")
                        species_name = line[0]
                        if species_name in dict.keys():
                            replace_name = dict[species_name]
                            w.write(replace_name + "\t" + "\t".join(line[1:]) + "\n")

    def set_output(self, dict):
        self.logger.info("开始设置结果文件目录")
        msums = glob.glob(self.work_dir + '/Maaslin/maaslin-*.txt')
        linksites_data = os.path.join(self.output_dir, 'Maaslin.xls')
        linksum_data = os.path.join(self.output_dir, 'Maaslin.txt')
        linkline_data = os.path.join(self.output_dir,'Maaslin.message.xls')
        for i in linksites_data,linksum_data, linkline_data:
            if os.path.exists(i):
               os.remove(i)
        os.link(self.work_dir + '/Maaslin/maaslin.tsv', linksites_data)
        os.link(msums[0], linksum_data)
        os.link(self.work_dir + '/Maaslin.message.xls', linkline_data)
        all_files = os.listdir(self.output_dir)
        for file in all_files:
            file_path = os.path.join(self.output_dir, file)
            new_file_path = os.path.join(self.output_dir, file + "_new")
            self.replace_back(dict, file_path, new_file_path)
            os.remove(file_path)
            os.rename(new_file_path, file_path)
        self.logger.info("设置结果文件目录完成！")
        self.end()



