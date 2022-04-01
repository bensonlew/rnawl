# -*- coding: utf-8 -*-
# __author__ = 'wentianliu'
# last modify 20190104

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CstacksAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(CstacksAgent, self).__init__(parent)
        options = [
            {"name": "group_list", "type": "string"},
            # {"name": "group_list", "type": "infile", "format": "noref_wgs.group_list"},  # 是否有分组方案
            {"name": "ustacks_output", "type": "string"}  # 03.ustacks结果路径 ./Demo/01.stacks/03.ustacks/
        ]
        self.add_option(options)
        self.step.add_steps('Cstacks')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Cstacks.start()
        self.step.update()

    def step_end(self):
        self.step.Cstacks.finish()
        self.step.update()

    def check_options(self):
        if not self.option('ustacks_output'):
            raise OptionError('必须输入:ustacks_output', code="35500203")

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 16
        self._memory = "230G"

    def end(self):
        super(CstacksAgent, self).end()


class CstacksTool(Tool):
    def __init__(self, config):
        super(CstacksTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.cstacks_st_sh_path = "bioinfo/noRefWGS/cstacks_st.sh"
        self.cstacks_sh_path = "bioinfo/noRefWGS/cstacks.sh"
        self.cstacks_path = self.config.SOFTWARE_DIR + "/bioinfo/noRefWGS/cstacks"

    def mk_list(self):
        """
        生成group.list和sample.list。
        :return:
        """
        if self.option("group_list"):
            os.symlink(self.option("group_list"), os.path.join(self.work_dir, 'group.list'))
            group_dict = {}
            with open(self.option("group_list"), "r")as fr:
                lines = fr.readlines()
                if len(lines) == 0:
                    self.set_error("group.list为空！", code="35500207")
                for line in lines:
                    tmp = line.strip().split("\t")
                    if len(tmp) != 2:
                        self.set_error("group.list应该为两列！", code="35500208")
                    sample = tmp[0]
                    group = tmp[1]
                    if group_dict in group_dict.keys():
                        group_dict[group].append(sample)
                    else:
                        group_dict[group] = [sample]
                write_lines = ""
                for value in group_dict.values():
                    max_sample = 0
                    max_size = 0
                    for i in value:
                        alleles_size = round(os.path.getsize(os.path.join(self.option("ustacks_output"),
                                                                          (i + '.alleles.tsv.gz')))/1024/1024, 2)
                        snps_size = round(os.path.getsize(os.path.join(self.option("ustacks_output"),
                                                                       (i + '.snps.tsv.gz')))/1024/1024, 2)
                        tags_size = round(os.path.getsize(os.path.join(self.option("ustacks_output"),
                                                                       (i + '.tags.tsv.gz')))/1024/1024, 2)
                        size = alleles_size + snps_size + tags_size
                        if size > max_size:
                            max_sample = i
                            max_size = size
                    if max_sample == 0:
                        pass
                    else:
                        write_lines = write_lines + str(max_sample) + '\n'
                sample_path = os.path.join(self.work_dir, "sample.list")
                with open(sample_path, "w")as fw:
                    fw.write(write_lines)
        else:
            ustacks_output = os.path.join(self.option("ustacks_output"), "ustacks.list")
            sample_list = []
            samples = []
            with open(ustacks_output, "r")as fr:
                lines = fr.readlines()
                for line in lines:
                    temp = line.strip().split("\t")
                    sample = temp[0]
                    alleles_size = round(os.path.getsize(os.path.join(self.option("ustacks_output"),
                                                                      (sample + '.alleles.tsv.gz'))) / 1024 / 1024, 2)
                    snps_size = round(os.path.getsize(os.path.join(self.option("ustacks_output"),
                                                                   (sample + '.snps.tsv.gz'))) / 1024 / 1024, 2)
                    tags_size = round(os.path.getsize(os.path.join(self.option("ustacks_output"),
                                                                   (sample + '.tags.tsv.gz'))) / 1024 / 1024, 2)
                    size = alleles_size + snps_size + tags_size
                    sample_list.append([sample, size])
                    samples.append(sample)
            sample_lines = ""
            group_lines = ""
            if len(sample_list) <= 20:
                for i in samples:
                    sample_lines = sample_lines + str(i) + "\n"
                    group_lines = group_lines + str(i) + "\t" + str(1) + "\n"
            else:
                sample_list.sort(key=lambda x: x[1], reverse=True)
                for i in sample_list[:20]:
                    sample_lines = sample_lines + str(i[0]) + "\n"
                    group_lines = group_lines + str(i[0]) + "\t" + str(1) + "\n"
            sample_path = os.path.join(self.work_dir, "sample.list")
            group_path = os.path.join(self.work_dir, "group.list")
            with open(sample_path, "w")as fw:
                fw.write(sample_lines)
            with open(group_path, "w")as fw:
                fw.write(group_lines)

    def run_cstacks(self, sample, is_head):
        """
        第一个样本执行cstacks_st.sh，其他执行cstacks.sh。
        is_head 是否为第一个，1为第一个，0为不是第一个。
        """
        input_path = os.path.join(self.option('ustacks_output'), sample)
        log_path = os.path.join(self.output_dir, ("cstacks." + sample + ".log"))
        catalog_path = os.path.join(self.output_dir, "catalog")
        if int(is_head) == 1:
            allfiles = os.listdir(self.output_dir)
            oldfiles = [os.path.join(self.output_dir, i) for i in allfiles]
            for oldfile in oldfiles:
                os.remove(oldfile)
            cmd = "{} {} {} {} {}".format(self.cstacks_st_sh_path, self.cstacks_path, input_path, self.output_dir,
                                          log_path)
        else:
            cmd = "{} {} {} {} {} {}".format(self.cstacks_sh_path, self.cstacks_path, input_path, self.output_dir,
                                             catalog_path, log_path)
        self.logger.info("开始进行cstacks" + str(sample))
        command = self.add_command(("cstacks" + str(sample)).lower(), cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info(str(sample) + "cstacks完成！")
        else:
            self.set_error(str(sample) + "cstacks出错！", code="35500209")

    def run(self):
        super(CstacksTool, self).run()
        self.mk_list()
        sample_path = os.path.join(self.work_dir, "sample.list")
        new_sample_path = os.path.join(self.output_dir, "sample.list")
        os.link(sample_path, new_sample_path)
        sample_list = []
        with open(sample_path, "r")as fr:
            lines = fr.readlines()
            if len(lines) == 0:
                pass
            else:
                for line in lines:
                    tmp = line.strip().split()
                    sample_list.append(tmp[0])
                self.run_cstacks(str(sample_list[0]), 1)
                for i in sample_list[1:]:
                    self.run_cstacks(str(i), 0)
        self.end()
