# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os


class VcftoolsStatAgent(Agent):
    """
    工具：vcftools
    对vcf进行window-pi、TajimaD、weir-fst-pop统计
    """
    def __init__(self, parent):
        super(VcftoolsStatAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf", "required": True},  # pop.recode.vcf或者比较分析的结果
            {"name": "group_file", "type": "infile", "format": 'ref_rna_v2.common'},  # 分组文件,输出文件以分组文件名称.分割命名
            {"name": "window_size", "type": "int", "default": 10000000},  # --window-pi,窗口大小
            {"name": "window_step", "type": "int", "default": 2000000},  # --window-pi-step,窗口步长
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('group_file').is_set:
            raise OptionError('必须设置输入分组文件。')
        if not self.option('window_size') > 0:
            raise OptionError('必须设置输入大于0的window size。')
        if not self.option('window_step') > 0:
            raise OptionError('必须设置输入大于0的window step。')
        if not self.option('window_size') > self.option('window_step'):
            raise OptionError('window size必须大于window step。')
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(VcftoolsStatAgent, self).end()


class VcftoolsStatTool(Tool):
    def __init__(self, config):
        super(VcftoolsStatTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.program = {
            'vcftools': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/dna_evolution/vcftools'),
            'parafly': 'program/parafly-r2013-01-21/bin/bin/ParaFly',
        }
        self.samples = list()
        self.cmd_list = list()

    def split_group(self):
        group_dict = dict()
        with open(self.option('group_file').prop['path'], 'r') as g:
            for line in g:
                if line.startswith('#'):
                    continue
                items = line.strip().split('\t')
                if len(items) == 2:
                    if items[1] not in group_dict.keys():
                        group_dict[items[1]] = list()
                        self.samples.append(items[1])
                    group_dict[items[1]].append(items[0])
        if len(self.samples) > 2:
            self.set_error('More than two groups have been found in your group file.')
        for each in self.samples:
            s_file = os.path.join(self.work_dir, each+'.list')
            with open(s_file, 'w') as s:
                s.write('\n'.join(group_dict.get(each)) + '\n')

    def run_vcftools_window_pi(self):
        """
        window-pi方法
        """
        for each in self.samples:
            group_list = os.path.join(self.work_dir, each + '.list')
            out = os.path.join(self.output_dir, each)
            cmd = "{} --vcf {} --remove-indels --keep {}".format(self.program['vcftools'],
                                                                 self.option("vcf_file").prop["path"],
                                                                 group_list)
            cmd += " --window-pi {} --window-pi-step {} --out {}".format(self.option("window_size"),
                                                                         self.option("window_step"), out)
            self.cmd_list.append(cmd)

    def run_vcftools_tajima_d(self):
        """
        TajimaD方法
        """
        for each in self.samples:
            group_list = os.path.join(self.work_dir, each + '.list')
            out = os.path.join(self.output_dir, each)
            cmd = "{} --vcf {} --remove-indels --keep {}".format(self.program['vcftools'],
                                                                 self.option("vcf_file").prop["path"],
                                                                 group_list)
            cmd += " --TajimaD {} --out {}".format(self.option("window_step"), out)
            self.cmd_list.append(cmd)

    def run_vcftools_weir_fst_pop(self):
        """
        weir-fst-pop方法
        """
        group1 = os.path.join(self.work_dir, self.samples[0] + '.list')
        group2 = os.path.join(self.work_dir, self.samples[1] + '.list')
        out = os.path.join(self.output_dir, self.samples[0] + '_vs_' + self.samples[1])
        cmd = "{} --vcf {} --remove-indels --weir-fst-pop {}".format(self.program['vcftools'],
                                                                     self.option("vcf_file").prop["path"],
                                                                     group1)
        cmd += " --weir-fst-pop {} --fst-window-size {}".format(group2, self.option("window_size"))
        cmd += " --fst-window-step {} --out {}".format(self.option("window_step"), out)
        self.cmd_list.append(cmd)

    def run_multi_cmds(self):
        """
        将多个cmd命令并行执行
        """
        cmd_name = 'run_vcftools_sweep'
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            f.write('\n'.join(self.cmd_list) + '\n')
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.program['parafly'], cmd_file, 2, wrong_cmd)
        command = self.add_command("more_" + cmd_name, cmd_more, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}完成！".format(cmd_name))
        else:
            self.set_error("{}出错！".format(cmd_name), code="34503602")
            raise Exception("{}出错！".format(cmd_name))

    def run(self):
        super(VcftoolsStatTool, self).run()
        self.split_group()
        self.run_vcftools_window_pi()
        self.run_vcftools_weir_fst_pop()
        self.run_vcftools_tajima_d()
        self.run_multi_cmds()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "selective_sweep_stats_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.selective_sweep.vcftools_stat",
            "instant": False,
            "options": dict(
                vcf_file='/mnt/ilustre/users/sanger-dev/workspace/20210513/Single_selective_sweep_filter_4504_9628/VcftoolsFilter/output/pop.recode.vcf',
                group_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/test/tool_052021/selective_sweep/group.txt',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)