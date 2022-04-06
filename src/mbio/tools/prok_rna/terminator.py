# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
import glob
import pandas as pd


class TerminatorAgent(Agent):
    """
    Terminator prediction
    """
    def __init__(self, parent):
        super(TerminatorAgent, self).__init__(parent)
        options = [
            {"name": "rock_index", "type": "infile", "format": "prok_rna.common_dir"},  # rockhopper index dir
        ]
        self.add_option(options)
        self.step.add_steps("terminator")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.terminator.start()
        self.step.update()

    def stepfinish(self):
        self.step.terminator.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('rock_index').is_set:
            raise OptionError("必须提供rockhopper index文件夹")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(TerminatorAgent, self).end()


class TerminatorTool(Tool):
    def __init__(self, config):
        super(TerminatorTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))
        self.program = {
            'transterm': 'bioinfo/prok_rna/transterm_hp_v2.09/transterm',
            'python': 'miniconda2/bin/python',
            'parafly': 'program/parafly-r2013-01-21/bin/bin/ParaFly',
        }
        self.file = {
            'dat': os.path.join(self.config.SOFTWARE_DIR, "bioinfo/prok_rna/transterm_hp_v2.09/expterm.dat"),
        }
        self.cmd_list = list()

    def run(self):
        """
        运行
        :return:
        """
        super(TerminatorTool, self).run()
        self.predict_terminator()
        self.run_multi_cmds()
        self.run_stats()
        self.set_output()
        self.end()

    def process_ptt(self, each):
        modified = os.path.join(self.work_dir, os.path.basename(each))
        with open(each, 'r') as i, open(modified, 'w') as o:
            o.write(i.readline())
            o.write(i.readline())
        each_df = pd.read_table(each, header=2)
        each_df = each_df[each_df["Length"] != 0]
        # filter off positive strands extended from the start position of positive strand and negative strand
        if not each_df.empty:
            end_threshold = each_df.iloc[-1, 0].split('..')[1]
            start_threshold = each_df.iloc[0, 0].split('..')[0]
            if not each_df[each_df['Strand'] == '-'].empty:
                end_threshold = each_df[each_df['Strand'] == '-'].iloc[-1, 0].split('..')[1]
            for line in each_df[each_df['Strand'] == '+'].iterrows():
                start, end = line[1]['Location'].split('..')
                if int(start) < int(start_threshold) and int(end) > int(end_threshold):
                    each_df.drop(line[0], inplace=True)
        each_df['Gene'] = each_df['Gene'].apply(lambda x: '-')
        each_df.to_csv(modified, mode='a', index=False, sep='\t')
        return modified

    def predict_terminator(self):
        self.logger.info('开始进行终止子预测')
        transterm = os.path.join(self.config.SOFTWARE_DIR, self.program['transterm'])
        all_dir = os.listdir(self.option('rock_index').prop['path'])
        for each in all_dir:
            abs_path = os.path.join(self.option('rock_index').prop['path'], each)
            if os.path.isdir(abs_path):
                fna_path = os.path.join(abs_path, '{}.fna'.format(each))
                ptt_path = self.process_ptt(os.path.join(abs_path, '{}.ptt'.format(each)))
                cmd = "{} -p {} ".format(transterm, self.file['dat'])
                cmd += "{} {} ".format(fna_path, ptt_path)
                cmd += '--bag-output {} '.format(os.path.join(self.work_dir, each + '_detail.xls'))
                cmd += "> {}".format(os.path.join(self.work_dir, each + '_result.xls'))
                self.cmd_list.append(cmd)

    def run_multi_cmds(self):
        """
        将多个cmd命令并行执行
        """
        cmd_name = 'run_transterm'
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            f.write('\n'.join(self.cmd_list) + '\n')
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(os.path.join(self.config.SOFTWARE_DIR, self.program['parafly']), cmd_file, 2, wrong_cmd)
        command = self.add_command("more_" + cmd_name, cmd_more, ignore_error=True, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("{}完成！".format(cmd_name))
        else:
            self.set_error("{}出错！".format(cmd_name), code="34503602")
            raise Exception("{}出错！".format(cmd_name))

    def run_stats(self):
        results = glob.glob(os.path.join(self.work_dir, '*_detail.xls'))
        merge = os.path.join(self.output_dir, 'terminator_prediction.xls')
        pattern = re.compile(r'\S+')
        with open(merge, 'w') as m:
            m.write('gene\tlocation\tstart\tend\tstrand\t5tail\t5stem\tloop\t3stem\t3tail\tscore\tpos\n')
        for each in results:
            location = os.path.basename(each).split('_detail')[0]
            with open(each, 'r') as f, open(merge, 'a') as m:
                for line in f:
                    items = pattern.findall(line)
                    if items[1] == 'NONE':
                        continue
                    elif len(items) > 12:
                        m.write(items[0] + '\t' + str(location) + '\t' + items[1] + '\t' + items[3] + '\t' +
                                items[4] + '\t' + items[7] + '\t' + items[8] + '\t' + items[9] + '\t' + items[10] +
                                '\t' + items[11] + '\t' + items[-2] + '\t' + items[-1] + '\n')

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "terminator_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna_v3.terminator",
            "instant": False,
            "options": dict(
                rock_index="/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210412/Prokrna_g0uk_3f5rs2hu3jkoaf37fqm5gd/RockhopperIndex/rock_index",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
