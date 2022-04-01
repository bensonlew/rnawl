## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "hongyu.chen"
import re
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class DedupMultiAgent(Agent):
    """
    新的查重模块
    用于全库排查守护进程：将2-10的查重放在一个tool里完成，减少同时投递的tool数目，以免任务队列超限造成的无谓等待
    包括脚本：pt_dup_new.R
    version v1.0
    author: hongyu.chen
    last_modify: 20180426
    """
    def __init__(self, parent):
        super(DedupMultiAgent, self).__init__(parent)
        options = [
            {"name": "mom_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "preg_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "ref_point", "type": "infile", "format": "paternity_test.rda"},
            {"name": "err_min", "type": "int", "default": 2},
            {"name": "dad_id", "type": "string"},
            {"name": "father_path", "type": "string"},  # 输入父本tab文件的所在路径
            {"name": "dad_list", "type": "string"},  # 要进行分析的父本列表
            {"name": "mem", "type": "int", "default": 4}
        ]
        self.add_option(options)
        self.step.add_steps("dedup_analysis")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.dedup_analysis.start()
        self.step.update()

    def stepfinish(self):
        self.step.dedup_analysis.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('ref_point'):
            raise OptionError("必须提供参考位点文件")
        if not self.option('mom_tab'):
            raise OptionError("必须提供母本tab")
        if not self.option('preg_tab'):
            raise OptionError("必须提供胎儿tab")
        if not self.option('father_path'):
            raise OptionError("必须提供查重部分父本tab")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '{}G'.format(self.option("mem"))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(DedupMultiAgent, self).end()


class DedupMultiTool(Tool):
    """
    查重tool
    # use: pt_dup_new.R WQ123M.tab WQ123S.tab 2
    /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/pt/targets.bed.rda WQ123F WQ123-F1,WQ123-F2
    """
    def __init__(self, config):
        super(DedupMultiTool, self).__init__(config)
        self._version = '1.0.1'
        self.R_path = 'program/R-3.3.1/bin/'
        self.script_path = self.config.SOFTWARE_DIR + '/bioinfo/medical/scripts/'
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def get_father_list(self):
        """
        获取dad_id前后200个父本tab文件，暂时不用
        :return:
        """
        dad_list = []
        results = os.listdir(self.option("father_path"))
        m = re.match(r'WQ([0-9]*)-F.*', self.option("dad_id"))
        if m:
            for dad in results:
                if int(m.group(1)) - 200 <= int(re.match(r'WQ([0-9]*)-F.*', str(dad)).group(1)) \
                        <= int(m.group(1)) + 201:
                    dad_list.append(dad)
        else:
            raise Exception("父本{}命名不规范！".format(self.option("dad_id")))
        return dad_list

    def run_tf(self, err_min):
        # dad_list = self.get_father_list()
        # self.logger.info(self.option("dad_list"))
        dedup_cmd = "{}Rscript {}pt_dup_new.R {} {} {} {} {} {} {} {}".format(self.R_path, self.script_path,
                                                                              self.option("mom_tab").prop['path'],
                                                                              self.option("preg_tab").prop['path'],
                                                                              err_min,
                                                                              self.option("ref_point").prop['path'],
                                                                              "pt_result_{}".format(err_min),
                                                                              self.option("father_path"),
                                                                              self.option("dad_id"),
                                                                              self.option("dad_list"))
        self.logger.info(dedup_cmd)
        self.logger.info("开始进行查重分析，err_min: {}".format(err_min))
        cmd = self.add_command("dedup_cmd_" + str(err_min), dedup_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行查重成功，err_min: {}".format(err_min))
        else:
            self.set_error('运行查重出错，err_min: {}'.format(err_min))
            raise Exception("运行查重出错，err_min: {}".format(err_min))

    def set_output(self, err_min):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        dir_name = "pt_result_" + str(err_min)
        self.linkdir(os.path.join(self.work_dir, dir_name), os.path.join(self.output_dir, dir_name))
        self.logger.info('设置文件夹{}路径成功'.format(dir_name))

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到指定目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                next_dirname = dirname + '/' + os.path.basename(oldfiles[i])
                self.linkdir(oldfiles[i], next_dirname)

    def run(self):
        super(DedupMultiTool, self).run()
        for i in range(2, self.option("err_min")):
            self.run_tf(i)
            self.set_output(i)
        self.end()
