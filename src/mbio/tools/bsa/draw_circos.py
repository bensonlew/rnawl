# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class DrawCircosAgent(Agent):
    """
    该tool用于整理出画circos的数据
    perl /mnt/ilustre/users/qingmei.cui/newmdt/sanger/BSA/BSA_pl_cui/bin/draw-circos.pl --windows 100000
    --gwindows 10000 -vcf /mnt/ilustre/users/qingmei.cui/newmdt/sanger/BSA/Data/NiNanJie_Mutmap/
    BSA_Test_TAIR10/pop.final.vcf --index sliding-win.result --column 4
    --chrlist /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.2.24.var/02.ref-config/ref.chrlist
    --gff /mnt/ilustre/users/qingmei.cui/newmdt/sanger/BSA/Data/NiNanJie_Mutmap/BSA_Test_TAIR10/pop.summary
    version 1.0
    author: HONGDONG
    last_modified:2018.03.06
    """
    
    def __init__(self, parent):
        super(DrawCircosAgent, self).__init__(parent)
        options = [
            {"name": "win", "type": "int", "default": 100000},
            {"name": "gwin", "type": "int", "default": 10000},
            {"name": "p_v_vcf", "type": "string"},  # pop.final.vcf
            {"name": "s_w_result", "type": "string"},  # sliding-win.result
            {"name": "column", "type": "int", "default": 4},
            {"name": "chrlist", "type": "string"},
            {"name": "pop_summary", "type": "string"}  # pop.summary
        ]
        self.add_option(options)
        self.step.add_steps('DrawCircos')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.DrawCircos.start()
        self.step.update()
        
    def step_end(self):
        self.step.DrawCircos.finish()
        self.step.update()
        
    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('pop_summary'):
            raise OptionError('必须提供pop_summary结果表', code="31500201")
        if not self.option('p_v_vcf'):
            raise OptionError('必须提供pop.final.vcf结果表', code="31500202")
        if not self.option("s_w_result"):
            raise OptionError("必须提供sliding-win.result结果文件", code="31500203")
        if not self.option("chrlist"):
            raise OptionError("必须提供chrlist文件！", code="31500204")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '10G'
        
    def end(self):
        super(DrawCircosAgent, self).end()


class DrawCircosTool(Tool):

    def __init__(self, config):
        super(DrawCircosTool, self).__init__(config)
        self._version = "1.0.1"
        self.script_path = self.config.PACKAGE_DIR + '/bsa/draw-circos.pl'
        self.perl_path = 'miniconda2/bin/perl '

    def run_drawcircos(self):
        one_cmd = self.perl_path + self.script_path + \
                  " --windows {} --gwindows {} -vcf {} --index {} --column {} --chrlist {} --gff {}".format(
                      self.option("win"), self.option("gwin"), self.option("p_v_vcf"), self.option("s_w_result"),
                      self.option("column"), self.option("chrlist"), self.option("pop_summary"))
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.set_error("运行one_cmd出错！", code="31500201")
            self.set_error("运行one_cmd出错！", code="31500204")

    def set_output(self):
        """
        circos.chrlist
        gene.num.csv
        snp.win.csv
        indel.win.csv
        sliding.win.csv
        为插件所需文件
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info('开始设置文件夹路径')
        file_path = os.path.join(self.work_dir, "draw.circos/")
        results = os.listdir(file_path)
        for f in results:
            if f in ["circos.chrlist", "gene.num.csv", "snp.win.csv", "indel.win.csv", "sliding.win.csv"]:
                os.link(os.path.join(file_path, f), self.output_dir + "/" + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        """
        运行
        """
        super(DrawCircosTool, self).run()
        self.run_drawcircos()
        self.set_output()
        self.end()


