# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class DrawCircosAgent(Agent):
    """
    该tool用于整理出画circos的数据
    version 1.0
    author: HONGDONG
    last_modified:20180420
    """
    
    def __init__(self, parent):
        super(DrawCircosAgent, self).__init__(parent)
        options = [
            {"name": "windows", "type": "int", "default": 100000},
            {"name": "snp", "type": "string"},  # snp.anno.primary.vcf
            {"name": "indel", "type": "string"},  # indel.anno.primary.vcf
            {"name": "gff", "type": "string"},   # 这个是gene.gff文件，不是ref.gff文件
            {"name": "sv", "type": "string"},
            {"name": "chrlist", "type": "string"},
            {"name": "cnv", "type": "string"}
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
        if not self.option('snp'):
            raise OptionError('必须提供snp结果表', code="34502201")
        if not self.option('indel'):
            raise OptionError('必须提供indel结果表', code="34502202")
        if not self.option("gff"):
            raise OptionError("必须提供gff结果文件", code="34502203")
        if not self.option("chrlist"):
            raise OptionError("必须提供chrlist文件！", code="34502204")
        if not self.option("sv"):
            raise OptionError("必须提供sv文件！", code="34502205")
        if not self.option("cnv"):
            raise OptionError("必须提供cnv文件！", code="34502206")
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
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/WGS/circos-0.69-6/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin')
        self.script_path = self.config.PACKAGE_DIR + '/wgs/draw.circos.pl'
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        self.image_magick = "program/ImageMagick/bin/convert"

    def run_drawcircos(self):
        """
        perl draw.circos.pl --windows 100000 --snp snp.anno.primary.vcf --indel indel.anno.primary.vcf
         --chrlist ref.chrlist --gff genes.gff --outdir 13.variant-stat --sv A8-10.sv.anno --cnv A8-10.cnv.anno
        :return:
        """
        one_cmd = "{}{} --windows {} --snp {} --indel {} --chrlist {} --gff {} --outdir {} --sv {} --cnv {}"\
            .format(self.perl_path, self.script_path, self.option("windows"), self.option("snp"), self.option("indel"),
                    self.option("chrlist"), self.option("gff"), self.work_dir, self.option("sv"), self.option("cnv"))
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.set_error("运行one_cmd出错！", code="34502201")
            self.set_error("运行one_cmd出错！", code="34502205")

    def convert_png_pdf(self):
        cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + self.work_dir + \
              "/circos.png" + ' ' + self.work_dir + "/circos.pdf"
        cmd_obj = self.add_command("convert", cmd).run()
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info("cmd执行成功")
        else:
            self.set_error("cmd运行出错!", code="34502202")

    def set_output(self):
        """
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info('开始设置文件夹路径')
        sample_id = os.path.basename(self.option("cnv")).split('.')[0]
        results = os.listdir(self.work_dir)
        for f in results:
            if f in ["circos.png", 'circos.svg', 'circos.pdf']:
                os.link(os.path.join(self.work_dir, f), self.output_dir + "/{}.{}".format(sample_id, f.split('.')[1]))
        self.logger.info('设置文件夹路径成功')

    def run(self):
        """
        运行
        """
        super(DrawCircosTool, self).run()
        self.run_drawcircos()
        self.convert_png_pdf()
        self.set_output()
        self.end()


