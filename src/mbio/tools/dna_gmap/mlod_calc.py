# -*- coding: utf-8 -*-
# __author__ = 'Zhao Binbin'
# modified 2018.06.05

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class MlodCalcAgent(Agent):
    """
    自动分群中lod计算
    lasted modified by hongdong@21080627
    1)解决pep8错误
    2)完善infile与outfile
    3)完善当是最后一个genotype的时候不要-p参数
    4)解决tool运行命令的组建逻辑有问题
    """
    def __init__(self, parent):
        super(MlodCalcAgent, self).__init__(parent)
        options = [
            {"name": "sub_genotype", "type": "infile", "format": "dna_gmap.marker"},
            {"name": "only", "type": 'int', "default": 200},
            {"name": "sub_mlod", "type": "outfile", "format": "dna_gmap.mlod"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sub_genotype").is_set:
            raise OptionError("必须输入sub.genotype文件", code="34800901")

    def set_resource(self):
        self._cpu = 2
        self._memory = "3G"

    def end(self):
        super(MlodCalcAgent, self).end()


class MlodCalcTool(Tool):
    def __init__(self, config):
        super(MlodCalcTool, self).__init__(config)
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.calculateMLOD_path = self.config.PACKAGE_DIR + "/dna_gmap/calculateMLOD.pl"

    def run_mlod_calc(self):
        """
        mlod calculation
        /mnt/ilustre/users/sanger-dev/app/program/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/calculateMLOD.pl
        -i /mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/yichuantupu/sub/sub.0.genotype
        -o /mnt/ilustre/users/sanger-dev/workspace/20180612/Single_test_mlod_calc_0607/MlodCalc/output/sub.0.mlod
        -p 200 -s
        """
        list_ = self.option("sub_genotype").prop['path'].split("/")
        name = list_[-1]
        newname = name.split(".")[1]
        self.logger.info(newname)
        cmd = "{} {} -i {} -o {} -s".format(self.perl_path, self.calculateMLOD_path,
                                            self.option("sub_genotype").prop['path'],
                                            self.output_dir + "/sub." + newname + ".mlod")
        if self.option("only"):
            cmd += " -p {}".format(self.option("only"))
        self.logger.info(cmd)
        command = self.add_command("mlod_calc", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("perl运行完成")
        elif command.return_code == 25:
            self.set_error("different number of individuals between $marker1 and $marker2", code="34800901")
        else:
            self.set_error("perl运行失败", code="34800902")
        self.option("sub_mlod").set_path(self.output_dir + "/sub." + newname + ".mlod")

    def run(self):
        super(MlodCalcTool, self).run()
        self.run_mlod_calc()
        self.end()


