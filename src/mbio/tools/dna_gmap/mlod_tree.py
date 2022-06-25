# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# modified 2018.07.27

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class MlodTreeAgent(Agent):
    """
    只用于生成tree不进行Linkage grouping
    lasted modified by hongdong @ 20180727
    """
    def __init__(self, parent):
        super(MlodTreeAgent, self).__init__(parent)
        options = [
            {"name": "total_mlod", "type": "infile", "format": "dna_gmap.mlod"},  # 输入文件
            {"name": "marker", "type": "infile", "format": "dna_gmap.marker"},  # 输入文件
            {"name": "key_file", "type": "string", "default": "Total"},  # 输入Key值
            {"name": "chr_num", "type": "int"},
            {"name": "start_lod", "type": "int", "default": 3},
            {"name": "end_lod", "type": "int", "default": 20},
            {"name": "stepsize_lod", "type": "int", "default": 1}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("total_mlod").is_set:
            raise OptionError("必须输入total_mlod文件", code="34801001")
        if not self.option("marker").is_set:
            raise OptionError("必须要输入marker文件", code="34801002")
        if not self.option("chr_num"):
            raise OptionError("必须要输入chr_num文件", code="34801003")

    def set_resource(self):
        """
        :return:
        """
        self._cpu = 2
        self._memory = "30G"

    def end(self):
        super(MlodTreeAgent, self).end()


class MlodTreeTool(Tool):
    def __init__(self, config):
        super(MlodTreeTool, self).__init__(config)
        self.perl_path = 'miniconda2/bin/perl'
        self.linkage_tree_path = self.config.PACKAGE_DIR + "/dna_gmap/linkage-tree.pl"
        self.total_mlod = ''

    def run_groupping(self):
        """
        /mnt/ilustre/users/sanger-dev/app/program/perl-5.24.0/bin/perl
        /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/dna_gmap/linkage-tree.pl
         -i /mnt/ilustre/users/sanger-dev/workspace/20180726/Gmap_tsg_31211/Groupping/Total.mlod
         -k Total -d /mnt/ilustre/users/sanger-dev/workspace/20180726/Gmap_tsg_31211/Groupping/output_test
         -g /mnt/ilustre/users/sanger-dev/workspace/20180726/Gmap_tsg_31211/BinnerCalculate/output/Total.bin.marker
         -n 15 -b 2 -e 50 -s 1 -redo -dumper
        """
        cmd = "{} {} -i {} -k {} -d {} -g {} -n {} -b {} -e {} -s {} -redo -dumper"\
            .format(self.perl_path, self.linkage_tree_path, self.total_mlod, self.option("key_file"), self.output_dir,
                    self.option("marker").prop['path'], self.option("chr_num"), self.option("start_lod"),
                    self.option("end_lod"), self.option("stepsize_lod"))
        command = self.add_command("groupping", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("perl运行完成")
        else:
            self.set_error("perl运行失败", code="34801001")

    def mlod_filter(self):
        """
        对Mlod总文件的过滤，默认过滤掉相互作用值小于2的，数据，大数据，要放到tool中
        add by hongdong@20180726
        :return:
        """
        out_file = os.path.join(self.work_dir, "Total.mlod")
        with open(self.option("total_mlod").prop['path'], "r") as r, open(out_file, 'w') as w:
            for line in r:
                if line.startswith("#"):
                    w.write(line)
                else:
                    temp = line.strip().split('\t')
                    if int(temp[2]) >= 2:
                        w.write(line)
        self.total_mlod = out_file

    def run(self):
        super(MlodTreeTool, self).run()
        self.mlod_filter()
        self.run_groupping()
        self.end()


