# -*- coding: utf-8 -*-
# __author__ = 'Zhao Binbin'
# modified 2018.06.05

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GrouppingAgent(Agent):
    """
    自动分群中lod计算
    lasted modified by hongdong @ 20180627
    1）修改infile与outfile
    2）修改pep8报错
    3）tool调用的参数与命令与产品线对应不上
    4）修改根据mlod文件动态设定内存大小
    """
    def __init__(self, parent):
        super(GrouppingAgent, self).__init__(parent)
        options = [
            {"name": "group_type", "type": "string", "default": 'ref'},  # mlod or ref
            {"name": "total_mlod", "type": "infile", "format": "dna_gmap.mlod"},  # 输入文件
            {"name": "marker", "type": "infile", "format": "dna_gmap.marker"},  # 输入文件
            {"name": "key_file", "type": "string", "default": "Total"},  # 输入Key值
            {"name": "scaf_marker", "type": "string", "default": "no"},  # yes就添加sca 否则不添加
            {"name": "chr_num", "type": "int"},
            {"name": "start_lod", "type": "int", "default": 3},
            {"name": "end_lod", "type": "int", "default": 20},
            {"name": "stepsize_lod", "type": "int", "default": 1},
            {"name": "min_group", "type": "int", "default": 20},
            {"name": "max_group", "type": "int", "default": 500},
            {"name": "total_lg", "type": "outfile", "format": "dna_gmap.lg"}
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("group_type"):
            if self.option("group_type") not in ['ref', 'mlod']:
                raise OptionError("group_type必须为ref 或者 mlod", code="34800401")
            if self.option("group_type") == 'mlod':
                if not self.option("chr_num"):
                    raise OptionError("当分群方式为mlod的时候，必须输入chr_num参数", code="34800402")
            if self.option("group_type") == 'ref':
                if not self.option("marker").is_set:
                    raise OptionError("当分群方式为ref的时候，必须输入marker参数", code="34800403")
        else:
            raise OptionError("必须输入group_type参数", code="34800404")
        if not self.option("total_mlod").is_set:
            raise OptionError("必须输入total_mlod文件", code="34800405")
        if self.option("scaf_marker"):
            if self.option("scaf_marker") not in ['no', 'yes']:
                raise OptionError("scaf_marker必须是yes或者no", code="34800406")

    def set_resource(self):
        """
        根据mlod文件的大小去动态设置内存
        :return:
        """
        if self.option("group_type") == 'ref':
            mem = 10
        else:
            size = os.path.getsize(self.option('total_mlod').prop['path'])
            mem = (size/1000/1000/1000) * 10
            if mem == 0:
                mem = 100
        self._cpu = 2
        self._memory = "{}G".format(mem)

    def end(self):
        super(GrouppingAgent, self).end()


class GrouppingTool(Tool):
    def __init__(self, config):
        super(GrouppingTool, self).__init__(config)
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.link_by_ref_path = self.config.PACKAGE_DIR + "/dna_gmap/linkage_by_ref.pl"
        self.link_by_mlod_path = self.config.PACKAGE_DIR + "/dna_gmap/linkage_by_mlod.pl"
        self.total_mlod = ''

    def run_groupping(self):
        """
        分群
        """
        if self.option("group_type") == "ref":  # ref存在的条件下
            cmd = "{} {} -i {} -2 {} -o {} -k {}".format(self.perl_path, self.link_by_ref_path, self.total_mlod,
                                                         self.option("marker").prop['path'], self.output_dir,
                                                         self.option("key_file"))
            if self.option("scaf_marker") == "yes":  # 默认不添加sca，如果要添加的话，增加-add参数
                cmd += " -add"
            command = self.add_command("groupping", cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("perl运行完成")
            else:
                self.set_error("perl运行失败", code="34800401")
        elif self.option("group_type") == "mlod":
            cmd = "{} {} -i {} -k {} -d {} -g {} -n {} -b {} -e {} -s {} -minGroup" \
                  " {} -maxGroup {} -redo -dumper".format(self.perl_path, self.link_by_mlod_path, self.total_mlod,
                                                          self.option("key_file"), self.output_dir,
                                                          self.option("marker").prop['path'], self.option("chr_num"),
                                                          self.option("start_lod"), self.option("end_lod"),
                                                          self.option("stepsize_lod"), self.option("min_group"),
                                                          self.option("max_group"))
            command = self.add_command("groupping", cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("perl运行完成")
            else:
                self.set_error("perl运行失败", code="34800402")
        if os.path.exists(self.output_dir + "/Total.lg"):
            self.option("total_lg").set_path(self.output_dir + "/Total.lg")

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
        super(GrouppingTool, self).run()
        self.mlod_filter()
        self.run_groupping()
        self.end()


