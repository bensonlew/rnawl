# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 2018.06.27

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class HandGrouppingAgent(Agent):
    """
    进行手动分群
    lasted modified by hongdong @ 20180627
    """
    def __init__(self, parent):
        super(HandGrouppingAgent, self).__init__(parent)
        options = [
            {"name": "gtree_hash", "type": "infile", "format": "dna_gmap.gtree_hash"},
            {"name": "tree_nodes", "type": "string"},  # "19/32;19/41;5/16677"
            {"name": "total_lg", "type": "outfile", "format": "dna_gmap.lg"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("gtree_hash").is_set:
            raise OptionError("请设置gtree_hash文件", code="34800501")
        if not self.option("tree_nodes"):
            raise OptionError("请设置tree_nodes文件", code="34800502")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(HandGrouppingAgent, self).end()


class HandGrouppingTool(Tool):
    def __init__(self, config):
        super(HandGrouppingTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.linkage_by_user = self.config.PACKAGE_DIR + "/dna_gmap/linkage_by_user.pl"
        self.tree_nodes_file = ''

    def run_linkage_by_user(self):
        """
        /mnt/ilustre/users/sanger-dev/app/program/perl/perls/perl-5.24.0/bin/perl
        ~/biocluster/src/mbio/packages/dna_gmap/linkage_by_user.pl -i Total.gTree.delFragment.mergeNode.hash
        -s test_node -o Total.lg
        """
        cmd = "{} {} -i {} -o {} -s {}".format(self.perl_path, self.linkage_by_user,
                                               self.option("gtree_hash").prop['path'], self.output_dir + "/Total.lg",
                                               self.tree_nodes_file)
        self.run_cmd(cmd, "linkage_by_user")

    def set_tree_nodes(self):
        tree_nodes_file = self.work_dir + "/tree_node"
        if os.path.exists(tree_nodes_file):
            os.remove(tree_nodes_file)
        nodes = self.option("tree_nodes").split(";")
        with open(tree_nodes_file, "w") as w:
            for m in nodes:
                w.write("{}\n".format(m))
        self.tree_nodes_file = tree_nodes_file

    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", variables=(cmd_name), code="34800501")

    def run(self):
        super(HandGrouppingTool, self).run()
        self.set_tree_nodes()
        self.run_linkage_by_user()
        self.option("total_lg").set_path(self.output_dir + "/Total.lg")
        self.end()
