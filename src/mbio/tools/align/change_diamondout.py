# -*- coding: utf-8 -*-
# __author__ = 'shijin'
from biocluster.agent import Agent
from biocluster.tool import Tool
import re
import os

class ChangeDiamondoutAgent(Agent):
    """
    修改diamond软件的输出
    """
    def __init__(self, parent):
        super(ChangeDiamondoutAgent, self).__init__(parent)
        options = [
            {"name": "nr_out", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "kegg_out", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "string_out", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "blast_nr_xml", "type": "outfile", "format": "align.blast.blast_xml"},
            {"name": "blast_string_xml", "type": "outfile", "format": "align.blast.blast_xml"},
            {"name": "blast_kegg_xml", "type": "outfile", "format": "align.blast.blast_xml"}
            ]
        self.add_option(options)
        self.step.add_steps('change')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.change.start()
        self.step.update()

    def step_end(self):
        self.step.change.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '2G'

    def end(self):
        super(ChangeDiamondoutAgent, self).end()


class ChangeDiamondoutTool(Tool):
    def __init__(self, config):
        super(ChangeDiamondoutTool, self).__init__(config)
        self.path_dict = {}
        self.ori = []
        self.repl = []
        if self.option("nr_out").is_set:
            self.path_dict["blast_nr_xml"] = self.option("nr_out").prop["path"]
        if self.option("kegg_out").is_set:
            self.path_dict["blast_kegg_xml"] = self.option("kegg_out").prop["path"]
        if self.option("string_out").is_set:
            self.path_dict["blast_string_xml"] = self.option("string_out").prop["path"]

    def xml_deal(self):
        for key in self.path_dict.keys():
            path = self.path_dict[key]
            with open(path,"r") as file:
                for line in file:
                    line = line.strip()
                    if line.lstrip().startswith("<Hit_id>"):
                        m = re.match("<Hit_id>(.+)</Hit_id>", line.lstrip())
                        if m:
                            self.ori.append(m.group(1))
                            line = file.next()
                            n = re.match("<Hit_def>(.+)</Hit_def>", line.lstrip())
                            try:
                                self.repl.append(n.group(1))
                            except:
                                print line
            with open(path,"r") as file, open(path + "_new", "w") as fw:
                i = 0
                for line in file:
                    if line.lstrip().startswith("<BlastOutput_version>"):
                        line = line.replace("diamond 0.8.35", "BLASTX 2.3.0+")
                    if line.lstrip().startswith("<Hit_id>"):
                        m = re.match("<Hit_id>(.+)</Hit_id>", line.lstrip())
                        if m:
                            line = line.replace(self.ori[i],self.repl[i])
                            i += 1
                    fw.write(line)
            name = os.path.basename(path)
            os.system("mv {} {}".format(path+"_new",self.output_dir + "/" + name))
            self.option(key).set_path(self.output_dir + "/" + name)


    def run(self):
        """
        运行
        :return:
        """
        super(ChangeDiamondoutTool, self).run()
        self.xml_deal()

