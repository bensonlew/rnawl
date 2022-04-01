# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2018/11/16'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,glob
from mbio.packages.metabolome.common import link_file
from biocluster.core.exceptions import OptionError


class TtssAgent(Agent):
    def __init__(self, parent):
        super(TtssAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "query_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "module", "type": "string", "default": "all"},
            {"name": "threshold", "type": "string", "default": "sensitive"},
            # "sensitive" (equal cutoff=0.95) "selective" (equal cutoff=0.9999) or "cutoff=float" or "float"
            {"name": "mem", "type": "int", "default": 10}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("query").is_set and not self.option("query_dir").is_set:
            raise OptionError("must type into query option", code="33301001") # 必须设置输入序列文件

        if self.option("module") not in ["all", "plant", "animal"]:
            raise OptionError("module not exits: %s" , variables=( self.option("module")), code="33301002")
        if self.option("threshold") not in ["sensitive", "selective"]:
            if self.option("threshold").startswith("cutoff="):
                test_str = self.ooption("threshold").lstrip("cutoff=")
                try:
                    float(test_str)
                except:
                    raise OptionError("threshold option error:%s" , variables=( self.option("threshold")), code="33301003")
            else:
                try:
                    if not 0 < float(self.option("threshold")) < 1:
                        raise OptionError("threshold option error:%s" , variables=( self.option("threshold")), code="33301004")
                except:
                    raise OptionError("threshold option error:%s" , variables=( self.option("threshold")), code="33301005")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        if self.option("query").is_set:
            number = os.path.getsize(self.option('query').prop['path']) / float(1024*1024*1024)
            if number > 1:
                total_memory = int(45 * number)
                self._memory = str(total_memory) + "G"
            else:
                self._memory = str(self.option('mem')) + "G"
        else:

            self._memory = str(self.option('mem')) + "G"


class TtssTool(Tool):
    def __init__(self, config):
        super(TtssTool, self).__init__(config)
        self._version = "1.0"
        self.ttss_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/TTSS/TTSS_GUI-1.0.1.jar'
        self.java_path = self.config.SOFTWARE_DIR +  "/program/sun_jdk1.8.0/bin/java"
        self.ttss_sh = 'bioinfo/metaGenomic/TTSS/ttss.sh'
        self.module_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/TTSS/module'
        if self.option("query").is_set:
            number = os.path.getsize(self.option('query').prop['path']) / float(1024*1024*1024)
            if number > 1:
                total_memory = int(45 * number)
                self._memory = str(total_memory) + "G"
            else:
                self._memory = str(self.option('mem')) + "G"
        else:

            self._memory = str(self.option('mem')) + "G"


    def run_ttss(self):
        """
        description
        :return:
        """
        if self.option("module") == "all":
            module = "TTSS_STD-2.0.2.jar"
        elif self.option("module") == "plant":
            module = "TTSS_PLANT-1.0.1.jar"
        elif self.option("module") == "animal":
            module = "TTSS_ANIMAL-1.0.1.jar"
        from_path = os.path.join(self.module_path, module)
        to_path = os.path.join(self.work_dir, 'module', module)
        if not os.path.isfile(to_path):
            os.makedirs(os.path.dirname(to_path))
            os.link(from_path, to_path)
        processfile = os.path.join(self.work_dir, "ttss_predict_tmp")
        self.java_memory = int(15 * (float(self._memory.strip("G")) / 45)) #add by qingchen.zhang 20191112
        if self.option("query").is_set:
            # cmd = "{} {} -d64 -Xmx{}G -jar {} -f {} -m {} -t {} -o {} -q".\
            cmd = "{} {} {} {} {} {} {} {}".\
            format(self.ttss_sh, self.java_path, self.java_memory, self.ttss_path, self.option("query").path, module, self.option('threshold'), processfile)
            # param2 is xms memory
            command = self.add_command('predict', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("effectiveT3 sucess")
            elif command.return_code in [1]:
                self.java_memory += 15
                self.add_state('memory_limit', 'memory is low!') #add by qingchen.zhang 20191112
            else:
                self.set_error("effectiveT3 error", code="33301001")
        else:
            count = 0
            command_list = []
            file_list = os.listdir(self.option("query_dir").path)
            for file in file_list:
                command_max = 5
                if command_max * 7 > self.option("mem") + 10 * self.version:
                    # self.memory(20 * 7)
                    self.add_state("memory_limit", "memory is low!")
                count += 1
                file_path = os.path.join(self.option("query_dir").path, file)
                cmd = "{} java {} {} {} {} {} {}".\
                format(self.ttss_sh, 3, self.ttss_path, file_path, module, self.option('threshold'), processfile + str(count))
                command = self.add_command('predict' + str(count), cmd).run()
                command_list.append(command)
                if len(command_list) == command_max or count == len(file_list):
                    self.wait(*command_list)
                    for one in command_list:
                        if one.return_code == 0:
                            self.logger.info("ttss %s success" % one.name)
                        else:
                            self.set_error("ttss %s error" , variables=( one.name), code="33301002")
                    command_list = []


    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        result = os.path.join(self.output_dir, "ttss_predict.txt")
        self.get_result(result)
        self.option('result', result)

    def get_result(self, outfile):
        with open(outfile, "w") as ot:
            ot.write("gene_id\tscore\n")
        for file in glob.glob(self.work_dir + "/ttss_predict_tmp*"):
            with open(file, "r") as f, open(outfile, 'a') as o:
                lines = f.readlines()
                for line in lines[1:-1]:
                    line = line.strip().split(';')
                    query = line[0].rsplit("_1", 1)[0]
                    if line[3] == 'true':
                        o.write("%s\t%s\n" % (query, line[2]))

    def run(self):
        super(TtssTool, self).run()
        self.run_ttss()
        self.set_output()
        self.end()
