#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import glob
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile
import unittest


class MethylationModule(Module):
    """
    三代数据甲基化位点预测
    last_modify:
    """

    def __init__(self, work_id):
        super(MethylationModule, self).__init__(work_id)
        options = [
            {"name": "methy_dir", "type": "infile", "format": "bacgenome.methy_dir"}, # 包含多个样本的h5文件以及list 样本对用关系
            {"name": "ref_fa", "type": "string", "default": ""},
            #{"name": "ref_list", "type": "infile", "format": "bacgenome.methy_file"},
            #{"name": "ref_input", "type": "infile", "format": "bacgenome.methy_file","check":"check_ref_file"},
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        #self.methy = self.add_tool('bacgenome.methylation')
        self.tool_list=[]
        self.sam_path = {}


    def check_options(self):
        if not self.option("methy_dir"):
            raise OptionError("必须输入三代数据文件夹")
        if not self.option("ref_fa"):
            raise OptionError("必须设置输入参考序列")
        return True

    def run_each_methy(self):
        """
        run methylation of each sample
        :return:
        """
        self.logger.info("start run methylation")
        #self.get_sample_ref()
        file_dict = self.option("methy_dir").prop["filedict"]
        methy_samples = self.option("methy_dir").prop["methy_sample"]
        methy_files = self.option("methy_dir").prop["methy_files"]
        for each in methy_samples:
            self.logger.info(each)
            methy_file_type = self.option("methy_dir").get_suffix(file_dict[each])
            if methy_file_type == 'fq':
                continue
            self.methy = self.add_tool('bacgenome.methylation')
            # if self.ref_file.has_key(each):
            #     eachref = self.ref_file[each]
            # else:
            #     self.set_error("sample is diff in ref_list and methy_dir")
            options = {
                'input': methy_files[each],
                'sample': each,
                'input_type': methy_file_type,
                "ref_input" : self.option('ref_fa')
            }
            self.methy.set_options(options)
            self.tool_list.append(self.methy)
        if len(self.tool_list) > 1:
            self.on_rely(self.tool_list, self.set_output)
        elif len(self.tool_list) == 1:
            self.tool_list[0].on('end', self.set_output)
        else:
            self.logger.info("没有甲基化的原始数据，跳过甲基化tool的运行")
            self.end()
        for tool in self.tool_list:
            tool.run()

    def get_sample_ref(self):
        ref_list = self.option("ref_list").prop["path"]
        self.ref_sample = self.option("ref_list").prop["ref_sample"]
        self.ref_file = self.option("ref_list").prop["ref_dict"]

    def run(self):
        super(MethylationModule, self).run()
        self.run_each_methy()
        #self.set_output()


    def set_output(self):
        self.logger.info("设置结果目录")
        for tool in self.tool_list:  #zouguanqing
            files = os.listdir(tool.output_dir)
            for dir in files:
                if os.path.exists(self.output_dir + '/' + dir):
                    shutil.rmtree(self.output_dir + '/' + dir)
                for file in os.listdir(tool.output_dir + '/' + dir):
                    if os.path.exists(self.output_dir + '/' + file):
                        os.remove(self.output_dir + '/' + file)
                    os.link(tool.output_dir + '/'+dir+'/'+file, self.output_dir + '/' + file)
        self.end()



    def end(self):
        super(MethylationModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "methylation_module",
            "type": "module",
            "name": "bacgenome.methylation",
            "instant": False,
            "options": dict(
                #methy_dir="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/smrt/pacbio_dir", ##情况1
                methy_dir="/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/bac/bac_update/methylation", ##情况2
                #ref_list="/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/smrt/fasta_file/ref.list",  #yuan
                ref_fa = "/mnt/ilustre/users/sanger-dev/sg-users/yuanshaohua/measurement/smrt/fasta_file/reference.fasta"  ##zou
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()