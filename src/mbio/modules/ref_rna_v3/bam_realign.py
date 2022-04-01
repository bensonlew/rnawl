# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20190220

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import json
import unittest

class BamRealignModule(Module):
    """
     bam_realign module，用于将bam文件分别以tool的形式投出去。
    """
    def __init__(self, work_id):
        super(BamRealignModule, self).__init__(work_id)
        options = [
            #{"name": "bam_list", "type": 'infile', "format": "wgs_v2.bam_list"},sequence.fasta,wgs_v2.bam_list,align.bwa.bam_dir
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": 'bam_realign_list', "type": "outfile", 'format': "ref_rna_v2.common"},
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam_dir"},  # bam格式文件
        ]
        self.add_option(options)
        self.bam_realign_tools = []
        self.samples = []
        self.picards = []

    def check_options(self):
        if not self.option("in_bam"):
            raise OptionError("请输入bam文件夹", code="25600702")
        return True

    def picard_run(self):
        for i in os.listdir(self.option("in_bam").prop["path"]):
            self.samples.append(i[:-4])
            f_path = os.path.join(self.option("in_bam").prop["path"], i)
            self.logger.info(f_path)  # 打印出f_path的信息，是上一步输出文件的路径
            picard = self.add_tool('ref_rna_v2.picard_rna')
            # picard = self.add_tool('whole_transcriptome.snp.samtools_presnp')
            picard.set_options({
                "in_bam": f_path
            })
            self.picards.append(picard)
            #picard.on("end", self.gatk_run, i[:-4])  # event["data"]传为样本名称
        for j in range(len(self.picards)):
            self.picards[j].on("end", self.set_output, 'picards')
        if len(self.picards) == 1:
            self.picards[0].on("end", self.run_bam_realign)
        else:
            self.on_rely(self.picards, self.run_bam_realign)
        for picard in self.picards:
            picard.run()

    def run_bam_realign(self):
        with open(os.path.join(self.work_dir,"bam_list"),"w") as blst:
          for i in os.listdir(os.path.join(self.work_dir,"picards")):
            if i.endswith("bam"):
                bam_path=os.path.join(os.path.join(self.work_dir,"picards"),i)
                blst.write(i.split(".")[0]+"\t"+bam_path+"\n")
        self.bam_list=os.path.join(self.work_dir,"bam_list")
        with open(self.bam_list) as f:
            lines = f.readlines()
            for line in lines:
                sample_name = line.strip().split("\t")[0]
                sample_path = line.strip().split("\t")[1]
                bam_realign = self.add_tool("ref_rna_v2.bam_realign")
                options = ({
                    "bam_file": sample_path,
                    "fa_file": self.option("fa_file").prop["path"],
                    "name": sample_name + ".realign.bam"
                })
                bam_realign.set_options(options)
                self.bam_realign_tools.append(bam_realign)
            for j in range(len(self.bam_realign_tools)):
                self.bam_realign_tools[j].on("end", self.set_output, 'bam_realign_dir')
            if self.bam_realign_tools:
                if len(self.bam_realign_tools) > 1:
                    self.on_rely(self.bam_realign_tools, self.end)
                elif len(self.bam_realign_tools) == 1:
                    self.bam_realign_tools[0].on('end', self.end)
            else:
                self.set_error("bam_realign_tools列表为空！", code="25600702")
            for tool in self.bam_realign_tools:
                gevent.sleep(1)
                tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'bam_realign_dir':
            self.linkdir(obj.output_dir, 'bam_realign_dir')
        if event['data'] == 'picards':
            self.linkdirfors(obj.output_dir, 'picards')
        else:
            pass

    def linkdirfors(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
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
                os.link(oldfiles[i], newdir)


    def linkdir(self, dirpath, dirname):
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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(BamRealignModule, self).run()
        self.picard_run()

    def end(self):
        self.set_realign_list()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(BamRealignModule, self).end()

    def set_realign_list(self):
        if len(self.bam_realign_tools) > 1:
            dir_path = self.output_dir + "/bam_realign_dir"
        else:
            dir_path = self.output_dir + "/bam_realign_dir"
        with open(self.output_dir + "/bam.list", 'w') as w:
            for m in os.listdir(dir_path):
                n = re.match('(.*)\.realign\.bam$', m)
                if n:
                    w.write('{}\t{}\n'.format(n.group(1), dir_path + "/{}".format(m)))
        self.option("bam_realign_list").set_path(self.output_dir + "/bam.list")

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "bam_realign" + str(random.randint(1, 10000))+"samtools",
            "type": "module",
            "name": "ref_rna_v3.bam_realign",
            "instant": False,
            "options": dict(
                #ref_dict=test_dir + "/" + "Mus_musculus.GRCm38.dna_rm.toplevel.clean.dict",
                #in_bam="/mnt/ilustre/users/sanger-dev/workspace/20190704/Snp_tsg_34708_8507_5413/bam_folder",
                in_bam="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test1128/testdata/bam",
                #bam_list="/mnt/ilustre/users/sanger-dev/workspace/20190618/Single_snp4914/PicardRna/output/bam_list",
                fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.fa",
                #call_type="sentieon",
                #ref_fasta=test_dir+"/"+"Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                # scm="complete",
                # scd="correlation",
                # corr_method='pearson',
                #output=None,
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

