# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200411

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import gevent
import re
import os
import unittest
import subprocess
from biocluster.config import Config


class TransCoordModule(Module):
    """
    该module为了基因组坐标转化功能，针对没有已经chain文件的分析项目开发的module,其目的功能是通过给定的
    已知fasta和目标fasta，通过互相比对生成chain文件，后通过liftover进行分析，来生成最终的结果文件
    """
    def __init__(self, work_id):
        super(TransCoordModule, self).__init__(work_id)
        options = [
            {"name": "raw_fasta", "type": "infile", "format": "ref_rna_v2.common"}, #当前基因组
            {"name": "target_fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 目标基因组
        ]
        self.add_option(options)
        self.file_list = []
        self.fa2chain_tools = []
        self.chain_nets = []


    def check_options(self):
        if not self.option("raw_fasta").is_set:
            raise OptionError("缺少raw_fasta参数",)
        if not self.option("target_fasta"):
            raise OptionError("缺少target_fasta参数")
        return True



    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def linkdir(self, dirpath, dirname):
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
                # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'fasta_split':
            self.linkdir(obj.output_dir, 'fasta_split')
        elif event['data'] == 'info':
            self.linkdir(obj.output_dir, 'info')
        elif event['data'] == 'fa2chain':
            self.linkdir(os.path.join(obj.output_dir,"chain"), 'chain')
            self.linkdir(os.path.join(obj.output_dir, "lastz"), 'lastz')
            self.linkdir(os.path.join(obj.output_dir, "psl"), 'psl')
        elif event['data'] == 'chain_merge':
            self.linkdir(obj.output_dir, 'new_chain')
        elif event['data'] == 'chain_net':
            self.linkdir(obj.output_dir, 'net')
        elif event['data'] == "netchainsubset":
            final_file=os.path.join(self.output_dir,"lift.chain")
            if os.path.exists(final_file):
                    os.remove(final_file)
            result_file = os.path.join(self.netchainsubset.output_dir,"lift.chain")
            os.link(result_file, final_file)
            self.end()
        else:
            pass

    def run_fasta_split(self):
        self.fasta_split=self.add_tool("tool_lab.trans_coord.fasta_split_sgene")
        options = {
            "fasta_file": self.option("raw_fasta").prop["path"],
        }
        self.fasta_split.set_options(options)
        self.fasta_split.on("end", self.run_prepare_tools)
        self.fasta_split.on("end", self.set_output, "fasta_split")
        self.logger.info("split_fasta is running!")
        self.fasta_split.run()

    def run_prepare_tools(self):
        self.prepare_tool=self.add_tool("tool_lab.trans_coord.prepare4trans")
        options = {
            "raw_fasta_file": self.option("raw_fasta").prop["path"],
            "target_fasta_file": self.option("target_fasta").prop["path"]
        }
        self.prepare_tool.set_options(options)
        self.prepare_tool.on("end", self.run_fa2chain)
        self.prepare_tool.on("end", self.set_output, "info")
        self.logger.info("prepare_tool is running!")
        self.prepare_tool.run()


    def run_fa2chain(self):
        singlefas=os.listdir(os.path.join(self.fasta_split.output_dir))
        for fasta_file in singlefas:
            fa2chain=self.add_tool("tool_lab.trans_coord.fa2chain")
            fasta_path=os.path.join(os.path.join(self.work_dir,"fasta_split"),fasta_file)
            options = {
                "fasta_file": fasta_path,
                "info_dir":os.path.join(self.prepare_tool.output_dir)
            }
            fa2chain.set_options(options)
            self.fa2chain_tools.append(fa2chain)
        for j in range(len(self.fa2chain_tools)):
            self.fa2chain_tools[j].on('end', self.set_output, 'fa2chain')
        if self.fa2chain_tools:
            if len(self.fa2chain_tools) > 1:
                self.on_rely(self.fa2chain_tools, self.run_chain_merge)
            elif len(self.fa2chain_tools) == 1:
                self.fa2chain_tools[0].on('end', self.run_chain_merge)
            for tool in self.fa2chain_tools:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("fa2chain_tools列表为空！")


    def run_chain_merge(self):
        self.chain_merge = self.add_tool("tool_lab.trans_coord.chain_merge")
        options = {
            "chain_dir": os.path.join(self.work_dir,"chain")
        }
        self.chain_merge.set_options(options)
        self.chain_merge.on("end", self.set_output, "chain_merge")
        self.chain_merge.on("end", self.run_chain_net)
        self.logger.info("prepare_tool is running!")
        self.chain_merge.run()

    def run_chain_net(self):
        new_chains = os.listdir(os.path.join(self.chain_merge.output_dir,"chain"))
        for chain_file in new_chains:
            chain_net = self.add_tool("tool_lab.trans_coord.chain_net")
            chain_path = os.path.join(os.path.join(self.chain_merge.output_dir,"chain"), chain_file)
            options = {
                "chain_file": chain_path,
                "info_dir": os.path.join(self.work_dir, "info")
            }
            chain_net.set_options(options)
            self.chain_nets.append(chain_net)
        for j in range(len(self.chain_nets)):
            self.chain_nets[j].on('end', self.set_output, 'chain_net')
        if self.chain_nets:
            if len(self.chain_nets) > 1:
                self.on_rely(self.chain_nets, self.run_netchainsubset)
            elif len(self.chain_nets) == 1:
                self.chain_nets[0].on('end', self.run_netchainsubset)
            for tool in self.chain_nets:
                gevent.sleep(1)
                tool.run()
        else:
            self.set_error("fa2chain_tools列表为空！")

    def run_netchainsubset(self):
        self.netchainsubset=self.add_tool("tool_lab.trans_coord.netchainsubset")
        options = {
            "net_dir": os.path.join(self.work_dir,"net"),
            "target_chain": os.path.join(self.work_dir, "new_chain","target.chain")
        }
        self.netchainsubset.set_options(options)
        self.netchainsubset.on("end", self.set_output, "netchainsubset")
        self.logger.info("prepare_tool is running!")
        self.netchainsubset.run()


    def run(self):
        super(TransCoordModule, self).run()
        self.logger.info("开始切割")
        self.run_fasta_split()

    def set_bam_list(self):
        """
        samtools和freebayes与gatk的bamlist格式不一致，这里进行处理下
        """
        self.logger.info("开始设置bamlist！")
        with open(self.option("bam_list"), "r") as r, open("{}/bam.list".format(self.work_dir), "w") as w:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                w.write(temp[1] + "\n")
        self.logger.info("设置bamlist成功！")

    def end(self):
        super(TransCoordModule, self).end()



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
            "id": "trans_coord" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.trans_coord",
            "instant": False,
            "options": dict(
                # raw_fasta="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM14920v2.fasta",
                # target_fasta="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ASM18445v3.fasta"
                raw_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                target_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.fa"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()




if __name__ == '__main__':
    unittest.main()