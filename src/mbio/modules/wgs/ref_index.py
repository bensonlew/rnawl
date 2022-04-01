# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20180414

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os


class RefIndexModule(Module):
    """
    用于构建转基因的参考基因组文件
    last modified by hongdong@20180521

    """
    def __init__(self, work_id):
        super(RefIndexModule, self).__init__(work_id)
        options = [
            {"name": "insert_fa", "type": "infile", "format": "sequence.fasta"},  # 用于wgs.bwa_index
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},    # 用于wgs.bwa_index
            {"name": "dbtype", "type": 'string'}  # 用于wgs.makeblastdb
        ]
        self.add_option(options)
        self.bwa_index = self.add_tool("wgs.bwa_index")
        self.samtools_faidx = self.add_tool("wgs.samtools_faidx")
        self.makeblastdb = self.add_tool("wgs.makeblastdb")

    def check_options(self):
        if not self.option("insert_fa").is_set:
            raise OptionError("缺少insert_fa参数", code="24501201")
        if not self.option("ref_fa"):
            raise OptionError("缺少ref_fa参数", code="24501202")
        return True

    def bwa_index_run(self):
        self.bwa_index.set_options({
            "insert_fa": self.option("insert_fa"),
            "ref_fa": self.option("ref_fa")
        })
        # self.bwa_index.on('end', self.set_output, 'bwa_index')
        self.bwa_index.on('end', self.samtools_faidx_run)
        self.bwa_index.on('end', self.makeblastdb_run)
        self.bwa_index.run()

    def samtools_faidx_run(self):
        self.samtools_faidx.set_options({
            "pop_fa": self.bwa_index.output_dir + "/pop.fa",
        })
        # self.samtools_faidx.on('end', self.end)
        self.samtools_faidx.on('end', self.set_output, "samtools_faidx")
        self.samtools_faidx.run()

    def makeblastdb_run(self):
        self.makeblastdb.set_options({
            "pop_fa": self.bwa_index.output_dir + "/pop.fa",
            "dbtype": self.option("dbtype")
        })
        # self.makeblastdb.on('end', self.end)
        self.makeblastdb.on('end', self.set_output, 'makeblastdb')
        self.makeblastdb.run()

    def set_output(self, event):
        obj = event['bind_object']
        # if event['data'] == 'bwa_index':
        #     self.linkdir(obj.output_dir, self.output_dir)
        if event['data'] == 'samtools_faidx':
            self.linkdir(self.bwa_index.output_dir, self.output_dir)
        if event['data'] == 'bwa_index':
            self.linkdir(obj.output_dir, self.output_dir)

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
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def run(self):
        super(RefIndexModule, self).run()
        self.on_rely([self.makeblastdb, self.samtools_faidx], self.end)
        self.bwa_index_run()

    def end(self):
        super(RefIndexModule, self).end()
