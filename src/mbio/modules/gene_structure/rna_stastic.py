# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
# last_modify:2017.2.17

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class RnaStasticModule(Module):
    def __init__(self,work_id):
        super(RnaStasticModule, self).__init__(work_id)
        self.step.add_steps('rna_editing')
        options = [
            {"name": "rna_bam_dir", "type": "infile", "format": "align.bwa.bam_dir"},
            {"name": "ref_hg19.fa", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)
        self.editing = self.add_tool("ref_rna.gene_structure.rnaediting")
        self.tools =[]
        self._end_info = 0

    def check_options(self):
        if not self.option("rna_bam_dir").is_set:
            raise OptionError("必须提供bam文件！")
        if not self.option("ref_hg19.fa").is_set:
            raise OptionError("请提供参考基因组文件")

        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start() 
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def editing_run(self):
        self.editing.set_options({
            "rna_bam_dir": self.option("rna_bam_dir"),
            "ref_hg19.fa": self.option("ref_hg19.fa"),
        })
        self.editing.on('start', self.set_step, {'start': self.step.rna_editing})
        self.editing.on('end', self.set_step, {'end': self.step.rna_editing})
        self.editing.on('end', self.set_output, 'editing')
        self.editing.run()

    def finish_update(self,event):
        step = getattr(self.step,event['data'])
        step.finish()
        self.step.update()

    def stastic_run(self):
        n=0
        data = os.path.join(self.work_dir, "Rnaediting/output/edit_reassin.txt")
        for i in [1,2]:
            self.step.add_steps('stastic_{}'.format(n))
            stastic = self.add_tool("ref_rna.stastic")
            stastic.set_options({
                "data": data,
                "row": i,
            })
            step = getattr(self.step, 'stastic_{}'.format(n))
            step.start()
            stastic.on('end',self.finish_update,'stastic_{}'.format(n))
            self.tools.append(stastic)
            self.on_rely(self.tools, self.end)
            n = n+1
        for tool in self.tools:
            tool.run()



    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
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

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'editing' or event['data'] == 'stastic':
            self.linkdir(obj.output_dir, self.output_dir)


    def run(self):
        super(RnaStasticModule, self).run()
        self.editing.on('end', self.stastic_run)
        self.editing_run()


    def end(self):
        repaths = [
            [".", "", "转录因子分析结果输出目录"],
            ["TF_result.txt", "txt", "分析结果文件信息"],
            ["stastic_result.txt", "txt", "统计结果信息"],
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        # sdir.add_regexp_rules(regexps)
        super(RnaStasticModule, self).end()
