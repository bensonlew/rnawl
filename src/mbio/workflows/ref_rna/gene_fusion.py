# -*- coding: utf-8 -*-
# __author__ = 'shijin'

"""有参转录组基础分析，无拼接部分"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import shutil
from mbio.files.meta.otu.group_table import GroupTableFile


class GeneFusionWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        """
        self._sheet = wsheet_object
        super(GeneFusionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fq_type", "type": "string"},  # PE OR SE

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.qc = self.add_module("ref_rna.qc.quality_control")
        self.gene_fusion = self.add_module("ref_rna.gene_fusion.gene_fusion")
        self.step.add_steps("qcstat","genefusion")
        
    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("fq_type") in ["PE","SE"]:
            raise OptionError("测序类型不正确")
        return True
        
    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()    
        
    def run_qc(self):
        self.qc.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type')
        })
        # self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.qcstat})
        self.qc.on('end', self.set_step, {'end': self.step.qcstat})
        self.qc.run()
        

    def run_gene_fusion(self):
        self.gene_fusion.set_options({
            "sickle_dir":self.qc.option("sickle_dir"),
            "seqprep_dir":self.qc.option("seqprep_dir"),
            "fastq_dir":self.option('fastq_dir'),
            'fq_type': self.option('fq_type')
        })
        self.gene_fusion.on('start', self.set_step, {'start': self.step.genefusion})
        self.gene_fusion.on('end', self.set_step, {'end': self.step.genefusion})
        self.gene_fusion.run()
        # self.gene_fusion.on("end",self.set_output,"gene_fusion")
        
        
    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            if mode == 'link':
                shutil.copytree(olddir, newdir, symlinks=True)
            elif mode == 'copy':
                shutil.copytree(olddir, newdir)
            else:
                raise Exception('错误的移动文件方式，必须是\'copy\'或者\'link\'')
        else:
            allfiles = os.listdir(olddir)
            oldfiles = [os.path.join(olddir, i) for i in allfiles]
            newfiles = [os.path.join(newdir, i) for i in allfiles]
            self.logger.info(newfiles)
            for newfile in newfiles:
                if os.path.isfile(newfile) and os.path.exists(newfile):
                    os.remove(newfile)
                elif os.path.isdir(newfile) and os.path.exists(newfile):
                    shutil.rmtree(newfile)
            for i in range(len(allfiles)):
                if os.path.isfile(oldfiles[i]):
                    os.system('cp {} {}'.format(oldfiles[i], newfiles[i]))
                else:
                    os.system('cp -r {} {}'.format(oldfiles[i], newdir))
                    
    def set_output(self, event):
        obj = event["bind_object"]
        # 设置qc报告文件
        if event['data'] == 'qc':
            self.move2outputdir(obj.output_dir, 'QC_stat')
        
            
    def run(self):
        self.qc.on("end",self.run_gene_fusion)
        self.gene_fusion.on("end",self.end)
        self.run_qc()
        super(GeneFusionWorkflow, self).run()
        
    def end(self):
        super(GeneFusionWorkflow, self).end()
