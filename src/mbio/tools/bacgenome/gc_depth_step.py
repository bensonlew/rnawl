# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou hao.gao'
# version 1.0
# last_modify: 2018.8.14
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import shutil

class GcDepthStepAgent(Agent):
    """
    细菌基因组depth_gc评估
    """
    def __init__(self, parent):
        super(GcDepthStepAgent, self).__init__(parent)
        options = [
            #{"name": "sam", "type": "infile", "format": "align.bwa.sam_dir"},  # sam
            {"name": "sam", "type": "string"},
            {"name": "windl", "type": "string", "default": "1000,3000,5000,8000,10000"},  # 滑动窗口大小
            {"name":"ref","type":"infile","format":"sequence.fasta"}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sam"):
            raise OptionError("必须添加sam文件！", code="31403201")


    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GcDepthStepAgent, self).end()

class GcDepthStepTool(Tool):
    def __init__(self, config):
        super(GcDepthStepTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.R_path = 'program/R-3.3.1_gcc5.1/bin/Rscript'
        self.depth = self.work_dir + "/" + 'cov.depth.txt'
        self.samtools = "bioinfo/align/samtools-1.6/bin/samtools"
        self.bac_depth = self.config.PACKAGE_DIR + '/bacgenome/bac_depth.sh'

    def run_pick(self):
        #self.out_sam = self.option('sam').prop['path']
        self.out_sam = self.option('sam')
        self.logger.info('开始生成bam文件')
        self.out_bam = self.work_dir + '/map.bam'
        self.logger.info(self.out_bam)
        cmd_bam = "{} view -F 12 -S {} -b -o {}".format(self.samtools, self.out_sam, self.out_bam)
        self.logger.info(cmd_bam)
        command1 = self.add_command('to-bam', cmd_bam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成"%cmd_bam)
        else:
            self.set_error("%s运行失败", variables=(cmd_bam), code="31403201")

    def run_depth(self):
        self.logger.info('对bam文件进行排序')
        self.sort_bam = self.work_dir + '/sort'
        cmd = '{} sort {} -o {}'.format(self.samtools, self.out_bam, self.sort_bam)
        command = self.add_command('sort-bam',cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("%s运行完成"%cmd)
        else:
            self.set_error("%s运行失败", variables=(cmd), code="31403202")

        self.logger.info('开始计算覆盖度')
        self.tmp_sh = self.work_dir + '/tmp.sh'
        self.depth = self.work_dir + '/depth.txt'
        #os.system('echo {} depth {} \> {} > {}'.format(self.config.SOFTWARE_DIR + '/' + self.samtools, self.sort_bam, self.depth, self.tmp_sh))
        self.tmp_sh = "{} {} {} {}".format(self.config.PACKAGE_DIR + '/bacgenome/samtools_depth.sh', self.config.SOFTWARE_DIR ,self.sort_bam, self.depth)
        cmd1 = '/program/sh {}'.format(self.tmp_sh)
        self.logger.info('1111111111')
        command1 = self.add_command('depth', cmd1).run()
        self.logger.info('bbbbbbbbb')
        self.wait(command1)
        self.logger.info('ccccccccc')
        if command1.return_code == 0:
            self.logger.info("%s运行成功"%cmd1)
        else:
            self.set_error("%s运行失败", variables=(cmd1), code="31403203")


    def change_depth_format(self):
        self.logger.info('开始转化depth文件格式')
        self.depth_final = self.work_dir + '/depth.final.out'
        fw = open(self.depth_final, 'w')
        with open(self.depth) as f:
            line1 = f.next()
            line1 = line1.strip()
            spline1= line1.split('\t')
            k = spline1[0]
            fw.write('>'+k+'\n'+spline1[2])
            for i in f:
                i = i.strip()
                spi = i.split('\t')
                if spi[0] == k:
                    fw.write(' '+spi[2])
                else:
                    k = spi[0]
                    fw.write('\n>'+k+"\n"+spi[2])
        fw.close()
        self.logger.info('转化depth文件格式完成，生成%s'%self.depth_final)


    def run_depth_gc(self):
        self.ref=self.option('ref').prop['path']
        wind_arrry = (self.option("windl")).split(',')
        self.logger.info("run_depth_gc运行开始")

        for wind in wind_arrry:
            out_dir = self.work_dir + "/" + 'depth_gc_' + str(wind) + "/"
            cmd = '{} {}gc_depth_v2.0.pl {} {} -windl {} -step 500 -outdir {}'.format(self.perl_path,self.perl_script,self.ref,self.depth_final,wind,out_dir)
            gc_depth =  'depth_gc_' + str(wind)
            command = self.add_command(gc_depth, cmd, ignore_error=True).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info(gc_depth + "运行完成")
            elif command.return_code == 255:
                self.logger.error("windl %s 结果为空，跳过此后的运行" % wind)
                break
            else:
                self.set_error("%s运行出错!", variables=(gc_depth), code="31403204")
            cmd_1 = '%s %s' % (self.R_path, out_dir + "depth_gc.cmd.r")
            self.logger.info(cmd_1)
            gc_depth_r = 'depth_gc_' + str(wind) +'_r'
            command1 = self.add_command(gc_depth_r, cmd_1, ignore_error=True).run()
            #command1.run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("R程序计算depth_gc成功")
            else:
                shutil.rmtree(out_dir)
                break
        self.logger.info("run_depth_gc运行结束")

    def set_output(self):
        #for j in 1000, 3000, 5000, 8000, 10000:
        for j in map(int, self.option('windl').split(',')):
            g1 = self.output_dir + "/" + "depth_gc_" + str(j) + "/"
            l1 = self.work_dir + "/" + "depth_gc_" + str(j) + "/"
            if not os.path.exists(l1) or not os.path.getsize(l1+'/gc_depth.wind'):
                break  # 没有结果不考贝
            if not os.listdir(l1):
                break  # 结果为空不考贝
            if os.path.exists(g1):
                pass
            else:
                os.mkdir(g1)
            for type in ['png', 'svg']:
                if os.path.exists(g1 + 'gc_depth.' + type):
                    os.remove(g1 + 'gc_depth.' + type)
                os.link(l1 + 'gc_depth.' + type, g1 + 'gc_depth.' + type)


    def run(self):
        super(GcDepthStepTool, self).run()
        if os.path.exists("sort.tmp.0000.bam"):
            os.system("rm sort.tmp.*.bam")
        self.run_pick()
        self.run_depth()
        self.change_depth_format()
        self.run_depth_gc()
        self.set_output()
        self.end()


