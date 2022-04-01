#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import re
from mbio.packages.alpha_diversity.estimator_size import est_size
##################################
from mbio.files.meta.otu.otu_table import OtuTableFile
#######################################


class EstimatorsAgent(Agent):
    """
    estimators:用于生成所有样本的指数表
    version 1.0
    author: qindanhua
    last_modify: 2016.05.06
    """
    #ESTIMATORS = ['ace', 'bergerparker', 'boneh', 'bootstrap', 'bstick', 'chao', 'coverage', 'default', 'efron','geometric', 'goodscoverage', 'heip', 'invsimpson', 'jack', 'logseries', 'npshannon', 'nseqs','qstat', 'shannon', 'shannoneven', 'shen', 'simpson', 'simpsoneven', 'smithwilson', 'sobs', 'solow']
    ESTIMATORS = ['ace', 'bergerparker', 'boneh', 'bootstrap', 'bstick', 'chao', 'coverage', 'default', 'efron','geometric', 'goodscoverage', 'heip', 'invsimpson', 'jack', 'logseries', 'npshannon', 'nseqs','qstat', 'shannon', 'shannoneven', 'shen', 'simpson', 'simpsoneven', 'smithwilson', 'sobs', 'solow','pd']

    def __init__(self, parent):
        super(EstimatorsAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table,meta.otu.tax_summary_dir"},  # 输入文件
            {"name": "indices", "type": "string", "default": "ace,chao,shannon,simpson"},  # 指数类型
            {"name": "level", "type": "string", "default": "otu"},  # level水平
            # {"name": "estimators", "type": "outfile", "format": "meta.alpha_diversity.estimators"}  # 输出结果
            #########################################added 1 line by yiru 20170421
            {"name": "newicktree","type": "infile", "format": "meta.beta_diversity.newick_tree"} # pd指数计算时需要的树文件
            ########################################
        ]
        self.add_option(options)
        self.step.add_steps('estimators')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.estimators.start()
        self.step.update()

    def step_end(self):
        self.step.estimators.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if self.option("level") not in ['otu', 'domain', 'kindom', 'phylum', 'class', 'order',
                                        'family', 'genus', 'species']:
            raise OptionError("请选择正确的分类水平", code="32700101")
        for estimators in self.option('indices').split(','):
            if estimators not in self.ESTIMATORS:
                raise OptionError("error:%s,请选择正确的指数类型", variables=(estimators), code="32700102")
        ########################## added 8 lines by yiru 20170425 
        IdxList = self.option("indices").split(",")
        if "pd" in IdxList:
            if not self.option('newicktree').is_set:
                raise OptionError('选择pd指数时必须提供进化树文件', code="32700103")
        if not self.option('otu_table').is_set:
            raise OptionError('必须提供输入文件', code="32700104")
        else:
            otulist = [line.split('\t')[0] for line in open(self.gettable().prop['path'])][1:]  # 获取所有OTU/物种名
        ##########################

    ##################################################added 7 lines by yiru 20170421
    def gettable(self):
        """
        根据level返回进行计算的otu表-对象
        :return tablepath:
        """
        if self.option('otu_table').format == "meta.otu.tax_summary_dir":
            newtable = OtuTableFile()
            newtable.set_path(self.option('otu_table').get_table(self.option('level')))
            return newtable
        else:
            return self.option('otu_table')
    ###############################################

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 11
        self._memory = '30G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        # print self.get_upload_files()
        super(EstimatorsAgent, self).end()


class EstimatorsTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(EstimatorsTool, self).__init__(config)
        self.indices = '-'.join(self.option('indices').split(','))
        self.special_est = ['boneh', 'efron', 'shen', 'solow']
        #################################################added 2 lines by yiru 20170425
        self.real_otu = self.gettable()  # 获取真实的OTU表路劲
        self.biom = self.biom_otu_table()  # 传入otu表需要转化为biom格式
        ######################################################
        self.otu_table = ''


    def shared(self):
        """
        执行生成shared文件的perl脚本命令
        """
        self.otu_table = self.option("otu_table").prop['path']
        if self.option("otu_table").format == "meta.otu.tax_summary_dir":
            self.otu_table = self.option("otu_table").get_table(self.option("level"))
        self.logger.info("otutable format:{}".format(self.option("otu_table").format))
        self.logger.info("转化otu_table({})为shared文件({})".format(self.otu_table, "otu.shared"))
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR+"/bioinfo/meta/scripts/otu2shared.pl "+" -i "+self.otu_table +
                                    " -l 0.97 -o " + self.option("level")+".shared", shell=True)
            self.logger.info("OK")
            return True
        except subprocess.CalledProcessError:
            self.logger.info("转化otu_table到shared文件出错")
            return False

    def mothur_check(self, command, line):
        if re.match(r"\[ERROR\]:", line):
            command.kill()
            self.set_error("mothur命令报错", code="32700101")

    def mothur(self):
        """
        运行mothur软件生成各样本指数表
        """
        indices = self.indices.split('-')
        cmd = ""   #guanqing.zou 20180515 解决pd单独运行问题
        for index in indices:
            if index != "pd": 
                cmd = 'bioinfo/meta/mothur-1.30/mothur.1.30 "#summary.single(shared=%s.shared,groupmode=f,calc=%s)"' % (self.option('level'),self.indices)
                if index in self.special_est:
                    size = est_size(self.otu_table)
                    cmd = 'bioinfo/meta/mothur-1.30/mothur.1.30 "#summary.single(shared=%s.shared,groupmode=f,calc=%s,size=%s)"' % \
                    (self.option('level'), self.indices, size)
            else:
                self.run_pd()

        ###guanqing.zou 20180515 解决pd单独运行问题
        if cmd =="":         
            self.get_est_table()

        else:
            self.logger.info(cmd)
            self.logger.info("开始运行mothur")
            command = self.add_command("mothur", cmd)
            command.run()
            self.wait(command)
            if command.return_code == None:
                command.rerun()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行mothur完成")
                    self.get_est_table()
                else:
                    self.set_error("运行mothur运行出错!", code="32700102")
                    return False
            elif command.return_code == 0:
                self.logger.info("运行mothur完成")
                self.get_est_table()           
            else:
                self.set_error("运行mothur运行出错!", code="32700103")
                return False
    
    #####################################################added 4 functions by yiru 20170425
    def gettable(self):
        """
        根据level返回进行计算的otu表路径
        :return tablepath:
        """
        if self.option('otu_table').format == "meta.otu.tax_summary_dir":
            return self.option('otu_table').get_table(self.option('level'))
        else:
            return self.option('otu_table').path

    def run_pd(self):
        cmd = 'program/Python/bin/alpha_diversity.py'
        cmd += ' -i %s -m PD_whole_tree -o adiv_chao1_pd.txt -t %s' % (self.biom,self.option("newicktree").prop['path'])
        self.logger.info("开始pd运算")
        self.logger.info(cmd)
        command = self.add_command("run_pd", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("pd运算完成")         
        else:
            self.set_error("运行mothur运行出错!", code="32700104")
            return False

    def get_est_table(self):
        cmd = 'program/Python/bin/python {}/bioinfo/meta/scripts/make_estimators_table.py'.format(self.config.SOFTWARE_DIR)
        command = self.add_command("get_est_table", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("OK")
            self.set_output()
        else:
            self.set_error("生成estimate文件出错!", code="32700105")
            return False

    def biom_otu_table(self):
        """
        将otutable转化成biom格式
        :return biom_path:返回生成的biom文件路径
        """
        if self.option('otu_table').format == "meta.otu.tax_summary_dir":
            # 如果提供的是otu分类文件夹，需要重新创建类，再使用类的convert_to_biom方法
            newtable = OtuTableFile()
            newtable.set_path(self.real_otu)
            newtable.check()
        else:
            newtable = self.option('otu_table')
        newtable.get_info()
        biom_path = os.path.join(self.work_dir, 'temp.biom')
        if os.path.isfile(biom_path):
            os.remove(biom_path)
        newtable.convert_to_biom(biom_path)
        return biom_path
    ##############################################################
    def set_output(self):
        """
        将结果文件链接至output
        """
        self.logger.info("set output")
        if len(os.listdir(self.output_dir)) != 0:
            os.remove(self.output_dir+'/estimators.xls')
            os.link(self.work_dir+'/estimators.xls', self.output_dir+'/estimators.xls')
        else:
            os.link(self.work_dir+'/estimators.xls', self.output_dir+'/estimators.xls')
            # self.option('estimators').set_path(self.output_dir+'/estimators')
        self.logger.info("done")
        self.end()

    ##############################################################run function modified by yiru 20170421
    def run(self):
        super(EstimatorsTool,self).run()
        if self.shared():
            self.mothur()
        else:
            self.set_error("shared运行出错", code="32700106")
    ###############################################################
    #def run(self):
        #"""
        #运行
        #"""
        #super(EstimatorsTool, self).run()
        #if self.shared():
            #if self.mothur():
                #self.end()
        #else:
            #self.set_error("shared运行出错!")

