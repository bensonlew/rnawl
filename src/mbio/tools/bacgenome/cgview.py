#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
import xml.etree.ElementTree as et


class CgviewAgent(Agent):
    """
    用于扫描图的gbk文件生成
    version 1.0
    author: gaohao
    last_modify: 2018.04.08
    """

    def __init__(self, parent):
        super(CgviewAgent, self).__init__(parent)
        options = [
            {"name": "gen_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "rrna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "trna_gff", "type": "infile", "format": "gene_structure.gff3"},  #
            {"name": "pro_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta"},  #
            {"name": "anno", "type": "infile", "format": "sequence.profile_table"},  #
            {"name": "sample_name", "type": "string"},
            {"name": "analysis", "type": "string", "default": "uncomplete"},  ###流程分析模式complete，uncomplete
            {"name": "xml_file", "type": "infile", "format": "bacgenome.cgview_xml"},  #
            {"name": "seq_type_circle", "type": "string","default": "Circular"},  # 开闭环
            {"name": "seq_type", "type": "string"}, # 物种名称
            {"name": "species_name", "type": "string"}
        ]
        self.add_option(options)
        self.list =[]


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('xml_file').is_set:
            if not self.option('gen_gff').is_set:
                raise OptionError("请设置基因组基因预测基因gff文件！", code="31400801")
            if not self.option('rrna_gff').is_set:
                raise OptionError("请设置基因组rRNA预测基因gff文件！", code="31400802")
            if not self.option('trna_gff').is_set:
                raise OptionError("请设置基因组tRNA预测基因gff文件！", code="31400803")
            if not self.option('pro_fa').is_set:
                raise OptionError("请设置基因组基因蛋白序列文件！", code="31400804")
            if not self.option('genome_fa').is_set:
                raise OptionError("请设置基因组组装序列文件！", code="31400805")
            if not self.option('anno').is_set:
                raise OptionError("请设置基因组基因注释汇总表！", code="31400806")
            if not self.option('analysis'):
                raise OptionError("请提供分析流程类型！", code="31400807")
        #else:
        #    if not self.option('species_name'):
        #       raise OptionError("请设置基因圈图的物种名称！", code="31400808")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 5
        self._memory = '20G'

    def end(self):
        super(CgviewAgent, self).end()


class CgviewTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(CgviewTool, self).__init__(config)
        if not self.option('xml_file').is_set:
            self.genegff = self.option('gen_gff').prop['path']
            self.rrnagff = self.option('rrna_gff').prop['path']
            self.trnagff = self.option('trna_gff').prop['path']
            self.genomefa = self.option('genome_fa').prop['path']
            self.pro = self.option('pro_fa').prop['path']
            self.anno = self.option('anno').prop['path']
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        #self.xml_builder = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/cgview/cgview_xml_builder/cgview_xml_builder.pl"
        self.xml_builder = self.config.PACKAGE_DIR + "/bacgenome/cgview_xml_builder.pl"   #20190523 比起原来perl，增加了个性化的cog颜色添加功能
        self.cgview_jar = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/cgview/cgview.jar"
        self.java_path = "/program/sun_jdk1.8.0/bin/java"
        self.cog_color = self.config.SOFTWARE_DIR + '/database/COG/cog_color.list'
        self.list = []

        self.desc_map = {
            "A":"RNA processing and modification",
            "C":"Energy production and conversion",
            "B":"Chromatin structure and dynamics",
            "E":"Amino acid transport and metabolism",
            "D":"Cell cycle control, cell division, chromosome partitioning",
            "G":"Carbohydrate transport and metabolism",
            "F":"Nucleotide transport and metabolism",
            "I":"Lipid transport and metabolism",
            "H":"Coenzyme transport and metabolism",
            "K":"Transcription",
            "J":"Translation, ribosomal structure and biogenesis",
            "M":"Cell wall/membrane/envelope biogenesis",
            "L":"Replication, recombination and repair",
            "O":"Posttranslational modification, protein turnover, chaperones",
            "N":"Cell motility",
            "Q":"Secondary metabolites biosynthesis, transport and catabolism",
            "P":"Inorganic ion transport and metabolism",
            "S":"Function unknown",
            "U":"Intracellular trafficking, secretion, and vesicular transport",
            "T":"Signal transduction mechanisms",
            "W":"Extracellular structures",
            "V":"Defense mechanisms",
            "Y":"Nuclear structure",
            "Z":"Cytoskeleton",
            "R":"General function prediction only"
        }

    def run_gbk_un(self):
        cmd = "{} {}GBK_generation_uncomplete.pl {} {} {} {} {} {} {}".format(self.perl_path, self.perl_script,self.genegff,self.trnagff,self.rrnagff,self.pro ,self.genomefa,self.anno,self.option('sample_name'))
        self.logger.info(cmd)
        self.logger.info("开始运行run_gbk_uncomplete")
        command = self.add_command("run_gbk_uncomplete", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            os.system('sed \'s/linear/circular/g\' %s -i '% (self.work_dir + '/' + self.option('sample_name') + '.gbk'))
            self.logger.info("运行run_gbk_uncomplete完成")
        else:
            self.set_error("运行run_gbk_uncomplete运行出错!", code="31400801")

    def run_gbk(self):
        cmd = "{} {}GBK_generation_complete.pl {} {} {} {} {} {} ".format(self.perl_path, self.perl_script,self.genegff,self.trnagff,self.rrnagff,self.pro ,self.genomefa,self.anno)
        self.logger.info(cmd)
        self.logger.info("开始运行run_gbk_complete")
        command = self.add_command("run_gbk_complete", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_gbk_complete完成")
        else:
            self.set_error("运行run_gbk_complete运行出错!", code="31400802")

    def run_cgview_un(self):
        cmd = "{} {} -sequence {} -output {} -genes {}".format(self.perl_path, self.xml_builder,self.work_dir + '/' + self.option('sample_name') + '.gbk',self.work_dir + '/cgview' + '.xml',self.work_dir + '/' + self.option('sample_name') + '.cgview.cog')
        self.logger.info(cmd)
        self.logger.info("开始运行run_xml_uncomplete")
        command = self.add_command("run_xml_uncomplete", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_xml_uncomplete完成")
        else:
            self.set_error("运行run_xml_uncomplete运行出错!", code="31400803")

    def run_cgview(self):
        for file in self.list:
            os.system('sed \'s/linear/circular/g\' %s -i ' % (self.work_dir + '/' + file +  '.gbk'))
            cmd = "{} {} -sequence {} -output {} -genes {}".format(self.perl_path, self.xml_builder,
                                                         self.work_dir + '/' + file + '.gbk',
                                                         self.work_dir + '/' + file + '.xml',self.work_dir + '/' + file + '.cgview.cog')
            self.logger.info(cmd)
            self.logger.info("开始运行run_xml_complete")
            path = file.lower() + ' run_xml_complete'
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s完成" %path)
            else:
                self.set_error("运行%s运行出错!" , variables=(path), code="31400804")

    def change_cog_desc(self,xml_file):   #
        tree = et.parse(xml_file)
        root = tree.getroot()
        legends = root.findall('legend')
        legend2 = legends[1]
        legend2.attrib['font'] = 'SansSerif, plain, 33'
        type_list = []
        for legend_item in legend2.findall('legendItem'):
            text = legend_item.get('text')
            m = re.match('(\w)\sCOG', text)
            if m:
                type = m.group(1)
                if type in self.desc_map:
                    type_list.append(type)
                    desc = self.desc_map[type]
                    legend_item.set('text','{} {}'.format(type,desc))

        ## 补可能缺少的 lable
        type_loc = {'A':0,'B':1,'U':20, 'V':21,'W':22,'Y':23,'Z':24}
        color_map = {'A':"rgb(255,225,255)",'B':"rgb(255,193,193)",
                     'U':"rgb(121,205,205)",'V':"rgb(0,191,255)","W":"rgb(67,110,238)",
                     "Y":"rgb(16,78,139)","Z":"rgb(0,0,238)"}

        for type in ['A','B','U','V','W','Y','Z']:
            if type not in type_list:
                l=et.Element('legendItem')
                l.set('drawSwatch','true')
                l.set('swatchColor', color_map[type])
                l.set('swatchOpacity','1.0')
                l.set('text', type + ' '+ self.desc_map[type])
                legend2.insert(type_loc[type],l)

        tree.write(xml_file)


    def run_modify_xml(self):
        if os.path.exists(self.work_dir + '/all.xml'):
            os.remove(self.work_dir + '/all.xml')
        os.link(self.option('xml_file').prop['path'],self.work_dir + '/all.xml')
        os.system('sed \'s/upper-left/middle-center/g\' %s -i'%(self.work_dir + '/all.xml'))

        os.system('sed \'s/Length://g\' %s -i' % (self.work_dir + '/all.xml'))
        if self.option('species_name'):
            os.system('sed \'s/Accession: unknown/%s/g\' %s -i' % (self.option('species_name'),self.work_dir + '/all.xml'))
        else:
            os.system('sed \'s/<legendItem text="Accession: unknown" \\/>//g\' %s -i' % (self.work_dir + '/all.xml'))

        self.change_cog_desc(self.work_dir+"/all.xml")

        if self.option('seq_type_circle') == "Circular":
            pass
        else:
            os.system("sed -i 's/isLinear=\"false\"/isLinear=\"true\"/g' {}".format(self.work_dir + '/all.xml'))
        for type in ['png','svg']:
            cmd = "{} -jar -Xmx1500m -d64 -Djava.awt.headless=true {} -i {} -o {} -f {}  -W 4500 ".format(self.java_path,self.cgview_jar,self.work_dir + '/all.xml', self.work_dir + '/cgview' +  '.' + type, type)
            self.logger.info(cmd)
            path ='gra_cgview_un' + type
            self.logger.info("开始运行run_gra_uncomplete")
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s 完成" %path)
            else:
                self.set_error("运行%s 运行出错!" , variables=(path), code="31400805")

    def run_gra_cgview_un(self):
        if self.option('seq_type_circle') == "Circular":
            pass
        else:
            os.system("sed -i 's/isLinear=\"false\"/isLinear=\"true\"/g' {}".format(self.work_dir + '/cgview' + '.xml'))
        for type in ['png','svg']:
            cmd = "{} -jar -Xmx1500m -d64 -Djava.awt.headless=true {} -i {} -o {} -f {} -W 3700 ".format(self.java_path,self.cgview_jar,self.work_dir + '/cgview' + '.xml', 'cgview.' + type, type)
            self.logger.info(cmd)
            path ='gra_cgview_un' + type
            self.logger.info("开始运行run_gra_uncomplete")
            command = self.add_command(path, cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行%s 完成" %path)
            else:
                self.set_error("运行%s 运行出错!" , variables=(path), code="31400806")

    def run_gra_cgview(self):
        for file in self.list:
            if self.option('seq_type_circle') == "Circular":
                pass
            else:
                os.system(
                    "sed -i 's/isLinear=\"false\"/isLinear=\"true\"/g' {}".format(self.work_dir + '/' + file + '.xml'))
            for type in ['png', 'svg']:
                cmd = "{} -jar -Xmx1500m -d64 -Djava.awt.headless=true {} -i {} -o {} -f {} -W 3700 ".format(self.java_path, self.cgview_jar,
                                                                      self.work_dir + '/' + file + '.xml', file + '.' + type, type)
                self.logger.info(cmd)
                path = file.lower() + '_cgview_un' + type
                self.logger.info("开始运行run_gra_complete")
                command = self.add_command(path, cmd)
                command.run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行%s 完成" % path)
                else:
                    self.set_error("运行%s 运行出错!" , variables=( path), code="31400807")

    def set_output(self):
        self.logger.info("set output")
        if self.option('xml_file').is_set:
            for type in ['png', 'svg']:
                if os.path.exists(self.output_dir + '/cgview' + '.' + type):
                    os.remove(self.output_dir + '/cgview' + '.' + type)
                os.link(self.work_dir + '/cgview' + '.' + type,self.output_dir + '/cgview' + '.' + type)
        else:
            if self.option('analysis') in ['uncomplete']:
                if not os.path.exists(self.output_dir + '/xml'):
                    os.mkdir(self.output_dir + '/xml')
                if os.path.exists(self.output_dir + '/xml/cgview' + '.' + 'xml'):
                    os.remove(self.output_dir + '/xml/cgview' + '.' + 'xml')
                os.link(self.work_dir + '/cgview' + '.' + 'xml',
                            self.output_dir + '/xml/cgview' + '.' + 'xml')
                for type in ['png', 'svg']:
                    if os.path.exists(self.output_dir + '/cgview' + '.' + type):
                        os.remove(self.output_dir + '/cgview' + '.' + type)
                    os.link(self.work_dir + '/cgview' + '.' + type,self.output_dir + '/cgview' + '.' + type)
            elif self.option('analysis') in ['complete']:
                for file in self.list:
                    path =self.output_dir + '/' + file
                    if not os.path.exists(path):
                        os.mkdir(path)
                    if not os.path.exists(path + '/xml'):
                        os.mkdir(path + '/xml')
                    if os.path.exists(self.output_dir + '/' + file + '/xml/cgview' + '.' + 'xml'):
                        os.remove(self.output_dir + '/' + file + '/xml/cgview' + '.' + 'xml')

                    os.link(self.work_dir + '/' + file  + '.' + 'xml',
                            self.output_dir + '/' + file + '/xml/cgview' + '.' + 'xml')
                    for type in ['png', 'svg']:
                        if os.path.exists(self.output_dir + '/' + file + '/cgview' + '.' + type):
                            os.remove(self.output_dir + '/' + file + '/cgview' + '.' + type)
                        else:
                            os.link(self.work_dir + '/' + file + '.' + type,
                                    self.output_dir + '/' + file + '/cgview' +'.' + type)

    def get_list(self):
        if self.option('analysis') in ['complete']:
            with open(self.genomefa, 'r') as f:
                lines = f.readlines()
                k = ''
                klen = 0
                for line in lines:
                    if re.search(r'^>', line):
                        line = line.replace('>', '')
                        tmp = line.rstrip('\n').split(' ')[0]
                        if klen > 1000 and klen < 20000000:             #guanqing.zou 20181017 完成图，长度小于1kb的序列不做cgview  20190213 : 增加肠道大于20M 的不做cgview
                             self.list.append(k)
                        else:
                            self.logger.info(k + ' 小于1000或者大于20M')
                        k = tmp
                        klen = 0
                    else:
                        line = line.strip()
                        klen += len(line)
                if klen > 1000 and klen < 20000000:
                    self.list.append(k)
                else:
                    self.logger.info(k + ' 小于1000或者大于20M')


    def run(self):
        """
        运行
        """
        super(CgviewTool, self).run()
        current_cog_color = self.work_dir + '/cog_color.list'  #20190523
        if os.path.exists(current_cog_color):
            os.remove(current_cog_color)
        os.link(self.cog_color,current_cog_color)

        if self.option('xml_file').is_set:
            self.run_modify_xml()
            self.set_output()
            self.end()
        else:
            if self.option('analysis') in ['uncomplete']:
                if os.path.getsize(self.genomefa) > 20000000:    #zouguanqing 20190213
                    self.logger.info('序列超过20M，不运行')
                    self.end()
                else:
                    self.run_gbk_un()
                    self.run_cgview_un()
                    self.run_gra_cgview_un()
                    self.set_output()
                    self.end()
            elif self.option('analysis') in ['complete']:
                self.get_list()  # 这步控制小于1k或大于20m的序列不做cgview
                self.run_gbk()
                self.run_cgview()
                self.run_gra_cgview()
                self.set_output()
                self.end()

