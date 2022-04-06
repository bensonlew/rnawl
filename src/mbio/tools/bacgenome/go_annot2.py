# -*- coding: utf-8 -*-
# __author__ = 'qignchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import re
import subprocess
import pandas as pd
import math
import shutil
import xml.etree.ElementTree as ET
from bson.objectid import ObjectId


class GoAnnotAgent(Agent):
    """
    改换go注释的方法
    """
    def __init__(self, parent):
        super(GoAnnotAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml"},
            {"name": "go2level_out", "type": "outfile", "format": "annotation.go.level2"},
            {"name": "golist_out", "type": "outfile", "format": "annotation.go.go_list"},
            {"name": "blast2go_annot", "type": "outfile", "format": "annotation.go.blast2go_annot"},
        ]
        self.add_option(options)
        self.step.add_steps('go_annotation')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = 'BLAST2GO'  # 投递到指定的队列BLAST2GO

    def step_start(self):
        self.step.go_annotation.start()
        self.step.update()

    def step_end(self):
        self.step.go_annotation.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = '60G'

    def end(self):
        super(GoAnnotAgent, self).end()


class GoAnnotTool(Tool):

    def __init__(self, config):
        super(GoAnnotTool, self).__init__(config)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0")
        self.set_environ(JRE_HOME=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/lib")
        self.set_environ(CLASSPATH=self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/jre/lib")
        self.annot_go = self.config.PACKAGE_DIR + "/prok_rna/goAnnot.py"
        self.mongodb_nr = Config().get_mongo_client(mtype="ref_rna", ref=True)[Config().get_mongo_dbname("ref_rna", ref=True)]['NR_sequence_20200604']
        self._version = "1.0"
        self.b2g_user = "biocluster102"
        self.b2g_password = "sanger-dev-123"

    def run(self):
        super(GoAnnotTool, self).run()
        self.run_b2g()

    def get_nrxml_gi_description(self, xml_path):
        """
        目的是将accession号换成gi号
        """
        xml = ET.parse(xml_path)
        root = xml.getroot()
        self.logger.info("开始改注释信息内容!")
        BlastOutput_iterations = root.find('BlastOutput_iterations')
        for one_query in BlastOutput_iterations.findall('Iteration'):
            iteration_hits = one_query.findall('Iteration_hits')
            for iteration_hit in iteration_hits:
                hits = iteration_hit.findall('Hit')
                for hit in hits:
                    seq_id = hit.find('Hit_id')
                    seq_id2 = seq_id.text
                    result = self.mongodb_nr.find_one({"origin_sequence_id": seq_id2})
                    if result:
                        hit_id = result["_id"]
                        self.logger.info(hit_id)
                        if not isinstance(hit_id, str):
                            if not isinstance(hit_id, ObjectId):
                                seq_id.text = hit_id.strip()
                        else:
                            hits.remove(hit)
                    else:
                        hits.remove(hit)
                        print "没找到gid:{}".format(seq_id2)
        xml.write('tmp.txt')
        with open('tmp.txt', 'rb') as f, open('tmp.xml', 'wb') as w:
            lines = f.readlines()
            a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
            w.write(a)
            w.writelines(lines)
        os.remove('tmp.txt')
        shutil.move('tmp.xml', xml_path)

    def run_b2g(self):
        self.blast_nr_out = self.work_dir + '/temp_blast_nr.xml'
        self.option("blastout").change_blast_version(self.blast_nr_out)
        self.get_nrxml_gi_description(self.blast_nr_out)
        cmd = '/program/sun_jdk1.8.0/bin/java -Xmx30g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
        cmd += self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/ext/*: es.blast2go.prog.B2GAnnotPipe'
        cmd += ' -in {} -prop {}/bioinfo/annotation/b2g4pipe_v2.5/b2gPipe.properties -annot -out {}'.format(self.blast_nr_out, self.config.SOFTWARE_DIR, self.work_dir + '/blast2go')
        self.logger.info('运行b2g程序')
        self.logger.info(cmd)
        b2g = self.add_command('b2g', cmd)
        b2g.run()
        self.wait('b2g')
        if b2g.return_code == 0:
            self.logger.info('运行b2g完成')
            linkfile = self.output_dir + '/blast2go.annot'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(self.work_dir + '/blast2go.annot', linkfile)
            self.option('blast2go_annot', linkfile)
            self.logger.debug("b2g end")
            if os.path.exists(self.work_dir + '/blast2go.annot'):
                if os.path.getsize(self.work_dir + '/blast2go.annot') != 0:
                    self.run_gomerge()
                else:
                    self.end()
            else:
                self.end()
        else:
            self.set_error('运行b2g出错', code="31201501")

    def run_gomerge(self):
        cmd1 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goMerge.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd1 += ' %s %s' % (
            self.work_dir + '/blast2go.annot', 'GO.list')
        # Config.DB_HOST,Config.DB_USER,Config.DB_PASSWD
        self.logger.info("运行mergeGO.py")
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            if os.path.exists(self.output_dir + '/GO.list'):
                os.remove(self.output_dir + '/GO.list')
            if os.path.exists(self.output_dir + '/query_gos.list'):
                os.remove(self.output_dir + '/query_gos.list')
            os.link(self.work_dir + '/GO.list',
                    self.output_dir + '/query_gos.list')
            self.option('golist_out', self.output_dir + '/query_gos.list')
        except subprocess.CalledProcessError:
            self.set_error('运行mergeGO.py出错', code="31201502")
        self.run_annotation()

    def run_annotation(self):
        cmd2 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goAnnot.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd2 += ' %s %s %s %s' % (
            self.work_dir + '/GO.list', 'localhost', self.b2g_user, self.b2g_password)  # 10.100.203.193
        self.logger.info("运行goAnnot.py")
        self.logger.info(cmd2)
        try:
            subprocess.check_output(cmd2, shell=True)
            self.logger.info("运行goAnnot.py完成")
            outfiles = ['go1234level_statistics.xls', 'go123level_statistics.xls', 'go12level_statistics.xls']
            for item in outfiles:
                linkfile = self.output_dir + '/' + item
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(self.work_dir + '/' + item, linkfile)
        except subprocess.CalledProcessError:
            self.set_error("运行goAnnot.py出错", code="31201503")
        # self.run_gosplit()
        self.run_not_level()

    def run_gosplit(self):
        cmd3 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goSplit.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        cmd3 += ' %s' % self.work_dir + '/go_detail.xls'
        self.logger.info("运行goSplit.py")
        self.logger.info(cmd3)
        try:
            subprocess.check_output(cmd3, shell=True)
            self.logger.info("运行goSplit.py完成")
            outfiles = ['go2level.xls', 'go3level.xls', 'go4level.xls']
            for item in outfiles:
                linkfile = self.output_dir + '/' + item
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(self.work_dir + '/' + item, linkfile)
            self.option('go2level_out', self.output_dir + '/go2level.xls')
        except subprocess.CalledProcessError:
            self.set_error("运行goSplit.py出错", code="31201504")
        self.end()

    def run_not_level(self):
        cmd = '{}/miniconda2/bin/python {}/annotation/go/go_desc.py'.format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
        cmd += ' {} {}'.format(self.work_dir + '/blast2go.annot', self.work_dir + '/go_statistics.xls' )
        self.logger.info('运行go_desc.py')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd,shell=True)
            self.logger.info('运行go_desc.py 完成')
            if os.path.exists(self.output_dir + '/go_statistics.xls'):
                os.remove(self.output_dir + '/go_statistics.xls')
            os.link(self.work_dir + '/go_statistics.xls', self.output_dir + '/go_statistics.xls')
        except:
            self.set_error("运行go_desc出错")
        self.end()



    # def run_b2g(self):
    #     self.blast_nr_out = self.work_dir + '/temp_blast_nr.xml'
    #     self.option("blastout").change_blast_version(self.blast_nr_out)
    #     cmd = '/program/sun_jdk1.8.0/bin/java -Xmx30g -cp ' + self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/*:'
    #     cmd += self.config.SOFTWARE_DIR + '/bioinfo/annotation/b2g4pipe_v2.5/ext/*: es.blast2go.prog.B2GAnnotPipe'
    #     cmd += ' -in {} -prop {}/bioinfo/annotation/b2g4pipe_v2.5/b2gPipe.properties -annot -out {}'.format(self.blast_nr_out, self.config.SOFTWARE_DIR, self.work_dir + '/blast2go')
    #     self.logger.info('运行b2g程序')
    #     self.logger.info(cmd)
    #     b2g = self.add_command('b2g', cmd)
    #     b2g.run()
    #     self.wait('b2g')
    #     if b2g.return_code == 0:
    #         self.logger.info('运行b2g完成')
    #         linkfile = self.output_dir + '/blast2go.annot'
    #         if os.path.exists(linkfile):
    #             os.remove(linkfile)
    #         os.link(self.work_dir + '/blast2go.annot', linkfile)
    #         self.option('blast2go_annot', linkfile)
    #         self.logger.debug("b2g end")
    #         if os.path.exists(self.work_dir + '/blast2go.annot'):
    #             if os.path.getsize(self.work_dir + '/blast2go.annot') != 0:
    #                 self.run_gomerge()
    #             else:
    #                 self.end()
    #         else:
    #             self.end()
    #     else:
    #         self.set_error('运行b2g出错', code="31201501")
    #
    # def run_gomerge(self):
    #     cmd1 = 'miniconda2/bin/python {}/annotation/goMerge.py'.format(
    #         self.config.PACKAGE_DIR)
    #     cmd1 += ' %s %s' % (self.option("blast2go_annot").prop['path'], 'GO.list')
    #     self.logger.info("运行mergeGO.py")
    #     self.logger.info(cmd1)
    #     command = self.add_command("run_gomerge", cmd1)
    #     command.run()
    #     self.wait()
    #     if command.return_code == 0:
    #         self.logger.info("运行mergeGO.py")
    #     else:
    #         self.set_error('running mergeGO.py error', code="31205501")
    #
    # # def run_annotation(self):
    # #     cmd2 = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/goAnnot.py'.format(
    # #         self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
    # #     cmd2 += ' %s %s %s %s' % (
    # #         self.work_dir + '/GO.list', 'localhost', self.b2g_user, self.b2g_password)  # 10.100.203.193
    # #     self.logger.info("运行goAnnot.py")
    # #     self.logger.info(cmd2)
    # #     try:
    # #         subprocess.check_output(cmd2, shell=True)
    # #         self.logger.info("运行goAnnot.py完成")
    # #         outfiles = ['go1234level_statistics.xls', 'go123level_statistics.xls', 'go12level_statistics.xls']
    # #         for item in outfiles:
    # #             linkfile = self.output_dir + '/' + item
    # #             if os.path.exists(linkfile):
    # #                 os.remove(linkfile)
    # #             os.link(self.work_dir + '/' + item, linkfile)
    # #     except subprocess.CalledProcessError:
    # #         self.set_error("运行goAnnot.py出错", code="31201503")
    #
    # # def run_gomerge(self):
    # #     cmd1 = 'miniconda2/bin/python {}/bioinfo/annotation/scripts/goMerge.py'.format(
    # #         self.config.SOFTWARE_DIR)
    # #     cmd1 += ' %s %s' % (self.option("blast2go_annot").prop['path'], 'GO.list')
    # #     self.logger.info("运行mergeGO.py")
    # #     self.logger.info(cmd1)
    # #     #try:
    # #         #subprocess.check_output(cmd1, shell=True)
    # #     #except subprocess.CalledProcessError:
    # #         #self.set_error('running mergeGO.py error', code="31205501")
    # #     command = self.add_command("run_gomerge", cmd1)
    # #     command.run()
    # #     self.wait()
    # #     if command.return_code == 0:
    # #         self.logger.info("运行mergeGO.py")
    # #     else:
    # #         self.set_error('running mergeGO.py error', code="31205501")
    #
    # def run_annotation(self):
    #     cmd2 = 'miniconda2/bin/python {} {}'.format(self.annot_go, self.work_dir + '/GO.list')
    #     self.logger.info("运行goAnnot.py")
    #     self.logger.info(cmd2)
    #     command = self.add_command("run_annotation", cmd2, ignore_error=True)
    #     command.run()
    #     self.wait()
    #     if command.return_code == 0:
    #         self.logger.info("运行run_annotation结束")
    #     elif command.return_code in [-9, 1]:  # change return_code by qingchen.zhang @ 20190610
    #         self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @ 20190610
    #     else:
    #         self.set_error("运行run_annotation出错")
    #
    # def run_go_gene_anno(self):
    #     cmd3 = 'miniconda2/bin/python {}/annotation/go/go_anno.py'.format(
    #         self.config.PACKAGE_DIR)
    #     cmd3 += ' -i %s -o %s' % (self.work_dir + '/go1234level_statistics.xls',self.work_dir + '/all.go.annotation.xls')
    #     self.logger.info("运行go_anno.py")
    #     command = self.add_command("run_go_gene_anno", cmd3, ignore_error=True)
    #     command.run()
    #     self.wait()
    #     if command.return_code == 0:
    #         self.logger.info("运行run_go_gene_anno结束")
    #     elif command.return_code in [-9, 1]:  # change return_code by qingchen.zhang @ 20190610
    #         self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @ 20190610
    #     else:
    #         self.set_error("运行run_go_gene_anno出错")
    #
    # def set_output(self):
    #     for i in ["go1234level_statistics.xls", 'go123level_statistics.xls', 'go12level_statistics.xls', 'GO.list', "go_statistics.xls"]:
    #         if os.path.exists(self.output_dir + '/' + i):
    #             os.remove(self.output_dir + '/' + i)
    #         os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)
    #
    # def format_go(self):
    #     """
    #     根据输入文件整理一下
    #     :return:
    #     """
    #     go_file = self.work_dir  + "/go_statistics2.xls"
    #     with open(self.option("blast2go_annot").prop['path'], 'r') as f, open(go_file, 'w') as w:
    #         for line in f:
    #             spline=line.strip().split("\t")
    #             gene = spline[0]
    #             go_name = list(set(spline[1].split(";")))
    #             go_des = spline[3]
    #             go_name.sort()
    #             for go in go_name:
    #                 w.write("{}\t{}\t{}\n".format(gene, go, go_des))
    #
    # def run_not_level(self):
    #     cmd = '{}/miniconda2/bin/python {}/annotation/go/go_desc.py'.format(self.config.SOFTWARE_DIR, self.config.PACKAGE_DIR)
    #     cmd += ' {} {}'.format(self.work_dir  + "/go_statistics2.xls", self.work_dir + '/go_statistics.xls' )
    #     self.logger.info('运行go_desc.py')
    #     self.logger.info(cmd)
    #     try:
    #         subprocess.check_output(cmd,shell=True)
    #         self.logger.info('运行go_desc.py 完成')
    #     except:
    #         self.set_error("运行go_desc出错")

