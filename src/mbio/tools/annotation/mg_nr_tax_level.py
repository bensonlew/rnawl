# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last modify by shaohua.yuan
# last modify date: 20170913

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from mbio.packages.annotation.nr_stat import nr_stat
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon
import os
from pymongo import MongoClient
from mbio.packages.align.blast.xml2table import xml2table
from biocluster.config import Config

class MgNrTaxLevelAgent(Agent):
    """
    nr注释的level统计
    """

    def __init__(self, parent):
        super(MgNrTaxLevelAgent, self).__init__(parent)
        options = [
            {"name": "nr_taxon_profile_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "nr_gene_anno_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "nr_align_dir", "type": "infile", "format": "annotation.mg_anno_dir"}
        ]
        self.add_option(options)
        self._memory_increase_step = 50  # 每次重运行增加内存70G by qingchen.zhang @ 20200410

    def check_options(self):
        if not self.option("nr_taxon_profile_dir").is_set:
            raise OptionError("必须设置输入taxon丰度文件夹", code="31203201")
        #if not self.option("nr_taxon_anno_dir").is_set:
            #raise OptionError("必须设置输入taxon注释文件夹", code="31203202")
        if not self.option("nr_gene_anno_dir").is_set:
            raise OptionError("必须设置输入nr注释结果table文件夹", code="31203203")
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '15G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
        ])
        super(MgNrTaxLevelAgent, self).end()


class MgNrTaxLevelTool(Tool):
    def __init__(self, config):
        super(MgNrTaxLevelTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "miniconda2/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/taxon/scripts/mg_nr_taxlevel.py'
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/cat_seq.sh'
        self.result_name = ''
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        #self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.gi_tax = self.mongodb.NR_sequence

    def run(self):
        """
        运行
        :return:
        """
        super(MgNrTaxLevelTool, self).run()
        #self.merge_anno_table()
        #self.merge_align_table()
        self.merge_gene_table()
        self.merge_profile_table()
        self.merge_nr_align_table()
        self.tax_level()
        os.system("sed -i '1i#Query\tTaxid\tDomain\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tAlign_len\tIdentity(%)'" + " %s" %(self.work_dir +"/gene_nr_anno.xls"))
        self.logger.info("完成基因注释表的表头插入")
        os.system("sed -i '1iScore\tE-Value\tHSP-Len\tIdentity-%\tSimilarity-%\tQuery-Name\tQ-Len\tQ-Begin\tQ-End\tQ-Frame\tHit-Name\tHit-Len\tHsp-Begin\tHsp-End\tHsp-Frame\tHit-Description'" + " %s" %(self.work_dir + "/tmp_nr_align.xls"))
        self.logger.info("完成比对结果表的表头插入")
        self.set_output()
        self.end()

    def merge_anno_table(self):
        nr_anno = 0
        anno_file = os.listdir(self.option('nr_taxon_anno_dir').prop['path'])
        self.anno_name = os.path.join(self.work_dir, "tmp_taxons_anno.xls")
        if os.path.exists(self.anno_name):
            os.remove(self.anno_name)
        for i in anno_file:
            nr_anno += 1
            file_path = os.path.join(self.option('nr_taxon_anno_dir').prop['path'], i)
            cmd = '{} {} {}'.format(self.sh_path, file_path, self.anno_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat anno" + str(nr_anno)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31203201")
                raise Exception("cat {} error".format(i))

        if self.option('out_type') == 1:  ##zouguanqing 20190314
            nr_anno = 0
            tmp_ncbi_result_deal = '/'.join(self.option('nr_taxon_anno_dir').prop['path'].split('/')[:-1]) + '/tmp_ncbi_result_deal'
            anno_file_new = os.listdir(tmp_ncbi_result_deal)
            self.anno_name_new = os.path.join(self.work_dir, "query_taxons.xls")
            if os.path.exists(self.anno_name_new):
                os.remove(self.anno_name_new)
            for i in anno_file_new:
                nr_anno += 1
                file_path = os.path.join(tmp_ncbi_result_deal, i)
                cmd = '{} {} {}'.format(self.sh_path, file_path, self.anno_name_new)
                self.logger.info("start cat {}".format(i))
                command_name = "cat anno new " + str(nr_anno)
                command = self.add_command(command_name, cmd).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("cat {} done".format(i))
                else:
                    self.set_error("cat %s error", variables=(i), code="31203201")
                    raise Exception("cat {} error".format(i))

        else:   ####### zouguanqing
            nr = mg_taxon()
            self.logger.info("start nr_stat(detail_to_level)")
            nr.detail_to_level(self.anno_name, self.work_dir)

        gene_tax = {}
        with open(self.anno_name, "r") as f1:
            for line in f1:
                line = line.strip().split("\t")
                gene = line[0]
                genes = gene.split("_")
                newgene = "_".join(genes[0:len(genes) - 1])
                gi = int(line[1])
                detail = self.gi_tax.find_one({"_id": gi})
                if detail:
                    taxid = detail["taxid"]
                    gene_tax[newgene] = taxid
        #nr = nr_stat()


        with open(self.work_dir + "/query_taxons.xls", "r") as f, open(self.work_dir + "/tmp_gene_nr_anno.xls", \
                                                                       "w") as outfile:
            outfile.write("#Query\tTaxid\tDomain\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
            for line in f:
                line = line.strip().split("\t")
                gene = line[0]
                genes = gene.split("_")
                newgene = "_".join(genes[0:len(genes) - 1])
                #print "query_taxons:",newgene
                tax = "\t".join(line[1].replace(" ","_").split(";"))
                if gene_tax.has_key(newgene):
                    #print "map:",newgene
                    outfile.write(newgene + "\t" + str(gene_tax[newgene]) + "\t" + tax + "\n")

    def merge_align_table(self):
        nr_align = 0
        align_file = os.listdir(self.option('nr_align_dir').prop['path'])
        align_name = os.path.join(self.work_dir, "tmp_nr_align.xls")
        if os.path.exists(align_name):
            os.remove(align_name)
        for i in align_file:
            nr_align += 1
            file_path = os.path.join(self.option('nr_align_dir').prop['path'], i)
            cmd = '{} {} {}'.format(self.sh_path, file_path, align_name)
            self.logger.info("start cat {}".format(i))
            command_name = "cat align" + str(nr_align)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("cat {} done".format(i))
            else:
                self.set_error("cat %s error", variables=(i), code="31203202")
        with open(self.work_dir + "/tmp_nr_align.xls", "r") as f, open(self.output_dir + "/nr_align_table.xls","w") as outf :
            data_length = {}
            data_identity = {}
            head = f.next().strip()
            outf.write( head + "\n")
            for line in f:
                line = line.strip()
                line1= line.split("\t")
                if line1[0] != "Score":
                    identity = line1[3]
                    length = line1[2]
                    gene = line1[5]
                    new = gene.split("_")
                    newgene = "_".join(new[0:len(new) - 1])
                    data_length[newgene] = length
                    data_identity[newgene] = identity
                    outf.write( line + "\n" )
        with open(self.work_dir + "/tmp_gene_nr_anno.xls","r") as f2, open(self.output_dir + "/gene_nr_anno.xls", \
                "w") as outfile:
            for line in f2:
                line = line.strip()
                line1 = line.split("\t")
                if line1[0] == "#Query":
                    head = line
                    outfile.write(head + "\tIdentity(%)\tAlign_len\n")
                else:
                    gene = line1[0]
                    outfile.write(line + "\t" + data_identity[gene] + "\t" + data_length[gene] + "\n" )

    def merge_gene_table(self):
        """
        合并所有注释结果
        :return:
        """
        anno_file = os.listdir(self.option('nr_gene_anno_dir').prop['path'])
        anno_name = os.path.join(self.work_dir, "gene_nr_anno.xls")
        if os.path.exists(anno_name):
            os.remove(anno_name)
        cmd = '{}'.format(self.sh_path)
        for i in anno_file:
            file_path = os.path.join(self.option('nr_gene_anno_dir').prop['path'], i)
            cmd += ' ' + file_path
        cmd += ' ' + anno_name
        self.logger.info(cmd)
        command_name = "cat_gene"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("合并注释结果{}成功".format(command_name))
        else:
            self.set_error("合并注释结果失败", code="31203202")

    def merge_profile_table(self):
        """
        合并所有丰度结果
        :return:
        """
        profile_file = os.listdir(self.option('nr_taxon_profile_dir').prop['path'])
        self.result_name = os.path.join(self.work_dir, "tmp_taxons_profile.xls")
        if os.path.exists(self.result_name):
            os.remove(self.result_name)
        cmd = '{}'.format(self.sh_path)
        for i in profile_file:
            file_path = os.path.join(self.option('nr_taxon_profile_dir').prop['path'], i)
            cmd += ' ' + file_path
        cmd += ' ' + self.result_name
        self.logger.info(cmd)
        command_name = "cat_profile"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("合并丰度结果{}成功".format(command_name))
        else:
            self.set_error("合并丰度结果失败", code="31203203")

    def merge_nr_align_table(self):
        """
        合并所有比对结果
        :return:
        """
        align_file = os.listdir(self.option('nr_align_dir').prop['path'])
        align_name = os.path.join(self.work_dir, "tmp_nr_align.xls")
        if os.path.exists(align_name):
            os.remove(align_name)
        cmd = '{}'.format(self.sh_path)
        for i in align_file:
            file_path = os.path.join(self.option('nr_align_dir').prop['path'], i)
            cmd += ' ' + file_path
        cmd += ' ' + align_name
        self.logger.info(cmd)
        command_name = "cat_align"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            os.system("sed -i '/Score/d' %s" %(self.work_dir + "/tmp_nr_align.xls"))
            self.logger.info("合并注释结果{}成功".format(command_name))
        else:
            self.set_error("合并注释结果失败", code="31203205")

    def tax_level(self):
        self.logger.info("start nr_tax_level")
        cmd2 = self.python_path + ' {} -i {} -l 1,2,3,4,5,6,7,8 -o {}'. \
            format(self.python_script, self.result_name, self.output_dir)
        command2 = self.add_command('nr_tax_level', cmd2).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("nr_tax_level succeed")
        elif command2.return_code in [1]:  # 内存超出是返回值为1 add by qingchen.zhang @ 20200410
            self.logger.info("return code: %s" % command2.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by qingchen.zhang @ 20200410
        else:
            self.set_error("nr_tax_level failed", code="31203204")
            raise Exception("nr_tax_level failed")

    def set_output(self):
        self.logger.info("start set_output")
        if os.path.exists(self.output_dir + "/gene_nr_anno.xls"):
            os.remove(self.output_dir + "/gene_nr_anno.xls")
        os.link(self.work_dir + "/gene_nr_anno.xls", self.output_dir + "/gene_nr_anno.xls")
        if os.path.exists(self.output_dir + "/nr_align_table.xls"):
            os.remove(self.output_dir + "/nr_align_table.xls")
        os.link(self.work_dir + "/tmp_nr_align.xls", self.output_dir + "/nr_align_table.xls")
        self.logger.info("Have got correct results")
