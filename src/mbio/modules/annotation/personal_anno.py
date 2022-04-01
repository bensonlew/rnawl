# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:2018.10.23

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import unittest


class PersonalAnnoModule(Module):
    '''
    宏基因组个性化数据库注释
    go, phi, mvirdb, tcdb, qs, pfam, sec, t3ss,probio, p450
    '''
    def __init__(self, work_id):
        super(PersonalAnnoModule, self).__init__(work_id)
        options = [
            {"name": "database", "type": "string"},  # 输入注释数据库，分号分隔
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table","required":True},  # 基因丰度表
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "blastout", "type": "infile", "format": "sequence.profile_table"},  # nr blast table表
            {"name": "nr_gene_anno", "type": "infile", "format": "sequence.profile_table"},  # nr gene anno file
            {"name": "nr_gene_anno_lca", "type": "infile", "format": "sequence.profile_table"},  # lca注释结果
            {"name": "nr_gene_anno_de", "type": "infile", "format": "sequence.profile_table"},  # deunclassified结果
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"}
        ]
        self.add_option(options)  ##### 检查option是否list格式，其中每个opt是否字典格式
        #self.split_fasta = self.add_tool("sequence.split_fasta")
        self.anno_tools = []
        self.name_tool_dict = {}
        #self.go_module = self.add_module("annotation.mg_go_anno")

    def check_options(self):
        database_list = ["go", "phi", "mvirdb", "tcdb", "qs", "pfam", "sec", "t3ss", "probio", "p450"]
        input_list = self.option("database").split(";")
        if not set(input_list).issubset(database_list):
            raise OptionError("database must be in [go, phi, mvirdb, tcdb, qs, pfam, sec, t3ss, probio, p450]", code="21202301")
        for each in input_list:
            use_fasta = ["phi", "mvirdb", "tcdb", "qs", "pfam", "sec", "t3ss", "p450"]
            if each in use_fasta and not self.option("query").is_set:
                raise OptionError("annotation must input fasta file!", code="21202302")
            if each == "go" and not self.option("blastout").is_set:
                raise OptionError("go annotation must input blastout file!", code="21202303")
            if each == "probio" or each == "ttss":
                if not self.option("nr_gene_anno").is_set and not self.option("nr_gene_anno_lca") and not self.option("nr_gene_anno_de"):
                    raise OptionError("probio annotation must input nr_gene_anno file!", code="21202304")
        return True

    def add_anno_module(self):
        databases = self.option("database").split(";")
        if "go" in databases:
            self.go_anno()
        if "qs" in databases:
            self.qs_anno()
        if "probio" in databases:
            self.probio_anno()
        if "phi" in databases:
            self.phi_anno()
        if "sec" in databases:
            self.sec_anno()
        if "t3ss" in databases:
            self.t3ss_anno()
        if "p450" in databases:
            self.p450_anno()
        if "pfam" in databases:
            self.pfam_anno()
        if "tcdb" in databases:
            self.tcdb_anno()
        if "mvirdb" in databases:
            self.mvirdb_anno()

    def go_anno(self):
        self.logger.info("add go annotation")
        self.go_module = self.add_module("annotation.mg_go_anno")
        self.anno_tools.append(self.go_module)
        #self.name_tool_dict[self.go_module] = "go"
        self.go_module.set_options({
            "blastout": self.nr_blastout,
            "reads_profile_table": self.reads_profile
        })

    def qs_anno(self):
        self.logger.info("add qs annotation")
        self.qs_module = self.add_module("annotation.qs_annotation")
        self.anno_tools.append(self.qs_module)
        #self.name_tool_dict[self.qs_module] = "qs"
        self.qs_module.set_options({
            "query": self.query,
            "lines": 100000,
            "reads_profile_table": self.reads_profile,
            "group_table": self.option('group_table')
        })

    def probio_anno(self):  ############### 考虑内部拆分
        self.logger.info("add probio annotation")
        #self.name_tool_dict["probio"] = self.probio_module
        if self.option("nr_gene_anno").is_set:
            self.probio_module = self.add_module("annotation.probio_anno")
            self.anno_tools.append(self.probio_module)
            self.probio_module.set_options({
                "nr_gene_anno": self.option("nr_gene_anno"),
                "reads_profile_table": self.reads_profile
            })
        if self.option("nr_gene_anno_lca").is_set:
            self.probio_module_lca = self.add_module("annotation.probio_anno")
            self.anno_tools.append(self.probio_module_lca)
            #self.name_tool_dict["probio"] = self.probio_module
            self.probio_module_lca.set_options({
                "nr_gene_anno": self.option("nr_gene_anno_lca"),
                "reads_profile_table": self.reads_profile
            })
        if self.option("nr_gene_anno_de").is_set:
            self.probio_module_de = self.add_module("annotation.probio_anno")
            self.anno_tools.append(self.probio_module_de)
            #self.name_tool_dict["probio"] = self.probio_module
            self.probio_module_de.set_options({
                "nr_gene_anno": self.option("nr_gene_anno_de"),
                "reads_profile_table": self.reads_profile
            })


    def p450_anno(self):
        self.logger.info("add p450 annotation")
        self.p450_module = self.add_module("annotation.cyps_annotation")
        self.anno_tools.append(self.p450_module)
        #self.name_tool_dict["p450"] = p450_module
        self.p450_module.set_options({
            "query": self.query,
            "lines": 100000,
            "reads_profile_table": self.reads_profile
        })

    def pfam_anno(self):
        self.logger.info("add pfam annotation")
        self.pfam_module = self.add_module("annotation.pfam_annotation")
        self.anno_tools.append(self.pfam_module)
        #self.name_tool_dict["pfam"] = pfam_module
        self.pfam_module.set_options({
            "query": self.query,
            #"lines": 200000,
            "lines": 20000, ### 测试用参数值
            "reads_profile_table": self.reads_profile,
            "database":"pfam",
        })

    def tcdb_anno(self):
        self.logger.info("add tcdb annotation")
        self.tcdb_module = self.add_module("annotation.tcdb_annotation")
        self.anno_tools.append(self.tcdb_module)
        #self.name_tool_dict["tcdb"] = tcdb_module
        self.tcdb_module.set_options({
            "query": self.query,
            "lines": 100000,
            "reads_profile_table": self.reads_profile,
        })

    def phi_anno(self):
        self.logger.info("add phi annotation")
        self.phi_module = self.add_module("annotation.phi_annotation")
        self.anno_tools.append(self.phi_module)
        #self.name_tool_dict["phi"] = phi_module
        self.phi_module.set_options({
            "query": self.query,
            "lines": 100000,
            "reads_profile_table": self.reads_profile,
        })

    def sec_anno(self):
        self.logger.info("add sec annotation")
        self.sec_module = self.add_module("annotation.sec_anno")
        self.anno_tools.append(self.sec_module)
        #self.name_tool_dict["sec"] = sec_module
        self.sec_module.set_options({
            "query": self.query,
            "lines": 100000,
            "reads_profile_table": self.reads_profile,
            "signalp_out_format": "short",
            "signalp_type": "gram-,gram+,euk",
            "anno_tool": "signalp",
            "tax_table": self.nr_gene_anno
        })

    def t3ss_anno(self):
        self.logger.info("add t3ss annotation")
        self.t3ss_module = self.add_module("annotation.sec_anno")
        self.anno_tools.append(self.t3ss_module)
        tax_table_list = []
        options ={
            "query": self.query,
            "reads_profile_table": self.reads_profile,
            "anno_tool": "EffectiveT3",
            #"tax_table_list": self.nr_gene_anno.prop["path"]
        }
        if self.nr_gene_anno:
            tax_table_list.append(self.nr_gene_anno.prop["path"])
        if self.nr_gene_anno_lca:
            tax_table_list.append(self.nr_gene_anno_lca.prop["path"])
        if self.nr_gene_anno_de:
            tax_table_list.append(self.nr_gene_anno_de.prop["path"])
        options["tax_table_list"] = ",".join(tax_table_list)
        self.t3ss_module.set_options(options)


    def mvirdb_anno(self):
        self.logger.info("add mvirdb annotation")
        self.mvirdb_module = self.add_module("annotation.mvirdb_annotation")
        self.anno_tools.append(self.mvirdb_module)
        #self.name_tool_dict["mvirdb"] = mvirdb_module
        self.mvirdb_module.set_options({
            "query": self.query,
            "lines": 100000,
            "reads_profile_table": self.reads_profile,
        })

    def set_output(self):
        '''
        '''
        for eachtool in self.anno_tools:
            self.link(eachtool)
        self.end()

    def link(self, eachtool):
        """
        link文件到本module的output目录
        """
        if hasattr(self,"go_module") and eachtool == self.go_module:
            out_dir = "Go/"
        elif hasattr(self,"probio_module") and eachtool == self.probio_module:
            out_dir = "Probio/"
        elif hasattr(self,"probio_module_lca") and eachtool == self.probio_module_lca:
            out_dir = "Probio_LCA/"
        elif hasattr(self,"probio_module_de") and eachtool == self.probio_module_de:
            out_dir = "Probio_Deunclassifed/"
        elif hasattr(self,"qs_module") and eachtool == self.qs_module:
            out_dir = "Qs/"
        elif hasattr(self,"p450_module") and eachtool == self.p450_module:
            out_dir = "P450/"
        elif hasattr(self,"pfam_module") and eachtool == self.pfam_module:
            out_dir = "Pfam/"
        elif hasattr(self,"tcdb_module") and eachtool == self.tcdb_module:
            out_dir = "Tcdb/"
        elif hasattr(self,"phi_module") and eachtool == self.phi_module:
            out_dir = "Phi/"
        elif hasattr(self,"sec_module") and eachtool == self.sec_module:
            out_dir = "Sec/"
        elif hasattr(self,"t3ss_module") and eachtool == self.t3ss_module:
            out_dir = "Ttss/"
        elif hasattr(self,"mvirdb_module") and eachtool == self.mvirdb_module:
            out_dir = "Mvirdb/"
        path = os.path.join(self.output_dir, out_dir)
        if not os.path.exists(path):
            os.mkdir(path)
        allfiles = os.listdir(eachtool.output_dir)
        oldfiles = [os.path.join(eachtool.output_dir, i) for i in allfiles]
        newfiles = [os.path.join(path, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        self.nr_gene_anno = ""
        self.nr_gene_anno_lca = ""
        self.nr_gene_anno_de = ""
        super(PersonalAnnoModule, self).run()
        self.logger.info("start personal annotation")
        self.reads_profile = self.option("reads_profile_table")
        if self.option("query").is_set:
            self.query = self.option("query")
        if self.option("blastout").is_set:
            self.nr_blastout = self.option("blastout")
        if self.option("nr_gene_anno").is_set:
            self.nr_gene_anno = self.option("nr_gene_anno")
        if self.option("nr_gene_anno_lca").is_set:
            self.nr_gene_anno_lca = self.option("nr_gene_anno_lca")
        if self.option("nr_gene_anno_de").is_set:
            self.nr_gene_anno_de = self.option("nr_gene_anno_de")
        self.add_anno_module()
        '''
        need_link = []
        for eachtool in self.anno_tools:
            mydict = {"CypsAnnotation":"P450", "MvirdbAnnotation": "Mvirdb","PhiAnnotation":"phi","ProbioAnno":"Probio",
                      "QsAnnotation":"Qs","TcdbAnnotation":"Tcdb"}
            out_dir = mydict[eachtool._name] + "/"
            self.logger.info("test>>>>>>>>>>>>>>>>>>>>>>")
            self.logger.info(out_dir)
            need_link.append((eachtool,"*",out_dir))
            path = os.path.join(self.output_dir, out_dir)
            self.logger.info("--------")
            self.logger.info(path)
            self.logger.info(os.path.exists(path))
            if not os.path.exists(path):
                os.mkdir(path)
            self.add_result_mapping(need_link)
        #self.add_result_mapping([(self.qs_module,"*xls","qs")])
        '''
        self.on_rely(self.anno_tools, self.set_output)
        for eachmodule in self.anno_tools:
            eachmodule.run()
            self.logger.info(eachmodule._name)

    def end(self):
        self.logger.info(">>>>>>>>>>>>>>>>>>end<<<<<<<<<<<<")
        super(PersonalAnnoModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "personal_anno" + str(random.randint(1, 10000)),
            #"id": "personal_anno",
            "type": "module",
            "name": "annotation.personal_anno",
            "instant": True,
            "options": dict(
                #query="/mnt/ilustre/tsanger-data/rerewrweset/files/m_195/195_5a4f29cb4cf76/tsanger_27029/workflow_results/geneset/uniGeneset/gene.uniGeneset.faa",
                query="/mnt/ilustre/users/sanger-dev/workspace/20181220/MetaGenomic_tsg_33075/MetaDiamond/SplitFasta/output/fasta_4",
                reads_profile_table="/mnt/ilustre/users/sanger-dev/workspace/20181220/MetaGenomic_tsg_33075/output/geneset/gene_profile/reads_number.xls",
                #reads_profile_table="/mnt/ilustre/tsanger-data/rerewrweset/files/m_195/195_5a4f29cb4cf76/tsanger_27029/workflow_results/geneset/gene_profile/reads_number.xls",
                #nr_gene_anno="/mnt/ilustre/tsanger-data/rerewrweset/files/m_195/195_5a4f29cb4cf76/tsanger_27029/workflow_results/nr/gene_nr_anno.xls",
                #blastout="/mnt/ilustre/users/sanger-dev/workspace/20181113/MetaGenomic_tsg_32885/upload_results/nr/test_nr_align_table.xls",
                database="qs",
                #nr_gene_anno_lca="/mnt/ilustre/tsanger-data/rerewrweset/files/m_195/195_5a4f29cb4cf76/tsanger_27029/workflow_results/nr_lca/gene_nr_anno.xls",
                #nr_gene_anno_de="/mnt/ilustre/tsanger-data/rerewrweset/files/m_195/195_5a4f29cb4cf76/tsanger_27029/workflow_results/nr_deunclassify/gene_nr_anno.xls",
                #database="go;qs;probio;p450;phi;tcdb;mvirdb;pfam;t3ss;sec",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
