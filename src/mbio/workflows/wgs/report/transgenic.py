# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180516

import re
import os
import gevent
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class TransgenicWorkflow(Workflow):
    """
    交互分析：转基因的接口
    lasted modified by hongdong @ 20180516
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TransgenicWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "samples", "type": "string"},  # 页面传进来的所有的样本id，样本为逗号分隔
            {"name": "insert_seq", "type": "string"},  # 转基因序列
            {"name": "genome_version_id", "type": "string"},  # 查找基因组信息
            {"name": "clean_path", "type": 'string'},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.json_path = self.config.SOFTWARE_DIR + "/database/dna_wgs_geneome/"  # 参考组配置文件
        self.ref_index = self.add_module("wgs.ref_index")
        self.mapping = self.add_module("wgs.transgene_mapping")
        self.make_config = self.add_tool("wgs.make_config")
        self.target_dir = ""
        self.soap_denovo_tools = []
        self.assembly_stat = []
        self.insert_bam = []
        self.blastn = []
        self.ref_db = ""

    def check_options(self):
        if not self.option("samples"):
            raise OptionError("缺少samples参数", code="14501201")
        if not self.option("insert_seq"):
            raise OptionError("缺少insert_seq参数", code="14501202")
        if not self.option("genome_version_id"):
            raise OptionError("缺少genome_version_id参数", code="14501203")
        return True

    def get_ref_path(self):
        """
        根据genome_version_id去sg_species_version中找到ref文件
        :return:
        """
        result = self.api.api("wgs.api_base").col_find_one("sg_species_version",
                                                           {"_id": ObjectId(self.option("genome_version_id"))})
        if result:
            self.ref_db = self.json_path + result['ref'].lstrip("/")
        else:
            self.set_error("在sg_species_version中没有找到%s对应信息！", variables=(self.option("genome_version_id")), code="14501213")

    def run_ref_index(self):
        self.ref_index.set_options({
            "insert_fa": self.option("insert_seq"),
            "ref_fa": self.ref_db,
            "dbtype": "nucl"
        })
        self.ref_index.on("end", self.mapping_run)
        self.ref_index.run()

    def mapping_run(self):
        self.mapping.set_options({
            "samples": self.option("samples"),
            "clean_path": self.option('clean_path'),
            "ref_fa": self.ref_index.output_dir + "/pop.fa"
        })
        self.mapping.on('end', self.retrive_insert_bam_run)
        # self.mapping.on("end", self.set_output, 'merged_bam')
        self.mapping.run()

    def retrive_insert_bam_run(self):
        for m in os.listdir(self.mapping.output_dir + "/samtools_merge"):
            retrive_insert_bam = self.add_tool("wgs.retrive_insert_bam")
            retrive_insert_bam.set_options({
                "bam": self.mapping.output_dir + "/samtools_merge/" + m,
                "id": m.split('.')[0],
                'fai': self.ref_index.output_dir + '/pop.fa.fai'
            })
            self.insert_bam.append(retrive_insert_bam)
        for j in range(len(self.insert_bam)):
            self.insert_bam[j].on("end", self.set_output, 'insert_bam_fq')
        if self.insert_bam:
            if len(self.insert_bam) > 1:
                self.on_rely(self.insert_bam, self.make_config_run)
            elif len(self.insert_bam) == 1:
                self.insert_bam[0].on('end', self.make_config_run)
        else:
            self.set_error("insert_bam列表为空！", code="14501214")
        for tool in self.insert_bam:
            gevent.sleep(1)
            tool.run()

    def make_config_run(self):
        """
        要将产品线的fq名字30211907.sca5-2283205-2284356.fq， 改成JY102.chr1:10000-20000.fastq.gz
        :return:
        """
        self.make_config.set_options({
            "fastq_dir":  self.insert_bam[0].output_dir if len(self.insert_bam) == 1
            else self.output_dir + "/insert_bam"
        })
        self.make_config.on("end", self.set_output, 'config')
        self.make_config.on("end", self.soap_denovo_run)
        self.make_config.run()

    def soap_denovo_run(self):
        for file_ in os.listdir(self.make_config.output_dir):
            if file_ == "fastq_list.txt":
                continue
            else:
                soap_denovo = self.add_tool("wgs.soap_denovo")
                soap_denovo.set_options({
                    "config_file": os.path.join(self.make_config.output_dir, file_)
                })
                self.soap_denovo_tools.append(soap_denovo)
        for j in range(len(self.soap_denovo_tools)):
            self.soap_denovo_tools[j].on("end", self.set_output, 'soap_denovo')
        if self.soap_denovo_tools:
            if len(self.soap_denovo_tools) > 1:
                self.on_rely(self.soap_denovo_tools, self.assembly_stat_run)
            elif len(self.soap_denovo_tools) == 1:
                self.soap_denovo_tools[0].on('end', self.assembly_stat_run)
        else:
            self.set_error("soap_denovo_tools列表为空！", code="14501215")
        for tool in self.soap_denovo_tools:
            gevent.sleep(1)
            tool.run()

    def assembly_stat_run(self):
        soap_data = self.set_soap_file()
        for m in soap_data:
            assembly_stat = self.add_tool("wgs.assembly_stat")
            assembly_stat.set_options({
                "denovo_scafseq": m.values()[0].split(';')[0],
                "denovo_scafstatistics": m.values()[0].split(';')[1]
            })
            self.assembly_stat.append(assembly_stat)
        for j in range(len(self.assembly_stat)):
            self.assembly_stat[j].on("end", self.set_output, 'assembly_stat')
        if self.assembly_stat:
            if len(self.assembly_stat) > 1:
                self.on_rely(self.assembly_stat, self.end)
            elif len(self.assembly_stat) == 1:
                self.assembly_stat[0].on('end', self.end)
        else:
            self.set_error("assembly_stat列表为空！", code="14501216")
        for tool in self.assembly_stat:
            gevent.sleep(1)
            tool.run()

    def set_soap_file(self):
        soap_data = []
        if len(self.soap_denovo_tools) == 1:
            path = self.soap_denovo_tools[0].output_dir
        else:
            path = os.path.join(self.output_dir, "soap_denovo")
        result = os.listdir(path)
        for m in result:
            n = re.match(r"(.*)\.denovo\.scafSeq$", m)
            if n:
                if os.path.getsize(os.path.join(path, m)) == 0:
                    continue
                if not os.path.exists(os.path.join(path, "{}.denovo.scafStatistics".format(n.group(1)))):
                    self.set_error("缺少%s文件！".format("%s.denovo.scafStatistics", variables=(n.group(1))), code="14501217")
                soap_data.append({n.group(1): ';'.join([os.path.join(path, m),
                                                        os.path.join(path, "{}.denovo.scafStatistics"
                                                                     .format(n.group(1)))])})
        return soap_data

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'insert_bam_fq':
            self.linkdir(obj.output_dir, self.output_dir + '/insert_bam_fq')
        if event['data'] == 'soap_denovo':
            self.linkdir(obj.output_dir, self.output_dir + '/soap_denovo')
        if event['data'] == 'config':
            self.linkdir(obj.output_dir, self.output_dir + '/config')
        if event['data'] == 'blastn':
            self.linkdir(obj.output_dir, self.output_dir + "/blast")
        if event['data'] == 'merged_bam':
            self.linkdir(obj.output_dir, self.output_dir + "/merged_bam")
        if event['data'] == 'assembly_stat':
            self.linkdir(obj.output_dir, self.output_dir + "/assembly_stat")

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

    def set_db(self):
        self.logger.info("设置assembly的导表！")
        # reads 数量统计表
        file_reads = os.path.join(self.output_dir, "qc_stat")
        self.api.api('wgs.assembly').add_qc_file(file_reads, self.option("main_id"))
        assembly_stat = os.path.join(self.output_dir, "assembly_stat")
        self.api.api('wgs.assembly').add_pc_file(assembly_stat, self.option("main_id"))
        self.api.api('wgs.api_base').update_db_record("sg_assembly",
                                                      {"_id": ObjectId(self.option("main_id"))},
                                                      {"seq_path": self.target_dir + "/soap_denovo"})
        self.api.api('wgs.assembly').add_anno_file(os.path.join(self.output_dir, "final_anno"), self.option("main_id"))
        self.logger.info("设置assembly的导表成功！")

    def run(self):
        self.get_ref_path()
        self.logger.info("初始化ref成功！")
        self.run_ref_index()
        # self.mapping_run()
        # self.retrive_insert_bam_run()
        super(TransgenicWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        # self.set_db()
        # self.set_output_file()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(TransgenicWorkflow, self).end()

    def set_output_file(self):
        """
        设定output的目录中文件结构，删除不需要的文件, 进行自由添加
        :return:
        """
        for files in os.listdir(self.output_dir + "/soap_denovo"):
            m = re.match(r'(.*)\.denovo\.scafSeq$', files)
            if not m:
                os.remove(self.output_dir + "/soap_denovo/{}".format(files))
        # rm_dir = list()
        # rm_dir.append(self.output_dir + "/diamond_go")
        # rm_dir.append(self.output_dir + "/diamond_nog")
        # rm_dir.append(self.output_dir + "/diamond_uniport")
        # rm_dir.append(self.output_dir + "/diamond_kegg")
        # rm_dir.append(self.output_dir + "/diamond_nr")
        # for dirs in rm_dir:
        #     if os.path.exists(dirs):
        #         code = os.system("rm -r {}".format(dirs))
        #         if code == 0:
        #             self.logger.info("删除文件夹{}成功！".format(dirs))

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        if self._sheet.client not in ['client01', 'client03']:
            self.set_error("client%s类型不正确！", variables=(self._sheet.client), code="14501218")
        if self._sheet.client == 'client01':
            self.target_dir = os.path.join("/mnt/ilustre/data", self._sheet.output.strip().split(':')[1])
        else:
            self.target_dir = os.path.join("/mnt/ilustre/tsanger-data", self._sheet.output.split(':')[1])
