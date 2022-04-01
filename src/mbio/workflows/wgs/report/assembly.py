# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180425

import re
import os
import gevent
from biocluster.config import Config
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.file import getsize, exists, list_dir


class AssemblyWorkflow(Workflow):
    """
    交互分析：数据组装
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AssemblyWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "samples", "type": "string"},  # 页面传进来的所有的样本id，样本为逗号分隔
            {"name": "poss", "type": "string"},  # 页面传进来的所有的选择区域，每个局域以逗号分隔，多个区域的时候以|分隔
            {"name": "unmapping", "type": 'string'},  # 页面选择添加的时候为true，否则为false
            {"name": "bam_path", "type": "string"},  # bam文件的列表
            {"name": "species_version_id", "type": "string"},  # 参考基因组配置文件
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"},
            {"name": "is_bucket", "type": "string", "default": "false"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.assembly_modules = []
        self.target_dir = ""
        self.analysis_num = []
        self.anno_blast = []
        self.gene_anno = []
        self.blastn = []
        self.summary_set = []
        self.ref_db = ""
        self.s3transfer = None
        self.assembly_api = None
        self.path_ = ''

    def check_options(self):
        if not self.option("samples"):
            raise OptionError("缺少samples参数", code="14500101")
        if not self.option("poss"):
            raise OptionError("缺少poss参数", code="14500102")
        if not self.option("unmapping"):
            raise OptionError("缺少unmapping参数", code="14500103")
        return True

    def set_which_analysis(self):
        """
        根据samples与poss去设定要运行多少模块
        :return:
        """
        chrlist = self.get_chr_list()
        chr_len = self.get_chr_max_len()
        samples = self.option("samples").split(",")
        poss = self.option("poss").split("|")
        poss = list(set(poss))  # 该值最坏的情况就是空，页面传进来是,,
        for pos in poss:
            chr_ = pos.split(",")[0]
            if chr_ and chr_ not in chrlist:
                self.set_error("染色体编号%s不在该物种基因组中存在！正常结束，请修改参数", variables=(chr_), code="14500125")
        for sample in samples:
            for pos in poss:
                pos_ = self.set_chr_len(pos, chr_len)
                analysis_dict = dict(sample=sample, pos=pos_)
                self.analysis_num.append(analysis_dict)
        self.logger.info(self.analysis_num)

    def download_from_s3_(self):
        """
        从对象存储中下载文件到指定路径
        :return:
        """
        if not os.path.exists(self.work_dir + "/temp"):
            os.mkdir(self.work_dir + "/temp")
        self.logger.info("开始下载对象存储中的文件！")
        transfer = MultiFileTransfer()
        for sample in self.option("samples").split(","):
            source = os.path.join(self.option("bam_path"), "{}.sort.bam".format(sample))
            source1 = os.path.join(self.option("bam_path"), "{}.sort.bam.bai".format(sample))
            if not exists(source) or not exists(source1):
                self.set_error("文件%s不存在！", variables=(source), code="14500126")
            transfer.add_download(source, '{}/temp/'.format(self.work_dir))
            transfer.add_download(source1, '{}/temp/'.format(self.work_dir))
        transfer.perform()
        self.logger.info("下载对象存储中的文件成功！")

    def get_chr_list(self):
        chrlist_path = Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/" + \
                       os.path.join(os.path.dirname(os.path.dirname(self.ref_db)).lstrip('/'), "total.chrlist")
        chrlist = []
        self.logger.info(chrlist_path)
        with open(chrlist_path, "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                chrlist.append(temp[0])
        return chrlist

    def get_chr_max_len(self):
        """
        获取每个染色体的最大值
        :return:
        """
        ref_dict = Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/" + os.path.join(
            os.path.dirname(os.path.dirname(self.ref_db)).lstrip('/'), "ref.dict")
        if not os.path.exists(ref_dict):
            self.set_error("文件%s不存在！", variables=(ref_dict), code="14500127")
        chr_len = {}
        with open(ref_dict, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split("\t")
                chr_name = temp[1].strip().split(":")[1]
                length = temp[2].strip().split(':')[1]
                chr_len[chr_name] = int(length)
        return chr_len

    def set_chr_len(self, pos, chr_len):
        # self.logger.info(chr_len)
        temp = pos.split(",")
        if temp[1] or temp[2]:
            if not temp[0]:
                self.set_error("输入参数%s不合法,有起始碱基位置与终止碱基位置，必须要有染色体编号！", variables=(pos), code="14500128")
        if temp[1]:
            try:
                start_ = int(temp[1])
            except Exception:
                self.set_error("输入参数不合格，起始碱基位置必须为整型", code="14500129")
            else:
                if start_ < 0 or start_ > chr_len[temp[0]]:
                    start = "1"
                else:
                    start = temp[1]
        else:
            start = ''
        if temp[2]:
            try:
                end_ = int(temp[2])
            except Exception:
                self.set_error("输入参数不合格，终止碱基位置必须为整型", code="14500130")
            else:
                if end_ < 0:
                    end = ""
                elif end_ > chr_len[temp[0]]:
                    end = str(chr_len[temp[0]])
                else:
                    end = temp[2]
        else:
            end = ''
        self.logger.info("pos:{}".format(','.join([temp[0], start, end])))
        return ','.join([temp[0], start, end])

    def check_bam_is_ok(self):
        """
        检查对应的bam文件存在
        :return:
        """
        if self.option("is_bucket") == "true":
            self.path_ = os.path.join(self.work_dir, "temp/")
        else:
            self.path_ = self.option("bam_path")
        self.logger.info("self.path_:{}".format(self.path_))
        self.logger.info("samples:{}".format(self.option("samples").split(",")))
        for m in self.option("samples").split(","):
            if not os.path.exists(os.path.join(self.path_, "{}.sort.bam".format(m))):
                self.set_error("文件%s不存在！".format(os.path.join(self.path_, "%s.sort.bam", variables=(m))), code="14500131")

    def run_assembly(self):
        for item in self.analysis_num:
            assembly = self.add_module("wgs.assembly")
            options = {
                "bam_file": os.path.join(self.path_, "{}.sort.bam".format(item['sample'])),
                "sample_id": item['sample'],
                "pos": item['pos'],
                "unmapping": self.option("unmapping")
            }
            assembly.set_options(options)
            self.assembly_modules.append(assembly)
        for j in range(len(self.assembly_modules)):
            self.assembly_modules[j].on("end", self.set_output, 'assembly')
        if self.assembly_modules:
            if len(self.assembly_modules) > 1:
                self.on_rely(self.assembly_modules, self.gene_anno_run)
            elif len(self.assembly_modules) == 1:
                self.assembly_modules[0].on('end', self.gene_anno_run)
        else:
            self.set_error("assembly_modules列表为空！", code="14500132")
        for module in self.assembly_modules:
            gevent.sleep(1)
            module.run()

    def gene_anno_run(self):
        soap_data = self.set_soap_file()
        for m in soap_data:
            anno = self.add_module("wgs.gene_anno")
            options = {
                "fasta": m.values()[0],
                "sample_id": m.keys()[0]
            }
            anno.set_options(options)
            self.gene_anno.append(anno)
        for j in range(len(self.gene_anno)):
            self.gene_anno[j].on("end", self.set_output, 'gene_anno')
        if self.gene_anno:
            if len(self.gene_anno) > 1:
                self.on_rely(self.gene_anno, self.blast_run)
            elif len(self.gene_anno) == 1:
                self.gene_anno[0].on('end', self.blast_run)
        else:
            self.set_error("gene_anno列表为空！", code="14500133")
        for module in self.gene_anno:
            gevent.sleep(1)
            module.run()

    def set_soap_file(self):
        soap_data = []
        path = os.path.join(self.output_dir, "soap_denovo")
        result = os.listdir(path)
        for m in result:
            n = re.match(r"(.*)\.denovo\.scafSeq$", m)
            if n:
                if os.path.getsize(os.path.join(path, m)) == 0:
                    continue
                soap_data.append({n.group(1): os.path.join(path, m)})
        return soap_data

    def blast_run(self):
        soap_data = self.set_soap_file()
        for m in soap_data:
            blastn = self.add_tool("wgs.blast_n")
            options = {
                "query_fa": m.values()[0],
                "sample_id": m.keys()[0],
                "dbname_nsq": Config().SOFTWARE_DIR + "/database/dna_wgs_geneome/" + self.ref_db,
                "outfmt": "6"
            }
            blastn.set_options(options)
            self.blastn.append(blastn)
        for j in range(len(self.gene_anno)):
            self.blastn[j].on("end", self.set_output, 'blastn')
        if self.blastn:
            if len(self.blastn) > 1:
                self.on_rely(self.blastn, self.summary_set_run)
            elif len(self.blastn) == 1:
                self.blastn[0].on('end', self.summary_set_run)
        else:
            self.set_error("blastn列表为空！", code="14500134")
        for module in self.blastn:
            gevent.sleep(1)
            module.run()

    def summary_set_run(self):
        for m in os.listdir(self.output_dir + "/gene_anno"):
            if not os.path.exists(self.output_dir + "/blast/{}.blast".format(m.split('.')[0])):
                self.set_error("文件%s不存在！".format(self.output_dir + "/blast/%s.blast", variables=(m.split('.')[0])), code="14500135")
            summary_set_ = self.add_tool("wgs.summary_set")
            options = {
                "anno_summary": self.output_dir + "/gene_anno/{}".format(m),
                "sample_id": m.split('.')[0],
                "blast_result": self.output_dir + "/blast/{}.blast".format(m.split('.')[0])
            }
            summary_set_.set_options(options)
            self.summary_set.append(summary_set_)
        for j in range(len(self.summary_set)):
            self.summary_set[j].on("end", self.set_output, 'final_anno')
        if self.summary_set:
            if len(self.summary_set) > 1:
                self.on_rely(self.summary_set, self.end)
            elif len(self.summary_set) == 1:
                self.summary_set[0].on('end', self.end)
        else:
            self.set_error("summary_set列表为空！", code="14500136")
        for tool in self.summary_set:
            gevent.sleep(1)
            tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'assembly':
            self.linkdir(obj.output_dir, self.output_dir)
        if event['data'] == 'gene_anno':
            self.linkdir(obj.output_dir, self.output_dir)
        if event['data'] == 'final_anno':
            self.linkdir(obj.output_dir, self.output_dir + '/final_anno')
        if event['data'] == 'blastn':
            self.linkdir(obj.output_dir, self.output_dir + "/blast")

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
                # else:
                #     os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                for file_ in os.listdir(oldfiles[i]):
                    if not os.path.exists(os.path.join(newdir, os.path.basename(oldfiles[i]))):
                        os.mkdir(os.path.join(newdir, os.path.basename(oldfiles[i])))
                    temp_path = os.path.join(newdir, os.path.basename(oldfiles[i]))
                    if os.path.exists(os.path.join(temp_path, file_)):
                        os.remove(os.path.join(temp_path, file_))
                    os.link(os.path.join(oldfiles[i], file_), os.path.join(temp_path, file_))

    def set_db(self):
        self.logger.info("设置assembly的导表！")
        # reads 数量统计表
        file_reads = os.path.join(self.output_dir, "qc_stat")
        self.assembly_api.add_qc_file(file_reads, self.option("main_id"))
        assembly_stat = os.path.join(self.output_dir, "assembly_stat")
        self.assembly_api.add_pc_file(assembly_stat, self.option("main_id"))
        base_api = self.api.api('wgs.api_base')
        if self.option("project_type"):
            base_api._project_type = self.option("project_type")
        base_api.update_db_record("sg_assembly", {"_id": ObjectId(self.option("main_id"))},
                                  {"seq_path": self._sheet.output.rstrip('/') + "/soap_denovo/"})
        if len(self.summary_set) == 1:
            self.assembly_api.add_anno_file(self.summary_set[0].output_dir,  self.option("main_id"))
        else:
            self.assembly_api.add_anno_file(os.path.join(self.output_dir, "final_anno"), self.option("main_id"))
        self.logger.info("设置assembly的导表成功！")

    def run(self):
        self.assembly_api = self.api.api('wgs.assembly')
        if self.option("is_bucket") == "true":
            self.download_from_s3_()
        if self.option("project_type"):
            self.assembly_api._project_type = self.option("project_type")
        self.ref_db = self.assembly_api.get_ref_db(self.option("species_version_id"))
        self.check_bam_is_ok()
        self.set_which_analysis()
        self.run_assembly()
        super(AssemblyWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        self.set_db()
        self.set_output_file()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssemblyWorkflow, self).end()

    def set_output_file(self):
        """
        设定output的目录中文件结构，删除不需要的文件, 进行自由添加
        :return:
        """
        for files in os.listdir(self.output_dir + "/soap_denovo"):
            m = re.match(r'(.*)\.denovo\.scafSeq$', files)
            if not m:
                os.remove(self.output_dir + "/soap_denovo/{}".format(files))
        rm_dir = list()
        rm_dir.append(self.output_dir + "/diamond_go")
        rm_dir.append(self.output_dir + "/diamond_nog")
        rm_dir.append(self.output_dir + "/diamond_uniport")
        rm_dir.append(self.output_dir + "/diamond_kegg")
        rm_dir.append(self.output_dir + "/diamond_nr")
        for dirs in rm_dir:
            if os.path.exists(dirs):
                code = os.system("rm -r {}".format(dirs))
                if code == 0:
                    self.logger.info("删除文件夹{}成功！".format(dirs))

    # def get_target_dir(self):
    #     """
    #     获取远程磁盘的路径
    #     :return:
    #     """
    #     if self._sheet.client not in ['client01', 'client03']:
    #         raise Exception("client{}类型不正确！".format(self._sheet.client))
    #     if self._sheet.client == 'client01':
    #         self.target_dir = os.path.join("/mnt/ilustre/data", self._sheet.output.strip().split(':')[1])
    #     else:
    #         self.target_dir = os.path.join("/mnt/ilustre/tsanger-data", self._sheet.output.split(':')[1])
