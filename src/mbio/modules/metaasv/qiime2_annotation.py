# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
import os
import gevent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir


class Qiime2AnnotationModule(Module):
    """
    metaasv qiime2 比对+注释
    核酸类型
    """
    def __init__(self, work_id):
        self.DATABASE = ['unite7.2/its_fungi',"unite8.0/its_fungi",
                         'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA', 'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                         'maarjam081/AM', 'Protist_PR2_v4.5',
                         'silva132/16s_archaea', 'silva132/16s_bacteria','silva132/18s_eukaryota', 'silva132/16s',
                         'silva138/16s_archaea', 'silva138/16s_bacteria','silva138/18s_eukaryota', 'silva138/16s',
                         'greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria',
                         'nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria', 'nt_v20210917/16s',
                         'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi', 'nt_v20210917',
                         'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea', 'nt', 'nt/16s','nt/18S','nt/its','nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s','Human_HOMD_v15.2', 'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi',"nt_v20200604",
                         'fgr/amoA_archaea_202012','fgr/amoA_bacteria_202012','fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012','fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012','fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012', 'fgr/mmoX_202012']
        super(Qiime2AnnotationModule, self).__init__(work_id)
        options = [
            {"name": "qza_fasta", "type": "infile", "format": "metaasv.qza"},  # 输入文件
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 输入fasta文件
            {"name": "database", "type": "string"}, ##输入的数据库
            {"name": "confidence", "type": "float", "default": 0.7},  # 置信度
            {"name": "identity", "type": "float", "default": 0.8},  # blast参数 identity
            {"name": "coverage", "type": "float", "default": 0.8},  # blast参数 coverage
            {"name": "ref_fasta", "type": "infile", 'format': 'sequence.fasta'}, # 输入fasta文件  database 为custom
            {"name": "ref_taxon", "type": "infile", 'format': 'taxon.seq_taxon'},  # 输入taxonomy文件  database为custom
            {'name': 'anno_method', 'type': 'string'},  # 注释方法
            {"name": "top_num", "type": "int", "default": 1},  # 序列比对最大输出条数，默认1条
            {"name": "num_threads", "type": "int", "default": 6},  # cpu数
        ]
        self.add_option(options)
        self.software_dir = Config().SOFTWARE_DIR
        self.convert_file = []

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("qza_fasta").is_set and (not self.option("fasta").is_set):
            raise OptionError("必须设置输入文件，请检查！")
        if self.option("database") == 'custom_mode':
            if not self.option("ref_fasta").is_set:
                raise OptionError("使用自定义数据库模式时必须设置ref_fasta")
            if not self.option("ref_taxon").is_set:
                raise OptionError("使用自定义数据库模式时必须设置ref_taxon")
        else:
            if self.option("database") not in self.DATABASE:
                raise OptionError("数据库%s不被支持", variables=(self.option("database")))
        if not 1.0 > self.option('identity') >= 0:
            raise OptionError('identity值设定必须为[0-1)之间：%s', variables=(self.option('identity')))
        if not 1.0 > self.option('coverage') >= 0:
            raise OptionError('coverage值设定必须为[0-1)之间：%s', variables=(self.option('coverage')))
        if not 1 >= self.option('confidence') >= 0:
            raise OptionError('confidence值设定必须为[0-1)之间：%s', variables=(self.option('confidence')))

        return True

    def run_convert_format(self):
        """
        对上传的taxon和ref_fasta序列进行转格式；
        fasta和taxonomy
        :return:
        """
        self.convert = self.add_tool("metaasv.file_to_qza")
        fasta_options = {
            "input_file": self.option("ref_fasta").prop['path'],
            "type": "fasta",
            "prefix": "ref_fasta"
        }
        self.convert.set_options(fasta_options)
        self.convert_file.append(self.convert)
        self.convert_taxon = self.add_tool("metaasv.file_to_qza")
        taxon_options = {
            "input_file": self.option("ref_taxon").prop['path'],
            "type": "taxonomy",
            "prefix": "ref_taxon"
        }
        self.convert_taxon.set_options(taxon_options)
        self.convert_file.append(self.convert_taxon)
        if self.option("anno_method") in ['blast']:
            self.on_rely(self.convert_file, self.run_blast)
        elif self.option("anno_method") in ['vsearch']:
            self.on_rely(self.convert_file, self.run_vsearch)
        elif self.option("anno_method") in ['bayes']:
            self.on_rely(self.convert_file, self.run_bayes)
        for tool in self.convert_file:
            tool.run()
            gevent.sleep(0)

    def run_convert_fasta(self):
        """
        流程2对fasta文件转格式
        :return:
        """
        self.convert_fasta = self.add_tool("metaasv.file_to_qza")
        fasta_options = {
            "input_file": self.option("fasta").prop['path'],
            "type": "fasta",
            "prefix": "ref_fasta"
        }
        self.convert_fasta.set_options(fasta_options)
        if self.option("anno_method") in ['blast']:
            self.convert_fasta.on("end", self.run_blast)
        elif self.option("anno_method") in ['vsearch']:
            self.convert_fasta.on("end", self.run_vsearch)
        elif self.option("anno_method") in ['bayes']:
            self.convert_fasta.on("end", self.run_bayes)
        self.convert_fasta.run()

    def run_blast(self):
        """
        blast 方法比对和注释
        :return:
        """
        self.blast = self.add_tool("metaasv.qiime2_blast")
        options = {
            "input_qza": self.option("qza_fasta"),
            "database": self.db_path,
            "database_type": self.option("database"),
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "top_num": self.option("top_num")
        }
        if self.option("qza_fasta").is_set:
            options["input_qza"] = self.option("qza_fasta")
        else:
            options["input_qza"] = os.path.join(self.convert_fasta.output_dir, "ref_fasta.qza")
        if self.option("database") in ["custom_mode"]:
            options["ref_fasta"] = os.path.join(self.convert.output_dir, "ref_fasta.qza")
            options["ref_taxon"] = os.path.join(self.convert_taxon.output_dir, "ref_taxon.qza")
        self.blast.set_options(options)
        self.blast.on("end", self.set_output)
        self.blast.run()

    def run_vsearch(self):
        """
        vsearch 方法比对和注释
        :return:
        """
        self.vsearch = self.add_tool("metaasv.qiime2_vsearch")
        options = {
            "database": self.db_path,
            "database_type": self.option("database"),
            "identity": float(self.option("identity")) * 100,
            "coverage": float(self.option("coverage")) * 100,
            "top_num": self.option("top_num"),
            "num_threads": self.option("num_threads"),
        }
        if self.option("qza_fasta").is_set:
            options["input_qza"] = self.option("qza_fasta")
        else:
            options["input_qza"] = os.path.join(self.convert_fasta.output_dir, "ref_fasta.qza")
        if self.option("database") in ["custom_mode"]:
            options["ref_fasta"] = os.path.join(self.convert.output_dir, "ref_fasta.qza")
            options["ref_taxon"] = os.path.join(self.convert_taxon.output_dir, "ref_taxon.qza")
        self.vsearch.set_options(options)
        self.vsearch.on("end", self.set_output)
        self.vsearch.run()

    def run_bayes(self):
        """
        bayes 方法比对和注释
        :return:
        """
        self.bayes = self.add_tool("metaasv.qiime2_bayes")
        options = {
            "database": self.db_path,
            "database_type": self.option("database"),
            "confidence": self.option("confidence"),
            "num_threads": self.option("num_threads")
        }
        if self.option("qza_fasta").is_set:
            options["input_qza"] = self.option("qza_fasta")
        else:
            options["input_qza"] = os.path.join(self.convert_fasta.output_dir, "ref_fasta.qza")
        if self.option("database") in ["custom_mode"]:
            options["ref_fasta"] = os.path.join(self.convert.output_dir, "ref_fasta.qza")
            options["ref_taxon"] = os.path.join(self.convert_taxon.output_dir, "ref_taxon.qza")
        self.bayes.set_options(options)
        self.bayes.on("end", self.set_output)
        self.bayes.run()

    def set_output(self):
        """
        设置结果文件
        :return:
        """
        if self.option("anno_method") in ['blast']:
            link_dir(self.blast.output_dir, self.output_dir)
        elif self.option("anno_method") in ['vsearch']:
            link_dir(self.vsearch.output_dir, self.output_dir)
        elif self.option("anno_method") in ['bayes']:
            link_dir(self.bayes.output_dir, self.output_dir)
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(Qiime2AnnotationModule, self).run()
        if self.option("database") in ['custom_mode']:
            if self.option("anno_method") in ['blast']:
                self.db_path = os.path.join(self.software_dir, "database/taxon_db/qiime2_qza")
            elif self.option("anno_method") in ['vsearch']:
                self.db_path = os.path.join(self.software_dir, "database/taxon_db/qiime2_qza")
            elif self.option("anno_method") in ['bayes']:
                self.db_path = os.path.join(self.software_dir, "database/taxon_db/qiime2_naive_bayes")
            self.run_convert_format()
        else:
            if self.option("anno_method") in ['blast']:
                self.db_path = os.path.join(self.software_dir, "database/taxon_db/qiime2_qza")
                if self.option("fasta").is_set:
                    self.run_convert_fasta()
                else:
                    self.run_blast()
            elif self.option("anno_method") in ['vsearch']:
                self.db_path = os.path.join(self.software_dir, "database/taxon_db/qiime2_qza")
                if self.option("fasta").is_set:
                    self.run_convert_fasta()
                else:
                    self.run_vsearch()
            elif self.option("anno_method") in ['bayes']:
                self.db_path = os.path.join(self.software_dir, "database/taxon_db/qiime2_naive_bayes")
                if self.option("fasta").is_set:
                    self.run_convert_fasta()
                else:
                    self.run_bayes()

    def end(self):
        """
        结束
        :return:
        """
        super(Qiime2AnnotationModule, self).end()
