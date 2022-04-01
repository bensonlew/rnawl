# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20210901

import os
import re
import json
import glob
import shutil
import zipfile
import pandas as pd
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class MetabolomePredealAgent(Agent):
    """
    代谢预处理
    """
    def __init__(self, parent=None):
        super(MetabolomePredealAgent, self).__init__(parent)
        options = [
            {"name": "lib_search_dir", "type": "infile", "format": "datasplit.path_dir"},  # 搜库文件夹
            {"name": "qcs", "type": "string"},         # 保留的qc样本,逗号分隔
            {"name": "filter_free", "type": "string"}, # 目标代谢物
            {"name": "samples", "type": "string"},     # 样本名称,逗号分隔
            # {"name": "olds", "type": "string"},      # 需修改的样本名称,逗号分隔,同samples
            {"name": "news", "type": "string"},        # 修改后的样本名称,逗号分隔,和samples一一对应
            {"name": "groups", "type": "string"},      # 分组名,逗号分隔,和samples一一对应
            {"name": "controls", "type": "string"},    # 对照组,逗号分隔
            {"name": "others", "type": "string"},      # 实验组,逗号分隔,和controls一一对应
            {"name":"kegg_amino_acid","type":"string","default":"yes"}, #是否保留
            {"name":"sure_score","type":"int","default":30},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
        ]
        self.add_option(options)
        self.queue = "chaifen"

    def check_options(self):
        if not self.option("lib_search_dir").is_set:
            raise OptionError("请设置lib_search_dir")
        if self.option("samples"):
            # raise OptionError("请设置samples")
            if not self.option("news"):
                raise OptionError("设置了samples必须设置news")
            if not self.option("groups"):
                raise OptionError("设置了samples必须设置groups")
            samples = self.option("samples").split(",")
            news = self.option("news").split(",")
            groups = self.option("groups").split(",")
            if len(list(set(samples))) != len(samples):
                raise OptionError("样本名samples:%s里有样本重复,请检查" % self.option("samples"))
            if len(list(set(news))) != len(news):
                raise OptionError("修改样本名news:%s里有样本重复,请检查" % self.option("news"))
            if len(samples) != len(news):
                raise OptionError("修改样本名news长度和样本名samples长度不同,请检查")
            if len(samples) != len(groups):
                raise OptionError("分组名groups长度和样本名samples长度不同,请检查")
        if self.option("controls"):
            if not self.option("others"):
                raise OptionError("设置了controls必须设置others")
            controls = self.option("controls").split(",")
            others = self.option("others").split(",")
            if len(controls) != len(controls):
                raise OptionError("对照组controls长度和实验组others长度不同,请检查")
            cos = [c + '_' + o if c > o else o + '_' + c for c, o in zip(controls, others)]
            if len(list(set(cos))) != len(cos):
                raise OptionError("对照组:%s和实验组:%s有重复", (self.option("controls"), self.option("others")))
            # for i in controls:
            #     if i not in groups:
            #         raise OptionError("对照组里%s不在分组里", i)
            # for i in others:
            #     if i not in groups:
            #         raise OptionError("实验组里%s不在分组里", i)

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(MetabolomePredealAgent, self).end()


class MetabolomePredealTool(Tool):
    def __init__(self, config):
        super(MetabolomePredealTool, self).__init__(config)
        self._version = 1.0
        self.python = "program/Python/bin/python"
        self.python3 = "program/Python35/bin/python"
        python_lib = os.path.join(self.config.SOFTWARE_DIR, "program/Python35/lib")
        python_path = os.path.join(self.config.SOFTWARE_DIR, "program/Python35/bin")
        gcc_path = os.path.join(self.config.SOFTWARE_DIR, 'gcc/5.1.0/bin')
        lib_path = os.path.join(self.config.SOFTWARE_DIR, 'gcc/5.1.0/lib64')
        r_path = os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.5.1/bin')
        self.set_environ(PATH=python_path, LD_LIBRARY_PATH=python_lib)
        self.set_environ(PATH=gcc_path, LD_LIBRARY_PATH=lib_path)
        self.set_environ(PATH=r_path)
        self.del_qc_path = os.path.join(self.config.PACKAGE_DIR, "datasplit/LC_measure_drop_qc.py")
        self.filter_free_path = os.path.join(self.config.PACKAGE_DIR, "datasplit/get_filter_free_list_mohu.py")
        self.predeal_path = os.path.join(self.config.PACKAGE_DIR, "datasplit/predeal_LC_matrix_easy.py")
        # self.norm_path = os.path.join(self.config.PACKAGE_DIR, "tool_lab/meta_Norm.py")
        self.norm_path = os.path.join(self.config.PACKAGE_DIR, "datasplit/meta_Norm.py")
        self.database_dir = os.path.join(self.config.SOFTWARE_DIR, "database/datasplit/metabolome")

    def get_all_files(self):
        """
        从文件夹里获取文件路径
        """
        self.logger.info("获取文件路径")
        self.pos_i, self.pos_q, self.neg_i, self.neg_q, self.group, self.control, self.changename, self.order = [''] * 8
        sk_dir = self.option("lib_search_dir").prop["path"]
        for file in os.listdir(sk_dir):
            if file in ['group.txt', 'group_auto.txt']:
                self.group = os.path.join(sk_dir, file)
            if file in ['control.txt', 'control_auto.txt', 'compare.txt']:
                self.control = os.path.join(sk_dir, file)
            if file in ['rename.txt', 'chname.txt', 'changename.txt', 'changename_auto.txt']:
                self.changename = os.path.join(sk_dir, file)
            if ('neg' in file or 'NEG' in file) and 'tions.csv' in file:
                self.neg_i = os.path.join(sk_dir, file)
            if ('pos' in file or 'POS' in file) and 'tions.csv' in file:
                self.pos_i = os.path.join(sk_dir, file)
            if ('neg' in file or 'NEG' in file) and 'ments.csv' in file:
                self.neg_q = os.path.join(sk_dir, file)
            if ('pos' in file or 'POS' in file) and 'ments.csv' in file:
                self.pos_q = os.path.join(sk_dir, file)
            if file in ['order.txt']:
                self.order = os.path.join(sk_dir, file)
        fen_zip = glob.glob(os.path.join(sk_dir, '*zip'))
        if not self.group and fen_zip:
            zip_file = zipfile.ZipFile(fen_zip[0])
            for fn in zip_file.namelist():
                if fn.endswith('txt'):
                    zip_file.extract(fn, sk_dir)
            self.group = os.path.join(sk_dir, 'group.txt')
            self.control = os.path.join(sk_dir, 'compare.txt')
            self.changename = os.path.join(sk_dir, 'rename.txt')
        if not all([self.neg_i, self.neg_q, self.pos_i, self.pos_q]):
            self.set_db("failed", "下机文件不全,请检查")
            self.set_error('下机文件不全,请检查')
        if not self.group and not self.option("groups"):
            self.set_db("failed", "缺少分组信息")
            self.set_error('缺少分组信息')
        if not self.control and not self.option("controls"):
            self.set_db("failed", "缺少对照组信息")
            self.set_error('缺少对照组信息')

    def run_del_measure_qc(self):
        """
        删除QC,重新生成1.*NEG_Compound Measurements.csv和*POS_Compound Measurements.csv
        """
        with open(self.pos_q) as pr:
            t_header = pr.readline()
            _ = pr.readline()
            header = pr.readline()
            qcs = re.findall('(QC\d+)', header)
            if not qcs:
                qcs = re.findall('(QC\d+)', t_header)
        qcs = list(set(qcs))
        del_qcs, qcs1 = [], []
        if self.option("qcs"):
            qcs1 = self.option("qcs").split(",")
        for qc in qcs:
            if qc not in qcs1:
                del_qcs.append(qc)
        if len(del_qcs) == 0:
            self.logger.info("不用进行删除QC")
            return
        self.logger.info("删除QC:%s" % ",".join(del_qcs))
        neg_q = os.path.join(self.work_dir, os.path.basename(self.neg_q) + "_dropqc.csv")
        pos_q = os.path.join(self.work_dir, os.path.basename(self.pos_q) + "_dropqc.csv")
        del_qcs1 = ",".join(del_qcs)
        cmd = "{} {} -quant ".format(self.python3, self.del_qc_path)
        cmd1 = cmd + "{} -del_qcs {} -sort_qc 0 -drop_qc {}".format(self.neg_q, del_qcs1, neg_q)
        command1 = self.add_command("del_measure_qc_neg", cmd1, ignore_error=True).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("del_measure_qc_neg运行成功")
        else:
            self.set_db("failed", "del_measure_qc_neg运行失败")
            self.set_error("del_measure_qc_neg运行失败")
        cmd2 = cmd + "{} -del_qcs {} -sort_qc 0 -drop_qc {}".format(self.pos_q, del_qcs1, pos_q)
        command2 = self.add_command("del_measure_qc_pos", cmd2, ignore_error=True).run()
        self.wait()
        if command2.return_code == 0:
            self.logger.info("del_measure_qc_pos运行成功")
        else:
            self.set_db("failed", "del_measure_qc_pos运行失败")
            self.set_error("del_measure_qc_pos运行失败")
        self.neg_q = neg_q
        self.pos_q = pos_q

    def get_filter_free_list(self):
        self.logger.info("开始筛选保留物质")
        inter_file = os.path.join(self.work_dir, 'inter.list')
        with open(inter_file, 'wb') as iw:
            for m in self.option("filter_free").split('\n'):
                m = m.strip()
                if m:
                    iw.write(m + '\n')
        ff_file = os.path.join(self.work_dir, 'filter_free.list')
        cmd = "{} {} -pos {} ".format(self.python3, self.filter_free_path, self.pos_i)
        cmd1 = cmd + "-neg {} -inter {} -out {} -outn {}".format(self.neg_i, inter_file, ff_file, self.out_f_n)
        command1 = self.add_command("filter_free_list", cmd1, ignore_error=True).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("filter_free_list运行成功")
        else:
            self.set_db("failed", "filter_free_list运行失败")
            self.set_error("filter_free_list运行失败")

    def get_new_file(self):
        self.logger.info("根据传入样本信息生成文件")
        if self.option("news"):
            old = self.option("samples").split(",")
            news = self.option("news").split(",")
            self.changename = os.path.join(self.work_dir, 'changename_auto.txt')
            with open(self.changename, 'w') as cw:
                for o, n in zip(old, news):
                    cw.write('%s\t%s\n' % (o, n))
        if self.option("samples"):
            # samples = self.option("samples").split(",")
            samples = self.option("news").split(",")
            groups = self.option("groups").split(",")
            have_qc = 'QC' in samples
            self.group = os.path.join(self.work_dir, 'group_auto.txt')
            with open(self.pos_q) as pr:
                t_header = pr.readline()
                _ = pr.readline()
                header = pr.readline()
                qcs = re.findall('(QC\d+)', header)
                if not qcs:
                    qcs = re.findall('(QC\d+)', t_header)
            with open(self.group, 'w') as gw:
                ss_ = []
                gw.write('#sample\tgroup\n')
                for s, g in zip(samples, groups):
                    if s not in ss_:
                        ss_.append(s)
                        gw.write(s + '\t' + g + '\n')
                if not have_qc:
                    for qc in qcs:
                        if qc not in ss_:
                            ss_.append(qc)
                            gw.write('%s\tQC\n' % qc)
        if self.option("controls"):
            controls = self.option("controls").split(",")
            others = self.option("others").split(",")
            # for i in controls:
            #     if i not in groups:
            #         raise OptionError("对照组里%s不在分组里", i)
            # for i in others:
            #     if i not in groups:
            #         raise OptionError("实验组里%s不在分组里", i)
            self.control = os.path.join(self.work_dir, 'control_auto.txt')
            control_other = list(zip(controls, others))
            with open(self.control, 'w') as cw:
                cw.write('#control\tother\n')
                for c, o in control_other:
                    cw.write(c + '\t' + o + '\n')

    def run_predealLc_matrix(self):
        self.logger.info("开始进行预处理")
        if not os.path.exists(self.group):
            self.set_db("failed", "缺少分组信息")
            self.set_error("缺少分组信息")
        sure_score = self.option("sure_score")
        min_row, max_row, min_score, max_score = 500, 800, 30, 50
        max_len, cutoff_percent = 80, "0"
        cmd = "{} {} -pos_i {} -pos_q {} ".format(self.python3, self.predeal_path, self.pos_i, self.pos_q)
        cmd += "-neg_i {} -neg_q {} -group_table {} ".format(self.neg_i, self.neg_q, self.group)
        cmd += "-chname {} -filter_free {} -database_dir {} ".format(self.changename, self.out_f_n, self.database_dir)
        cmd += "-cutoff_percent {} -min_row {} -max_row {} ".format(cutoff_percent, min_row, max_row)
        cmd += "-min_score {} -max_score {} -sure_score {} ".format(min_score, max_score, sure_score)
        cmd += "-max_len {} -output {} ".format(max_len, self.work_dir)
        cmd += "-kegg_amino_acid {} ".format(self.option("kegg_amino_acid"))
        command = self.add_command("predeal_lc_matrix", cmd, ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("predeal_lc_matrix运行成功")
        else:
            self.logger.info("predeal_lc_matrix运行失败")
            self.set_db("failed", "predeal_lc_matrix运行失败")
            self.set_error("predeal_lc_matrix运行失败")

    def run_order_norm(self):
        self.logger.info("存在order,开始进行矫正")
        o_df = pd.read_csv(self.order, sep='\t', dtype={0:str, 1:str})
        g_df = pd.read_csv(self.group, sep='\t', dtype={0:str, 1:str})
        c_d = dict()
        if self.changename and os.path.exists(self.changename):
            c_df = pd.read_csv(self.changename, sep='\t', dtype={0:str, 1:str})
            c_d = dict(zip(c_df.iloc[:, 0], c_df.iloc[:, 1]))
        samples = o_df.iloc[:, 0]
        def m(s):
            try:
                return c_d[s]
            except KeyError:
                return s
        samples = samples.apply(m)
        samples = samples[samples.isin(g_df['#sample'].tolist())]
        self.order = os.path.join(self.work_dir, 'order.txt')
        with open(self.order, 'w') as ow:
            ow.write('#sample\torder\n')
            for n, s in enumerate(samples.tolist()):
                ow.write('%s\t%i\n' % (s, (n + 1)))
        pos = os.path.join(self.work_dir, 'pos_measurement.xls')
        neg = os.path.join(self.work_dir, 'neg_measurement.xls')
        mark, qc_f, smp_f, method, top  = "ID", "0", "0", "svr", "5"
        cmd = "{} {} -pos {} -neg {} -order {}".format(self.python, self.norm_path, pos, neg, self.order)
        cmd += " -mark {} -qc_f {} -smp_f {} -method {} -top {}".format(mark, qc_f, smp_f, method, top)
        cmd += " -rsrc_path {}".format(self.config.PACKAGE_DIR)
        command1 = self.add_command("mult_run_norm", cmd, ignore_error=True).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("mult_run_norm运行成功")
        else:
            self.set_db("failed", "mult_run_norm运行失败")
            self.set_error("mult_run_norm运行失败")

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        group_path = os.path.join(self.output_dir, os.path.basename(self.group))
        control_path = os.path.join(self.output_dir, os.path.basename(self.control))
        os.system("cp {} {}".format(self.group, group_path))
        os.system("cp {} {}".format(self.control, control_path))
        pos_path_ = os.path.join(self.work_dir, "pos_measurement.xls")
        neg_path_ = os.path.join(self.work_dir, "neg_measurement.xls")
        pos_path = os.path.join(self.output_dir, "pos_measurement.xls")
        neg_path = os.path.join(self.output_dir, "neg_measurement.xls")
        self.link_path(pos_path_, pos_path)
        self.link_path(neg_path_, neg_path)
        if self.order:
            deal_pos_path_ = os.path.join(self.work_dir, "pos/Dealed_pos_measurement.xls")
            deal_neg_path_ = os.path.join(self.work_dir, "neg/Dealed_neg_measurement.xls")
            deal_pos_path = os.path.join(self.output_dir, "Dealed_pos_measurement.xls")
            deal_neg_path = os.path.join(self.output_dir, "Dealed_neg_measurement.xls")
            self.link_path(deal_pos_path_, deal_pos_path)
            self.link_path(deal_neg_path_, deal_neg_path)

    def run_md5sum(self):
        """
        生成md5校验码
        """
        fastq_dir = os.path.join(self.output_dir)
        self.fq_list = []
        md5_info = {}
        for f in os.listdir(fastq_dir):
            fq_path = os.path.join(fastq_dir, f)
            if os.path.isdir(fq_path):
                continue
            md5 = os.popen("md5sum {}".format(fq_path)).readlines()[0].split(" ")[0]
            md5_info[f] = md5
        with open(os.path.join(fastq_dir, "md5sum.txt"), "wb") as w:
            for f in md5_info.keys():
                w.write(md5_info[f] + "  " + f + "\n")

    def link_path(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_db(self, status, desc):
        if self.option("main_id"):
            api_db = self.api.api("datasplit.metabolome_predeal")
            data_json = os.path.dirname(self.work_dir) + "/data.json"
            s3_output_dir = json.loads(open(data_json).read())["output"]
            # pos_path = os.path.join(self.s3_output_dir, "pos_measurement.xls")
            # neg_path = os.path.join(self.s3_output_dir, "neg_measurement.xls")
            # group_path = os.path.join(self.s3_output_dir, os.path.basename(self.group))
            # compare_path = os.path.join(self.s3_output_dir, os.path.basename(self.control))
            # if self.order:
            #     deal_pos_path = os.path.join(self.s3_output_dir, "Dealed_pos_measurement.xls")
            #     deal_neg_path = os.path.join(self.s3_output_dir, "Dealed_neg_measurement.xls")
            predeal_id = self.option("main_id")
            # status = "end"
            # desc = "Job has been finished"
            # log = ""
            log = desc
            log_path = os.path.join(self.work_dir, "predeal_lc_matrix.o")
            if os.path.exists(log_path):
                with open(log_path, "rb") as f:
                    log = f.readlines()
            api_db.update_sg_predeal(predeal_id, self.output_dir, s3_output_dir, status, desc, log)

    def run(self):
        super(MetabolomePredealTool, self).run()
        self.get_all_files()
        self.run_del_measure_qc()
        self.out_f_n = os.path.join(self.work_dir, 'filter_free_no_dup.list')
        if self.option("filter_free"):
            self.get_filter_free_list()
        if not os.path.exists(self.group):
            self.get_new_file()
        self.run_predealLc_matrix()
        if self.order:
            self.run_order_norm()
        self.set_output()
        self.run_md5sum()
        self.set_db("end", "Job has been finished")
        self.end()
