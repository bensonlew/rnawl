# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180514

import os
import re
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class BlastSeqWorkflow(Workflow):
    """
    交互分析：局部组装与转基因部分的blast接口与序列下载接口--首先获取到需要比对的序列，然后进行比对
    接口中需要，sg_assembly的主表id，task_id， sample_id, scaffold_id
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BlastSeqWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sample_id", "type": "string"},  # 样本id
            {"name": "scaffold_id", "type": "string"},  # 序列id
            {"name": "seq_path", "type": "infile", 'format': 'bsa.dir'},  # 序列的路径
            {"name": "ref_db", "type": "string"},  # 构建的参考组的db文件，makeblastdb生成, 也可能是插入片段的序列
            {"name": "types", "type": "string", "default": "blast"},
            {"name": "blast_type", "type": "string", "default": "ref"},  # ref/insert
            {"name": "insert_seq", "type": "infile", "format": "sequence.fasta"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.blast_n = self.add_tool("wgs.blast_n")
        self.makeblastdb = self.add_tool("wgs.makeblastdb")
        self.target_dir = ""

    def check_options(self):
        if not self.option('sample_id'):
            raise OptionError('必须输入样本id', code="14500201")
        if not self.option('scaffold_id'):
            raise OptionError('必须输入序列id', code="14500202")
        if not self.option('seq_path'):
            raise OptionError('必须输入seq_path所在路径', code="14500203")
        if not self.option('ref_db'):
            raise OptionError('必须输入db路径', code="14500204")
        if self.option("types") not in ['blast', "seq"]:
            # raise OptionError("分析类型{}必须为blast或者seq".format(self.option("types")))
            raise OptionError("分析类型%s必须为blast或者seq", variables=("anno_analysis"), code="14500205")
        return True

    def blast_n_run(self):
        options = {
            "query_fa": os.path.join(self.output_dir, "seq.fa"),
            "dbname_nsq": self.option("ref_db") if self.option("blast_type") == "ref" else
            os.path.join(self.work_dir, "insert_blastdb/insert_seq")
        }
        self.blast_n.set_options(options)
        self.blast_n.on("end", self.set_output, "blast")
        # self.blast_n.on("end", self.set_db)
        self.blast_n.run()

    def get_fasta(self):
        """
        更具scanfload去获取对应序列
        /mnt/ilustre/users/sanger-dev/workspace/20180510/Assembly_wgs_test_0510094133_9523_7275/output/soap_denovo
        :return:
        """
        self.logger.info("开始获取比对序列")
        seq_path = os.path.join(self.output_dir, "seq.fa")
        file_path = os.path.join(self.target_dir, "{}.denovo.scafSeq".format(self.option("sample_id")))
        if not os.path.exists(file_path):
            self.set_error("文件%s不存在！", variables=(file_path), code="14500207")
        with open(file_path, "r") as r, open(seq_path, "w") as w:
            data = r.readlines()
            self.logger.info("dd:{}".format(len(data)))
            n = 0
            for line in data:
                # self.logger.info("1:{}".format(line))
                m = re.match(r'^>{}'.format(self.option('scaffold_id')), line)
                if m:
                    # self.logger.info("2:{}".format('^>{}'.format(self.option('scaffold_id'))))
                    w.write(line)
                    break
                n += 1
            # self.logger.info("3:{}".format(n))
            while True:
                if len(data) == n + 1:
                    self.logger.info("cc:{}".format(n))
                    break
                if not re.match(r'^>.*', data[n+1]):
                    # self.logger.info("4:{}".format(data[n+1]))
                    w.write(data[n+1])
                    n += 1
                else:
                    break
        if os.path.getsize(seq_path) == 0:
            self.set_error("文件%s大小为0，不能进行后面的计算！", variables=(seq_path), code="14500208")
        self.logger.info("获取比对序列成功")
        if self.option("types") != "blast":
            self.set_db()

    def makeblastdb_run(self):
        self.makeblastdb.set_options({
            "pop_fa": os.path.join(self.work_dir, "insert_blastdb/insert_seq.fa")
        })
        self.makeblastdb.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'blast':
            self.linkdir(obj.output_dir, 'blast')
        self.set_db()

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
        self.logger.info("保存结果到mongo！")
        api = self.api.api("wgs.api_base")
        if self.option("project_type"):
            api._project_type = self.option("project_type")
        if self.option("types") == "blast":
            update_dict = {"file_path": "{}/blast/blastn.blast".format(self._sheet.output.rstrip('/'))}
        else:
            update_dict = {"file_path": "{}/seq.fa".format(self._sheet.output.rstrip('/'))}
        api.update_db_record("sg_blast_seq", {"_id": ObjectId(self.option("main_id"))}, update_dict)
        self.logger.info("保存结果到mongo成功！")
        self.end()

    def run(self):
        self.get_target_dir()
        if self.option("types") == "blast":
            self.get_fasta()
            if self.option("blast_type") == "ref":
                self.blast_n_run()
            else:
                self.cp_insert_fa_work_dir()
                self.makeblastdb.on("end", self.blast_n_run)
                self.makeblastdb_run()
            super(BlastSeqWorkflow, self).run()
        else:
            self.start_listener()
            self.fire("start")
            self.get_fasta()
            # self.end()

    def cp_insert_fa_work_dir(self):
        if not os.path.exists(os.path.join(self.work_dir, "insert_blastdb")):
            os.mkdir(os.path.join(self.work_dir, "insert_blastdb"))
        file_ = os.path.join(self.work_dir, "insert_blastdb/insert_seq.fa")
        if os.path.exists(file_):
            os.remove(file_)
        code = os.system("cp {} {}".format(self.option("insert_seq").prop['path'], file_))
        if code != 0:
            self.set_error("复制%s到%s失败！", variables=(self.option("insert_seq").prop['path'], file_), code="14500209")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(BlastSeqWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        self.target_dir = self.option("seq_path").prop['path']
