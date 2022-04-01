# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
import fitz
import datetime
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config
from bson.objectid import ObjectId
import time

class VcfDistributionWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VcfDistributionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "group_file","type": "infile","format":"denovo_rna_v2.common"},
            {"name":"wsize","type":"float","default":1},
            {"name": "max_missing", "type": "float", "default": 30.0},  # 缺失率
            {"name": "maf", "type": "float", "default": 0.05},  # 次要等位基因频率min
            {"name":"advanced_params","type":"bool"},#用于前端设置高级参数，无实际作用
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.filter = self.add_tool("tool_lab.vcf_filter")
        # self.plink = self.add_tool("tool_lab.vcf_plink")
        self.distribution = self.add_tool('tool_lab.vcf_distribution_v2')

    def check_option(self):
        """
        参数检查
        """
        if self.option('source') == 'tool_lab' and not self.option("vcf_file"):
            raise OptionError("必须输入vcf文件")
        # if not self.option("group_file"):
        #     raise OptionError("必须输入分组文件")
        if not self.option("max_missing"):
            raise OptionError("必须设置最大缺失率")
        if not self.option("maf"):
            raise OptionError("必须设置maf")
    
    def run_filter(self):
        option = {
            'vcf_file':self.vcf_file,   # modified by zhangyitong on 20210819
            'max_missing':self.option("max_missing")/100,
            "maf":self.option('maf')
        }
        # if self.option("group_file").is_set:
        if self.group_file:
            chr_list = []
            with open(self.group_file,'r') as cf:
                while 1:
                    line = cf.readline()
                    if not line:
                        break
                    chr_list.append(line.rstrip())
            if len(chr_list) > 0:
                option["chr"] = ";".join(chr_list)
            else:
                self.set_error("染色体列表为空")
        else:
            option["chr"] = "all"
        self.filter.set_options(option)
        self.filter.on('end',self.run_distribution)
        self.filter.run()

    # def run_plink(self):
    #     self.filter_vcf = os.path.join(self.filter.output_dir,"pop.recode.vcf")
    #     self.plink.set_options({
    #         'vcf_file': self.filter_vcf,
    #     })
    #     self.plink.on('end',self.run_distribution)
    #     self.plink.run()

    def run_distribution(self):
        self.distribution.set_options({
            'vcf_file':os.path.join(self.filter.output_dir,"pop.recode.vcf"),
            'wsize':self.option('wsize')*1000000
        })
        self.distribution.on('end',self.set_output)
        self.distribution.run()
    
    def set_output(self):
        # self.linkdir(self.distribution.output_dir, 'distribution')
        pictrue_path = self.distribution.output_dir + "/Col1.Col0.SNP-Density..pdf"
        pdf = fitz.open(pictrue_path)
        page = pdf[0]
        pm = page.getPixmap(alpha=False)
        pm.writePNG(os.path.join(self.output_dir,"variants.density.png"))
        try:
            os.link(pictrue_path,os.path.join(self.output_dir,"variants.density.pdf"))
        except Exception as e:
            self.set_error("设置结果失败：{}".format(e))
        self.set_db()
    
    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
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
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("开始导表")
        pictrue_path = self._sheet.output + "/variants.density.png"
        self.logger.info(pictrue_path)
        api_vcf_distri = self.api.api("tool_lab.vcf_distribution")
        api_vcf_distri.add_detail(self.option("main_id"),pictrue_path)
        self.end()

    def run(self):
        self.group_file = ''
        if self.option('source') == 'project':
            file_path, chrom = self.check_file_path()
            if exists(file_path):
                self.file_path = file_path
            elif exists(file_path.replace('08SNP', 'SNP')):
                self.file_path = file_path.replace('08SNP', 'SNP')
            else:
                self.set_error('未能找到该分析记录的vcf文件')
            self.vcf_file = self.download_s3_file(self.file_path, 'final.vcf')
            if chrom and exists(chrom):
                self.chrom_path = chrom
                self.group_file = self.download_s3_file(self.chrom_path, 'chrom.txt')
        if self.option('source') == 'tool_lab':
            self.vcf_file = self.option('vcf_file').prop['path']
            if self.option('group_file').is_set:
                self.group_file = self.option('group_file').prop['path']
        self.run_filter()
        super(VcfDistributionWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VcfDistributionWorkflow, self).end()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        file_path = upset['file_path']
        if 'chrom_path' in upset:
            chrom_path = upset['chrom_path']
        else:
            chrom_path = ''
        return file_path, chrom_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path,), code='13700502')
        return to_path


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    add_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f") +\
        str(random.randint(1000,10000))
    data = {
        'name': 'test_distribution',
        'id': 'vcf_distri_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "vcf_file": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/test_WGS/2.distribution/filter.vcf",
        "wsize":1,
        "main_id" : add_time
        }
    }
    wsheet = Sheet(data=data)
    wf = VcfDistributionWorkflow(wsheet)
    wf.run()