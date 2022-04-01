# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
from Bio import SeqIO
import unittest
import time
import re
import os
import glob
import sys
import shutil

class BacDatapretonbAgent(Agent):
    def __init__(self, parent):
        super(BacDatapretonbAgent, self).__init__(parent)
        options = [
            {"name": "list_url", "type": "string"},
            {"name": "sample_info_url", "type": "string"},
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_option(self):
        """
        参数检查
        """
        if not self.option("list_url"):
            raise OptionError("没有找到fasta_dir")
        if not self.option("sample_info_url"):
            raise OptionError("没有找到list.xls")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(BacDatapretonbAgent, self).end()

class BacDatapretonbTool(Tool):
    def __init__(self, config):
        super(BacDatapretonbTool, self).__init__(config)
        self.all_fasta_dir = "/mnt/ilustre/users/sanger-dev/junjian/fasta"
        self.pictrue_dir = "/mnt/ilustre/users/sanger-dev/junjian/picture"
        # self.all_fasta_dir = "/mnt/lustre/users/sanger-dev/sg-users/wjwang/bac/fasta"
        # self.pictrue_dir = "/mnt/lustre/users/sanger-dev/sg-users/wjwang/bac/picture"

    def run(self):
        super(BacDatapretonbTool,self).run()
        self.download_file()
        self.get_fasta_dir()
        self.set_output()
        self.end()

    def download_file(self):
        """
        从前端给的下载链接下载
        """
        list_link = self.option("list_url")
        sample_info_link = self.option("sample_info_url")
        os.chdir(self.work_dir)
        download_list = os.system("wget {} -o list_download.log".format(list_link))
        if download_list == 0:
            self.logger.info("下载list.xls成功")
        else:
            self.set_error("下载list.xls失败")
        download_sample_info = os.system("wget {} -o sample_info_link".format(sample_info_link))
        if download_sample_info == 0:
            self.logger.info("下载sample_info表成功")
        else:
            self.set_error("下载sample_info失败")

    def get_fasta_dir(self):
        list_file = os.path.join(self.work_dir, "list.xls")
        sn_primer = {}
        fasta_dir = os.path.join(self.work_dir,"fasta_dir")
        os.mkdir(fasta_dir)
        with open(list_file,"r") as lf:
            while 1:
                file_list = []
                file_name = ""
                n = 0
                line = lf.readline()
                if not line:
                    break
                line = line.strip().split('\t')
                if len(line) < 5: 
                    continue
                mjNo = line[1]
                if sn_primer.has_key(mjNo):
                    sn_primer[mjNo].append(line[2])
                else:
                    sn_primer[mjNo] = [line[2]]
                for ff in glob.glob(self.all_fasta_dir+'/*'+mjNo+'*.ab1')+glob.glob(self.all_fasta_dir+'/*'+mjNo+'*.seq'):
                    fn = os.path.basename(ff)
                    if re.search(line[2],fn):
                        if not fn.startswith(mjNo):    
                            file_name = mjNo + '-' + line[2] + "."+ fn.split(".")[-1]
                        else:
                            file_name = mjNo + '-' + line[2] + "."+ fn.split(".")[-1]
                        file_list.append(file_name)
                    else:
                        continue     
                    try:
                        shutil.copyfile(ff,os.path.join(fasta_dir, file_name))
                    except Exception as e:
                        self.logger.info("设置结果目录失败{}".format(e))
                    if fn.endswith(".seq"):
                        os.system('if [ `head -1 \"'+os.path.join(fasta_dir,file_name)+'\"|grep "^>"` ];then true;else sed -i "1i>'+mjNo+ "-"+ line[2]+ '" \"'+os.path.join(fasta_dir,file_name)+'\";fi')
                if len(file_list) <2:
                    self.logger.info(file_list)
                    time.sleep(60)
                    self.logger.info("1 min")
                    # 找不到样本， 可能是横杠的问题，将下划线改为横杠
                    mjNo1 = mjNo.replace('_','-')
                    file_list = []
                    for ff in glob.glob(self.all_fasta_dir+'/*'+mjNo1+'*.ab1')+glob.glob(self.all_fasta_dir+'/*'+mjNo1+'*.seq'):
                        fn = os.path.basename(ff)
                        if re.search(line[2],fn):
                            if not fn.startswith(mjNo1):    
                                file_name = mjNo + '-' + line[2] + "."+ fn.split(".")[-1]
                            else:
                                file_name = mjNo + '-' + line[2] + "."+ fn.split(".")[-1]
                            file_list.append(file_name)
                        else:
                            continue     
                        try:
                            shutil.copyfile(ff,os.path.join(fasta_dir, file_name))
                        except Exception as e:
                            self.logger.info("设置结果目录失败{}".format(e))
                        if fn.endswith(".seq"):
                            os.system('if [ `head -1 \"'+os.path.join(fasta_dir,file_name)+'\"|grep "^>"` ];then true;else sed -i "1i>'+mjNo+ "-"+ line[2]+ '" \"'+os.path.join(fasta_dir,file_name)+'\";fi')
                    if len(file_list) <2:
                        self.logger.info(file_list)
                        self.logger.info("5 min")
                        time.sleep(300)
                        file_list = []
                        for ff in glob.glob(self.all_fasta_dir+'/*'+mjNo+'*.ab1')+glob.glob(self.all_fasta_dir+'/*'+mjNo+'*.seq'):
                            fn = os.path.basename(ff)
                            if re.search(line[2],fn):
                                if not fn.startswith(mjNo):    
                                    file_name = mjNo + '-' + line[2] + "."+ fn.split(".")[-1]
                                else:
                                    file_name =  mjNo + '-' + line[2] + "."+ fn.split(".")[-1]
                                file_list.append(file_name)
                            else:
                                continue     
                            try:
                                shutil.copyfile(ff,os.path.join(fasta_dir, file_name))
                            except Exception as e:
                                self.logger.info("设置结果目录失败{}".format(e))
                            if fn.endswith(".seq"):
                                os.system('if [ `head -1 \"'+os.path.join(fasta_dir,file_name)+'\"|grep "^>"` ];then true;else sed -i "1i>'+mjNo+ "-"+ line[2]+ '" \"'+os.path.join(fasta_dir,file_name)+'\";fi')
                        if len(file_list) <2:
                            self.set_error("样本编号{}的{}引物的原始数据查找不到,请检查数据，如数据没有问题，请在十分钟后重新投递任务".format(line[1],line[2])) 
        for i in sn_primer.keys():
            if len(sn_primer[i]) < 2:
                self.set_error("{}样本序列少于2，请检查list".format(i))


    def set_output(self):
        with open(os.path.join(self.work_dir, "sample.detail.xls"),"r") as sinfo:
            sinfo.readline()
            line = sinfo.readline()
            temp= line.split("\t")
            temp= line.split("\t")
            if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                new_pictrue_name = "{}.jpg".format(temp[4])
            else:
                time.sleep(60)
                self.logger.info("1 min")
                if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                    pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                    new_pictrue_name = "{}.jpg".format(temp[4])
                else:
                    time.sleep(180)
                    self.logger.info("3 min")
                    if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                        pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                        new_pictrue_name = "{}.jpg".format(temp[4])
                    else:
                        time.sleep(300)
                        self.logger.info("5 min")
                        if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                            pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                            new_pictrue_name = "{}.jpg".format(temp[4])
                        else:
                            time.sleep(600)
                            self.logger.info("10 min")
                            if len(glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))) > 0:
                                pictrue_path = glob.glob(os.path.join(self.pictrue_dir+'/'+temp[4]+'='+'*.jpg'))[0]
                                new_pictrue_name = "{}.jpg".format(temp[4])
                            else:
                                self.set_error("找不到订单号{}的胶图,请检查胶图是否上传，如已上传，请10分钟后重新运行".format(temp[4]))     
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "list.xls"),
                        os.path.join(self.output_dir, "list.xls"))
            os.link(os.path.join(self.work_dir, "sample.detail.xls"),
                        os.path.join(self.output_dir, "sample.detail.xls"))
            shutil.copytree(os.path.join(self.work_dir,"fasta_dir"),
                        os.path.join(self.output_dir, "fasta_dir"))
            shutil.copyfile(pictrue_path,
                        os.path.join(self.output_dir, new_pictrue_name))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))

                        

class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bac_datapre_' + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.bac_datapre",
            "options": {
                 "list_url":"http://ngs.tngs.com/xlsx/list.xls",
                 "sample_info_url":"http://ngs.tngs.com/xlsx/sample.detail.xls"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
