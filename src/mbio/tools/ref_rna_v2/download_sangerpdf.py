# -*- coding: utf-8 -*-
# __author__ = 'xuxi'
# last_modifiy:2021.01.26

import os,re,shutil,glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import urllib
import json
import subprocess
import datetime

class DownloadSangerpdfAgent(Agent):
    """
    docker内的sanger_pdf.py的help文档如下：
        usage: sanger_pdf.py [-h] -l LOGIN [-m MODE] url

        this script was used to download pdf from interaction page in http://www.sanger.com

        positional arguments:
          url                   a url link of sanger website. Example: http://report.sanger.com/drna/tf_family/task_id/sanger_281030.html

        optional arguments:
          -h, --help            show this help message and exit
          -l LOGIN, --login LOGIN
                                the username and password to login in www.sanger.com 
                                This argument is NECESSARY REQUIRED!!!
                                !!! Note:Username and Password must be separated by "____username_password____".
                                Example: myname____username_password____mypassword
          -m MODE, --mode MODE  [single|all|part[n:m]]
                                the mode of downloading pdf. 
                                ("single") mode will only download pdf from appointed url of one webpage. 
                                ("all") mode will download pdf from all webpages of whole project which related to appointed url webpage.
                                ("part[n:m]") mode will download part pdfs from all webpages of whole project which related to appointed url webpage.
                                -----------------
                                assume links are ['u1','u2','u3','u4','u5','u6','u7']
                                example_index_1 = [0:4]   ---> ['u1', 'u2', 'u3', 'u4']
                                example_index_2 = [2:6]   ---> ['u3', 'u4', 'u5', 'u6']
                                example_index_3 = [3:]    ---> ['u4', 'u5', 'u6', 'u7']  
                                For axample : if you want download fisrt 4 links, please input 'part[0:4]'
                                -----------------
                                
                                Default: single
    """
    def __init__(self, parent):
        super(DownloadSangerpdfAgent, self).__init__(parent)
        options = [
            {"name": "url", "type": "string", 'default': ""},
            {"name": "userpass", "type": "string", 'default': ""},
            {"name": "mode", "type": "string", 'default': "single"},
            {"name": "imageversion", "type": "string", 'default': "latest"}, 
        ]
        self.add_option(options)
        self.step.add_steps('download_sangerpdf')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        
    def step_start(self):
        self.step.download_sangerpdf.start()
        self.step.update()

    def step_end(self):
        self.step.download_sangerpdf.finish()
        self.step.update()    
        
    def check_options(self):
        if not self.option("url"):
            raise OptionError("请传入下载sanger交互页面pdf的url ！")
        if not self.option("userpass"):
            raise OptionError("请传入下载sanger交互页面pdf所需的账户和密码 ！")

    def set_resource(self):
        self._cpu = 4
        self._memory = '10G'
        
    def end(self):
        super(DownloadSangerpdfAgent, self).end()


class DownloadSangerpdfTool(Tool):

    def __init__(self, config):
        super(DownloadSangerpdfTool, self).__init__(config)

    def check_image_local_exit_or_pull(self, image_version_):
        inspect_cmd = 'docker inspect hub.majorbio.com/rnagroup/puppeteer_sanger:{}'.format(image_version_)
        image_version_inspect = subprocess.Popen(inspect_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        image_version_status = image_version_inspect.communicate()
        image_version_inspect.wait()
        if "Error: No such object" in image_version_status[0]:
            self.logger.info("本计算节点不存在{}版本的docker镜像，将尝试从库中拉取".format(image_version_))
            pull_cmd = 'docker pull hub.majorbio.com/rnagroup/puppeteer_sanger:{}'.format(image_version_)
            image_version_pull = subprocess.Popen(pull_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            image_version_status_ofpull = image_version_pull.communicate()
            image_version_pull.wait()
            if "Downloaded newer image for" in image_version_status_ofpull[0]:
                self.logger.info("成功拉取了版本为{}的hub.majorbio.com/rnagroup/puppeteer_sanger镜像".format(image_version_))
            else:
                raise OptionError("拉取版本为{}的hub.majorbio.com/rnagroup/puppeteer_sanger镜像失败,拉取镜像的返回值为{}".format(image_version_, image_version_status_ofpull), code="33707407")
        elif "RepoTags" in image_version_status[0]:
            self.logger.info("本计算节点已存在{}版本的docker镜像,将不会从库中拉取".format(image_version_))
        else:
            raise OptionError("无法判断{}版本的docker镜像是否存在！".format(image_version_))

    def prepare_image(self):
        url = 'https://hub.majorbio.com/api/v2.0/projects/rnagroup/repositories/puppeteer_sanger/artifacts?with_tag=true&with_scan_overview=true&with_label=true&page_size=15&page=1'
        data = urllib.urlopen(url).read()
        json_datas = json.loads(data)
        version_createtime_and_name = {}
        for json_data in json_datas:
            version_createtime_and_name[json_data['tags'][0]['push_time'].split('.')[0]] = json_data['tags'][0]['name']
        
        if self.option("imageversion") == "latest":
            # 获取最新的版本号，检查该版本在本地是否存在 若不存在的话就从远程库里拉取
            sorted_createtime = sorted((datetime.datetime.strptime(d, "%Y-%m-%dT%H:%M:%S") for d in version_createtime_and_name.keys()), reverse=True)
            latest_version = version_createtime_and_name[sorted_createtime[0].strftime("%Y-%m-%dT%H:%M:%S")]
            self.logger.info("找到远程仓库中最新的docker镜像版本号为：{}".format(latest_version))
            self.check_image_local_exit_or_pull(latest_version)
            return latest_version
        else:
            # 检查该版本是否在远程库里是否存在，存在的话，再检查在本地是否存在 若不存在就从远程库里拉取
            remotehub_all_versions = version_createtime_and_name.values
            if self.option("imageversion") in remotehub_all_versions:
                self.check_image_local_exit_or_pull(self.option("imageversion"))
                return self.option("imageversion")
            else:
                raise OptionError("该版本的docker镜像在仓库里都不存在，可能该版本已被删除或者根本不存在！！")

    def downloadpdf(self):
        current_output_dir = self.output_dir
        path_list = current_output_dir.split('/')
        current_user =  path_list[int(path_list.index('users')+1)]
        if current_user in ['sanger','isanger']:
            good_current_user = current_user
        else:
            raise OptionError("判断当前用户失败", code="33707406")

        cmd = "/usr/bin/docker run --user {5} --memory 10g --oom-kill-disable --shm-size=8gb --cpus=4 --rm -i --privileged \
        -v {0}:/home/{5}/output \
        -w /home/{5}/output hub.majorbio.com/rnagroup/puppeteer_sanger:{1} \
        python3 -u /app/sanger_puppeteer/sanger_pdf.py \
        -l {2} -m {3} {4}".format(self.output_dir, self.image_version, self.option("userpass"), self.option("mode"), self.option("url"), good_current_user)
        self.logger.info("开始运行docker容器获取pdf")
        command = self.add_command("downloadpdf", cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("docker容器获取pdf完成!")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("重运行docker容器获取pdf完成!")
            else:
                self.set_error("重运行docker容器获取pdf 依然出错！")

    def run(self):
        """
        运行
        """
        super(DownloadSangerpdfTool, self).run()
       
        self.logger.info("开始运行download_sangerpdf")
        if self.option("url") and self.option("userpass"):
            self.image_version = self.prepare_image()
            self.downloadpdf()

        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet             
        data = {
            "id": "DownloadSangerpdf_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.download_sangerpdf",
            "instant": False,
            "options": dict(
                url="http://report.sanger.com/wholerna/expcorr/task_id/sanger_316281.html",
                userpass="sgtest@majorbio.com____username_password____test01",
                mode="single",
                imageversion="latest",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        # data['options']['method'] = 'edgeR'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'DESeq2'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '2'
        # data['options']['method'] = 'DEGseq'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()

if __name__ == '__main__':
    unittest.main()