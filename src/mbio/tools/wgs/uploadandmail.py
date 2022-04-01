# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last modify 20190603
import os
import json
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.dna_evolution.send_email import SendEmail
from biocluster.api.file.lib.transfer import MultiFileTransfer


class UploadandmailAgent(Agent):
    """
    用于将基因组配置文件进行上传到指定位置，然后通知产品线去检查是否合格
    """
    def __init__(self, parent):
        super(UploadandmailAgent, self).__init__(parent)
        options = [
            {"name": "file_dir", "type": "string"},
            {"name": "mail_info", "type": "string"},
            {"name": "target_dir", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('GOanno')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.GOanno.start()
        self.step.update()

    def step_end(self):
        self.step.GOanno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("file_dir"):
            raise OptionError("请设置file_dir参数")
        if not self.option("mail_info"):
            raise OptionError("请设置mail_info参数")
        if not self.option("target_dir"):
            raise OptionError('请设置target_dir参数')

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(UploadandmailAgent, self).end()


class UploadandmailTool(Tool):
    def __init__(self, config):
        super(UploadandmailTool, self).__init__(config)

    def sendmail(self):
        """
        发送邮件mail_info = {'project_id': 4567, 'email_id': "hongdong.xuan@majorbio.com"}
        :return:
        """
        task_id = 'none'
        project_id = 'none'
        email_id = 'hongdong.xuan@majorbio.com'
        info = json.loads(self.option('mail_info'))
        if 'task_id' in info.keys():
            task_id = info['task_id']
        if "project_id" in info.keys():
            project_id = info['project_id']
        if 'email_id' in info.keys():
            email_id = info['email_id']
        a = SendEmail("1274095891@qq.com", "smtp.qq.com", "ocfnjnhnsbqubaej", "1274095891@qq.com", email_id,
                      "任务id:{}的基因组配置信息-请及时检查！".format(task_id), 465)
        a.send_msg("{}".format("http://www.{}.com/task/project_tasks/project_id/{}.html"
                               .format(task_id.split('_')[0], project_id)))
        # a.attachment(self.mapping_stat.output_dir + "/result.stat/Total.mapped.detail.xls")
        a.send_email()
        self.logger.info("邮件发送完成")

    def uploadfile(self):
        """
        上传文件ref.gff文件+ref.new.mRNA.fa文件
        :return:
        """
        gff = self.option("file_dir").rstrip('/') + 'ref.gff'
        mrna = self.option("file_dir").rstrip('/') + 'ref.new.mRNA.fa'
        target = os.path.dirname(self.option('target_dir')).rstrip('/') + '/'
        transfer = MultiFileTransfer()
        if os.path.exists(gff):
            transfer.add_upload(gff, target, base_path=os.path.dirname(gff))
        if os.path.exists(mrna):
            transfer.add_upload(mrna, target, base_path=os.path.dirname(mrna))
        transfer.perform()
        self.logger.info("文件上传成功！")

    def run(self):
        super(UploadandmailTool, self).run()
        self.uploadfile()
        self.sendmail()
        self.end()
