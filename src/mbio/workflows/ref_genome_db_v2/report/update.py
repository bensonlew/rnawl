# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import os
import subprocess
from biocluster.config import Config


class UpdateWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(UpdateWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="name", type="string", default=None),
            dict(name="genome_id", type="string", default=None),
            dict(name="result_dir", type="string", default=None),
            dict(name="task_status", type="string", default=None),
            dict(name="update_info", type='string'),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.update_mongo_script = self.config.PACKAGE_DIR + "/ref_genome_db_v2/update_db_to_sanger.py"
        self.update_disk_script = self.config.PACKAGE_DIR + "/ref_genome_db_v2/syn2nb2.sh"
        self.db = Config().get_mongo_client(mtype='ref_genome_db')[Config().get_mongo_dbname('ref_genome_db')]
        self.collection = self.db['sg_task']
        # self.tool = self.add_tool("ref_genome_db_v2.update")

    def run(self):
        # self.tool.on("end", self.set_db)
        # self.run_tool()
        self.start_listener()
        self.fire("start")
        self.stop_timeout_check()
        if os.path.exists(os.path.join(self.work_dir, "disk_copy_finish")):
            self.logger.info("已完成文件拷贝，跳过")
        else:
            self.update_disk()
        self.update_mongo()

    def update_disk(self):
        cmd = "{} {} {}".format(self.update_disk_script, self.option("result_dir"), self.option("task_status"))
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            os.system('touch \"disk_copy_finish\"')
            self.logger.info("文件上传完成")
        except subprocess.CalledProcessError:
            self.collection.update({'genome_id': self.option("genome_id")}, {'$set': {'task_status': self.option("task_status")}}, upsert=True)
            self.set_error("文件上传失败!!")

    def update_disk_bak(self):
        if self.option("task_status") == "sanger_finish":
            if os.path.exists(self.option("result_dir")):
                self.logger.info("sanger for sanger") # sanger运行工作流，sanger运行更新脚本
                target_dir = self.option("result_dir").replace("lustre", "ilustre").replace("sanger", "isanger")
                cmd = "scp -r -i ~/.ssh/id_rsa {} isanger@10.2.0.115:{}".format(
                    self.option("result_dir"), target_dir)
                self.logger.info(cmd)
                try:
                    code = os.system(cmd)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd))
                        self.set_error("向isanger服务器传递数据失败")
                except:
                    self.set_error("向isanger服务器传递数据失败")
                target_dir1 = self.option("result_dir").replace("sanger", "sanger-dev")
                cmd1 = "scp -r -i ~/.ssh/id_rsa {} sanger-dev@10.2.3.173:{}".format(
                    self.option("result_dir"), target_dir1)
                self.logger.info(cmd1)
                try:
                    code = os.system(cmd1)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd1))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd1))
                        self.set_error("向sanger-dev服务器传递数据失败")
                except:
                    self.set_error("向sanger-dev服务器传递数据失败")
            else:
                self.logger.info("isanger for sanger") # sanger运行工作流，isanger运行更新脚本
                target_dir = self.option("result_dir").replace("lustre", "ilustre").replace("sanger", "isanger")
                cmd = "scp -r -i ~/app/database/Genome_DB_finish/.key/nb_id_rsa sanger@10.2.0.110:{} {}".format(
                    self.option("result_dir"), target_dir)
                self.logger.info(cmd)
                try:
                    code = os.system(cmd)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd))
                        self.set_error("向isanger服务器传递数据失败")
                except:
                    self.set_error("向isanger服务器传递数据失败")
                target_dir1 = self.option("result_dir").replace("sanger", "sanger-dev")
                cmd1 = "scp -r -i ~/.ssh/id_rsa {} sanger-dev@10.2.3.173:{}".format(
                    target_dir, target_dir1)
                self.logger.info(cmd1)
                try:
                    code = os.system(cmd1)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd1))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd1))
                        self.set_error("向sanger-dev服务器传递数据失败")
                except:
                    self.set_error("向sanger-dev服务器传递数据失败")
        elif self.option("task_status") == "isanger_finish":
            if os.path.exists(self.option("result_dir")):
                self.logger.info("isanger for isanger") # isanger运行工作流，isanger运行更新脚本
                target_dir = self.option("result_dir").replace("ilustre", "lustre").replace("isanger", "sanger")
                cmd = "scp -r -i ~/.ssh/id_rsa {} sanger@10.2.0.110:{}".format(
                    self.option("result_dir"), target_dir)
                self.logger.info(cmd)
                try:
                    code = os.system(cmd)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd))
                        self.set_error("向sanger服务器传递数据失败")
                except:
                    self.set_error("向sanger服务器传递数据失败")
                target_dir1 = self.option("result_dir").replace("ilustre", "lustre").replace("isanger", "sanger-dev")
                cmd1 = "scp -r -i ~/.ssh/id_rsa {} sanger-dev@10.2.3.173:{}".format(
                    self.option("result_dir"), target_dir1)
                self.logger.info(cmd1)
                try:
                    code = os.system(cmd1)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd1))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd1))
                        self.set_error("向sanger-dev服务器传递数据失败")
                except:
                    self.set_error("向sanger-dev服务器传递数据失败")
            else:
                self.logger.info("sanger for isanger")  # isanger运行工作流，sanger运行更新脚本
                target_dir = self.option("result_dir").replace("ilustre", "lustre").replace("isanger", "sanger")
                cmd = "scp -r -i ~/app/database/Genome_DB_finish/.key/nb2_rsa isanger@10.2.0.115:{} {}".format(
                    self.option("result_dir"), target_dir)
                self.logger.info(cmd)
                try:
                    code = os.system(cmd)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd))
                        self.set_error("向sanger服务器传递数据失败")
                except:
                    self.set_error("向sanger服务器传递数据失败")
                target_dir1 = self.option("result_dir").replace("ilustre", "lustre").replace("isanger", "sanger-dev")
                cmd1 = "scp -r -i ~/.ssh/id_rsa {} sanger-dev@10.2.3.173:{}".format(
                    target_dir, target_dir1)
                self.logger.info(cmd1)
                try:
                    code = os.system(cmd1)
                    if code == 0:
                        self.logger.info("命令{}执行成功！".format(cmd1))
                    else:
                        self.logger.info("命令{}执行失败！".format(cmd1))
                        self.set_error("向sanger-dev服务器传递数据失败")
                except:
                    self.set_error("向sanger-dev服务器传递数据失败")
        else:
            self.set_error("暂不支持该状态更新")

    def update_mongo(self):
        cmd = "{} {} --genome_ids {}".format(self.config.SOFTWARE_DIR + '/program/Python/bin/python', self.update_mongo_script, self.option("genome_id"))
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("数据库更新完成")
        except subprocess.CalledProcessError:
            self.collection.update({'genome_id': self.option("genome_id")}, {'$set': {'task_status': self.option("task_status")}}, upsert=True)
            self.set_error("数据库更新失败")
        try:
            self.collection.update({'genome_id': self.option("genome_id")}, {'$set': {'task_status': 'sanger_update'}}, upsert=True)
            self.logger.info("sg_task主表task_status状态更新完成")
        except:
            self.set_error("sg_task主表task_status状态更新失败")
        self.end()

    def run_tool(self):
        options = dict(
            genome_id=self.option('genome_id'),
            result_dir=self.option('result_dir'),
        )
        self.tool.set_options(options)
        self.tool.run()

    def end(self):
        super(UpdateWorkflow, self).end()