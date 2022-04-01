# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.03.01

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.bsa.bsa_copy_demo import BsaCopyDemo
from mbio.packages.dna_gmap.gmap_copy_demo import GmapCopyDemo
from mbio.packages.dna_evolution.evolution_copy_demo import EvolutionCopyDemo
from mbio.packages.noref_wgs.noref_wgs_copy_demo import NorefWgsCopyDemo
from mbio.packages.wgs_v2.wgs_v2_copy_demo import WgsV2CopyDemo



class BsaCopyDemoWorkflow(Workflow):
    """
    拉取BSA demo数据
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BsaCopyDemoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": 'string', "default": ''},
            {"name": "target_task_id", "type": 'string', "default": ''},
            {"name": "target_project_sn", "type": 'string', "default": ''},
            {"name": "target_member_id", "type": 'string', "default": ''},
            {"name": "project_type", "type": "string", "default": "dna_bsa"}  # 项目类型，dna_bsa/dna_gmap
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("task_id"):
            raise OptionError("缺少demo的task_id,请检查")
        if not self.option("target_task_id"):
            raise OptionError("缺少demo的新task_id,请检查")
        if not self.option("target_project_sn"):
            raise OptionError("缺少demo的新project_sn,请检查")
        if not self.option("target_member_id"):
            raise OptionError("缺少demo的新member_id,请检查")
        print self.option("project_type")

    def run(self):
        self.start_listener()
        self.fire("start")
        if self.option("project_type") == "dna_bsa":
            print "是我bsa"
            copy_task = BsaCopyDemo(old_task_id=self.option("task_id"), new_task_id=self.option("target_task_id"),
                                    new_project_sn=self.option("target_project_sn"), new_member_id=self.option("target_member_id"))
        elif self.option("project_type") == "dna_evolution":
            copy_task = EvolutionCopyDemo(old_task_id=self.option("task_id"), new_task_id=self.option("target_task_id"),
                                    new_project_sn=self.option("target_project_sn"), new_member_id=self.option("target_member_id"))
            print "应该是我进化"
        elif self.option("project_type") == "dna_gmap":
            copy_task = GmapCopyDemo(old_task_id=self.option("task_id"), new_task_id=self.option("target_task_id"),
                                     new_project_sn=self.option("target_project_sn"), new_member_id=self.option("target_member_id"))
        elif self.option("project_type") == "noref_wgs":
            copy_task = NorefWgsCopyDemo(old_task_id=self.option("task_id"), new_task_id=self.option("target_task_id"),
                                     new_project_sn=self.option("target_project_sn"), new_member_id=self.option("target_member_id"))
            print "我是无参变异检测的demo"
        elif self.option("project_type") == "dna_wgs_v2":
            copy_task = WgsV2CopyDemo(old_task_id=self.option("task_id"), new_task_id=self.option("target_task_id"),
                                     new_project_sn=self.option("target_project_sn"),
                                     new_member_id=self.option("target_member_id"))
            print "我是有参变异检测v2的demo"
        copy_task.run()
        self.end()

    def end(self):
        super(BsaCopyDemoWorkflow, self).end()
