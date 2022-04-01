# -*- coding: utf-8 -*-
# __author__ = 'hongyu'
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import os
import datetime
import gevent
import re


class PtDailyDedupWorkflow(Workflow):
    """
    用于全库排查守护进程
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PtDailyDedupWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "back_time", "type": "int", "default": 7},          # 默认对七天内的M-S运行查重
            {"name": "err_min", "type": "int", "default": 50},
            {"name": "update_info", "type": "string"},
            {"name": "mem", "type": "int"}
        ]
        self.add_option(options)
        self.api_pt = self.api.api("medical.paternity_test_v3.paternity_test_v3")
        self.set_options(self._sheet.options())
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_data/"
        self.ref_all = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_all/"
        self.ms_ids = []
        self.ref_set = set(self.api_pt.get_ref_list())
        self.tab_exists = set([i.split(".")[0] for i in os.listdir(self.ref_all)])
        self.modules_dedup = []

    def check_options(self):
        """
        参数检查
        :return:
        """
        if not self.option("back_time"):
            raise OptionError("必须输入back_time参数")
        if not self.option("err_min"):
            raise OptionError("必须输入err_min参数")
        if not self.option("mem"):
            raise OptionError("必须输入mem参数")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def find_recent_ms(self):
        """
        根据self.option("back_time")定义的时间，查询近期查重的母本-胎儿
        :return:
        """
        time_now = datetime.datetime.now().timetuple()                         # 当前的时间元组
        time_past = datetime.datetime(time_now.tm_year, time_now.tm_mon, time_now.tm_mday) - datetime.timedelta(
            days=self.option("back_time"))
        time_past_ = time_past.strftime('%Y-%m-%d %H:%M:%S')                   # 格式："2018-04-24 00:00:00"
        self.logger.info("查找{}之后的M-S".format(time_past_))
        self.ms_ids = self.api_pt.find_ms_id(time_past_)
        if not self.ms_ids:
            self.logger.info("没有需要查重的M-S")
            gevent.spawn_later(5, self.end)

    def get_dad_no_analysis(self, mom_id, preg_id):
        """
        获取没有进行过查重分析的父本
        :return:
        """
        father_dedup_id = self.api_pt.find_dedup_err_id(5, mom_id, preg_id)
        if father_dedup_id:
            analysed = set(self.api_pt.find_dedup_father_id(father_dedup_id))
            not_analysed = self.ref_set - analysed
            dad_list = not_analysed & self.tab_exists
        else:
            dad_list = self.ref_set & self.tab_exists
            self.logger.info("库中未查找到该M-S的查重结果，后面将进行全库查重分析！")
        dad_list_ = [i for i in dad_list if "-M" in i or "-F" in i]                 # 只对F、M样品进行查重
        return dad_list_

    def dedup_run(self):
        """
        循环self.ms_ids，启动dedup
        :return:
        """
        n = 0
        for m_s_id in self.ms_ids:
            mom_id = m_s_id.split('_')[0]
            preg_id = "_".join(m_s_id.split('_')[1:])
            if not self.api_pt.tab_exist(mom_id):
                self.logger.info("{}tab文件不存在，跳过{}查重".format(mom_id, m_s_id))
                continue
            if (not self.api_pt.check_existence("sg_sample_tab_ms", {"sample_id": preg_id})) and (
                    not self.api_pt.tab_exist(preg_id)):
                self.logger.info("{}tab文件不存在，跳过{}查重".format(preg_id, m_s_id))
                continue
            dad_list = self.get_dad_no_analysis(mom_id, preg_id)
            if not dad_list:
                self.logger.info("查重父本列表为空，将跳过{}查重分析模块！".format(m_s_id))
                continue
            pt_daily_dedup = self.add_module("medical.paternity_test_v3.daily_dedup")
            self.step.add_steps("pt_dedup_{}".format(n))
            pt_daily_dedup.set_options({
                    "dad_id": "daily_dedup",     # 随便命名一个dad_id
                    "mom_id": mom_id,
                    "preg_id": preg_id,
                    "ref_point": self.ref_point,
                    "ref_data": self.ref_data,
                    "ref_all": self.ref_all,
                    "err_min": self.option("err_min"),
                    "mem": self.option("mem"),
                    "dad_list": ','.join(dad_list)
            })
            step = getattr(self.step, "pt_dedup_{}".format(n))
            step.start()
            pt_daily_dedup.on("end", self.finish_update, "pt_dedup_{}".format(n))
            n += 1
            self.modules_dedup.append(pt_daily_dedup)
        for j in range(len(self.modules_dedup)):
            self.modules_dedup[j].on('end', self.set_output, 'pt_daily_dedup')
        if self.modules_dedup:
            if len(self.modules_dedup) > 1:
                self.on_rely(self.modules_dedup, self.end)
            else:
                self.modules_dedup[0].on('end', self.end)
        else:
            self.logger.info("没有需要查重的M-S,正常退出")
            gevent.spawn_later(5, self.end)

        for i in self.modules_dedup:
            # gevent.sleep(0)
            i.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到指定目录
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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                next_dirname = dirname + '/' + os.path.basename(oldfiles[i])
                self.linkdir(oldfiles[i], next_dirname)

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "pt_daily_dedup":
            self.linkdir(obj.output_dir, self.output_dir + '/{}_{}'.format(str(obj.option('mom_id')),
                                                                           str(obj.option('preg_id'))))
            # for i in range(2, self.option('err_min')):                     # tool运行结果导表
            #     dir_path = obj.output_dir + '/pt_result_' + str(i)
            #     if not os.path.exists(dir_path):
            #         continue
            #     result_name = obj.option("dad_id") + "_" + obj.option("mom_id") + "_" + obj.option("preg_id")
            #     results = os.listdir(dir_path)
            #     for f in results:
            #         if f == result_name + '.txt':
            #             father_dedup_id = self.api_pt.find_dedup_err_id(i, obj.option("mom_id"),
            #                                                             obj.option("preg_id"))
            #             if not father_dedup_id:
            #                 break
            #             self.api_pt.import_dedup_data(dir_path + '/' + f, father_dedup_id)
            # father_ids = self.api_pt.search_father_id(obj.option("mom_id"), obj.option("preg_id"))  # 找出所有相关家系的主表ID
            # for father_id in father_ids:                                   # 更新批次排查结果
            #     for j in range(2, self.option('err_min')):
            #         self.api_pt.update_dedup_father(j, obj.option("mom_id"), obj.option("preg_id"), father_id)

    def run(self):
        self.find_recent_ms()
        self.dedup_run()
        super(PtDailyDedupWorkflow, self).run()

    def end(self):
        super(PtDailyDedupWorkflow, self).end()
