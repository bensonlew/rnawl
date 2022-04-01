# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer
import json
import os
import re


class ErrSiteFiterHomohybridWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        用途及运行逻辑：
        用户过滤错配位点的信息，通过超找错配位点，然后过滤掉tab中该位点,重新画图，出报告，该程序适用与胎儿过滤标准是杂合与纯合位点
        __auther__: hongdong
        last modified by hongdong@20191111
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(ErrSiteFiterHomohybridWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "dad_id", "type": "string"},
            {"name": "mom_id", "type": "string"},
            {"name": "son_id", "type": "string"},
            {"name": "update_info", "type": "string"},  # 用于后面更新主表
            {"name": "member_id", "type": "string"},  # 会员ID
            {"name": "father_err_id_old", "type": "string"},  # 旧的错配位点主表id
            {"name": "err_min", 'type': "int"},
            {"name": "father_id", "type": "string"},
            {"name": "old_dad", "type": "string"},
            {"name": "old_mom", "type": "string"},
            {"name": "old_son", "type": "string"},
            # 自由交互发起的过滤分析需要把dp_correction状态传过来，控制导出的tab文件样品名
            {"name": "dp_correction", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.dp_correction = None
        self.family_analysis = self.add_module("medical.paternity_test_v3.family_analysis_homohybrid")
        self.path = None  # 用于存储报告中的图片
        self.origin_son = None
        self.api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        self.father_err_id_new = json.loads(self.option('update_info')).keys()[0]
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"  # 每个位点的基因型概率
        self.is_free = False
        self.bucket = self.config.get_project_region_bucket(project_type="pt_v3").rstrip("/")

    def check_options(self):
        """检查参数设置"""
        if not self.option('dad_id'):
            raise OptionError('缺少参数dad_id')
        if not self.option('mom_id'):
            raise OptionError('缺少参数mom_id')
        if not self.option('son_id'):
            raise OptionError('缺少参数son_id')
        if not self.option('father_err_id_old'):
            raise OptionError('缺少参数father_err_id_old')
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def family_analysis_run(self):
        """家系分析生成父权制与调试表"""
        self.family_analysis.set_options({
            "dad_tab": self.output_dir + '/' + self.option("dad_id") + '.tab',
            "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
            "preg_tab": self.output_dir + '/' + self.option("son_id") + '.tab',
            "ref_point": self.ref_point,
            "err_min": self.option("err_min"),
            "is_free_combine": self.is_free
        })
        self.family_analysis.on('end', self.set_output, 'pt_analysis')
        self.family_analysis.run()

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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "pt_analysis":
            self.linkdir(obj.output_dir + '/family_analysis', self.output_dir)
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir)
            self.linkdir(obj.output_dir + '/result_info', self.output_dir)
            self.end()

    def set_type(self):
        """
        设置分析类型，是自由交互的还是普通家系的
        :return:
        """
        if self.option('old_dad'):
            self.is_free = True

    def run(self):
        """
        运行逻辑：
        首先检测该家系的所有的2-8错配的结果是不是都是no，如果都是no说明不是真父，故不能进行过滤分析
        （因为过滤后可能会将假父变为真父），这里的卡控任然不能完全避免上述问题，所以会在页面上进行过滤功能的限定，
        只能是最高权限的人使用，并且改报告人员，明确知道这个父本是真父的前提下才能使用。
        :return:
        """
        self.set_type()
        max_err = self.api_pt.find_err_max(self.option('father_id'), postype='homohybrid')
        if not self.api_pt.is_father(self.option('father_id'), max_err, postype='homohybrid'):
            raise Exception("父本与胎儿的错配结果检测均是不匹配，不能进行过滤分析！")
        err_postions = self.api_pt.err_postions_list(self.option("father_err_id_old"), postype='homohybrid')
        self.logger.info(err_postions)
        if len(err_postions) == 0:
            raise Exception("查询得到的错配位点列表为空，流程终止--可以忽略！")
        if not self.is_free:
            self.api_pt.attempt_to_export_tab_file(self.option("dad_id"), self.output_dir)
            self.api_pt.attempt_to_export_tab_file(self.option("mom_id"), self.output_dir)
            self.api_pt.attempt_to_export_tab_file(self.option("son_id"), self.output_dir, None, err_postions)
        else:
            self.api_pt.attempt_to_export_tab_file(self.option("old_dad"), self.output_dir, self.option("dad_id"))
            self.api_pt.attempt_to_export_tab_file(self.option("old_mom"), self.output_dir, self.option("mom_id"))
            if self.option("dp_correction") == "true":
                dp_correction_new_id = self.option("old_son") + "_" + self.option("old_mom")
                self.api_pt.attempt_to_export_tab_file(dp_correction_new_id, self.output_dir, self.option("son_id"),
                                                       err_postions)
            else:
                self.api_pt.attempt_to_export_tab_file(self.option("old_son"), self.output_dir, self.option("son_id"),
                                                       err_postions)

        self.family_analysis_run()
        super(ErrSiteFiterHomohybridWorkflow, self).run()

    def set_db(self):
        """结果文件批量入库"""
        api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        results = os.listdir(self.output_dir)
        self.logger.info(results)
        result_name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.option("son_id")
        for f in results:
            if f == result_name + "_family_joined_tab.txt":
                api_pt.add_sg_pt_father_detail(self.output_dir + '/' + f, self.father_err_id_new)
            elif f == self.option("mom_id") + "_" + self.option("son_id") + '_info_show.txt':
                api_pt.add_info_detail(self.output_dir + '/' + f, self.father_err_id_new)
            elif f == result_name + "_family_match_update.png":
                if self.config.RGW_ENABLE:
                    pic_path = "{}/pt_result_pic/".format(self.bucket) + self.father_err_id_new + '/output/'
                else:
                    pic_path = "rerewrweset/MEDfiles/pt_result_pic/" + self.father_err_id_new + '/'
                api_pt.add_result_pic(pic_path, result_name, self.option("father_id"), self.father_err_id_new,
                                      postype='homohybrid')
            elif f == result_name + '_test_pos_update.txt':
                api_pt.add_test_pos(self.output_dir + '/' + f, self.father_err_id_new, postype='homohybrid')
            elif f == result_name + "_family_analysis.txt":
                api_pt.add_analysis_tab(self.output_dir + '/' + f, self.father_err_id_new, self.option("father_id"),
                                        int(self.option("err_min")), postype='homohybrid')
            elif f == result_name + "_num_pos_map.txt":
                api_pt.add_num_pos_map(self.output_dir + '/' + f, self.father_err_id_new, self.option("father_id"))

    def end(self):
        self.set_db()
        self.upload_target_path()
        super(ErrSiteFiterHomohybridWorkflow, self).end()

    def upload_target_path(self):
        """
        将结果link到指定磁盘路径下
        :return:
        """
        self.logger.info("开始设置结果目录！")
        if self.config.RGW_ENABLE:
            transfer = MultiFileTransfer()
            transfer.add_upload(self.output_dir + "/", "{}/pt_result_pic/{}/".format(self.bucket,
                                                                                     self.father_err_id_new),
                                base_path=os.path.dirname(self.output_dir))
            transfer.perform()
        else:
            if self._sheet.client == 'client01':
                self.path = '/mnt/ilustre/data/rerewrweset/MEDfiles/pt_result_pic/{}'.format(self.father_err_id_new)
            else:
                self.path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/pt_result_pic/{}'\
                    .format(self.father_err_id_new)
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            self.linkdir(self.output_dir, self.path)
        self.logger.info("设置结果目录成功！")