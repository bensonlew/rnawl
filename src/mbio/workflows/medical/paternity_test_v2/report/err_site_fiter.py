# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import json
import os
import re


class ErrSiteFiterWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        用途及运行逻辑：
        用户过滤错配位点的信息，通过超找错配位点，然后过滤掉tab中该位点,重新画图，出报告，该过程不需要call snp与查重
        __auther__: hongdong
        last modified by hongdong@20171204
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(ErrSiteFiterWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "dad_id", "type": "string"},
            {"name": "mom_id", "type": "string"},
            {"name": "son_id", "type": "string"},
            {"name": "update_info", "type": "string"},  # 用于后面更新主表
            {"name": "member_id", "type": "string"},  # 会员ID
            {"name": "father_err_id_old", "type": "string"},  # 旧的错配位点主表id
            {"name": "err_min", 'type': "int"},
            {"name": "father_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.family_analysis = self.add_module("medical.paternity_test_v2.family_analysis")
        self.result_info = self.add_tool("medical.paternity_test_v2.result_info")
        self.path = None  # 用于存储报告中的图片
        self.api_pt = self.api.api('medical.paternity_test_v2')
        self.father_err_id_new = json.loads(self.option('update_info')).keys()[0]
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"  # 每个位点的基因型概率

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
            "err_min": self.option("err_min")
        })
        self.family_analysis.on('end', self.set_output, 'pt_analysis')
        self.family_analysis.on('end', self.result_info_run)
        self.family_analysis.run()

    def result_info_run(self):
        """画出相关结果图SNP分型图等"""
        rdata = self.option("dad_id") + '_' + self.option("mom_id") + '_' + self.option('son_id') + "_family_" \
                                                                                                    "joined_tab.Rdata"
        self.result_info.set_options({
            "tab_merged": os.path.join(self.family_analysis.output_dir + '/family_merge', rdata)
        })
        self.result_info.on("end", self.set_output, 'result_info')
        self.result_info.on('end', self.end)
        self.result_info.run()

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
        if event['data'] == "result_info":
            if self._sheet.client == 'client01':
                self.path = '/mnt/ilustre/data/rerewrweset/MEDfiles/pt_result_pic/' + self.father_err_id_new
            else:
                self.path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/pt_result_pic/' + self.father_err_id_new
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            self.linkdir(obj.output_dir, self.output_dir)
            self.linkdir(obj.output_dir, self.path + '/')

    def run(self):
        """
        运行逻辑：
        首先检测该家系的所有的2-8错配的结果是不是都是no，如果都是no说明不是真父，故不能进行过滤分析
        （因为过滤后可能会将假父变为真父），这里的卡控任然不能完全避免上述问题，所以会在页面上进行过滤功能的限定，
        只能是最高权限的人使用，并且改报告人员，明确知道这个父本是真父的前提下才能使用。
        :return:
        """
        if not self.api_pt.is_father(self.option('father_id'), 7):
            raise Exception("父本与胎儿2-8的错配结果检测均是不匹配，不能进行过滤分析！")
        err_postions = self.api_pt.err_postions_list(self.option("father_err_id_old"))
        self.logger.info(err_postions)
        if len(err_postions) == 0:
            raise Exception("查询得到的错配位点列表为空，流程终止--可以忽略！")
        if self.api_pt.tab_exist(self.option("dad_id")):
            self.api_pt.export_tab_file(self.option("dad_id"), self.output_dir)
        else:
            raise Exception("样本{}的tab文件不存在！".format(self.option("dad_id")))
        if self.api_pt.tab_exist(self.option("mom_id")):
            self.api_pt.export_tab_file(self.option("mom_id"), self.output_dir)
        else:
            raise Exception("样本{}的tab文件不存在！".format(self.option("mom_id")))
        if self.api_pt.tab_exist(self.option("son_id")):
            self.api_pt.export_tab_file_select(self.option("son_id"), self.output_dir, err_postions)
        else:
            raise Exception("样本{}的tab文件不存在！".format(self.option("son_id")))
        self.family_analysis_run()
        super(ErrSiteFiterWorkflow, self).run()

    def set_db(self):
        """结果文件批量入库"""
        api_pt = self.api.api('medical.paternity_test_v2')
        results = os.listdir(self.output_dir)
        self.logger.info(results)
        result_name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.option("son_id")
        for f in results:
            if f == result_name + "_family_joined_tab.txt":
                api_pt.add_sg_pt_father_detail(self.output_dir + '/' + f, self.father_err_id_new)
            elif f == self.option("mom_id") + "_" + self.option("son_id") + '_info_show.txt':
                api_pt.add_info_detail(self.output_dir + '/' + f, self.father_err_id_new)
            elif f == result_name + "_family.png":
                api_pt.add_result_pic("pt_result_pic/" + self.father_err_id_new + '/', result_name,
                                      self.option("father_id"), self.father_err_id_new)
            elif f == result_name + '_test_pos.txt':
                api_pt.add_test_pos(self.output_dir + '/' + f, self.father_err_id_new)
            elif f == result_name + "_family_analysis.txt":
                api_pt.add_analysis_tab(self.output_dir + '/' + f, self.father_err_id_new, self.option("father_id"),
                                        int(self.option("err_min")))

    def end(self):
        self.set_db()
        super(ErrSiteFiterWorkflow, self).end()
