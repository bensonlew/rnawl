# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import os
import re


class PtFreeCombineWorkflow(Workflow):
    """
    last_modify 20171208 by hongdongxuan
    运行逻辑:
    1)is_report字段为false时：
    此时是不需要进行出报告的，只是进行模糊匹配的父本的查错，此时也要运行家系分析，然后运行差错的模块，结果命名是dedup开始的记录，
    计算结果也是算2-10。
    2)is_report字段为true时:
    此时是要进行出报告的，这时必须的条件是有新父本编号，新母本编号与新胎儿编号，三者缺一不可，该条件是会进行家系分析，
    画图以及查重分析，运行的结果是2-10，结果命名是report开始的记录。
    3）导表的逻辑是，生成一个家系主表，然后关联所有的不同的错配的结果记录表，主表的名字是sg_pt_free_combine，结果细节表的话存
    在原有的表格中。
    """
    def __init__(self, wsheet_object):
        """
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(PtFreeCombineWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "dad_id", "type": "string"},  # 父本过滤器
            {"name": "mom_id", "type": "string"},  # 母本编号
            {"name": "preg_id", "type": "string"},  # 胎儿编号
            {"name": "err_min", "type": "int", "default": 2},  # 允许错配数
            {"name": "case_id", "type": "string"},   # 父本组
            {"name": "update_info", "type": "string"},
            {"name": 'new_mom_id', "type": "string"},  # 新母本编号
            {"name": 'new_dad_id', "type": "string"},  # 新父本编号
            {"name": 'new_preg_id', "type": "string"},  # 新胎儿编号
            {"name": 'father_id', "type": "string"},   # 主表id，接口中生成的
            {"name": 'member_id', "type": "string"},  # 用户ID
            {"name": "is_report", "type": "string"}
        ]
        self.add_option(options)
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_data/"
        self.ref_all = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_all/"
        self.api_pt = self.api.api("medical.paternity_test_v2")
        self.set_options(self._sheet.options())
        self.modules_family_analysis = []
        self.father_err_id = None
        self.tools_result = []
        self.tools_dedup = []
        self.dad_list = []
        self.is_dedup = False
        self.family = []
        self.rdata = None
        self.path = None
        self.num = None

    def check_options(self):
        """
        参数检查
        :return:
        """
        if not self.option("is_report"):
            raise OptionError("必须输入is_report参数")
        if not self.option("case_id"):
            raise OptionError("必须输入父本组编号")
        if not self.option("mom_id"):
            raise OptionError("必须输入母本编号")
        if not self.option("preg_id"):
            raise OptionError("必须输入胎儿编号")
        if not self.option('err_min'):
            raise OptionError("必须输入err_min")
        if self.option("is_report") == "true":
            if not(self.option('new_dad_id') and self.option('new_mom_id') and self.option('new_preg_id')):
                raise OptionError("is_report为true, 是要进行出报告, 必须要输入新的父本, 母本,胎儿的编号, "
                                  "否则后面无法进行家系分析！")
            elif not self.option('dad_id'):
                raise OptionError("is_report为true，必须要有父本过滤器这一选项！")
        return True

    def family_analysis_run(self):
        for i in range(2, self.num):
            pt_analysis = self.add_module("medical.paternity_test_v2.family_analysis")
            pt_analysis.set_options({
                "dad_tab": self.output_dir + '/' + self.family[0] + '.tab',
                "mom_tab": self.output_dir + '/' + self.family[1] + '.tab',
                "preg_tab": self.output_dir + '/' + self.family[2] + '.tab',
                "ref_point": self.ref_point,
                "err_min": i
            })
            self.modules_family_analysis.append(pt_analysis)
        for j in range(len(self.modules_family_analysis)):
            self.modules_family_analysis[j].on('end', self.set_output, 'pt_analysis')
        if self.modules_family_analysis:
            if len(self.modules_family_analysis) > 1:
                self.on_rely(self.modules_family_analysis, self.result_info_run)
            elif len(self.modules_family_analysis) == 1:
                self.modules_family_analysis[0].on('end', self.result_info_run)
        else:
            raise Exception("modules_family_analysis列表为空！")
        for t in self.modules_family_analysis:
            t.run()

    def result_info_run(self):
        for l in range(2, self.num):
            result_dir = os.path.join(self.output_dir, 'pt_result_' + str(l))
            results = os.listdir(result_dir)
            for f in results:
                if re.match(r'.*family_joined_tab\.Rdata$', f):
                    result_info = self.add_tool("medical.paternity_test_v2.result_info")
                    self.rdata = os.path.join(result_dir, f)
                    result_info.set_options({
                        "tab_merged": self.rdata,
                        "is_free_combine": True
                    })
                    self.tools_result.append(result_info)
        for j in range(len(self.tools_result)):
            self.tools_result[j].on('end', self.set_output, 'result_info')
        if self.tools_result:
            if len(self.tools_result) > 1:
                self.on_rely(self.tools_result, self.dedup_run)
            elif len(self.tools_result) == 1:
                self.tools_result[0].on('end', self.dedup_run)
        else:
            raise Exception("tools_result列表为空！")
        for t in self.tools_result:
            t.run()

    def dedup_run(self):
        """
        查重部分新机制：如果是全库查重的话，直接使用参考库就可以了。
        如果是部分查重的话，更据模糊匹配的结果中的父本进行查重
        :return:
        """
        for p in range(2, self.option('err_min') + 1):
            pt_analysis_dedup = self.add_tool("medical.paternity_test_v2.dedup_v2")
            pt_analysis_dedup.set_options({
                "dad_id": self.family[0],
                "mom_tab": self.output_dir + '/' + self.family[1] + '.tab',
                "preg_tab": self.output_dir + '/' + self.family[2] + '.tab',
                "ref_point": self.ref_point,
                "err_min": p,
                "dad_list": ','.join(self.dad_list),
                "father_path": self.ref_all if self.is_dedup else self.ref_data
            })
            self.tools_dedup.append(pt_analysis_dedup)
        for j in range(len(self.tools_dedup)):
            self.tools_dedup[j].on('end', self.set_output, 'dedup')
        if self.tools_dedup:
            if len(self.tools_dedup) > 1:
                self.on_rely(self.tools_dedup, self.end)
            else:
                self.tools_dedup[0].on("end", self.end)
        else:
            raise Exception('tools_dedup列表为空！')
        for tool in self.tools_dedup:
            tool.run()

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
            self.linkdir(obj.output_dir + '/family_analysis', self.output_dir + '/pt_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir + '/pt_result_' + str(obj.option('err_min')))
        if event['data'] == "result_info":
            dir_path = os.path.dirname(obj.option("tab_merged").prop['path'])
            if self._sheet.client == 'client01':
                self.path = '/mnt/ilustre/data/rerewrweset/MEDfiles/pt_result_pic/' + self.option("father_id")
            else:
                self.path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/pt_result_pic/' + self.option("father_id")
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            self.linkdir(obj.output_dir, dir_path)
            target_dir = os.path.basename(dir_path.rstrip())
            self.linkdir(obj.output_dir, self.path + '/' + target_dir)
        if event['data'] == "dedup":
            self.linkdir(obj.output_dir, self.output_dir + '/pt_result_' + str(obj.option('err_min')))

    def get_ref_list(self):
        """
        用于获取查重参考库中父本的所有的样本id
        :return:
        """
        ref_list = []
        ref_dad = os.listdir(self.ref_data)
        for m in ref_dad:
            ref_list.append(m.strip().split('.')[0])
        return ref_list

    def check_tab(self, dad_list):
        """
        用于检测，模糊匹配的结果是不是在参考库中存在，如果有不存在的话，要手动添加进去，否则终止流程
        :return:
        """
        no_tab = []
        for m in dad_list:
            if not os.path.exists(os.path.join(self.ref_all, "{}.tab".format(m))):
                no_tab.append(m)
        if no_tab:
            raise Exception("{}：上述样本的tab文件不在tab_all这个库中，请添加！".format(no_tab))

    def run(self):
        """
        因为该分析是自由交互，样本编号之间差异比较大，可能不是同一case下的样本也会进行组合进行分析，所以化繁为简，
        直接放在self.family这个列表中，避免了一系列的判断
        在父本查重的地方，
        :return:
        """
        if self.option('is_report') == 'true':
            # self.num = self.option("err_min") + 1
            self.num = 10  # 进行出报告的时候，所有错配都做
            self.logger.info("开始进行家系分析出报告！")
            self.family.append(self.option('new_dad_id'))
            self.family.append(self.option('new_mom_id'))
            self.family.append(self.option('new_preg_id'))
            self.api_pt.export_tab_file(self.option('mom_id'), self.output_dir, self.option('new_mom_id'))
            self.api_pt.export_tab_file(self.option('case_id') + '-' + self.option('dad_id'), self.output_dir,
                                        self.option('new_dad_id'))
            self.api_pt.export_tab_file(self.option('preg_id'), self.output_dir, self.option('new_preg_id'))
            self.dad_list = self.get_ref_list()
            self.logger.info("self.family:{}".format(self.family))
            self.logger.info("dad_list:{}".format(self.dad_list))
            self.family_analysis_run()
        elif self.option('is_report') == 'false':
            self.is_dedup = True
            self.num = 4
            self.logger.info("开始进行父本的查错！")
            if self.option('dad_id'):
                dad_id = self.option('dad_id')
            else:
                dad_id = ""
            self.api_pt.export_tab_file(self.option('mom_id'), self.output_dir)
            self.api_pt.export_tab_file(self.option('preg_id'), self.output_dir)
            self.dad_list = self.api_pt.find_father_id(self.option('case_id'), dad_id)
            self.logger.info("dad_list:{}".format(self.dad_list))
            self.check_tab(self.dad_list)  # 用于检查模糊匹配的tab文件是否都存在
            if len(self.dad_list) == 0:
                raise Exception("模糊匹配的父本列表为空，无法进行后面的分析！")
            if self.option('dad_id'):
                self.logger.info("家系只进行错配为2的分析，查重进行2-10！")
                if self.api_pt.tab_exist(self.option('case_id') + "-" + self.option('dad_id')):
                    self.api_pt.export_tab_file(self.option('case_id') + "-" + self.option('dad_id'), self.output_dir)
                    self.family.append(self.option('case_id') + "-" + self.option('dad_id'))
                else:
                    self.api_pt.export_tab_file(self.dad_list[0], self.output_dir)
                    self.family.append(self.dad_list[0])
                self.family.append(self.option('mom_id'))
                self.family.append(self.option('preg_id'))
                self.family_analysis_run()
            else:
                self.logger.info("不进行家系分析直接进行查重分析，计算错配2-10！")
                self.family.append(self.dad_list[0])
                self.family.append(self.option('mom_id'))
                self.family.append(self.option('preg_id'))
                self.dedup_run()
            self.logger.info("self.family:{}".format(self.family))
        super(PtFreeCombineWorkflow, self).run()

    def set_db(self):
        """
        结果文件批量入库
        :return:
        """
        api_pt = self.api.api('medical.paternity_test_v2')
        for i in range(2, self.option('err_min') + 1):
            dir_path = self.output_dir + '/pt_result_' + str(i)
            if not os.path.exists(dir_path):
                continue
            name = self.family[0] + "_" + self.family[1] + "_" + self.family[2] + "_error-" + str(i)
            self.father_err_id = api_pt.add_sg_father_err('no', 'no', self.option("father_id"), i, "2", name,
                                                          self.option('member_id'))
            results = os.listdir(dir_path)
            result_name = self.family[0] + "_" + self.family[1] + "_" + self.family[2]
            for f in results:
                if f == result_name + "_family_joined_tab.txt":
                    api_pt.add_sg_pt_father_detail(dir_path + '/' + f, self.father_err_id)
                elif f == self.family[1] + "_" + self.family[2] + '_info_show.txt':
                    api_pt.add_info_detail(dir_path + '/' + f, self.father_err_id)
                elif f == result_name + '.txt':
                    api_pt.import_dedup_data(dir_path + '/' + f, 'no', self.father_err_id)
                elif f == result_name + "_family.png":
                    api_pt.add_result_pic("pt_result_pic/" + self.option("father_id") + '/pt_result_' + str(i) + '/',
                                          result_name, self.option("father_id"), self.father_err_id)
                elif f == result_name + '_test_pos.txt':
                    api_pt.add_test_pos(dir_path + '/' + f, self.father_err_id)
                elif f == result_name + "_family_analysis.txt":
                    api_pt.add_analysis_tab(dir_path + '/' + f, self.father_err_id, self.option("father_id"), i)
            api_pt.update_father_err_status(self.father_err_id)
        api_pt.update_father_tab_batch(self.option("father_id"), self.family[0], 'free')  # 更新父本的建库批次与抽提批次
        if self.option('is_report') == 'true':  # 不进行出报告的时候 不需要组建家系信息
            api_pt.add_sg_family_info(self.option("father_id"), self.option("new_dad_id"), self.option("new_mom_id"),
                                      'free')

    def end(self):
        self.set_db()
        super(PtFreeCombineWorkflow, self).end()
