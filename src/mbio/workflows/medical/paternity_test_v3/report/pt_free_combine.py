# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
from biocluster.core.exceptions import OptionError
from biocluster.api.file.lib.transfer import MultiFileTransfer
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
            {"name": "is_report", "type": "string"},
            {"name": "dp_correction", "type": "string"}
        ]
        self.add_option(options)
        self.ref_point = self.config.SOFTWARE_DIR + "/database/human/pt_ref/targets.bed.rda"
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_v3_fm/"
        self.ref_all = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_all/"
        self.api_pt = self.api.api("medical.paternity_test_v3.paternity_test_v3")
        self.set_options(self._sheet.options())
        self.father_err_id_homohybrid = None
        self.father_err_id = None
        self.tools_result = []
        self.tools_dedup = []
        self.dad_list = []
        self.is_dedup = False
        self.family = []
        self.rdata = None
        self.path = None
        self.preg_tab = None
        self.dp_correction = None
        self.skip_family_analysis = False
        self.dp_correction_new_id = None
        self.err_min = self.option("err_min")
        self.err_list = [5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        self.err_list_homohybrid = [5, 6, 7, 8, 9, 10, 15, 20]
        self.err_list_valid = [i for i in self.err_list if i <= self.err_min]
        self.err_list_valid_homohybrid = [i for i in self.err_list_homohybrid if i <= self.err_min]
        self.f_id = ''
        self.m_id = ""
        self.s_id = ''
        self.num = []
        self.num_homohybrid = []
        self.similarity_tools = []
        self.bucket = self.config.get_project_region_bucket(project_type="pt_v3").rstrip("/")

    def check_options(self):
        """
        参数检查
        :return:
        """
        if not self.option("is_report"):
            raise OptionError("必须输入is_report参数")
        if not self.option("dp_correction"):
            raise OptionError("必须输入dp_correction参数")
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

    def dp_correction_run(self):
        """
        降噪tool运行
        :return:
        """
        self.logger.info("开始MS降噪")

        mom_tab = self.output_dir + '/' + self.family[1] + '.tab'
        preg_tab = self.output_dir + '/' + self.family[2] + '.tab'

        self.dp_correction = self.add_tool("medical.paternity_test_v3.dp_correction")
        self.dp_correction.set_options({
            "mom_tab": mom_tab,
            "preg_tab": preg_tab,
            "new_id": self.dp_correction_new_id
        })
        self.dp_correction.on("end", self.set_output, "dp_correction")
        self.dp_correction.run()

    def family_analysis_run(self):
        tasks = []
        if self.option("dp_correction") == "true":
            if self.option('is_report') == 'true':
                preg_tab = os.path.join(self.output_dir, "new_tab/{}.tab".format(self.family[2]))
            else:
                preg_tab = os.path.join(self.output_dir, "dp_correction/{}.tab".format(self.family[2]))
        else:
            preg_tab = os.path.join(self.output_dir, "{}.tab".format(self.family[2]))
        for i in self.num:
            pt_analysis = self.add_module("medical.paternity_test_v3.family_analysis")
            pt_analysis.set_options({
                "dad_tab": self.output_dir + '/' + self.family[0] + '.tab',
                "mom_tab": self.output_dir + '/' + self.family[1] + '.tab',
                "preg_tab": preg_tab,
                "ref_point": self.ref_point,
                "err_min": i,
                "is_free_combine": True
            })
            pt_analysis.on('end', self.set_output, 'pt_analysis')
            tasks.append(pt_analysis)
        for i in self.num_homohybrid:
            pt_analysis_homohybrid = self.add_module("medical.paternity_test_v3.family_analysis_homohybrid")
            pt_analysis_homohybrid.set_options({
                "dad_tab": self.output_dir + '/' + self.family[0] + '.tab',
                "mom_tab": self.output_dir + '/' + self.family[1] + '.tab',
                "preg_tab": preg_tab,
                "ref_point": self.ref_point,
                "err_min": i,
                "is_free_combine": True
            })
            pt_analysis_homohybrid.on('end', self.set_output, 'pt_analysis_homohybrid')
            tasks.append(pt_analysis_homohybrid)
        self.on_rely(tasks, self.dedup_run)
        for t in tasks:
            t.run()

    def dedup_run(self):
        """
        查重部分新机制：如果是全库查重的话，直接使用参考库就可以了。
        如果是部分查重的话，更据模糊匹配的结果中的父本进行查重
        :return:
        """
        if self.option("dp_correction") == "true":
            if self.option('is_report') == 'true':
                preg_tab = os.path.join(self.output_dir, "new_tab/{}.tab".format(self.family[2]))
            else:
                preg_tab = os.path.join(self.output_dir, "dp_correction/{}.tab".format(self.family[2]))
        else:
            preg_tab = os.path.join(self.output_dir, "{}.tab".format(self.family[2]))

        for p in self.err_list_valid:
            pt_analysis_dedup = self.add_tool("medical.paternity_test_v3.dedup_v2")
            pt_analysis_dedup.set_options({
                "dad_id": self.family[0],
                "mom_tab": self.output_dir + '/' + self.family[1] + '.tab',
                "preg_tab": preg_tab,
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
                if self.option('is_report') == 'true':
                    self.on_rely(self.tools_dedup, self.sample_similarity_run)
                else:
                    self.on_rely(self.tools_dedup, self.end)
            else:
                if self.option('is_report') == 'true':
                    self.on_rely(self.tools_dedup, self.sample_similarity_run)
                else:
                    self.tools_dedup[0].on("end", self.end)
        else:
            raise Exception('tools_dedup列表为空！')
        for tool in self.tools_dedup:
            tool.run()

    def sample_similarity_run(self):
        """个体识别"""
        for id_ in [self.family[0], self.family[1]]:
            similarity = self.add_tool("medical.paternity_test_v3.sample_similarity")
            similarity.set_options({
                "sample_id": id_
            })
            self.similarity_tools.append(similarity)
        if self.similarity_tools:
            if len(self.similarity_tools) > 1:
                self.on_rely(self.similarity_tools, self.end)
            else:
                self.similarity_tools[0].on("end", self.end)
        else:
            raise Exception("similarity_tools列表为空！")
        for tool in self.similarity_tools:
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
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir + '/pt_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/result_info', self.output_dir + '/pt_result_' +
                         str(obj.option('err_min')))
        if event['data'] == "pt_analysis_homohybrid":
            self.linkdir(obj.output_dir + '/family_analysis', self.output_dir + '/pt_homohybrid_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir + '/pt_homohybrid_result_' +
                         str(obj.option('err_min')))
            self.linkdir(obj.output_dir + '/result_info', self.output_dir + '/pt_homohybrid_result_' +
                         str(obj.option('err_min')))
        if event['data'] == "dedup":
            self.linkdir(obj.output_dir, self.output_dir + '/pt_result_' + str(obj.option('err_min')))
        if event['data'] == "dp_correction":
            """
            added 20181108 by hongyuchen
            如果进行了dp_correction，家系分析和查重所用的胎儿tab文件要发生变化，分两种情况：
            1、report: report模式会输入胎儿的新ID，降噪后样品名由原始的胎儿、母本组合而成，需要重命名为新的ID，重命名后的文件
            放在结果目录下new_tab目录中
            2、dedup：dedup模式不需要重命名降噪后的tab文件，可以直接使用，降噪tab文件在dp_correction目录中；在此情况下
            self.family中的胎儿ID要换成降噪后的胎儿ID
            
            针对self.family中的F、M、S样品名:
            report：使用F、M、S的新ID，无论是否进行dp_correction
            dedup：使用F、M、S的原始输入ID，当进行dp_correction后，S ID要换为降噪后的ID
            """
            self.linkdir(obj.output_dir, "dp_correction")
            # report要将降噪后的胎儿文件重命名
            if self.option('is_report') == 'true':
                if not os.path.exists(os.path.join(self.output_dir, "new_tab")):
                    os.mkdir(os.path.join(self.output_dir, "new_tab"))
                file_1 = os.path.join(self.output_dir, "dp_correction/{}.tab".format(self.dp_correction_new_id))
                file_2 = os.path.join(self.output_dir, "new_tab/{}.tab".format(self.option('new_preg_id')))
                with open(file_1, "r") as f, open(file_2, "w") as out:
                    for line in f:
                        each = line.rstrip("\n").split("\t")
                        each[0] = self.option('new_preg_id')
                        content = "\t".join(each)
                        out.write("{}\n".format(content))
            else:
                self.family[2] = self.dp_correction_new_id               # S ID换为降噪后的ID

            if self.skip_family_analysis:
                self.dedup_run()
            else:
                self.family_analysis_run()

    def get_ref_list_(self):
        """
        用于获取查重参考库中父本的所有的样本id  WQ1802807-S-L1_M
        :return:
        """
        ref_list = []
        pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        refs = pt.get_ref_list(month_back=6)
        for m in refs:
            if os.path.exists(os.path.join(self.ref_data, "{}.tab".format(m))):
                if "-F" in m or "-M" in m:
                    ref_list.append(m)
                elif "-S" in m and "_" in m:
                    ref_list.append(m)
        return ref_list

    def get_ref_id(self, dad_list):
        """
        获取参考组中F M 以及降噪后的胎儿
        :return:
        """
        for m in dad_list:
            if "-S" in m and "_" not in m:
                pass
            else:
                self.dad_list.append(m)

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

    def set_real_sample_id(self):
        """
        WQ1234-2-F1, WQ1234-F1兼容查找这两种的样本
        :return:
        """
        dad_id = self.option('dad_id') if self.option('dad_id') else ""
        self.f_id = self.option('case_id') + '-' + dad_id
        self.m_id = self.option('mom_id')
        # self.s_id = self.option('preg_id')   # WQ12345-S-1_M
        self.s_id = self.option('preg_id')
        self.logger.info("f_id:{}; m_id:{}, s_id:{}".format(self.f_id, self.m_id, self.s_id))

    def run(self):
        """
        因为该分析是自由交互，样本编号之间差异比较大，可能不是同一case下的样本也会进行组合进行分析，所以化繁为简，
        直接放在self.family这个列表中，避免了一系列的判断
        在父本查重的地方，
        :return:
        """
        self.set_real_sample_id()

        # 使用preg_id和mom_id，保证降噪后的new_id唯一
        self.dp_correction_new_id = self.option("preg_id") + "_" + self.option("mom_id")

        if self.option('is_report') == 'true':
            self.logger.info("开始进行家系分析出报告！")
            self.skip_family_analysis = False
            self.num.extend(self.err_list_valid)
            self.num_homohybrid.extend(self.err_list_valid_homohybrid)
            self.family.append(self.option('new_dad_id'))
            self.family.append(self.option('new_mom_id'))
            self.family.append(self.option('new_preg_id'))
            self.api_pt.attempt_to_export_tab_file(self.m_id, self.output_dir, self.option('new_mom_id'))
            self.api_pt.attempt_to_export_tab_file(self.f_id, self.output_dir, self.option('new_dad_id'))
            self.api_pt.attempt_to_export_tab_file(self.s_id, self.output_dir, self.option('new_preg_id'))
            self.dad_list = self.get_ref_list_()
            self.logger.info("dad_list:{}".format(self.dad_list))
        elif self.option('is_report') == 'false':
            self.num = [5]                     # dedup模式家系分析只进行err_min为5的
            self.num_homohybrid = [5]
            self.is_dedup = True
            self.logger.info("开始进行父本的查错！")
            if self.option('dad_id'):
                dad_id = self.option('dad_id')
            else:
                dad_id = ""
            self.api_pt.attempt_to_export_tab_file(self.m_id, self.output_dir)
            self.api_pt.attempt_to_export_tab_file(self.s_id, self.output_dir)
            dad_list_ = self.api_pt.find_father_id(self.option('case_id'), dad_id)
            self.get_ref_id(dad_list_)  # 过滤没有降噪的S样本
            self.logger.info("dad_list:{}".format(self.dad_list))
            self.check_tab(self.dad_list)  # 用于检查模糊匹配的tab文件是否都存在
            if len(self.dad_list) == 0:
                raise Exception("模糊匹配的父本列表为空，无法进行后面的分析！")
            if self.option('dad_id'):
                self.logger.info("家系只进行错配为5的分析，所有错配查重！")
                self.skip_family_analysis = False
                if self.api_pt.tab_exist(self.f_id):
                    self.api_pt.attempt_to_export_tab_file(self.f_id, self.output_dir)
                    self.family.append(self.f_id)
                else:
                    self.api_pt.attempt_to_export_tab_file(self.dad_list[0], self.output_dir)
                    self.family.append(self.dad_list[0])
                self.family.append(self.m_id)
                self.family.append(self.s_id)
            else:
                self.logger.info("不进行家系分析直接进行查重分析，计算所有错配！")
                self.skip_family_analysis = True
                self.family.append(self.dad_list[0])
                self.family.append(self.m_id)
                self.family.append(self.s_id)

        self.logger.info("family:{}".format(self.family))

        if self.option("dp_correction") == "true":
            self.dp_correction_run()
        elif not self.skip_family_analysis:
            self.family_analysis_run()
        else:
            self.dedup_run()

        super(PtFreeCombineWorkflow, self).run()

    def set_db(self):
        """
        结果文件批量入库
        :return:
        """
        api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        # 胎儿杂合导表部分
        for i in self.err_list_valid:
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
                    api_pt.add_info_detail(dir_path + '/' + f, self.father_err_id, self.family[0], self.family[1],
                                           self.family[2])
                elif f == result_name + '.txt':
                    api_pt.import_dedup_data(dir_path + '/' + f, 'no', self.father_err_id)
                    # api_pt.import_dedup_data(dir_path + '/' + f, 'no', self.father_err_id)
                elif f == result_name + "_family.png":
                    if self.config.RGW_ENABLE:
                        pic_path = "{}/pt_result_pic/".format(self.bucket) + self.option("father_id") + \
                                   '/pt_result_' + str(i) + '/'
                    else:
                        pic_path = "rerewrweset/MEDfiles/pt_result_pic/" + self.option("father_id") + \
                                   '/pt_result_' + str(i) + '/'
                    api_pt.add_result_pic(pic_path, result_name, self.option("father_id"), self.father_err_id)
                elif f == result_name + '_test_pos.txt':
                    api_pt.add_test_pos(dir_path + '/' + f, self.father_err_id)
                elif f == result_name + "_family_analysis.txt":
                    api_pt.add_analysis_tab(dir_path + '/' + f, self.father_err_id, self.option("father_id"), i)
            api_pt.update_father_err_status(self.father_err_id)
        # 胎儿纯合与杂合导表部分
        for i in self.err_list_valid_homohybrid:
            dir_path = self.output_dir + '/pt_homohybrid_result_' + str(i)
            if not os.path.exists(dir_path):
                continue
            name = self.family[0] + "_" + self.family[1] + "_" + self.family[2] + "_error-" + str(i)
            self.father_err_id_homohybrid = api_pt.add_sg_father_err('no', 'no', self.option("father_id"), i, "2",
                                                                     name, self.option('member_id'),
                                                                     postype='homohybrid')
            results = os.listdir(dir_path)
            result_name = self.family[0] + "_" + self.family[1] + "_" + self.family[2]
            for f in results:
                if f == result_name + "_family_match_update.png":
                    if self.config.RGW_ENABLE:
                        pic_path = "{}/pt_result_pic/".format(self.bucket) + self.option("father_id") + \
                                   '/pt_homohybrid_result_' + str(i) + '/'
                    else:
                        pic_path = "rerewrweset/MEDfiles/pt_result_pic/" + self.option("father_id") + \
                                   '/pt_result_' + str(i) + '/'
                    api_pt.add_result_pic(pic_path, result_name, self.option("father_id"),
                                          self.father_err_id_homohybrid, postype='homohybrid')
                elif f == result_name + '_test_pos_update.txt':
                    api_pt.add_test_pos(dir_path + '/' + f, self.father_err_id_homohybrid, postype='homohybrid')
                elif f == result_name + "_family_analysis.txt":
                    api_pt.add_analysis_tab(dir_path + '/' + f, self.father_err_id_homohybrid,
                                            self.option("father_id"), i, postype='homohybrid')
                elif f == result_name + "_num_pos_map.txt":
                    api_pt.add_num_pos_map(dir_path + '/' + f, self.father_err_id_homohybrid, self.option("father_id"))
            api_pt.update_father_err_status(self.father_err_id_homohybrid, postype='homohybrid')
        api_pt.update_father_tab_batch(self.option("father_id"), self.family[0], 'free')  # 更新父本的建库批次与抽提批次
        # if self.option('is_report') == 'true':  # 不进行出报告的时候 不需要组建家系信息
        #     api_pt.add_sg_family_info(self.option("father_id"), self.option("new_dad_id"), self.option("new_mom_id"),
        #                               'free')

    def end(self):
        self.set_db()
        self.upload_target_path()
        super(PtFreeCombineWorkflow, self).end()

    def upload_target_path(self):
        """
        将结果link到指定磁盘路径下
        :return:
        """
        self.logger.info("开始设置结果目录！")
        if self.config.RGW_ENABLE:
            transfer = MultiFileTransfer()
            transfer.add_upload(self.output_dir + "/", "{}/pt_result_pic/{}/".format(self.bucket,
                                                                                     self.option("father_id")),
                                base_path=self.output_dir + "/")
            transfer.perform()
        else:
            if self._sheet.client == 'client01':
                self.path = '/mnt/ilustre/data/rerewrweset/MEDfiles/pt_result_pic/{}'.format(self.option("father_id"))
            else:
                self.path = '/mnt/ilustre/tsanger-data/rerewrweset/MEDfiles/pt_result_pic/{}'\
                    .format(self.option("father_id"))
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            self.linkdir(self.output_dir, self.path)
        self.logger.info("设置结果目录成功！")