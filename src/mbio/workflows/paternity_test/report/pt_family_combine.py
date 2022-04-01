# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

'''医学检验所-无创产前亲子鉴定家系自由组合模块'''
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import re
from bson import ObjectId
from biocluster.config import Config
import json
import shutil


class PtFamilyCombineWorkflow(Workflow):
    """
    last_modify 20170914 by hongdongxuan
    运行逻辑，当只有父本组，父本过滤器，母本，胎儿，错配的时候只会进行查重，按照过滤器将所有的样本取出然后进行查重分析，不进行家系分析
    当有新父本编号与新母本编号，新胎儿编号的时候，会进行家系分析，并进行生成新的图片，家系名字按照新的样本来命名。
    """
    def __init__(self, wsheet_object):
        '''
        :param wsheet_object:
        '''
        self._sheet = wsheet_object
        super(PtFamilyCombineWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "dad_id", "type": "string"},  # 父本过滤器
            {"name": "mom_id", "type": "string"},
            {"name": "preg_id", "type": "string"},
            {"name": "err_min", "type": "int", "default": 2},  # 允许错配数
            {"name": "dad_group", "type": "string"},   # 父本组
            {"name": "update_info", "type": "string"},
            {"name": 'new_mom_id', "type": "string"},
            {"name": 'new_dad_id', "type": "string"},
            {"name": 'new_preg_id', "type": "string"},
            {"name": 'dedup_start', "type": "string"},  # 区域查重的开始
            {"name": 'dedup_end', "type": "string"},   # 区域查重的结束
            {"name": 'dedup_all', "type": "string"},  # 是否全库查重
            {"name": 'main_id', "type": "string"},   # 主表id，接口中生成的
            {"name": 'member_id', "type": "string"},  # 做记录
        ]
        self.add_option(options)
        self.ref_dir = self.config.SOFTWARE_DIR
        self.tab_file = self.api.tab_file
        self.sg_pt = self.api.sg_paternity_test
        self.ref_data = self.config.SOFTWARE_DIR + "/database/human/pt_ref/tab_data/"
        self.set_options(self._sheet.options())
        self.tools_analysis = []
        self.tools_result = []
        self.tools_dedup = []
        self.sure_end = "false"  # 用于区分流程是否只是进行家系排查分析，不进行家系父权制，画图等分析
        self.dedup_done = ''
        self.no_dedup = ''  # 用于判断后面是否进行查重，不查重的时候直接在set_output中激发end
        self.dad_list = []

    def check_options(self):
        """
        参数检查
        :return:
        """
        if not self.option("dad_group"):
            raise OptionError("必须输入父本组编号")
        if not self.option("mom_id"):
            raise OptionError("必须输入母本编号")
        if not self.option("preg_id"):
            raise OptionError("必须输入胎儿编号")
        if self.option('new_dad_id') or self.option('new_mom_id') or self.option('new_preg_id'):
            if not self.option('new_dad_id') or not self.option('new_mom_id') or not self.option('new_preg_id'):
                raise OptionError("必须同时输入新的父本，母本，胎儿的编号，否则后面无法进行家系分析！")
            if not self.option('mom_id') or not self.option('preg_id') or not self.option('dad_group') or \
                    not self.option('dad_id'):
                raise OptionError("必须输入具体的父本组，过滤器，母本编号，胎儿编号，否则后面无法进行家系分析！")
        return True

    """
    def cat_sample_run(self):
        self.cat_sample.set_options({
            "sample_1": self.preg_list[0],
            "sample_2": self.preg_list[1],
        })
        self.cat_sample.run()

    def fastq2mongo_dc_run(self):
        i = self.preg_list[0] + "_" + self.preg_list[1]
        self.fastq2mongo_dc.set_options({
            "sample_id": i,
            "fastq_path": self.cat_sample.option("fastq_path"),
            "cpu_number": 4,
            "ref_fasta": self.ref_dir + "/database/human/hg38.chromosomal_assembly/ref.fa",
            "targets_bedfile": self.ref_dir + "/database/human/pt_ref/snp.chr.sort.3.bed",
            "type": 'pt'
        })
        self.fastq2mongo_dc.on('end', self.set_output, 'fastq2mongo')
        self.fastq2mongo_dc.run()
    """

    def pt_analysis_run(self):
        for i in range(len(self.family)):
            pt_analysis = self.add_module("paternity_test.pt_analysis")
            dad_id = self.family[i][0]
            mom_id = self.family[i][1]
            preg_id = self.family[i][2]
            pt_analysis.set_options({
                "dad_tab": self.output_dir + '/' + dad_id + '.tab',
                "mom_tab": self.output_dir + '/' + mom_id + '.tab',
                "preg_tab": self.output_dir + '/' + preg_id + '.tab',
                "ref_point": self.ref_dir + "/database/human/pt_ref/targets.bed.rda",
                "err_min": self.option("err_min")
            })
            self.tools_analysis.append(pt_analysis)
        for j in range(len(self.tools_analysis)):
            self.tools_analysis[j].on('end', self.set_output, 'pt_analysis')
        if len(self.tools_analysis) > 1:
            self.on_rely(self.tools_analysis, self.result_info_run)
        elif len(self.tools_analysis) == 1:
            self.tools_analysis[0].on('end', self.result_info_run)
        for t in self.tools_analysis:
            t.run()

    def result_info_run(self):
        if len(self.family) == 1:
            result_dir = self.work_dir + "/PtAnalysis/FamilyMerge/output"
        else:
            result_dir = self.output_dir
        results = os.listdir(result_dir)
        for f in results:
            if re.match(r'.*family_joined_tab\.Rdata$', f):
                result_info = self.add_tool("paternity_test.result_info")
                self.rdata = os.path.join(result_dir, f)
                result_info.set_options({
                    "tab_merged": self.rdata
                })
                self.tools_result.append(result_info)
            else:
                pass
        for j in range(len(self.tools_result)):
            self.tools_result[j].on('end', self.set_output, 'result_info')
        if len(self.tools_result) > 1:
            if not self.option('dedup_start') and self.option('dedup_all') == 'false':
                self.on_rely(self.tools_result, self.end)
            else:
                self.on_rely(self.tools_result, self.dedup_run)
        else:
            if not self.option('dedup_start') and self.option('dedup_all') == 'false':
                self.no_dedup = "true"
                pass
                # self.tools_result[0].on('end', self.end)
            else:
                self.on_rely(self.tools_result, self.dedup_run)
        for t in self.tools_result:
            t.run()

    def dedup_run(self):
        """
        查重部分新机制：如果是全库查重的话，直接使用参考库就可以了。
        如果是部分查重的话根据输入的dedup_start:WQ12345678和dedup_end:WQ12349999中WQ后面跟的数字来确定查重的父本
        :return:
        """
        father_data = self.output_dir + "/" + "father"  # 自定义查重父本
        if not os.path.exists(father_data):
            os.mkdir(father_data)
        if self.option('dedup_all') != "true":
            api_read_tab = self.api.tab_file
            # temp_s = re.match('WQ([1-9].*)', self.option('dedup_start'))
            # temp_e = re.match('WQ([1-9].*)', self.option('dedup_end'))
            # num_list = []
            # if str(temp_s.group(1)) == str(temp_e.group(1)):
            #     num_list.append(temp_s.group(1))
            # else:
            #     num_list.append(temp_s.group(1))
            #     num_list.append(temp_e.group(1))
            # num_list = range(int(temp_s.group(1)), int(temp_e.group(1)))
            # self.logger.info("查重区域{}-{}".format(temp_s.group(1), temp_e.group(1)))
            num_list = []
            if self.option('dedup_start') == self.option('dedup_end'):
                num_list.append(self.option('dedup_start'))
            else:
                num_list.append(self.option('dedup_start'))
                num_list.append(self.option('dedup_end'))
            name_list = []
            for m in num_list:
                x = api_read_tab.dedup_sample_report(m)
                if len(x):
                    for k in range(len(x)):
                        name_list.append(x[k])
            name_list = list(set(name_list))
            self.logger.info("区域查重所有的父本:%s" % name_list)
            if len(name_list) == 0:
                self.dedup_done = 'true'
                self.end()
            for m in name_list:
                api_read_tab.export_tab_file(str(m), father_data)
            self.logger.info("导出所有的区域查重的父本文件完成！")
            for i in range(len(self.family)):
                mom_id = self.family[i][1]
                preg_id = self.family[i][2]
                dad_id = self.family[i][0]
                pt_analysis_dedup = self.add_tool("paternity_test.dedup")
                pt_analysis_dedup.set_options({
                    "mom_tab": self.output_dir + '/' + mom_id + '.tab',
                    "preg_tab": self.output_dir + '/' + preg_id + '.tab',
                    "ref_point": self.ref_dir + "/database/human/pt_ref/targets.bed.rda",
                    "err_min": self.option("err_min"),
                    "father_path": father_data + "/",
                    "dad_id": dad_id
                })
                self.tools_dedup.append(pt_analysis_dedup)
        else:
            for i in range(len(self.family)):
                mom_id = self.family[i][1]
                preg_id = self.family[i][2]
                dad_id = self.family[i][0]
                pt_analysis_dedup = self.add_tool("paternity_test.dedup")
                pt_analysis_dedup.set_options({
                    "mom_tab": self.output_dir + '/' + mom_id + '.tab',
                    "preg_tab": self.output_dir + '/' + preg_id + '.tab',
                    "ref_point": self.ref_dir + "/database/human/pt_ref/targets.bed.rda",
                    "err_min": self.option("err_min"),
                    "father_path": self.ref_data,
                    "dad_id": dad_id
                })
                self.tools_dedup.append(pt_analysis_dedup)
        for j in range(len(self.tools_dedup)):
            self.tools_dedup[j].on('end', self.set_output, 'dedup')
        if len(self.tools_dedup) > 1:
            self.on_rely(self.tools_dedup, self.end)
        else:
            self.tools_dedup[0].on("end", self.end)
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
        # if event['data'] == "fastq2mongo":
        #     tab_name = self.option('preg_id') + '.mem.sort.hit.vcf.tab'
        #     self.linkdir(obj.output_dir + '/fastq2tab', self.output_dir)
        #     new_tab_name = self.option('preg_id') + '.tab'
        #     os.link(obj.output_dir + '/fastq2tab/' + tab_name, os.path.join(self.output_dir, new_tab_name))
        if event['data'] == "pt_analysis":
            self.linkdir(obj.output_dir + '/family_analysis', self.output_dir)
            self.linkdir(obj.output_dir + '/family_merge', self.output_dir)
        if event['data'] == "result_info":
            self.linkdir(obj.output_dir, self.output_dir)
        if event['data'] == "dedup":
            self.linkdir(obj.output_dir, self.output_dir)
        if self.no_dedup == "true" or self.sure_end == "true":
            self.end()
        else:
            pass

    def run(self):

        """
        可以组成的家系全部放在self.family这个列表中
        如果存在dad_id,那么构成的家系肯定只有一组（更换id的新家系，不更换id的老家系）
        在组成家系的时候将tab文件从库中拿出
        :return:
        """
        self.family = []
        if self.option('new_dad_id'):  # 如果存在这个参数，表示后面的流程中要修改相应的tab文件的名字，替换完成自由交互出报告的功能
            combine = []
            old_mom_id = self.option('mom_id')
            old_dad_id = self.option('dad_group') + '-' + self.option('dad_id')
            old_preg_id = self.option('preg_id')
            self.mom_id = self.option('new_mom_id')  # 新的家系的各个样本的id名称
            self.dad_id = self.option('new_dad_id')
            self.preg_id = self.option('new_preg_id')
            combine.append(self.dad_id)
            combine.append(self.mom_id)
            combine.append(self.preg_id)
            self.tab_file.export_tab_file(old_mom_id, self.output_dir, self.mom_id)
            self.tab_file.export_tab_file(old_dad_id, self.output_dir, self.dad_id)
            self.tab_file.export_tab_file(old_preg_id, self.output_dir, self.preg_id)
            self.family.append(combine)
            self.logger.info("self.family:{}".format(self.family))
            self.pt_analysis_run()
        else:
            self.dad_list = self.tab_file.find_father_id(self.option('dad_group'), self.option('dad_id'))
            self.logger.info("dad_list:{}".format(self.dad_list))
            self.tab_file.export_tab_file(self.option('mom_id'), self.output_dir)
            self.tab_file.export_tab_file(self.option('preg_id'), self.output_dir)
            self.sure_end = "true"
            self.free_analysis_run()

        """
        if self.option('dad_id'):
            combine = []
            old_mom_id = self.option('mom_id')
            old_dad_id = self.option('dad_group') + '-' + self.option('dad_id')
            old_preg_id = self.option('preg_id')
            if self.option('new_dad_id'):
                self.mom_id = self.option('new_mom_id')  # 新的家系的各个样本的id名称
                self.dad_id = self.option('new_dad_id')
                self.preg_id = self.option('new_preg_id')
                combine.append(self.dad_id)
                combine.append(self.mom_id)
                combine.append(self.preg_id)
                self.tab_file.export_tab_file(old_mom_id, self.output_dir, self.mom_id)
                self.tab_file.export_tab_file(old_dad_id, self.output_dir, self.dad_id)
                self.tab_file.export_tab_file(old_preg_id, self.output_dir, self.preg_id)
            else:
                combine.append(old_dad_id)
                combine.append(old_mom_id)
                combine.append(old_preg_id)
                self.tab_file.export_tab_file(old_mom_id, self.output_dir)
                self.tab_file.export_tab_file(old_dad_id, self.output_dir)
                self.tab_file.export_tab_file(old_preg_id, self.output_dir)
            self.family.append(combine)
        else:
            dad_list = self.tab_file.find_father_id(self.option('dad_group'))
            self.tab_file.export_tab_file(self.option('mom_id'), self.output_dir)
            self.tab_file.export_tab_file(self.option('preg_id'), self.output_dir)
            for i in dad_list:
                combine = []
                combine.append(i)
                combine.append(self.option('mom_id'))
                combine.append(self.option('preg_id'))
                self.family.append(combine)
                self.tab_file.export_tab_file(i, self.output_dir)
        self.logger.info(self.family)
        self.pt_analysis_run()
        """
        """
        else:
            self.cat_sample = self.add_tool('paternity_test.cat_sample')
            self.fastq2mongo_dc = self.add_module("paternity_test.fastq2mongo_dc")
            self.cat_sample.on('end', self.fastq2mongo_dc_run)
            self.fastq2mongo_dc.on('end', self.pt_analysis_run)
            self.pt_analysis.on('end', self.result_info_run)
            self.result_info.on('end', self.dedup_run)
            self.cat_sample_run()
        """
        super(PtFamilyCombineWorkflow, self).run()

    def free_analysis_run(self):
        """
        更据父本组以及过滤器进行模糊匹配后的所有样本进行家系分析
        :return:
        """
        self.free_analysis = self.add_tool("paternity_test.dedup")
        all_sample_data = self.output_dir + "/" + "sample"  # 自定义自由分析的文件路径
        if not os.path.exists(all_sample_data):
            os.mkdir(all_sample_data)
        self.logger.info("模糊匹配的样本个数{}".format(len(self.dad_list)))
        for m in self.dad_list:
            self.tab_file.export_tab_file(m, all_sample_data)
        self.logger.info("数据库中导出样本全部完成！")
        self.free_analysis.set_options({
            "mom_tab": self.output_dir + '/' + self.option('mom_id') + '.tab',
            "preg_tab": self.output_dir + '/' + self.option('preg_id') + '.tab',
            "ref_point": self.ref_dir + "/database/human/pt_ref/targets.bed.rda",
            "err_min": self.option("err_min"),
            "father_path": all_sample_data + "/",
            "dad_id": "free"
        })
        self.free_analysis.on('end', self.set_output, "dedup")
        self.free_analysis.run()

    def end(self):  # 和亲子主流程的导表基本一致
        self.logger.info("开始end函数")
        api_main = self.api.sg_paternity_test
        if self.option('dedup_all') == 'true':
            dedup_num = 'all'
        else:
            if self.option('dedup_start'):
                dedup_num = self.option('dedup_start') + '-' + self.option('dedup_end')
            else:
                dedup_num = 'no'
        results = os.listdir(self.output_dir)
        self.logger.info("output: {}".format(results))
        if self.sure_end != "true":
            for i in range(len(self.family)):
                self.dedup_done = 'false'
                dad_id = self.family[i][0]
                mom_id = self.family[i][1]
                preg_id = self.family[i][2]
                self.father_id = api_main.add_sg_father(dad_id, mom_id, preg_id, self.option('main_id'),
                                                        self.option("member_id"), type='free')
                # 此处的main_id相当于别处的batch_id为该自由交互的主表，和正式流程里的不一样

                self.pt_father_id = api_main.add_pt_father(father_id=self.father_id, err_min=self.option('err_min'), dedup=dedup_num)  # 交互表id
                dedup = '.*' + mom_id + '_' + preg_id + '_family_analysis.txt'
                dedup_new = dad_id + '_' + mom_id + '_' + preg_id + '.txt'
                self.logger.info(dedup_new)
                for f in results:
                    if re.search(dedup, f):  # 结果表
                        self.logger.info('结果表')
                        api_main.add_analysis_tab(self.output_dir + '/' + f, self.pt_father_id)
                    elif f == dad_id + '_' + mom_id + '_' + preg_id + '_family_joined_tab.txt':  # 调试表
                        self.logger.info('调试表')
                        api_main.add_sg_pt_father_detail(self.output_dir + '/' + f, self.pt_father_id)
                    elif f == mom_id + '_' + preg_id + '_info_show.txt':  # 母本和胎儿的浓度信息
                        self.logger.info('母本和胎儿的浓度信息')
                        api_main.add_info_detail(self.output_dir + '/' + f, self.pt_father_id)
                    elif f == dad_id + '_' + mom_id + '_' + preg_id + '_test_pos.txt':  # 报告位点
                        self.logger.info('报告位点')
                        api_main.add_test_pos(self.output_dir + '/' + f, self.pt_father_id)
                    elif f == dad_id + '_' + mom_id + '_' + preg_id + '_family.png':  # 插入图
                        self.logger.info('插入图')
                        file_dir = self.output_dir + '/' + dad_id + '_' + mom_id + '_' + preg_id
                        api_main.add_pt_father_figure(file_dir, self.pt_father_id)
                    elif str(f) == str(dedup_new):  # 查重表
                        self.logger.info('查重表')
                        self.dedup_done = 'true'
                        api_main.import_dedup_data(self.output_dir + '/' + f, self.pt_father_id)
                ############
                # 写这一段是因为如果只有一个家系存在查重的时候很有可能在end开始时还没有完成相应set_output,所以去直接指定路径获取查重文件
                # 防止数据漏掉
                dedup_path = self.work_dir + '/Dedup/output/' + dedup_new
                if len(self.family) == 1 and os.path.exists(dedup_path) and self.dedup_done == 'false':
                    self.logger.info('查重表')
                    api_main.import_dedup_data(dedup_path, self.pt_father_id)
                ############

                if dad_id + '_' + mom_id + '_' + preg_id + '_family.png' not in results:  # 20170707 zhouxuan
                    api_main.has_problem(self.pt_father_id, dad_id)
                elif dad_id + '_' + mom_id + '_' + preg_id + '_fig2.png' not in results:
                    api_main.has_problem(self.pt_father_id, dad_id)
                if mom_id + '_' + preg_id + '_info_show.txt' not in results:
                    api_main.update_infoshow(self.pt_father_id, mom_id, preg_id)

                api_main.add_father_result(self.father_id, self.pt_father_id, dad_id)
                api_main.add_father_qc(self.father_id, self.pt_father_id)
                api_main.update_sg_pt_father(self.pt_father_id)
        else:
            for i in results:
                if str(i) == 'free_' + self.option('mom_id') + '_' + self.option('preg_id') + '.txt':
                    self.logger.info("开始导入家系排查的查重结果！")
                    api_main.import_dedup_data(self.output_dir + '/' + i, ObjectId(self.option('main_id')))
        super(PtFamilyCombineWorkflow, self).end()
