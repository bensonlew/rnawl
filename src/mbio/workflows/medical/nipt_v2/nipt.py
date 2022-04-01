# -*- coding: utf-8 -*-
# __author__ = 'hongdongxuan'
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from biocluster.config import Config
import json
import os
import re
import datetime


class NiptWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        医学检验所-无创产前筛查流程
        lasted modified by HONGDONG 20171218
        :param wsheet_object:
        """
        self._sheet = wsheet_object
        super(NiptWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},
            # {"name": "batch_id", "type": "string"},
            {'name': 'member_id', 'type': 'string'},
            {"name": "bw", "type": "int", "default": 10},
            {"name": "bs", "type": "int", "default": 1},
            {"name": "ref_group", "type": "int", "default": 1},
            {"name": "update_info", "type": "string"},
            {"name": "single", "type": "string", "default": "false"},  # 单端与双端
            {"name": "sample_list", "type": "string"},  # 产筛样本id
            {"name": "sanger_type", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.api_nipt = self.api.api('medical.nipt_analysis_v2')
        self.analysis_id = json.loads(self.option('update_info')).keys()[0]
        self.samples = None
        self.main_id = None
        self.modules = []
        self.bin_step = [{'bin': 10, 'step': 1}, {'bin': 5, 'step': 5}, {'bin': 1, 'step': 1},
                         {'bin': 500, 'step': 500}]

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option('fastq_path').is_set:
            raise OptionError('必须提供fastq文件所在的路径')
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

    def analysis_run(self):
        n = 0
        small_sample = []
        for i in self.samples:
            file_ = os.path.join(self.option('fastq_path').prop['path'], i + "_R1.fastq.gz")
            if os.path.getsize(file_) < 1048576:
                small_sample.append(file_)
                self.samples.remove(i)
                self.api_nipt.add_main_(self.option('member_id'), i, self.analysis_id)
                self.api_nipt.update_analysis_status(self.analysis_id, "1")  # 样本fastq文件小的时候也算在运行成功+1
        self.logger.info('small_sample：%s' % small_sample)
        self.logger.info("samples: %s" % self.samples)
        for sample in self.samples:
            nipt_analysis = self.add_module("medical.nipt_v2.nipt_analysis")
            self.step.add_steps('nipt_analysis{}'.format(n))
            nipt_analysis.set_options({
                "sample_id": sample,
                "fastq_path": self.option("fastq_path"),
                "bw": self.option('bw'),
                'bs': self.option('bs'),
                'ref_group': self.option('ref_group'),
                "single": self.option("single")
            })
            step = getattr(self.step, 'nipt_analysis{}'.format(n))
            step.start()
            nipt_analysis.on('end', self.finish_update, 'nipt_analysis{}'.format(n))
            self.modules.append(nipt_analysis)
            n += 1
            self.main_id = self.api_nipt.add_main(self.option('member_id'), sample, self.analysis_id)
        if self.modules:
            if len(self.modules) != 1:
                for j in range(len(self.modules)):
                    self.modules[j].on('end', self.import_tab)
            if len(self.modules) > 1:
                self.on_rely(self.modules, self.end)
            elif len(self.modules) == 1:
                self.modules[0].on('end', self.end)
        else:
            raise Exception("产筛的module列表为空！")
        for module in self.modules:
            module.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        # newdir = os.path.join(self.output_dir, dirname)
        newdir = dirname
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

    def import_tab(self, event):
        """
        用于进行样本计算完成后的单独导表，不在一起运行完成之后在进行导表，考虑到阻塞的问题，通过wpm请求任务时候不能阻塞不能超过10mins
        :return:
        """
        name = ''
        obj = event["bind_object"]
        file_name = os.listdir(obj.output_dir)
        for file__ in file_name:
            if re.match(r'.*\.bed\.2$', file__):
                name = file__.split(".")[0]
                self.logger.info("接下来要进行导表的样本：%s" % name)
                break
        if name not in self.samples:
            raise Exception("该样本名不在上机表中，流程终止")
        main_id = self.api_nipt.get_id(name, self.analysis_id)
        for m in self.bin_step:
            interaction_id = self.api_nipt.add_interaction(main_id, m['bin'], m['step'], self.option('ref_group'), name)
            for i in file_name:
                if i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_z.xls':
                    self.api_nipt.add_z_result(obj.output_dir + '/' + i, interaction_id)  # 插入不同交互的z值
                elif i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_zz.xls':
                    self.api_nipt.add_zz_result(obj.output_dir + '/' + i, interaction_id)  # 插入不同交互的zz值
                # elif i == name + '_result.txt':
                #     self.api_nipt.report_result(interaction_id, obj.output_dir + '/' + i, self.analysis_id)
            self.api_nipt.update_interaction(name, interaction_id)  # 更新每个交互的状态
        for i in file_name:
            # if i == name + '.bed.2' and not self.api_nipt.check_exist_bed(str(name)):
            #     self.logger.info("数据库中不存在：{},即将进行导表！".format(name))
            #     self.api_nipt.add_bed_file(obj.output_dir + '/' + i)
            if re.search(name + '.*_fastqc.html$', i):
                sanger_path = Config().get_netdata_config(self.option('sanger_type'))
                path = sanger_path[self.option('sanger_type') +
                                   "_path"] + "/rerewrweset/MEDfiles/nipt_fastqc/{}".format(
                    datetime.datetime.now().strftime("%Y%m%d"))
                if not os.path.exists(path):
                    os.mkdir(path)
                if not os.path.exists(path + '/' + i):
                    os.link(obj.output_dir + '/' + i, path + '/' + i)
                file_dir = os.path.basename(path)
                file_path = '{}/{}'.format(file_dir, i)
                self.api_nipt.add_fastqc(file_path, self.analysis_id)  # fastqc入库
            elif i == name + '_result.txt':
                self.api_nipt.report_result(main_id, obj.output_dir + '/' + i, self.analysis_id)
                self.api_nipt.update_main(main_id, obj.output_dir + '/' + i)  # 更新zz值到主表中
        self.api_nipt.add_qc(obj.output_dir + '/' + name + '.qc', obj.output_dir + '/' + name + '.gc',
                             obj.output_dir + '/new_gc.txt', self.analysis_id)
        self.api_nipt.update_analysis_status(self.analysis_id, "1")

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == 'nipt_analysis':
            allfiles = os.listdir(obj.output_dir)
            oldfiles = [os.path.join(obj.output_dir, i) for i in allfiles]
            newfiles = [os.path.join(self.output_dir, i) for i in allfiles]
            for i in range(len(allfiles)):
                os.link(oldfiles[i], newfiles[i])

    def set_single_output(self):
        allfiles = os.listdir(self.modules[0].output_dir)
        self.logger.info("allfiles:%s" % allfiles)
        oldfiles = [os.path.join(self.modules[0].output_dir, i) for i in allfiles]
        newfiles = [os.path.join(self.output_dir, i) for i in allfiles]
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def check_file(self, samples_list):
        """
        用于检查样本是否已经分析过了
        :param samples_list:
        :return:
        """
        samples_ = []
        for m in samples_list.strip().split(','):
            m = m.strip()   # 移除样本名前面的空格
            self.logger.info(m)
            if not self.api_nipt.check_exist_bed(m):
                samples_.append(m)
                self.logger.info("样本{}添加到分析队列！".format(m))
            else:
                self.api_nipt.update_analysis_status(self.analysis_id, "1")   # 之前计算过了 也进行加一
                self.logger.info("样本{}已经分析过，不添加到分析队列！".format(m))
        return samples_

    def run(self):
        """
        运行逻辑：根据传入进来的样本编号列表，首先检查样本的bed文件是否已经存在，存在就不在进行分析了，不存的时候在进行bed文件分析
        :return:
        """
        self.samples = self.check_file(self.option("sample_list"))
        if len(self.samples) == 0:
            self.start_listener()
            self.fire("start")
            self.logger.error("该批次的样本都分析过了，没有样本可以进行分析！")
            self.api_nipt.update_analysis_status(self.analysis_id, '2')
            self.end()
            # self.exit(exitcode=1, data='该批次的样本都分析过了，没有样本可以进行分析！', terminated=False)
        else:
            # self.api_nipt.update_analysis_status(self.analysis_id, "3", len(self.samples)),  # 初始化下all_count
            self.analysis_run()
            super(NiptWorkflow, self).run()

    def end(self):
        super(NiptWorkflow, self).end()
        if len(self.modules) == 1:
            self.logger.info("只有一个样本")
            self.set_single_output()
            self.logger.info("移动结果文件成功")
            name = self.samples[0]
            main_id = self.api_nipt.get_id(name, self.analysis_id)
            for m in self.bin_step:
                interaction_id = self.api_nipt.add_interaction(main_id, m['bin'], m['step'],
                                                               self.option('ref_group'), name)
                for i in os.listdir(self.modules[0].output_dir):
                    if i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_z.xls':
                        self.api_nipt.add_z_result(self.modules[0].output_dir + '/' + i, interaction_id)  # 插入不同交互的z值
                    elif i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_zz.xls':
                        self.api_nipt.add_zz_result(self.modules[0].output_dir + '/' + i, interaction_id)  # 插入不同交互的zz值
                    # elif i == name + '_result.txt':
                    #     self.api_nipt.report_result(interaction_id, self.modules[0].output_dir + '/' + i,
                    #                                 self.analysis_id)
                self.api_nipt.update_interaction(name, interaction_id)  # 更新每个交互的状态
            for i in os.listdir(self.modules[0].output_dir):
                # if i == name + '.bed.2' and not self.api_nipt.check_exist_bed(str(name)):
                #     self.api_nipt.add_bed_file(self.modules[0].output_dir + '/' + i)
                if re.search(name + '.*_fastqc.html$', i):
                    sanger_path = Config().get_netdata_config(self.option('sanger_type'))
                    path = sanger_path[self.option('sanger_type') +
                                       "_path"] + "/rerewrweset/MEDfiles/nipt_fastqc/{}".format(
                        datetime.datetime.now().strftime("%Y%m%d"))
                    if not os.path.exists(path):
                        os.mkdir(path)
                    if not os.path.exists(path + '/' + i):
                        os.link(self.modules[0].output_dir + '/' + i, path + '/' + i)
                    file_dir = os.path.basename(path)
                    file_path = '{}/{}'.format(file_dir, i)
                    self.api_nipt.add_fastqc(file_path, self.analysis_id)  # fastqc入库
                elif i == name + '_result.txt':
                    self.api_nipt.report_result(main_id, self.modules[0].output_dir + '/' + i, self.analysis_id)
                    self.api_nipt.update_main(main_id, self.modules[0].output_dir + '/' + i)  # 更新zz值到主表中
            self.api_nipt.add_qc(self.modules[0].output_dir + '/' + name + '.qc',
                                 self.modules[0].output_dir + '/' + name + '.gc',
                                 self.modules[0].output_dir + '/new_gc.txt', self.analysis_id)
            self.api_nipt.update_analysis_status(self.analysis_id, "1")
        if len(self.samples) != 0:
            self.api_nipt.add_sample_summary(self.samples)
            self.api_nipt.add_report_summary(self.samples, self.analysis_id)

