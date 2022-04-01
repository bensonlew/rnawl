# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
# last modified:201705

"""医学检验所-无创产前筛查流程"""
import time
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import re
from biocluster.config import Config


class NiptWorkflow(Workflow):
    def __init__(self, wsheet_object):
        '''
        :param wsheet_object:
        '''
        self._sheet = wsheet_object
        super(NiptWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "batch_id", "type": "string"},
            {'name': 'member_id', 'type': 'string'},
            {"name": "bw", "type": "int", "default": 10},
            {"name": "bs", "type": "int", "default": 1},
            {"name": "ref_group", "type": "int", "default": 1},
            {"name": "update_info", "type": "string"},
            {"name": "single", "type": "string", "default": "false"},
            {"name": "sanger_type", "type": "string"},  # 判断sanger or tsanger
            {"name": "direct_get_path", "type": "string"}

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tools = []
        self.sample_id = []
        self.main_id = None
        self.api_nipt = self.api.nipt_analysis
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

    def urgency_set(self):
        """
        用于设定加急样本以及样本名字检测！
        :return:
        """
        fastq = []
        fastq1 = os.listdir(self.option('fastq_path').prop['path'])
        # self.logger.info("fastq1:{}".format(fastq1))
        sample_list_urgency, sample_list_normal = self.api.pt_customer.get_ws_urgency(self.option('batch_id'))
        self.logger.info("sample_list_urgency:{}".format(sample_list_urgency))
        if self.api.pt_customer.get_urgency_type(self.option('batch_id'), 'nipt'):
            if len(sample_list_urgency) == 0:
                fastq = fastq1
            else:
                for m in sample_list_urgency:
                    if str("{}_R1.fastq.gz".format(m)) in fastq1:
                        fastq.append("{}_R1.fastq.gz".format(m))
                        fastq.append("{}_R2.fastq.gz".format(m))
                    else:
                        raise Exception("加急样本{}不在ws_dir中，请进行确认！".format(m))
        else:
            fastq = fastq1

        for i in fastq:
            m = re.match('(.*)_R1.fastq.gz', i)
            if m:
                check = self.api_nipt.check_exist_bed(m.group(1))
                if self.option('direct_get_path'):
                    if check:
                        self.logger.info('样本{}已存在于数据库中'.format(m.group(1)))
                    else:
                        self.sample_id.append(m.group(1))
                        self.logger.info('将样本{}添加到待分析队列'.format(m.group(1)))
                else:
                    if check:
                        raise Exception('请检查样本{}是否重名'.format(m.group(1)))
                    else:
                        self.sample_id.append(m.group(1))
                        self.logger.info('将样本{}添加到待分析队列'.format(m.group(1)))

    def analysis_run(self):
        self.urgency_set()
        n = 0
        small_sample = []
        for i in self.sample_id:
            file = os.path.join(self.option('fastq_path').prop['path'], i + "_R1.fastq.gz")
            if os.path.getsize(file) < 1048576:
                small_sample.append(file)
                self.sample_id.remove(i)
                self.api_nipt.add_main_(self.option('member_id'), i, self.option('batch_id'))
        self.logger.info('small_sample：%s' % small_sample)
        self.logger.info("self.sample_id: %s" % self.sample_id)
        if len(self.sample_id) == 0:
            self.logger.error("该批次的样本都分析过了，没有样本可以进行分析！")
            self.api.pt_customer.update_urgency_info(self.option('batch_id'), "nipt")
            self.exit(exitcode=1, data='该批次的样本都分析过了，没有样本可以进行分析！', terminated=False)
        for sample in self.sample_id:
            self.logger.info("sample_id: %s" % sample)
            nipt_analysis = self.add_module("nipt.nipt_analysis")
            self.step.add_steps('nipt_analysis{}'.format(n))
            nipt_analysis.set_options({
                "sample_id": sample,
                "fastq_path": self.option("fastq_path"),
                "bw": self.option('bw'),
                'bs': self.option('bs'),
                'ref_group': self.option('ref_group'),
                "single": self.option("single")
            }
            )
            step = getattr(self.step, 'nipt_analysis{}'.format(n))
            step.start()
            nipt_analysis.on('end', self.finish_update, 'nipt_analysis{}'.format(n))
            self.tools.append(nipt_analysis)
            n += 1

            self.main_id = self.api_nipt.add_main(self.option('member_id'), sample, self.option('batch_id'))
            # for m in self.bin_step:
            #     self.api_nipt.add_interaction(self.main_id, m['bin'], m['step'], self.option('ref_group'), sample)
        if len(self.tools) == 1:
            pass
        else:
            for j in range(len(self.tools)):
                # self.tools[j].on('end', self.set_output, 'nipt_analysis')
                self.tools[j].on('end', self.import_tab)

        if len(self.tools) > 1:
            self.on_rely(self.tools, self.end)
        elif len(self.tools) == 1:
            self.tools[0].on('end', self.end)

        for tool in self.tools:
            tool.run()

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
        for file in file_name:
            if re.match(r'.*\.bed\.2$', file):
                name = file.split(".")[0]
                self.logger.info("接下来要进行导表的样本：%s" % name)
                break
        if name not in self.sample_id:
            raise Exception("该样本名不在上机表中，流程终止")
        main_id = self.api_nipt.get_id(name, self.option('batch_id'))
        for m in self.bin_step:
            interaction_id = self.api_nipt.add_interaction(main_id, m['bin'], m['step'], self.option('ref_group'), name)
            for i in file_name:
                # self.logger.info("输入导库文件{}".format(i))
                if i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_z.xls':
                    self.api_nipt.add_z_result(obj.output_dir + '/' + i, interaction_id)  # 插入不同交互的z值
                elif i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_zz.xls':
                    self.api_nipt.add_zz_result(obj.output_dir + '/' + i, interaction_id)  # 插入不同交互的zz值
                elif i == name + '_result.txt':
                    self.api_nipt.report_result(interaction_id, obj.output_dir + '/' + i, self.option('batch_id'))
            self.api_nipt.update_interaction(name, interaction_id)  # 更新每个交互的状态
        for i in file_name:
            if i == name + '.bed.2' and not self.api_nipt.check_exist_bed(str(name)):
                self.logger.info("数据库中不存在：{},即将进行导表！".format(name))
                self.api_nipt.add_bed_file(obj.output_dir + '/' + i)
            elif re.search(name + '.*_fastqc.html$', i):
                sanger_path = Config().get_netdata_config(self.option('sanger_type'))
                path = sanger_path[self.option('sanger_type') + "_path"] + "/rerewrweset/nipt_fastqc"
                if not os.path.exists(path + '/' + i):
                    os.link(obj.output_dir + '/' + i, path + '/' + i)
                self.api_nipt.add_fastqc(obj.output_dir + '/' + i, self.option('batch_id'))  # fastqc入库
            elif i == name + '_result.txt':
                self.api_nipt.update_main(main_id, obj.output_dir + '/' + i)  # 更新zz值到主表中
        self.api_nipt.add_qc(obj.output_dir + '/' + name + '.qc', obj.output_dir + '/' + name + '.gc',
                             self.option('batch_id'))

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == 'nipt_analysis':
            allfiles = os.listdir(obj.output_dir)
            oldfiles = [os.path.join(obj.output_dir, i) for i in allfiles]
            newfiles = [os.path.join(self.output_dir, i) for i in allfiles]
            for i in range(len(allfiles)):
                os.link(oldfiles[i], newfiles[i])

    def set_single_output(self):
        allfiles = os.listdir(self.work_dir + "/NiptAnalysis/output")
        self.logger.info("allfiles:%s" % allfiles)
        oldfiles = [os.path.join(self.work_dir + "/NiptAnalysis/output", i) for i in allfiles]
        newfiles = [os.path.join(self.output_dir, i) for i in allfiles]
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def run(self):
        self.analysis_run()
        super(NiptWorkflow, self).run()

    def end(self):
        super(NiptWorkflow, self).end()
        if len(self.tools) == 1:
            self.logger.info("只有一个样本")
            self.set_single_output()
            self.logger.info("移动结果文件成功")
            name = self.sample_id[0]
            main_id = self.api_nipt.get_id(name, self.option('batch_id'))
            for m in self.bin_step:
                interaction_id = self.api_nipt.add_interaction(main_id, m['bin'], m['step'],
                                                               self.option('ref_group'), name)
                for i in os.listdir(self.tools[0].output_dir):
                    if i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_z.xls':
                        self.api_nipt.add_z_result(self.tools[0].output_dir + '/' + i, interaction_id)  # 插入不同交互的z值
                    elif i == name + '_' + str(m['bin']) + '_' + str(m['step']) + '_zz.xls':
                        self.api_nipt.add_zz_result(self.tools[0].output_dir + '/' + i, interaction_id)  # 插入不同交互的zz值
                    elif i == name + '_result.txt':
                        self.api_nipt.report_result(interaction_id, self.tools[0].output_dir + '/' + i, self.option('batch_id'))
                        self.api_nipt.update_main(main_id, self.tools[0].output_dir + '/' + i)  # 更新main_id到主表中
                self.api_nipt.update_interaction(name, interaction_id)  # 更新每个交互的状态
            for i in os.listdir(self.tools[0].output_dir):
                if i == name + '.bed.2' and not self.api_nipt.check_exist_bed(str(name)):
                    self.api_nipt.add_bed_file(self.tools[0].output_dir + '/' + i)
                elif re.search(name + '.*_fastqc.html$', i):
                    sanger_path = Config().get_netdata_config(self.option('sanger_type'))
                    path = sanger_path[self.option('sanger_type') + "_path"] + "/rerewrweset/nipt_fastqc"
                    if not os.path.exists(path + '/' + i):
                        os.link(self.tools[0].output_dir + '/' + i, path + '/' + i)
                    self.api_nipt.add_fastqc(self.tools[0].output_dir + '/' + i, self.option('batch_id'))  # fastqc入库
                elif i == name + '_result.txt':
                    self.api_nipt.update_main(main_id, self.tools[0].output_dir + '/' + i)  # 更新zz值到主表中
            self.api_nipt.add_qc(self.tools[0].output_dir + '/' + name + '.qc',
                                 self.tools[0].output_dir + '/' + name + '.gc', self.option('batch_id'))
        self.api.pt_customer.update_urgency_info(self.option('batch_id'), "nipt")
