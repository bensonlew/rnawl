# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from bson.objectid import ObjectId
import json
import os
import gevent
from collections import defaultdict


class ArghubWorkflow(Workflow):
    def __init__(self, wsheet):
        self._sheet = wsheet
        super(ArghubWorkflow, self).__init__(wsheet)
        options = [
            {'name': 'rrna', 'type': 'infile', 'format': 'sequence.fasta'},  # rrna 输入文件
            {'name': 'cds', 'type': 'infile', 'format': 'sequence.fasta'},  # cds 输入文件
            {'name': 'prot', 'type': 'infile', 'format': 'sequence.fasta'},  # protin 输入文件
            {'name': 'read1', 'type': 'infile', 'format': 'sequence.fastq'},  # reads 输入文件 read1
            {'name': 'read2', 'type': 'infile', 'format': 'sequence.fastq'},  # reads 输入文件 read2
            {'name': 'single', 'type': 'infile', 'format': 'sequence.fasta'},  # genome输入文件
            {'name': 'meta', 'type': 'infile','format': 'sequence.fasta'},  # genome输入文件
            {'name': 'input_type', 'type': 'string'},  # rrna cds prot read single meta
            {'name': 'gff', 'type': 'infile', 'format': 'gene_structure.gff3'},  # gff 文件 genome输入的可选输入文件
            {'name': 'sequence', 'type': 'string', },
            {'name': 'aligner', 'type': 'string', 'default': 'blast'},  # 使用的比对工具 blast, diamond, hmmer
            {'name': 'evalue', 'type': 'float', 'default': 1e-5},  # 比对设置的 e 值
            {'name': 'identity', 'type': 'float', 'default': 50},  # 比对identity阈值, 输入工具为hmmer时无此项
            {'name': 'db_type', 'type': 'string', 'default': 'plus'},  # core plus
            {'name': 'main_id', 'type': 'string',},
            {'name': 'update_info', 'type': 'string',}
        ]
        self.add_option(options)
        self.need_predict = False
        self.arghubs = []  # 抗性基因分析对象
        self.pres = []  # 分析前数据处理对象
        self.mge_depends = []
        self.out_tasks = []
        self.outfiles = {}  # 输出文件字典，用于导表
        self.set_options(wsheet.options())

        #self.IMPORT_REPORT_DATA = True
        #self.IMPORT_REPORT_AFTER_END = False

    def check_options(self):
        if not self.option('input_type'):
            raise OptionError('设置输入数据类型')
        elif self.option('input_type') == 'read' and not (
                self.option('read1').is_set and self.option('read2').is_set):
            raise OptionError('input_type 为 read时必须设置 read1 和 read2 参数')
        elif self.option('input_type') != 'read' and not self.option(
                self.option('input_type')).is_set:
            if self.option('sequence'):
                self.logger.info(os.path.join(self.work_dir, "input.txt"))
                with open(os.path.join(self.work_dir, "input.txt"), 'w') as w:
                    w.write(">{}\n".format(self.option('input_type')))
                    w.write("{}\n".format(self.option("sequence")))
                self.option(self.option('input_type'), os.path.join(self.work_dir, "input.txt"))
            else:
                raise OptionError('input_type 为 {0} 时必须设置 {0} 参数'.format(
                    self.option('input_type')))
        if self.option('input_type') in ['single', 'meta'
                                         ] and not self.option('gff').is_set:
            self.need_predict = True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def all_steps(self):
        analysis_type = "analysis_" + self.option('input_type')
        if self.option('input_type') == 'read':
            self.step.add_steps(analysis_type, "set_output")
        elif self.option('input_type') in ['single', 'meta']:
            self.step.add_steps(analysis_type, "input_prepare", "arg_predict", "mge_predict", "set_output")
        else:
            self.step.add_steps(analysis_type, "input_prepare", "arg_predict", "set_output")
        self.on('start', self.set_step, {'start': self.step._steps[analysis_type]})
        self.on('end', self.set_step, {'end': self.step._steps[analysis_type]})

    def run(self):
        '''
        运行入口，设置运行逻辑
        '''
        #task_info = self.api.api('task_info.arghub_task_info')
        #task_info.add_task_info()
        # 非read输入进行拆分
        self.all_steps()
        self.outfmt = self.add_module('arghub.outfmt')
        if self.option('input_type') == 'read':
            arghub = self.add_module('arghub.arghub')
            self.arghubs.append(arghub)
            arghub.on('end', self.set_output)
            self.set_one(arghub)
            arghub.run()
        else:
            self.split = self.add_tool('sequence.split_fasta')
            self.split.on('start', self.set_step, {'start': self.step.input_prepare})
            self.split.on('end', self.set_step, {'end': self.step.input_prepare})
            self.split.on('end', self.multi_arghubs)
            gevent.spawn_later(5, self.split_job)
        super(ArghubWorkflow, self).run()

    def split_job(self):
        '''
        对输入进行拆分
        基因组文件: 每5000条scaffold or contig进行拆分
        其它: 每10000条序列进行拆分
        '''
        opts = {
            'fasta': self.option(self.option('input_type')).prop['path'],
            # 'lines': lines,
        }
        self.split.set_options(opts)
        self.logger.info('####开始拆分序列####')
        self.split.run()

    def multi_arghubs(self):
        '''
        拆分后的序列作为输入并行运行
        '''
        self.logger.info('####并行运行arghub####')
        self.step.arg_predict.start()
        self.step.update()
        for seq_f in os.listdir(self.split.output_dir):
            if seq_f.endswith('.fai'):
                continue
            arghub = self.add_module('arghub.arghub')
            self.arghubs.append(arghub)
            seq_f = os.path.join(self.split.output_dir, seq_f)
            self.set_one(arghub, seq_f)

        mge = []
        if self.mge_depends:
            self.mge = self.add_module('arghub.arg_mge')
            self.on_rely(self.mge_depends, self.run_mge)  # mge 模块 self.mge
            self.on_rely(self.arghubs + [ self.mge,], self.run_outfmt)
            self.on_rely(self.mge, self.set_step, {'end': self.step.mge_predict})
        else:
            self.on_rely(self.arghubs, self.run_outfmt)
        self.on_rely(self.arghubs, self.set_step, {'end': self.step.arg_predict})

        if self.pres:
            map(lambda task: task.run(), self.pres)
        else:
            self.logger.info('##### in ELSE #####')
            map(lambda task: task.run(), self.arghubs)

    def set_one(self, arghub, seq=None):
        '''
        每一个拆分单元运行抗性基因分析过程
        1. 分析前处理: a: 无gff文件且input_type = single or meta 进行基因预测 module: arghub.predict
                      b: 有gff文件input_type = single or meta, 或者 input_type = cds, 进行序列处理 tool: arghub.input_prepare
        input_type = prot, read, rrna 时，不进行步骤1
        2. 抗性基因分析: a: prot, rrna输入，通过module: arghub完成
                        b: reads输入，通过tool: ariba完成
        other: 如果input_type = single or meta, 步骤1 完成后进行文件格式调整， 作为可移动元件的输入
        '''

        # 分析前处理
        if self.need_predict:  # 基因预测
            pre = self.add_module('arghub.predict')
            pre.set_options({'genome': seq, 'mode': self.option('input_type')})
        elif self.option('gff').is_set or self.option('input_type') == 'cds':  # 从序列获取蛋白序列
            pre = self.add_tool('arghub.input_prepare')
            if self.option('gff').is_set:
                opts = {'genome': seq, 'gff': self.option('gff').prop['path']}
            else:
                opts = {'cds': seq}
            self.logger.info("input_prepare opts:{}".format(opts))
            pre.set_options(opts)

        # 配置分析参数
        # common arghub opts
        opts = {
            'aligner': self.option('aligner'),
            'evalue': self.option('evalue'),
            'identity': self.option('identity'),
            'db_type': self.option('db_type')
        }
        if self.option('input_type') == 'read':
            self.logger.info('### {} ###'.format(self.option('read1').prop))
            opts['read1'] = self.option('read1').prop['path']
            opts['read2'] = self.option('read2').prop['path']
            arghub.set_options(opts)
        elif self.option('input_type') in ['prot', 'rrna']:
            opts[self.option('input_type')] = seq
            arghub.set_options(opts)
        else:
            pre.on('end', self.run_arghub, data=(arghub, opts))
            self.pres.append(pre)
            self.logger.info('pre = {}\nself.pres = {}'.format(pre, self.pres))
        # 额外的，如果输入为基因组类型的，就进行可移动分析
        if self.option('input_type') in ['single', 'meta']:
            uni = self.add_tool('arghub.uni_format')
            self.mge_depends.append(uni)
            pre.on('end', self.run_uni, data=uni)

    def run_uni(self, event):
        bind_obj = event['bind_object']
        uni = event['data']
        opts = {}
        self.logger.info('run uni {}'.format(bind_obj))
        for k in [
                'gff', 'prot', 'cds', 'cds_o', 'trna', 'rrna', 't_gff', 'r_gff'
        ]:
            if k in bind_obj.get_option_object() and bind_obj.option(k).is_set:
                self.logger.info('run uni' + k)
                opt = k.split('_o')[0]
                opts[opt] = bind_obj.option(k).prop['path']
        if self.need_predict:
            opts['predict'] = 1
        uni.set_options(opts)
        uni.run()

    def run_arghub(self, event):
        '''
        分析前处理步骤完成后，触发此方法
        '''
        bind_obj = event['bind_object']
        arghub, opts = event['data']
        for f in ['prot', 'rrna']:
            if bind_obj.option(f).is_set:
                opts[f] = bind_obj.option(f).prop['path']
        arghub.set_options(opts)
        arghub.run()

    def run_argout(self):
        '''
        分析结果文件结合数据库注释文件，输出最终的结果文件表
        '''
        self.logger.info('### start outformat ###')
        raw_output_list = os.path.join(self.work_dir,
                                       'arghub_raw_output_list.txt')
        with open(raw_output_list, 'w') as o:
            [o.write(m.option('output') + '\n') for m in self.arghubs]
        opts = {
            'file_list': raw_output_list,
            'file_type': self.option('input_type'),
        }
        self.argout.set_options(opts)
        self.argout.run()

    def run_mge(self):
        self.step.mge_predict.start()
        self.step.update()
        inputs_for_mge = os.path.join(self.work_dir, 'inputs_for_mge.txt')
        name_file_list = os.path.join(self.work_dir,
                                      'name_file_list.txt')  # 改名文件
        opts = {}
        opts['genome'] = self.option(self.option('input_type')).prop['path']
        opts['file_for_mge'] = inputs_for_mge
        for_mges = defaultdict(set)
        for_names = []
        for uni in self.mge_depends:
            self.logger.info('run mge uni: {}'.format(uni))
            if self.need_predict:
                for_names.append(uni.option('name_file').path + '\n')
                for_mges['trna'].add(uni.option('trna_o').path)
            #for opt in ['gff_o', 'prot_o', 'cds_o', 'rnt']:
            for opt in ['gff_o', 'prot_o', 'cds_o', 'rnt', 'rrna_o']:
                self.logger.info('run mge {} {}'.format(
                    opt,
                    uni.option(opt).path))
                k = opt.split('_o')[0]
                for_mges[k].add(uni.option(opt).path)
        if self.need_predict:
            opts['name_file_list'] = name_file_list
            with open(name_file_list, 'w') as w2:
                [w2.write(l) for l in for_names]
        else:
            for pre in self.pres:
                for_mges['trna'].add(pre.option('trna').path)

        with open(inputs_for_mge, 'w') as w:
            [
                w.write('{}\t{}\n'.format(k, i)) for k, v in for_mges.items()
                for i in v if i
            ]
        self.mge.set_options(opts)
        self.logger.info("mge opts: {}".format(opts))
        self.mge.run()

    def run_outfmt(self):
        '''
        抗性基因分析和可移动元件结果整合
        '''
        self.step.set_output.start()
        self.step.update()
        arg_output_list = os.path.join(self.work_dir,
                                       'arg_raw_output_list.txt')
        with open(arg_output_list, 'w') as o:
            [o.write(m.option('output') + '\n') for m in self.arghubs]
        opts = {
            'argout_list': arg_output_list,
            'arg_type': self.option('input_type'),
        }
        if self.option('input_type') in ['single', 'meta']:
            genome = self.option(self.option('input_type')).path
            opts.update({
                'cds': self.mge.option('cds'),
                'prot': self.mge.option('prot'),
                'trna': self.mge.option('trna'),
                'rrna': self.mge.option('rrna'),
                'rnt': self.mge.option('rnt'),
                'mge_result': self.mge.option('mge_result'),
                'mge_elem': self.mge.option('mge_elem'),
                'sample': genome.split('/')[-1].rpartition('.')[0]
            })
        if self.option('input_type') == 'prot':
            opts['prot'] = self.option('prot')
        elif self.option('input_type') == 'cds':
            #prot_list = os.path.join(self.work_dir, "prot_list.txt")
            #with open(prot_list, 'w') as w:
            #    [w.write(m.option("prot").path + '\n') for m in self.pres]
            #opts['prot_list'] = prot_list
            all_prot = os.path.join(self.work_dir, "all_prot.faa")
            if os.path.exists(all_prot):
                os.remove(all_prot)
            for m in self.pres:
                cmd = "cat {} >> {}".format(m.option("prot").path, all_prot)
                os.system(cmd)
            opts['prot'] = all_prot
            opts['cds'] = self.option('cds')
        elif self.option('input_type') == 'rrna':
            opts['rrna'] = self.option('rrna')
        if self.need_predict:
            opts['rename'] = self.mge.option('rename')
        self.outfmt.on('end', self.set_output)
        self.logger.info(opts)
        self.outfmt.set_options(opts)
        self.logger.info("outfmt opts: {}".format(opts))
        self.outfmt.run()

    def set_output(self):
        self.logger.info('### copy arg results to OUTPUT_DIR###')
        in_type = self.option('input_type')
        if self.option('input_type') == 'read':
            self.link(self.arghubs[0].option('output'))
            self.outfiles['arg_n'] = self.arghubs[0].option('output')
        else:
            for name, option in self.outfmt.get_option_object().items():
                if option.type == 'outfile' and option.is_set:
                    self.link(option.value.path)
                    self.outfiles[name] = option.value.path
        self.end()

    def end(self):
        '''
        任务结束，结束前对结果进行导表和结果文件上传
        '''
        self.logger.info('### ending stuff ###')
        self.send_files()
        self.run_api()
        self.step.set_output.finish()
        self.step.update()
        super(ArghubWorkflow, self).end()

    def run_api(self):
        '''
        导表
        '''
        self.logger.info('### 开始导表 ###')
        api = self.api.api('arghub.arghub_analysis')
        #params = {}
        #for name, option in self.get_option_object().items():
        #    if option.type in ['infile']:
        #        if option.is_set:
        #            params[name] = option.value.path
        #    elif option.value:
        #        params[name] = option.value
        #self.logger.info(params)

        main_id = self.option("main_id")
        api.update_main(self.option("main_id"))

        for name, flpath in self.outfiles.items():
            if name == 'arg_n':
                api.add_result(flpath, main_id, 'analysis_arg_detail')
            elif name == 'stat':
                #api.add_result(flpath, main_id, 'analysis_stat')
                pass
            else:
                api.add_result(flpath, main_id, 'analysis_mge_detail')
        self.logger.info('### 导表结束 ###')

    def send_files(self):
        relpath = [
            ['.', '', 'arghub分析结果文件夹', 0, ''],
        ]
        regexps = []
        send_dir = self.add_upload_dir(self.output_dir)
        send_dir.add_relpath_rules(relpath)
        send_dir.add_regexp_rules(regexps)


if __name__ == '__main__':
    from biocluster.wpm.client import worker_client
    from arghub import ArghubWorkflow
    from biocluster.wsheet import Sheet
    import time
    data = {
        'type': 'workflow',
        'name': 'arghub.arghub',
        'client': 'client03',
        'member_type': '123',
        'cmd_id': '123',
        'rerun': True,
    }
    test_dict = {
        #'arghub_test1': {  # 1 预测后进行分析
        #    'single':
        #    '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCF_004119875.1.fna',
        #    'input_type': 'single'
        #},
        #'arghub_test1-1': {  # 1 预测后进行分析
        #    'meta':
        #    '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCF_004119875.1.fna',
        #    'input_type': 'meta'
        #},
        'arghub_test1-2': {  # 1 预测后进行分析
            'cds':
            '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCF_004119875.1.fna.cds.fna',
            'input_type': 'cds'
        },
        # 'arghub_test1.5': {  # 1 预测后进行分析
        #    'single': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCA_000279975.1_ASM27997v1_genomic.fna',
        #    'input_type': 'single'
        # },
        # 'arghub_test1.6': {  # 1 预测后进行分析
        #    'single': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCA_001510735.1_ASM151073v1_genomic.fna',
        #    'input_type': 'single'
        # },
        #'arghub_test2': {  # 2 处理获取蛋白文件后进行分析
        #    'single':
        #    '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCF_004119875.1.fna',
        #    'gff':
        #    '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/NZ_CP035414.1.gff',
        #    'input_type': 'single'
        #},
        # 'arghub_test3-1': {  # 1 不需要前面的处理，进入分析
        #    'input_type': 'prot',
        #    'prot': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCF_004119875.1.fna.prot.faa',
        # },
        # 'arghub_test4-1': {  # 4 预测后进行分析
        #    'meta': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/LD200518_0016.contig.fa',
        #    'input_type': 'meta',
        # },
        # 'arghub_test5-1': {
        #    'read1': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/two.1.fq',
        #    'read2': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/two.2.fq',
        #    'input_type': 'read'
        # },
        # 'arghub_test6': {
        #    'read1': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/r1.fq',
        #    'read2': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/toolapp/r2.fq',
        #    'input_type': 'read'
        # },
        # 'arghub_test7-1': {  # 1 不需要前面的处理，进入分析
        #    'input_type': 'prot',
        #    'aligner': 'hmmscan',
        #    'prot': '/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/arghub/new/coding/test/input/GCF_004119875.1.fna.prot.faa',
        # },
    }
    for i, opts in test_dict.items():
        data['id'] = str(i)
        data['options'] = opts
        # submit
        worker = worker_client()
        result = worker.add_task(data)
        print result
        time.sleep(10)
        # local
        # wsheet = Sheet(data=data)
        # ArghubWorkflow(wsheet).run()
