# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
# from biocluster.api.file.lib.s3 import S3TransferManager
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer


class DnaSplitByBarcodeModule(Module):
    """
    多样性对多个文库进行二次拆分及质控
    author: wangzhaoyue
    last_modify: 2017.12.11
    """

    def __init__(self, work_id):
        super(DnaSplitByBarcodeModule, self).__init__(work_id)
        options = [
            {"name": "lib_path", "type": "infile", "format": "datasplit.path"},  # list,存放文库文件夹对应的路径信息,第一列路径
            {'name': 'ziplevel', 'type': "string", "default": "6"},  # 压缩级别
            {'name': 'combinatorial', "type": "string", "default": "2"},  # Use combinatorial barcode matching
            {'name': "library_info", "type": "infile", "format": "sequence.barcode_info"},  # 文库信息，包含文库中样本酶的信息，文库类型等
        ]
        self.lib_name = {}  # 获取文库对应的路径
        self.R1 = ''
        self.R2 = ''
        self.tools = []
        self.lib_info = defaultdict(list)
        self.sample_info = {}  # PE 样本对应的路径
        self.lib_sample = {}  # PE {文库名：样本名}
        self.PE_sample_name = {}   # PE {样本名：laneID:libID_样本名_R1.fastq.gz}
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('lib_path'):
            raise OptionError('必须输入文库文件夹对应的路径信息')
        if not self.option('library_info'):
            raise OptionError('必须输入文库的样本信息表')
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def axe_demux_run(self):
        n = 0
        for lib in self.lib_name:
            self.axe_demux = self.add_tool('datasplit.axe_demux')
            self.step.add_steps('axe_demux{}'.format(n))
            config = self.work_dir + '/' + lib + '.barcode.config'
            self.R1 = self.lib_name[lib][0]
            self.R2 = self.lib_name[lib][1]
            if self.R1 == '' or self.R2 == '':
                raise Exception('文库路径中的文库{}对应的R1,R2序列不全，请核实！'.format(lib))
            opts = {
                "fq1": self.R1,
                "fq2": self.R2,
                "ziplevel": self.option('ziplevel'),
                "combinatorial": self.option('combinatorial'),
                "library_info": config,
            }
            self.axe_demux.set_options(opts)
            step = getattr(self.step, 'axe_demux{}'.format(n))
            step.start()
            self.step.update()
            self.axe_demux.on('end', self.finish_update, 'axe_demux{}'.format(n))
            self.tools.append(self.axe_demux)
            n += 1
        if n == 0:
            self.set_output()
        elif n == 1:
            self.tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.tools, self.set_output)
        self.step.update()
        for tool in self.tools:
            tool.run()

    def run(self):
        """
        运行
        :return:
        """
        super(DnaSplitByBarcodeModule, self).run()
        self.get_info()
        time.sleep(2)
        self.axe_demux_run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        allfiles = os.listdir(dirpath)
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
                os.link(oldfiles[i], newdir)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(self.sample_info) != 0:
            for sample in self.sample_info.keys():
                if self.sample_info[sample][0].startswith("s3:"):
                    from_path = self.sample_info[sample][0]
                    to_path = self.output_dir + '/' + self.PE_sample_name[sample] + '_R1.fastq.gz'
                    # self.download_from_s3(from_files, to_path, cover=True)
                    # transfer = S3TransferManager(base_path=self.work_dir)
                    # transfer.add(from_uri=from_path, to_uri=to_path)
                    # transfer.wait()
                    download(from_path, to_path)
                    from_path = self.sample_info[sample][1]
                    to_path = self.output_dir + '/' + self.PE_sample_name[sample] + '_R2.fastq.gz'
                    # self.download_from_s3(from_files, to_path, cover=True)
                    # transfer = S3TransferManager(base_path=self.work_dir)
                    # transfer.add(from_uri=from_path, to_uri=to_path)
                    # transfer.wait()
                    download(from_path, to_path)
                else:
                    os.link(self.sample_info[sample][0], self.output_dir + '/' + self.PE_sample_name[sample] + '_R1.fastq.gz')
                    os.link(self.sample_info[sample][1], self.output_dir + '/' + self.PE_sample_name[sample] + '_R2.fastq.gz')
        for tool in self.tools:
            self.linkdir(tool.output_dir, self.output_dir)
        self.logger.info("设置结果目录成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(DnaSplitByBarcodeModule, self).end()

    def get_info(self):
        """
        按照文库，将对应的文库序列及文库信息分开；
        文库类型为RAD及GBS的需要二次拆分，PE不需要跳过
        :return:
        """
        with open(self.option('library_info').prop['path'])as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split('\t')
                if re.search("PE", tmp[4]):
                    self.PE_sample_name[tmp[5]] = tmp[1] + ':' + tmp[3] + ":" + tmp[5]
                    self.sample_info[tmp[5]] = ''
                    self.lib_sample[tmp[3]] = tmp[5]
                else:
                    self.lib_info[tmp[3]].append(line)
        for key in self.lib_info.keys():
            with open(self.work_dir + '/' + key + '.barcode.config', 'w+')as fw:
                fw.write('#RunID\tLaneID\tProjectID\tLibID\tLibType\tSampleID\tSampleNeed\tEnzyme1\tEnzyme2\n')
                for i in self.lib_info[key]:
                    fw.write(i)    # 将文库信息按照文库拆开
        with open(self.option('lib_path').prop['path'])as f:  # 文库对应序列路径
            for line in f:
                tmp = line.strip().split('\t')
                lib = tmp[0]
                if lib in self.lib_sample.keys():  # 如果是PE，单独以样本为键存路径
                    self.sample_info[self.lib_sample[lib]] = tmp[1].strip().split(";")  # 换行符要去掉
                else:
                    if lib in self.lib_info.keys():
                        self.lib_name[lib] = tmp[1].strip().split(";")  # 换行符要去掉
                    else:
                        raise Exception('文库路径中的文库名{}在样本信息表中不存在，请核实！'.format(lib))
