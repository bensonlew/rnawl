# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
# last_modify:2017.2.3

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class TfAnalysisModule(Module):
    def __init__(self,work_id):
        super(TfAnalysisModule, self).__init__(work_id)
        self.step.add_steps('tf_analysis', 'family_stastic')
        options = [
            {"name": "query_amino", "type": "infile", "format": "ref_rna.protein_regulation.gene_in_fasta"},  # 上游输入的氨基酸文件（含与差异基因的对应）
            {"name": "diff_gene_id", "type": "string"},
            {"name": "database", "type": "string", "default": "AnimalTFDB"},  # 还有PlantTFDB和AnimalTFDB
            # {"name": "data", "type": "infile", "format": "ref_rna.protein_regulation.txt"},  # 待统计数据
            {"name": "row", "type": "int"}  # 对数据中的第几列进行统计
        ]
        self.add_option(options)
        self.analysis = self.add_tool("ref_rna.protein_regulation.TF_predict")
        self.stastic = self.add_tool("ref_rna.stastic")
        self._end_info = 0

    def check_options(self):
        database_list = ["PlantTFDB", "AnimalTFDB", "iTAK"]
        # if not self.option('query_amino').is_set:
        #     raise OptionError("必须输入氨基酸序列")
        if not self.option('diff_gene_id'):
            raise OptionError("请输入差异基因对应id关系")
        if self.option('database') not in database_list:  # species的判定有问题
            raise OptionError("database选择不正确")
        # if not self.option('data').is_set:
        #     raise OptionError("请确认输入")
        if not self.option('row'):
            raise OptionError('请输入要统计的是第几列')

        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def analysis_run(self):
        self.analysis.set_options({
            "query_amino": self.option("query_amino"),
            "diff_gene_id": self.option("diff_gene_id"),
            "database": self.option("database"),
        })
        self.analysis.on('start', self.set_step, {'start': self.step.tf_analysis})
        self.analysis.on('end', self.set_step, {'end': self.step.tf_analysis})
        self.analysis.on('end', self.set_output, 'analysis')
        self.analysis.run()

    def stastic_run(self):
        data = os.path.join(self.work_dir, "TfPredict/output/TF_result.txt")
        if os.path.isfile(data):
            self.stastic.set_options({
                "data": data,
                "row": self.option('row'),
            })
            self.stastic.on('start', self.set_step, {'start': self.step.family_stastic})
            self.stastic.on('end', self.set_step, {'end': self.step.family_stastic})
            self.stastic.on('end', self.set_output,'stastic')
            self.stastic.run()
        else:
            raise OptionError('未能预测对应转录因子家族')

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
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'analysis' or event['data'] == 'stastic':
            self.linkdir(obj.output_dir, self.output_dir)


    def run(self):
        super(TfAnalysisModule, self).run()
        self.analysis.on('end', self.stastic_run)
        self.stastic.on('end', self.end)
        self.analysis_run()


    def end(self):
        repaths = [
            [".", "", "转录因子分析结果输出目录"],
            ["TF_result.txt", "txt", "分析结果文件信息"],
            ["stastic_result.txt", "txt", "统计结果信息"],
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        # sdir.add_regexp_rules(regexps)
        super(TfAnalysisModule, self).end()
