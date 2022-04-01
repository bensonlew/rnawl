# -*- coding: utf-8 -*-
# __author__ = 'ghd'


from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class HmdbModelModule(Module):
    def __init__(self, work_id):
        super(HmdbModelModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fastq,sequence.fastq_dir,sequence.profile_table"},
            {"name": "platform", "type": "string", "default": "16s"},
            {"name": "model", "type": "string", "default": "crc"},
            {"name": "model_type", "type": "string", "default": "randomforest"}
        ]
        self.add_option(options)
        # self.step.add_steps("map", "tax", "model")
        self.map = self.add_tool("hmdb.map")
        self.tax = self.add_tool("hmdb.add_tax")
        self.model = self.add_tool("hmdb.model")

    def check_options(self):
        type_hash = {
            'crc': ["ada", "bayes", "gradient", "lr", "randomforest", "tree"],
            't2d': ["ada", "gradient", "randomforest", "tree"],
            'obesity': ["ada", "bayes", "gradient", "randomforest", "svm", "tree"]
        }
        avaiable_paltform = {
            "crc": ["16s"],
            "t2d": ["shotgun"],
            "obesity": ["16s"]
        }
        if self.option("platform") not in ["16s", "shotgun"]:
            raise OptionError("选择的平台不正确：%s" % self.option("platform"))
        if self.option("model") not in ["crc", "t2d", "obesity"]:
            raise OptionError("选择的模型不正确：%s" % self.option("model"))
        if self.option("platform") not in avaiable_paltform[self.option("model")]:
            raise OptionError("%s的预测模型不可做%s平台的数据，允许的平台：%s" % (
                self.option('model'), self.option('platform'), avaiable_paltform[self.option('model')]))
        if self.option("model_type") not in type_hash[self.option("model")]:
            raise OptionError("选择的模型方法不正确：%s\n可用的方法：%s" % (self.option("model_type"), type_hash[self.option("model")]))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_map(self):
        opts = {}
        if self.option('query').format == "sequence.fastq":
            opts['query'] = self.option('query')
        elif self.option('query').format == 'sequence.fastq_dir':
            opts['query_dir'] = self.option('query')
        if self.option('platform') == "16s":
            opts['method'] = "local"
        elif self.option("platform") == "shotgun":
            opts['method'] = "global"
        self.map.set_options(opts)
        self.map.on("end", self.run_tax)
        self.map.run()

    def run_tax(self):
        opts = {}
        if self.option("query").format == "sequence.profile_table":
            opts['abund'] = self.option("query")
        else:
            opts["abund"] = self.map.option('abund')
        self.tax.set_options(opts)
        self.tax.on("end", self.run_model)
        self.tax.run()

    def run_model(self):
        self.model.set_options({
            'model': self.option("model"),
            'model_type': self.option('model_type'),
            'train_file': self.tax.option('new_abund')
        })
        self.output_dir = self.model.output_dir
        self.model.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            # [".", "", "物种组成分析结果目录"],
            ['model_result.xls', 'xls', '模型预测结果表']
        ])
        super(HmdbModelModule, self).end()

    def run(self):
        super(HmdbModelModule, self).run()
        if self.option("query").format in ["sequence.fastq", "sequence.fastq_dir"]:
            self.run_map()
        elif self.option("query").format == "sequence.profile_table":
            self.run_tax()
        self.model.on("end", self.end)
