# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify:2016.09.22

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError
from biocluster.config import Config

class CorrNetworkAnalysisModule(Module):
    def __init__(self, work_id):
        super(CorrNetworkAnalysisModule, self).__init__(work_id)
        self.step.add_steps('otu_association', 'corr_network_calc')
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "lable", "type": "float", "default": 0.03},  # 设定lable（距离水平）的值
            {"name": "method", "type": "string", "default": "pearson"}, # 设定计算相似性的算法
            {"name": "coefficient", "type": "float", "default": 0.08},  # 相关性系数阈值
            {"name": "significance", "type": "float", "default": 0.05},
            {"name": "group_name", "type": "string"}, # 输出文件夹名称
        ]
        self.add_option(options)
        self.otu_association = self.add_tool("meta.otu.otu_association")
        self.corr_network_calc = self.add_tool("meta.otu.corr_network_calc")
        self.heatmap_tool = self.add_tool("meta.association_model.corr_tree")
        self._end_info = 0

    def check_options(self):
        if not self.option("otutable").is_set:
            raise OptionError("必须要提供otutable文件", code="22700601")
        self.option('otutable').get_info()
        if self.option('otutable').prop['sample_num'] <= 2:
            raise OptionError('otu表的样本数目少于3，不可进行网络分析', code="22700602")
        # if self.option('grouptable').is_set:
        #     if self.option('grouptable').prop['sample_number'] <= 2:
        #         raise OptionError('分组表中的样本数目少于3，不可进行网络分析')
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def otuassociation_run(self):
        self.otu_association.set_options({
            "otutable": self.option('otutable').prop["path"],
            "lable": self.option("lable"),
            "method": self.option("method")
        })
        self.otu_association.on('end', self.set_output, 'otu_association')
        self.otu_association.on('start', self.set_step, {'start': self.step.otu_association})
        self.otu_association.on('end', self.set_step, {'end': self.step.otu_association})
        self.otu_association.run()

    def corrnetworkcalc_run(self):
        lable = self.option("lable")
        method = self.option("method")
        # corr_name = "shared."+ str(lable) + "." + str(method) +".otu.corr"
        corr_name = "shared." + str(lable) + "." + str(method) + ".corr"  # modified by hongdongxuan 20170324
        corr_table = os.path.join(self.work_dir, "OtuAssociation/output/", corr_name)
        print corr_table
        self.corr_network_calc.set_options({
            "corr_table": corr_table,
            "coefficient": self.option("coefficient"),
            "significance": self.option("significance")
        })
        self.corr_network_calc.on('end', self.set_output, 'corr_network_calc')
        self.corr_network_calc.run()

    def heatmap_run(self):
        self.heatmap_tool.set_options({
            "exp" : self.option('otutable').prop["path"],
            "pvalue" : self.option("significance"),
            "corr_method": self.option("method"), ##add by qingchen.zhang@20200116
        })
        self.heatmap_tool.on('end',self.set_output,'heatmap')
        self.heatmap_tool.run()


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
        if self.option('group_name'):
            group_name = self.option('group_name')+"/"
            if not os.path.exists(self.output_dir+"/"+group_name):
                os.mkdir(self.output_dir+"/"+group_name)
        else:
            group_name = ""

        if event['data'] == 'otu_association':
            self.linkdir(obj.output_dir, group_name + 'otu_association')
        elif event['data'] == 'corr_network_calc':
            self.linkdir(obj.output_dir, group_name + 'corr_network_calc')
        elif event['data'] == 'heatmap':
            self.linkdir(obj.output_dir, group_name + 'heatmap')
        else:
            pass

    def run(self):
        super(CorrNetworkAnalysisModule, self).run()
        self.otu_association.on('end', self.corrnetworkcalc_run)
        self.on_rely([self.corr_network_calc, self.heatmap_tool], self.end)
        #self.corr_network_calc.on('end', self.end)
        self.otuassociation_run()
        self.heatmap_run()


    def end(self):
        repaths = [
            [".", "", "物种相关性网络结果输出目录"],
            ["shared.txt", "txt", "shared文件"],
            ["corr_network_attributes.txt", "txt", "网络的单值属性表"],
            ["corr_network_by_cut.txt", "txt", "相关系数筛选后网络边文件"],
            ["corr_network_centrality.txt", "txt", "网络节点的中心系数表"],
            ["corr_network_clustering.txt", "txt", "网络节点的聚类系数表"],
            ["corr_network_degree_distribution.txt", "txt", "网络节点的度分布表"],
            ["corr_network_node_degree.txt", "txt", "网络节点的度统计表"]
        ]
        regexps = [
            [r"*\.otu\.corr", "corr", "物种相似性网络边文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CorrNetworkAnalysisModule, self).end()
