# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:2016.11.25

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError
from mbio.files.rna.diff_stat_table import DiffStatTableFile
# import shutil


class DiffAnalysisModule(Module):
    def __init__(self, work_id):
        super(DiffAnalysisModule, self).__init__(work_id)
        self.step.add_steps('cluster', 'network', 'go_rich', 'kegg_rich', 'go_regulate', 'kegg_regulate', 'cog_class')
        options = [
            {"name": "analysis", "type": "string", "default": "cluster,network,kegg_rich,go_rich,go_regulate,kegg_regulate,cog_class"},  # 选择要做的分析
            {"name": "diff_fpkm", "type": "infile", "format": "rna.express_matrix"},  # 差异基因表达量表
            {"name": "all_list", "type": "infile", "format": "rna.gene_list"},  # 全部基因名称文件
            # {"name": "go_all_list", "type": "infile", "format": "rna.gene_list"},  # go全部注释上的基因文件
            # {"name": "kegg_all_list", "type": "infile", "format": "rna.gene_list"},  # kegg全部注释上的基因文件
            {"name": "diff_list", "type": "infile", "format": "rna.gene_list"},  # 差异基因名称文件
            {"name": "distance_method", "type": "string", "default": "euclidean"},  # 计算距离的算法
            {"name": "log", "type": "int", "default": 10},  # 画热图时对原始表进行取对数处理，底数为10或2
            {"name": "method", "type": "string", "default": "hclust"},  # 聚类方法选择
            {"name": "sub_num", "type": "int", "default": 5},  # 子聚类的数目
            {"name": "softpower", "type": "int", "default": 9},
            {"name": "dissimilarity", "type": "float", "default": 0.25},
            {"name": "module", "type": "float", "default": 0.1},
            {"name": "network", "type": "float", "default": 0.2},
            {"name": "diff_list_dir", "type": "infile", "format": "rna.gene_list_dir"},  # 两两样本/分组的差异基因文件
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "gene_kegg_table", "type": "infile", "format": "annotation.kegg.kegg_table"},  # KEGG的pathway文件
            {"name": "gene_go_list", "type": "infile", "format": "annotation.go.go_list"},  # test
            {"name": "gene_go_level_2", "type": "infile", "format": "annotation.go.level2"},
            {"name": "diff_stat_dir", "type": "infile", "format": "rna.diff_stat_dir"},
            {"name": "cog_table", "type": "infile", "format": "annotation.cog.cog_table"},
            {"name":"is_genelist","type":"bool","default":False}, #是否设置gene_list参数
        ]
        self.add_option(options)
        self.cluster = self.add_tool("rna.cluster")
        self.network = self.add_tool("denovo_rna.express.network")
        self._end_info = 0
        self.kegg_rich_tool = []
        self.cog_class_tool = []
        self.go_rich_tool = []
        self.go_regulate_tool = []
        self.kegg_regulate_tool = []
        self.tools = []

    def check_options(self):
        analysis = self.option('analysis').split(',')
        if ('cluster' or 'network') in analysis and not self.option('diff_fpkm'):
            raise OptionError('缺少输入文件：diff_fpkm差异基因表达量矩阵')
        if 'network' in analysis and not self.option('diff_list'):
            raise OptionError('缺少网络分析的输入文件：diff_list差异基因文件')
        if 'go_rich' in analysis:
            # opts = ['diff_list_dir', 'go_all_list', 'gene_go_list']
            opts = ['diff_list_dir','gene_go_list']
            for i in opts:
                if not self._options[i]:
                    raise OptionError('缺少go富集分析的输入文件:{}'.format(i))
        if 'kegg_rich' in analysis:
            # opts = ['diff_list_dir', 'kegg_all_list', 'gene_kegg_table']
            opts = ['diff_list_dir','gene_kegg_table']
            for i in opts:
                if not self._options[i]:
                    raise OptionError('缺少kegg富集分析的输入文件:{}'.format(i))
        if 'kegg_regulate' in analysis:
            opts = [ 'gene_kegg_table', 'diff_stat_dir']
            for i in opts:
                if not self._options[i]:
                    raise OptionError('缺少kegg调控分析的输入文件:{}'.format(i))
        if 'go_regulate' in analysis:
            opts = ['diff_list_dir', 'gene_go_level_2']
            for i in opts:
                if not self._options[i]:
                    raise OptionError('缺少go调控分析的输入文件:{}'.format(i))
        if self.option("distance_method") not in ("manhattan", "euclidean"):
            raise OptionError("所选距离算法不在提供的范围内")
        if self.option('log') not in (10, 2):
            raise OptionError("所选log底数不在提供的范围内")
        if self.option("method") not in ("hclust", "kmeans", "both"):
            raise OptionError("所选方法不在范围内")
        if not isinstance(self.option("sub_num"), int):
            raise OptionError("子聚类数目必须为整数")
        if self.option("softpower") > 20 or self.option("softpower") < 1:
            raise OptionError("softpower值超出范围")
        if self.option('dissimilarity') > 1 or self.option("dissimilarity") < 0:
            raise OptionError("模块dissimilarity相异值超出范围")
        if self.option('module') > 1 or self.option("module") < 0:
            raise OptionError("模块module相异值超出范围")
        if self.option('network') > 1 or self.option("network") < 0:
            raise OptionError("模块network相异值超出范围")
        if self.option('correct') not in ['BY', 'BH', 'None', 'QVALUE']:
            raise OptionError('多重检验校正的方法不在提供的范围内')
        if 'cluster' or 'network' or 'kegg_rich' or 'go_rich' or 'kegg_regulate' or 'go_regulate'or 'cog_class' in self.option('analysis'):
            pass
        else:
            raise OptionError('没有选择任何分析或者分析类型选择错误：%s' % self.option('analysis'))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def cluster_run(self):
        self.cluster.set_options({
            "diff_fpkm": self.option("diff_fpkm"),
            "log": self.option("log"),
            "method": self.option("method"),
            "sub_num": self.option("sub_num"),
            "is_genelist": True,
            "diff_list_dir": self.option("diff_list_dir")
        })
        self.cluster.on('end', self.set_output, 'cluster')
        self.cluster.on('start', self.set_step, {'start': self.step.cluster})
        self.cluster.on('end', self.set_step, {'end': self.step.cluster})
        self.tools.append(self.cluster)
        # self.cluster.run()

    def network_run(self):
        self.network.set_options({
            "diff_fpkm": self.option("diff_fpkm"),
            "gene_file": self.option("diff_list"),
            "softpower": self.option("softpower"),
            "dissimilarity": self.option("dissimilarity"),
            "module": self.option("module"),
            "network": self.option("network")
        })
        self.network.on('end', self.set_output, 'network')
        self.network.on('start', self.set_step, {'start': self.step.network})
        self.network.on('end', self.set_step, {'end': self.step.network})
        self.tools.append(self.network)
        # self.network.run()

    def kegg_rich_run(self):
        self.step.kegg_rich.start()
        opts = {
            "kegg_table": self.option("gene_kegg_table"),
            "correct": self.option("correct"),
            # "all_list": self.option("kegg_all_list"),
        }
        files = os.listdir(self.option('diff_list_dir').prop['path'])
        for f in files:
            opts.update({"diff_list": os.path.join(self.option('diff_list_dir').prop['path'], f)})
            self.kegg_rich = self.add_tool("denovo_rna.express.kegg_rich")
            self.kegg_rich.set_options(opts)
            self.kegg_rich_tool.append(self.kegg_rich)
        if len(self.kegg_rich_tool) == 1:
            self.kegg_rich.on('end', self.set_output, 'kegg_rich')
            self.kegg_rich.on('end', self.set_step, {'end': self.step.kegg_rich})
        else:
            self.on_rely(self.kegg_rich_tool, self.set_output, 'kegg_rich')
        self.tools += self.kegg_rich_tool
        # for tool in self.kegg_rich_tool:
        #     tool.run()

    def go_rich_run(self):
        self.step.go_rich.start()
        opts = {
            # "all_list": self.option("go_all_list"),
            "go_list": self.option("gene_go_list")
        }
        files = os.listdir(self.option('diff_list_dir').prop['path'])
        for f in files:
            opts.update({"diff_list": os.path.join(self.option('diff_list_dir').prop['path'], f)})
            self.go_rich = self.add_tool("denovo_rna.express.go_enrich")
            self.go_rich.set_options(opts)
            self.go_rich_tool.append(self.go_rich)
        if len(self.go_rich_tool) == 1:
            self.go_rich.on('end', self.set_output, 'go_rich')
            self.go_rich.on('end', self.set_step, {'end': self.step.go_rich})
        else:
            self.on_rely(self.go_rich_tool, self.set_output, 'go_rich')
        self.tools += self.go_rich_tool
        # for tool in self.go_rich_tool:
        #     tool.run()

    def go_regulate_run(self):
        self.step.go_regulate.start()
        files = os.listdir(self.option('diff_stat_dir').prop['path'])
        for f in files:
            statfile = os.path.join(self.option('diff_stat_dir').prop['path'], f)
            tmpfile = DiffStatTableFile()
            tmpfile.set_path(statfile)
            diff_gene, regulate_dict = tmpfile.get_table_info()
            go_genes = self.option("gene_go_level_2").get_gene()
            same_gene = list(set(diff_gene) & set(go_genes))
            self.logger.info(str(len(go_genes)))
            self.logger.info(str(len(diff_gene)))
            if same_gene:
                self.go_regulate = self.add_tool("denovo_rna.express.go_regulate")
                self.go_regulate.set_options({
                    "diff_stat": statfile,
                    "go_level_2": self.option("gene_go_level_2")
                })
                self.go_regulate_tool.append(self.go_regulate)
            else:
                self.logger.info("gene_go_level_2与stat_table中无相同基因")
        if len(self.go_regulate_tool) == 1:
            self.go_regulate.on('end', self.set_output, 'go_regulate')
        elif len(self.go_regulate_tool) == 0:
            pass
        else:
            self.on_rely(self.go_regulate_tool, self.set_output, 'go_regulate')
        self.tools += self.go_regulate_tool
        # for tool in self.go_regulate_tool:
        #     tool.run()

    def kegg_regulate_run(self):
        self.step.kegg_regulate.start()
        self.logger.info(vars(self.option('diff_stat_dir')))
        files = os.listdir(self.option('diff_stat_dir').prop['path'])
        for f in files:
            self.kegg_regulate = self.add_tool("denovo_rna.express.kegg_regulate")
            self.kegg_regulate.set_options({
                "diff_stat": os.path.join(self.option('diff_stat_dir').prop['path'], f),
                "kegg_table": self.option("gene_kegg_table")
            })
            self.kegg_regulate_tool.append(self.kegg_regulate)
        if len(self.kegg_regulate_tool) == 1:
            self.kegg_regulate.on('end', self.set_output, 'kegg_regulate')
        else:
            self.on_rely(self.kegg_regulate_tool, self.set_output, 'kegg_regulate')
        self.tools += self.kegg_regulate_tool
        # for tool in self.go_regulate_tool:
        #     tool.run()

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
            # shutil.rmtree(newdir)
            os.makedirs(newdir)
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
        if event['data'] == 'cluster':
            self.linkdir(obj.output_dir, 'cluster')
            print '..........clust end'
        elif event['data'] == 'network':
            self.linkdir(obj.output_dir, 'network')
            print '..........network end'
        elif event['data'] == 'kegg_rich':
            for tool in self.kegg_rich_tool:
                dirname = 'kegg_rich/' + os.path.splitext(os.path.basename(tool.option('diff_list').path))[0]
                self.linkdir(tool.output_dir, dirname)
            self.set_step(event={'data': {'end': self.step.kegg_rich}})
            print '......ke:%s....kegg_rich end' % len(self.kegg_rich_tool)
        elif event['data'] == 'go_rich':
            for tool in self.go_rich_tool:
                dirname = 'go_rich/' + os.path.splitext(os.path.basename(tool.option('diff_list').path))[0]
                self.linkdir(tool.output_dir, dirname)
            self.set_step(event={'data': {'end': self.step.go_rich}})
            print '.......go:%s...go_rich end' % len(self.go_rich_tool)
        elif event['data'] == 'go_regulate':
            for tool in self.go_regulate_tool:
                dirname = 'go_regulate/' + os.path.splitext(os.path.basename(tool.option('diff_stat').path))[0].strip('_edgr_stat')
                self.linkdir(tool.output_dir, dirname)
            self.set_step(event={'data': {'end': self.step.go_regulate}})
            print '......gr:%s....go_regulate end' % len(self.go_regulate_tool)
        elif event['data'] == 'kegg_regulate':
            for tool in self.kegg_regulate_tool:
                dirname = 'kegg_regulate/' + os.path.splitext(os.path.basename(tool.option('diff_stat').path))[0].strip('_edgr_stat')
                self.linkdir(tool.output_dir, dirname)
            self.set_step(event={'data': {'end': self.step.kegg_regulate}})
            print '.......kr:%s...kegg_regulate end' % len(self.kegg_regulate_tool)
        elif event['data'] == 'cog_class':
            for tool in self.cog_class_tool:
                dirname = 'cog_class/' + os.path.splitext(os.path.basename(tool.option('diff_list').path))[0]
                self.linkdir(tool.output_dir, dirname)
            self.set_step(event={'data': {'end': self.step.cog_class}})
            print '.......kr:%s...cog class end' % len(self.cog_class_tool)
        else:
            pass

    def run(self):
        super(DiffAnalysisModule, self).run()
        analysis = self.option('analysis').split(',')
        if 'cluster' in analysis:
            self.cluster_run()
        if 'network' in analysis:
            self.network_run()
        if 'kegg_rich' in analysis:
            self.kegg_rich_run()
        if 'go_rich' in analysis:
            self.go_rich_run()
        if 'go_regulate' in analysis:
            self.go_regulate_run()
        if 'kegg_regulate' in analysis:
            self.kegg_regulate_run()
        if 'cog_class' in analysis:
            self.cog_class_run()
        if len(self.tools) != 1:
            self.on_rely(self.tools, self.end)
            print '......on rely end'
        else:
            self.tools[0].on('end', self.end)
        self.run_tools()

    def cog_class_run(self):
        self.step.cog_class.start()
        opts = {
            "cog_table": self.option("cog_table")
        }
        files = os.listdir(self.option('diff_list_dir').prop['path'])
        for f in files:
            opts.update({"diff_list": os.path.join(self.option('diff_list_dir').prop['path'], f)})
            self.cog_class = self.add_tool("annotation.cog.cog_class")
            self.cog_class.set_options(opts)
            self.cog_class_tool.append(self.cog_class)
        if len(self.cog_class_tool) == 1:
            self.cog_class_tool[0].on('end', self.set_output, 'cog_class')
        else:
            self.on_rely(self.cog_class_tool, self.set_output, 'cog_class')
        self.tools += self.cog_class_tool

    def run_tools(self):
        analysis = self.option('analysis').split(',')
        if 'cluster' in analysis:
            self.cluster.run()
        if 'network' in analysis:
            self.network.run()
        if 'kegg_rich' in analysis:
            for tool in self.kegg_rich_tool:
                tool.run()
        if 'go_rich' in analysis:
            for tool in self.go_rich_tool:
                tool.run()
        if 'go_regulate' in analysis:
            for tool in self.go_regulate_tool:
                tool.run()
        if 'kegg_regulate' in analysis:
            for tool in self.kegg_regulate_tool:
                tool.run()
        if 'cog_class' in analysis:
            for tool in self.cog_class_tool:
                tool.run()

    def end(self):
        repaths = [
            [".", "", "差异表达模块结果输出目录"],
            ["cluster", "", "cluster分析结果输出目录"],
            ["network", "", "network分析结果输出目录"],
            ["kegg_rich", "", "kegg_rich分析结果输出目录"],
            ["go_rich", "", "go_rich分析结果输出目录"],
            ["go_regulate", "", "go_regulate分析结果输出目录"],
            ["kegg_regulate", "", "kegg_regulate分析结果输出目录"],
            ["network/all_edges.txt", "txt", "edges结果信息"],
            ["network/all_nodes.txt ", "txt", "nodes结果信息"],
            ["network/removeGene.xls ", "xls", "移除的基因信息"],
            ["network/removeSample.xls ", "xls", "移除的样本信息"],
            ["network/softPower.pdf", "pdf", "softpower相关信息"],
            ["network/ModuleTree.pdf", "pdf", "ModuleTree图"],
            ["network/eigengeneClustering.pdf", "pdf", "eigengeneClustering图"],
            ["network/eigenGeneHeatmap.pdf", "pdf", "eigenGeneHeatmap图"],
            ["network/networkHeatmap.pdf", "pdf", "networkHeatmap图"],
            ["network/sampleClustering.pdf", "pdf", "sampleClustering图"],
            ["cluster/hclust", "", "hclust分析结果输出目录"],
            ["cluster/hclust/hc_gene_order", "txt", "按基因聚类的基因排序列表"],
            ["cluster/hclust/hc_sample_order", "txt", "按样本聚类的样本排序列表"],
            ["cluster/hclust/hclust_heatmap.xls", "xls", "层级聚类热图数据"],
            ["cluster/kmeans/kmeans_heatmap.xls", "xls", "层级聚类热图数据"]
        ]
        regexps = [
            [r"subcluster_", "xls", "子聚类热图数据"],
            [r"network/CytoscapeInput.*", "txt", "Cytoscape作图数据"],
            [r"go_regulate/.*?go_regulate_stat.xls$", "xls", "基因上下调在GO2level层级分布情况表"],
            [r"go_rich/go_enrich_.*", "xls", "go富集结果文件"],
            [r"go_rich/go_lineage.*", "png", "go富集有向无环图"],
            [r"kegg_rich/.*?kegg_enrichment.xls$", "xls", "kegg富集分析结果"],
            [r"kegg_regulate/.*?kegg_regulate_stat.xls$", "xls", "kegg调控统计表"],
            [r"kegg_regulate/.*?pathways$", "dir", "kegg调控统计pathway结果图片"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(DiffAnalysisModule, self).end()
