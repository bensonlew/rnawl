# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last_modify:2017.06.13

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class CommonAnnoStatModule(Module):
    """
    宏基因组流程必做注释nr\kegg\cog注释的结果统计模块
    """
    def __init__(self, work_id):
        super(CommonAnnoStatModule, self).__init__(work_id)
        options = [
            {"name": "nr_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # nr比对结果文件夹
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 样本序列丰度表
            {"name": "tax_level_dir", "type": "outfile", "format": "meta_genomic.nr_dir"},  # nr注释输出结果
            {"name": "string_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # string比对结果文件夹
            {"name": "cog_result_dir", "type": "outfile", "format": "meta_genomic.cog_dir"},  # cog统计结果文件夹
            {"name": "kegg_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # kegg比对结果文件夹
            {"name": "kegg_result_dir", "type": "outfile", "format": "meta_genomic.kegg_dir"}  # kegg统计结果文件夹
        ]
        self.add_option(options)
        self.step.add_steps('nr_ncbi', 'nr_stat', 'cat', 'nr_level', 'string2cog')
        self.ncbi_tools = []
        self.nr_stat_tools = []
        self.string2cog_tools = []
        self.kegg_anno_tools = []
        self.cat_kegg_tools = []
        self.kegg_anno = 0
        self.ncbi_number = 0
        self.nrstat_number = 0
        self.string2cog = 0
        self.cat_nr = self.add_tool('align.cat_hmmscanout')
        self.cat_cog = self.add_tool('align.cat_hmmscanout')
        self.nr_level = self.add_tool('meta_genomic.nr_tax_level')
        self.cog_stat = self.add_tool('meta_genomic.metagen_cog_stat')
        self.kegg_stat = self.add_tool('meta_genomic.metagen_kegg_stat')
        self.sum_tool =[]

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table")
        if not self.option("nr_xml_dir").is_set:
            if not self.option("kegg_xml_dir").is_set:
                if not self.option("string_xml_dir").is_set:
                    raise OptionError("必须输入一种xml_dir文件夹(nr/string/kegg)")
        return True

    def run(self):
        super(CommonAnnoStatModule, self).run()
        if self.option("nr_xml_dir").is_set:
            self.run_nr_ncbi()
            self.sum_tool.append(self.nr_level)
        if self.option("string_xml_dir").is_set:
            self.run_string2cog()
            self.sum_tool.append(self.cog_stat)
        if self.option("kegg_xml_dir").is_set:
            self.run_kegg_anno()
            self.sum_tool.append(self.kegg_stat)
        if len(self.sum_tool) == 1:
            self.sum_tool[0].on('end', self.end)
        else:
            self.on_rely(self.sum_tool, self.end)

    def run_nr_ncbi(self):
        xml_file = os.listdir(self.option('nr_xml_dir').prop['path'])
        for i in xml_file:
            file_path = os.path.join(self.option('nr_xml_dir').prop['path'], i)
            ncbi = self.add_tool("taxon.ncbi_taxon")
            ncbi.set_options({
                "blastout": file_path,
                "blastdb": "nr",
            })
            ncbi.on('end', self.set_output, "nr_ncbi")
            self.ncbi_tools.append(ncbi)
        if len(self.ncbi_tools) == 1:
            self.ncbi_tools[0].on('end', self.run_nr_stat)
        else:
            self.on_rely(self.ncbi_tools, self.run_nr_stat)
        for tool in self.ncbi_tools:
            tool.run()

    def run_nr_stat(self):
        file_dir = os.path.join(self.work_dir, 'tmp_ncbi_result')
        file_name = os.listdir(file_dir)
        for i in file_name:
            file_path = os.path.join(file_dir, i)
            nr_stat = self.add_tool("meta_genomic.metagen_nr")
            nr_stat.set_options({
                "taxon_out": file_path,
                "reads_profile_table": self.option('reads_profile_table'),
            })
            nr_stat.on('end', self.set_output, "nr_stat")
            self.nr_stat_tools.append(nr_stat)
        if len(self.nr_stat_tools) == 1:
            self.nr_stat_tools[0].on('end', self.run_nr_tax_level)
        else:
            self.on_rely(self.nr_stat_tools, self.run_cat_nr_stat)
        for tool in self.nr_stat_tools:
            tool.run()

    def run_cat_nr_stat(self):
        self.cat_nr.set_options({
            'hmmscan_out': self.work_dir + '/tmp_nrstat_result',
        })
        self.cat_nr.on('end', self.run_nr_tax_level)
        self.cat_nr.run()

    def run_nr_tax_level(self):
        os.mkdir(self.output_dir + '/nr_tax_level')
        self.nr_level.set_options({
            'nr_taxon_profile': self.cat_nr.option('hmmscan_result').prop['path'],
        })
        self.nr_level.on('end', self.set_output, 'nr_level')
        # self.nr_level.on('end', self.end)  # 单独测试nr的时候使用
        self.nr_level.run()

    def run_string2cog(self):
        xml_file = os.listdir(self.option('string_xml_dir').prop['path'])
        for i in xml_file:
            file_path = os.path.join(self.option('string_xml_dir').prop['path'], i)
            string2cog = self.add_tool("meta_genomic.cog_anno")
            string2cog.set_options({
                "string_xml": file_path
            })
            string2cog.on('end', self.set_output, "string2cog")
            self.string2cog_tools.append(string2cog)
        if len(self.string2cog_tools) == 1:
            self.string2cog_tools[0].on('end', self.run_cog_stat)
        else:
            self.on_rely(self.string2cog_tools, self.run_cat_string2cog)
        for tool in self.string2cog_tools:
            tool.run()

    def run_cat_string2cog(self):
        self.cat_cog.set_options({
            'hmmscan_out': self.work_dir + '/tmp_string2cog'
        })
        self.cat_cog.on('end', self.run_cog_stat)
        self.cat_cog.run()

    def run_cog_stat(self):
        os.mkdir(self.output_dir + '/cog_result_dir')
        self.cog_stat.set_options({
            'cog_table': self.cat_cog.option('hmmscan_result').prop['path'],
            'reads_profile_table': self.option('reads_profile_table'),
        })
        self.cog_stat.on('end', self.set_output, 'cog_stat')
        # self.cog_stat.on('end', self.end)
        self.cog_stat.run()

    def run_kegg_anno(self):
        os.mkdir(self.work_dir + '/tmp_kegg_anno')
        os.mkdir(self.work_dir + '/tmp_kegg_enzyme')
        os.mkdir(self.work_dir + '/tmp_kegg_module')
        os.mkdir(self.work_dir + '/tmp_kegg_pathway')
        xml_file = os.listdir(self.option('kegg_xml_dir').prop['path'])
        for i in xml_file:
            file_path = os.path.join(self.option('kegg_xml_dir').prop['path'], i)
            kegg_anno = self.add_tool("meta_genomic.metagen_kegg_anno")
            kegg_anno.set_options({
                "kegg_xml": file_path
            })
            kegg_anno.on('end', self.set_output, "kegg_anno")
            self.kegg_anno_tools.append(kegg_anno)
        if len(self.kegg_anno_tools) == 1:
            self.kegg_anno_tools[0].on('end', self.run_kegg_stat)
        else:
            self.on_rely(self.kegg_anno_tools, self.run_cat_kegg)
        for tool in self.kegg_anno_tools:
            tool.run()

    def run_cat_kegg(self):
        os.mkdir(self.work_dir + '/kegg_anno')
        path_list = ['tmp_kegg_anno', 'tmp_kegg_enzyme', 'tmp_kegg_module', 'tmp_kegg_pathway']
        for i in path_list:
            cat_kegg = self.add_tool('align.cat_hmmscanout')
            cat_kegg.set_options({
                'hmmscan_out': os.path.join(self.work_dir, i)
            })
            cat_kegg.on('end', self.set_output, i)
            self.cat_kegg_tools.append(cat_kegg)
        self.on_rely(self.cat_kegg_tools, self.run_kegg_stat)
        for tool in self.cat_kegg_tools:
            tool.run()

    def run_kegg_stat(self):
        os.mkdir(self.output_dir + '/kegg_result_dir')
        self.kegg_stat.set_options({
            "kegg_result_dir": self.work_dir + '/kegg_anno',
            "reads_profile": self.option('reads_profile_table'),
        })
        self.kegg_stat.on('end', self.set_output, 'kegg_stat')
        # self.kegg_stat.on('end', self.end)
        self.kegg_stat.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'nr_ncbi':
            if os.path.exists(self.work_dir + '/tmp_ncbi_result'):
                pass
            else:
                os.mkdir(self.work_dir + '/tmp_ncbi_result')
            self.ncbi_number += 1
            os.link(obj.output_dir + '/query_taxons_detail.xls',
                    self.work_dir + '/tmp_ncbi_result/' + str(self.ncbi_number) + '_taxon.xls')
        elif event['data'] == 'nr_stat':
            if os.path.exists(self.work_dir + '/tmp_nrstat_result'):
                pass
            else:
                os.mkdir(self.work_dir + '/tmp_nrstat_result')
            self.nrstat_number += 1
            os.link(obj.output_dir + '/tax_profile.xls',
                    self.work_dir + '/tmp_nrstat_result/' + str(self.nrstat_number) + '_tax_profile.xls')
        elif event['data'] == 'nr_level':
            self.linkdir(obj.output_dir, "nr_tax_level")
        elif event['data'] == 'string2cog':
            if os.path.exists(self.work_dir + '/tmp_string2cog'):
                pass
            else:
                os.mkdir(self.work_dir + '/tmp_string2cog')
            self.string2cog += 1
            os.link(obj.output_dir + '/cog_anno_result.xls',
                    self.work_dir + '/tmp_string2cog/' + 'string2cog_{}.xls'.format(str(self.string2cog)))
        elif event['data'] == 'cog_stat':
            self.linkdir(obj.output_dir, "cog_result_dir")
        elif event['data'] == 'kegg_stat':
            self.linkdir(obj.output_dir, "kegg_result_dir")
        elif event['data'] == 'kegg_anno':
            self.kegg_anno += 1
            os.link(obj.output_dir + '/kegg_anno.xls',
                    self.work_dir + '/tmp_kegg_anno/' + 'kegg_anno_{}.xls'.format(str(self.kegg_anno)))
            os.link(obj.output_dir + '/kegg_enzyme_list.xls',
                    self.work_dir + '/tmp_kegg_enzyme/' + 'kegg_enzyme_{}.xls'.format(str(self.kegg_anno)))
            os.link(obj.output_dir + '/kegg_module_list.xls',
                    self.work_dir + '/tmp_kegg_module/' + 'kegg_module_{}.xls'.format(str(self.kegg_anno)))
            os.link(obj.output_dir + '/kegg_pathway_list.xls',
                    self.work_dir + '/tmp_kegg_pathway/' + 'kegg_pathway_{}.xls'.format(str(self.kegg_anno)))
        elif event['data'] == 'tmp_kegg_anno':
            os.link(obj.option('hmmscan_result').prop['path'],
                    self.work_dir + '/kegg_anno/' + 'kegg_anno.xls')
        elif event['data'] == 'tmp_kegg_enzyme':
            os.link(obj.option('hmmscan_result').prop['path'],
                    self.work_dir + '/kegg_anno/' + 'kegg_enzyme_list.xls')
        elif event['data'] == 'tmp_kegg_module':
            os.link(obj.option('hmmscan_result').prop['path'],
                    self.work_dir + '/kegg_anno/' + 'kegg_module_list.xls')
        elif event['data'] == 'tmp_kegg_pathway':
            os.link(obj.option('hmmscan_result').prop['path'],
                    self.work_dir + '/kegg_anno/' + 'kegg_pathway_list.xls')
        else:
            pass

    def end(self):
        if self.option('nr_xml_dir').is_set:
            self.option('tax_level_dir', self.output_dir + '/nr_tax_level')
        if self.option('string_xml_dir').is_set:
            self.option('cog_result_dir', self.output_dir + '/cog_result_dir')
        if self.option('kegg_xml_dir').is_set:
            self.option('kegg_result_dir', self.output_dir + '/kegg_result_dir')
        super(CommonAnnoStatModule, self).end()

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
                file_name = os.listdir(oldfiles[i])
                os.mkdir(newfiles[i])
                for file_name_ in file_name:
                    os.link(os.path.join(oldfiles[i], file_name_), os.path.join(newfiles[i], file_name_))