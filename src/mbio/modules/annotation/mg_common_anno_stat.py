# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last_modify:2017.09.25
# last modify by shaohau.yuan

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class MgCommonAnnoStatModule(Module):
    """
    宏基因组流程必做注释nr\kegg\cog注释的结果统计模块
    """

    def __init__(self, work_id):
        super(MgCommonAnnoStatModule, self).__init__(work_id)
        options = [
            {"name": "nr_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # nr比对结果文件夹
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # 样本序列丰度表
            {"name": "tax_level_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},  # nr注释输出结果
            {"name": "cog_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # string比对结果文件夹
            {"name": "cog_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},  # cog统计结果文件夹
            {"name": "kegg_xml_dir", "type": "infile", "format": "align.blast.blast_xml_dir"},  # kegg比对结果文件夹
            {"name": "kegg_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},  # kegg统计结果文件夹
            {"name": "nr_method", "type": "string", "default": "best_hit"},  # 新增NR不同筛选结果方法,best_hit,lca,deunclassified
            {"name": "out_type", "type": "int", "default": 0},  # # type为0或1,0为未整理的物种信息，1为整理后的物种信息
            {"name": "nr_db", "type": "string", "default": "nr_v20200604"}  # # type为0或1,0为未整理的物种信息，1为整理后的物种信息
        ]
        self.add_option(options)
        self.ncbi_tools = []
        self.nr_stat_tools = []
        self.eggnog_tools = []
        self.kegg_anno_tools = []
        self.kegg_anno = 0
        self.ncbi_number = 0
        self.nrstat_number = 0
        self.eggnog = 0
        self.cog_align = 0
        self.nr_level = self.add_tool('annotation.mg_nr_tax_level')
        self.cog_stat = self.add_tool('annotation.mg_cog_stat')
        self.kegg_stat = self.add_tool('annotation.mg_kegg_stat')
        self.kegg_level = self.add_tool('annotation.mg_kegg_level')
        self.kegg_st_le_tools = []
        self.sum_tool = []

    def check_options(self):
        if not self.option("reads_profile_table").is_set:
            raise OptionError("必须设置参数reads_profile_table", code="21201101")
        if not self.option("nr_xml_dir").is_set:
            if not self.option("kegg_xml_dir").is_set:
                if not self.option("cog_xml_dir").is_set:
                    raise OptionError("必须输入一种xml_dir文件夹(nr/string/kegg)", code="21201102")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(MgCommonAnnoStatModule, self).run()
        if self.option("nr_xml_dir").is_set:
            self.run_nr_ncbi()
            self.sum_tool.append(self.nr_level)
        if self.option("cog_xml_dir").is_set:
            self.run_eggnog()
            self.sum_tool.append(self.cog_stat)
        if self.option("kegg_xml_dir").is_set:
            self.run_kegg_anno()
            self.sum_tool.append(self.kegg_stat)
            self.sum_tool.append(self.kegg_level)
        if len(self.sum_tool) == 1:
            self.sum_tool[0].on('end', self.set_output)
        else:
            self.on_rely(self.sum_tool, self.set_output)

    def set_run(self, event):
        '''
        控制tool的串联运行
        '''
        tools = event['data']
        if len(tools) == 1:
            tools[0].run()
        else:
            t2 = tools[1:]
            tools[0].on('end', self.set_run, data=t2)
            tools[0].run()

    def run_tools(self, tools, thread=60):
        '''
        对于一个 tools的列表，最多并行运行60个
        '''
        step = len(tools) / thread
        if (len(tools) % thread) > 0:
            step += 1
        for s in range(0, len(tools), step):
            my_tools = tools[s + 1: s + step]
            if my_tools:
                tools[s].on('end', self.set_run, data=my_tools)
            tools[s].run()

    def run_nr_ncbi(self):
        xml_file = os.listdir(self.option('nr_xml_dir').prop['path'])
        for i in xml_file:
            file_path = os.path.join(self.option('nr_xml_dir').prop['path'], i)
            ncbi = self.add_tool("annotation.mg_ncbi_taxon")
            opts = {
                "blastout": file_path,
                "blastdb": self.option("nr_db")
            }
            if self.option("nr_method"):
                opts["nr_method"] = self.option("nr_method")
            ncbi.set_options(opts)
            # ncbi.on('end', self.set_output, "nr_ncbi")
            self.ncbi_tools.append(ncbi)
        if len(self.ncbi_tools) == 1:
            self.ncbi_tools[0].on('end', self.run_nr_stat)
        else:
            self.on_rely(self.ncbi_tools, self.run_nr_stat)

        self.run_tools(self.ncbi_tools)

    def run_nr_stat(self):
        file_dir = os.path.join(self.work_dir, 'tmp_ncbi_result_deal')
        gene_nr_dir = os.path.join(self.work_dir, 'tmp_gene_nr_anno')
        nr_align_dir = os.path.join(self.work_dir, 'tmp_ncbi_align')
        self.remkdir(gene_nr_dir)
        self.remkdir(file_dir)
        self.remkdir(nr_align_dir)

        for tool in self.ncbi_tools:
            self.ncbi_number += 1
            eachtaxon_new = os.path.join(tool.output_dir, "query_taxons.xls")
            eachtaxon_new_link = file_dir + '/query_taxons_{}.xls'.format(str(self.ncbi_number))
            self.relink(eachtaxon_new, eachtaxon_new_link)
            gene_nr_anno = os.path.join(tool.output_dir, "gene_nr_anno.xls")
            gene_nr_anno_link = gene_nr_dir + "/gene_nr_anno_{}.xls".format(str(self.ncbi_number))
            self.relink(gene_nr_anno, gene_nr_anno_link)
            eachnr_new = os.path.join(tool.work_dir, "temp_blastable.xls")
            eachnr_new_link = nr_align_dir + '/temp_blastable_{}.xls'.format(str(self.ncbi_number))
            self.relink(eachnr_new, eachnr_new_link)
        file_name = os.listdir(file_dir)
        for i in file_name:
            file_path = os.path.join(file_dir, i)
            nr_stat = self.add_tool("annotation.mg_nr_stat")
            opts = {
                "taxon_out": file_path,
                "reads_profile_table": self.option('reads_profile_table'),
            }
            if "nr_method" in self.get_option_object().keys():
                opts["nr_method"] = self.option("nr_method")
            nr_stat.set_options(opts)
            self.nr_stat_tools.append(nr_stat)
        if len(self.nr_stat_tools) == 1:
            self.nr_stat_tools[0].on('end', self.run_nr_tax_level)
        else:
            self.on_rely(self.nr_stat_tools, self.run_nr_tax_level)
        for tool in self.nr_stat_tools:
            tool.run()

    def run_nr_tax_level(self):
        file_dir = os.path.join(self.work_dir, 'tmp_nrstat_result')
        self.remkdir(file_dir)

        for tool in self.nr_stat_tools:
            self.nrstat_number += 1
            nr_stat_old = os.path.join(tool.output_dir, "tax_profile.xls")
            self.logger.info("nr注释统计结果的路径:{}".format(nr_stat_old))
            nr_linkfile = file_dir + '/tax_profile_{}.xls'.format(str(self.nrstat_number))
            self.relink(nr_stat_old, nr_linkfile)
        opts = {
            'nr_taxon_profile_dir': self.work_dir + '/tmp_nrstat_result',
            'nr_gene_anno_dir': self.work_dir + '/tmp_gene_nr_anno',
            'nr_align_dir': self.work_dir + '/tmp_ncbi_align',
            }
        self.nr_level.set_options(opts)
        self.nr_level.run()

    def run_eggnog(self):
        xml_file = os.listdir(self.option('cog_xml_dir').prop['path'])
        for i in xml_file:
            file_path = os.path.join(self.option('cog_xml_dir').prop['path'], i)
            eggnog = self.add_tool("annotation.mg_cog_anno")
            eggnog.set_options({
                "cog_xml": file_path
            })
            #eggnog.on('end', self.set_output, "eggnog")
            self.eggnog_tools.append(eggnog)
        if len(self.eggnog_tools) == 1:
            self.eggnog_tools[0].on('end', self.run_cog_stat)
        else:
            self.on_rely(self.eggnog_tools, self.run_cog_stat)
        for tool in self.eggnog_tools:
            tool.run()

    def run_cog_stat(self):
        # self.remkdir(self.output_dir + '/cog_result_dir')
        file_dir = os.path.join(self.work_dir, 'tmp_eggnog')
        align_dir = os.path.join(self.work_dir, 'tmp_cog_align')
        self.remkdir(file_dir)
        self.remkdir(align_dir)
        for tool in self.eggnog_tools:
            self.eggnog += 1
            cog_old = os.path.join(tool.output_dir, "gene_cog_anno.xls")
            cog_linkfile = file_dir + '/eggnog_{}.xls'.format(str(self.eggnog))
            align_old = tool.output_dir + "/tmp_cog_table.xls"
            align_link = align_dir + '/eggnog_align_{}.xls'.format(str(self.eggnog))
            self.relink(cog_old, cog_linkfile)
            self.relink(align_old, align_link)
        self.cog_stat.set_options({
            'cog_table_dir': self.work_dir + '/tmp_eggnog',
            'align_table_dir': self.work_dir + '/tmp_cog_align',
            'reads_profile_table': self.option('reads_profile_table'),
        })
        # self.cog_stat.on('end', self.set_output, 'cog_stat')
        # self.cog_stat.on('end', self.end)
        self.cog_stat.run()

    def run_kegg_anno(self):
        self.remkdir(self.work_dir + '/tmp_kegg_anno')
        xml_file = os.listdir(self.option('kegg_xml_dir').prop['path'])
        for i in xml_file:
            file_path = os.path.join(self.option('kegg_xml_dir').prop['path'], i)
            kegg_anno = self.add_tool("annotation.mg_kegg_anno")
            kegg_anno.set_options({
                "kegg_xml": file_path
            })
            # kegg_anno.on('end', self.set_output, "kegg_anno")
            self.kegg_anno_tools.append(kegg_anno)
        if len(self.kegg_anno_tools) == 1:
            self.kegg_anno_tools[0].on('end', self.run_kegg_stat)
        else:
            self.on_rely(self.kegg_anno_tools, self.run_kegg_stat)
        for tool in self.kegg_anno_tools:
            tool.run()

    def run_kegg_stat(self):
        # self.remkdir(self.output_dir + '/kegg_result_dir')
        file_dir = os.path.join(self.work_dir, 'tmp_kegg_anno')
        align_dir = os.path.join(self.work_dir, 'tmp_kegg_align')
        self.remkdir(file_dir)
        self.remkdir(align_dir)
        for tool in self.kegg_anno_tools:
            self.kegg_anno += 1
            kegg_old = os.path.join(tool.output_dir, "gene_kegg_anno.xls")
            kegg_linkfile = file_dir + '/kegg_anno_result_{}.xls'.format(str(self.kegg_anno))
            align_old = tool.output_dir + "/tmp_kegg_table.xls"
            align_link = align_dir + '/kegg_align_{}.xls'.format(str(self.kegg_anno))
            self.relink(kegg_old, kegg_linkfile)
            self.relink(align_old, align_link)
            self.relink(tool.output_dir + "/kegg_enzyme_list.xls",
                        file_dir + '/kegg_enzyme_list_{}.xls'.format(str(self.kegg_anno)))
            self.relink(tool.output_dir + "/kegg_module_list.xls",
                        file_dir + '/kegg_module_list_{}.xls'.format(str(self.kegg_anno)))
            self.relink(tool.output_dir + "/kegg_pathway_list.xls",
                        file_dir + '/kegg_pathway_list_{}.xls'.format(str(self.kegg_anno)))
        self.kegg_stat.set_options({
            "kegg_result_dir": self.work_dir + '/tmp_kegg_anno/',
            "align_table_dir": self.work_dir + '/tmp_kegg_align/',
            "reads_profile": self.option('reads_profile_table'),
            "xml_file_dir": self.option('kegg_xml_dir')
        })
        #self.kegg_stat.on('end', self.set_output, 'kegg_stat')
        # self.kegg_stat.on('end', self.end)
        self.kegg_st_le_tools.append(self.kegg_stat)
        self.kegg_level.set_options({
            "kegg_result_dir": self.work_dir + '/tmp_kegg_anno/',
            "reads_profile": self.option('reads_profile_table'),
        })
        # self.kegg_level.on('end', self.set_output, 'kegg_level')
        self.kegg_st_le_tools.append(self.kegg_level)
        if len(self.kegg_st_le_tools) < 2:
            self.logger.info("kegg stat or level error!")
        for tool in self.kegg_st_le_tools:
            tool.run()

    def set_output(self):
        if self.option("nr_xml_dir").is_set:
            self.linkdir(self.nr_level.output_dir, "nr_tax_level")
        if self.option("cog_xml_dir").is_set:
            self.linkdir(self.cog_stat.output_dir, "cog_result_dir")
        if self.option("kegg_xml_dir").is_set:
            self.linkdir(self.kegg_stat.output_dir, "kegg_result_dir")
            self.linkdir(self.kegg_level.output_dir, "kegg_result_dir")
        self.end()

    def relink(self, oldfile, newfile):
        if os.path.exists(newfile):
            os.remove(newfile)
            os.link(oldfile, newfile)
        else:
            os.link(oldfile, newfile)

    def remkdir(self, dir_name):
        if os.path.exists(dir_name):
            pass
        else:
            os.mkdir(dir_name)

    def end(self):
        super(MgCommonAnnoStatModule, self).end()
        if self.option('nr_xml_dir').is_set:
            self.option('tax_level_dir', self.output_dir + '/nr_tax_level')
        if self.option('cog_xml_dir').is_set:
            self.option('cog_result_dir', self.output_dir + '/cog_result_dir')
        if self.option('kegg_xml_dir').is_set:
            self.option('kegg_result_dir', self.output_dir + '/kegg_result_dir')

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
