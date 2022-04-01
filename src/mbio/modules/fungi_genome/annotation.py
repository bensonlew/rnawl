# -*- coding: utf-8 -*-
# __author__ = 'juan.zhu'
# last_modify: 2018.03.12

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class AnnotationModule(Module):
    """
    细菌基因组基础注释模块
    分析内容：NR、swiss-prot、pfam、cog、go和kegg注释, cazy
    """

    def __init__(self, work_id):
        super(AnnotationModule, self).__init__(work_id)
        option = [
            {"name": "gene_seq", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测gff文件
            {"name": "database_list", "type": "string", "default": "nr_v20200604,swissprot_v20200617,pfam_v33.1,eggnog,kegg_v94.2,cazy_v8"},
            # 数据库list:nr,swissprot,pfam,eggnog,kegg,cazy
            {"name": "sample", "type": "string", "default": ""},  # 样品名称
            {"name": "tidy_nr", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_go", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_swissprot", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_ref", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_pfam", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cog", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_kegg", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cazy", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_antismash", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_summary", "type": "outfile", "format": "sequence.profile_table"},
            # {"name": "kegg_level", "type": "outfile", "format": "sequence.profile_table"},
            # {"name": "kegg_xml", "type": "outfile", "format": "align.blast.blast_xml"}  # "kegg_level"和"kegg_xml" 用于做pathway通路图(完成图时，将质粒和染色体的合并pathway level统计）
        ]
        self.add_option(option)
        self.diamond = self.add_module('align.diamond')  # nr
        self.blast = self.add_module('align.blast')  # swssiport
        self.hmmscan = self.add_module('align.hmmscan')  # pfam
        self.top_diamond = self.add_module('align.meta_diamond')  # cog、kegg
        self.cog_anno = self.add_tool('annotation.mg_cog_anno')
        self.kegg_anno = self.add_tool('annotation.kegg_anno')
        self.go_anno = self.add_module('bacgenome.go_anno')
        self.cazy_anno = self.add_tool('annotation.cazy_anno')
        self.anno_tidy = self.add_tool('fungi_genome.anno_tidy')
        self.step.add_steps('go_', 'kegg', 'anno_tidy','anno_cazy')
        self.tools = []
        self.database = ""

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        """
    检查参数
        :return:
        """
        if not self.option("gene_seq").is_set:
            raise OptionError("必须提供基因序列文件", code="21400101")
        if self.option("database_list") == "":
            raise OptionError("必须提供比对的数据名称", code="21400102")
        if self.option("sample") == "":
            raise OptionError("必须提供样品名称", code="21400103")
        return True

    def run_diamond_blast(self):
        opts = {
            'query': self.option('gene_seq'),
            'query_type': 'prot',
            'database': self.database,
            "blast": 'blastp',
            'outfmt': 6
        }
        #self.tool = ""
        if self.database in ["nr", "nr_v20200604"]:
            #self.tool = self.diamond
            opts['sensitive'] = 0##qingchen.zhang
            opts['lines'] = 50000  ##qingchen.zhang
            self.tools.append(self.diamond)
            self.tools.append(self.go_anno)
            self.diamond.on("end", self.run_go_anno)
            self.diamond.set_options(opts)
            self.diamond.run()
        elif self.database in ["swissprot", "swissprot_v20200617"]:
            #self.tool = self.blast
            self.tools.append(self.blast)
            self.blast.set_options(opts)
            self.blast.run()
        #self.tool.set_options(opts)
        #self.tool.run()

    def run_top_dimond(self):
        opts = {
            'query': self.option('gene_seq'),
            'query_type': 'prot',
            'database': self.database,
            'outfmt': 5
        }
        if self.database == "eggnog":
            self.top_diamond.set_options(opts)
            self.tools.append(self.top_diamond)
            self.tools.append(self.cog_anno)
            self.top_diamond.on("end", self.run_cog_anno)
            self.top_diamond.run()
        elif self.database in ["kegg", 'kegg_v94.2']:
            self.top_diamond2 = self.add_module('align.meta_diamond')
            self.top_diamond2.set_options(opts)
            self.tools.append(self.top_diamond2)
            self.tools.append(self.kegg_anno)
            self.top_diamond2.on("end", self.run_kegg_anno)
            self.top_diamond2.run()

    def run_hmmscm(self):
        opts = {
            'query': self.option('gene_seq'),
            'database': self.database,
        }
        if self.database in ["cazy", "cazy_v8"]:
            self.hmmscan.set_options(opts)
            self.tools.append(self.hmmscan)
            self.tools.append(self.cazy_anno)
            self.hmmscan.on("end", self.run_cazy_anno)
            self.hmmscan.run()
        elif self.database in ["pfam", "pfam_v33.1"]:
            self.hmmscan2 = self.add_module('align.hmmscan')
            self.hmmscan2.set_options(opts)
            self.tools.append(self.hmmscan2)
            self.hmmscan2.run()

    def run_go_anno(self):
        opts = {
            'blastout': self.diamond.option('outxml')
        }
        self.go_anno.set_options(opts)
        self.tools.append(self.go_anno)
        self.go_anno.on('start', self.set_step, {'start': self.step.go_})
        self.go_anno.on('end', self.set_step, {'end': self.step.go_})
        self.go_anno.on('end', self.set_output, "go")
        self.go_anno.run()

    def run_cog_anno(self):
        opts = {
            'cog_xml': self.top_diamond.option('outxml')
        }
        self.cog_anno.set_options(opts)
        self.tools.append(self.cog_anno)
        self.cog_anno.run()

    def run_kegg_anno(self):
        opts = {
            'kegg_xml': self.top_diamond2.option('outxml')
        }
        self.kegg_anno.set_options(opts)
        self.tools.append(self.kegg_anno)
        self.kegg_anno.on('start', self.set_step, {'start': self.step.kegg})
        self.kegg_anno.on('end', self.set_step, {'end': self.step.kegg})
        self.kegg_anno.on('end', self.set_output, "kegg")
        self.kegg_anno.run()

    def run_cazy_anno(self):
        opts = {
            'hmmscan_result': self.hmmscan.option('align_result'),
            'version': self.hmmscan.option("database"), ## fix by qingchen.zhang 兼容新老版本
        }
        self.cazy_anno.set_options(opts)
        self.tools.append(self.cazy_anno)
        self.cazy_anno.on('start', self.set_step, {'start': self.step.anno_cazy})
        self.cazy_anno.on('end', self.set_step, {'end': self.step.anno_cazy})
        self.cazy_anno.on('end', self.set_output, "anno_cazy")
        self.cazy_anno.run()

    def run_anno_tidy(self):
        opts = {
            'gene_gff': self.option('gene_gff'),
            'anno_nr': self.diamond.option('outtable'),
            'anno_go': self.go_anno.output_dir + "/go12level_statistics.xls",
            'anno_swissprot': self.blast.option('outtable'),
            'anno_pfam': self.hmmscan2.option('align_result'),
            'anno_cog': self.cog_anno.output_dir + "/gene_cog_anno.xls",
            'kegg_xml': self.top_diamond2.option('outxml'),
            'anno_kegg': self.kegg_anno.option('anno_kegg'),
            'anno_cazy': self.cazy_anno.output_dir + "/gene_cazy_parse_anno.xls",
            'summary': 'T',
            'sample': self.option("sample")
        }
        self.anno_tidy.set_options(opts)
        self.anno_tidy.on('start', self.set_step, {'start': self.step.anno_tidy})
        self.anno_tidy.on('end', self.set_step, {'end': self.step.anno_tidy})
        self.anno_tidy.on('end', self.set_output, "anno_tidy")
        self.anno_tidy.run()

    def run(self):
        super(AnnotationModule, self).run()
        if self.option('database_list') != "":
            base = self.option('database_list').split(",")
            for ba in base:
                self.database = ba
                if self.database in ["nr","nr_v20200604", "swissprot", "swissprot_v20200617"]:
                    self.run_diamond_blast()
                if self.database in ["eggnog", "kegg", 'kegg_v94.2']:
                    self.run_top_dimond()
                if self.database in ["pfam","pfam_v33.1", "cazy", "cazy_v8"]:
                    self.run_hmmscm()
        self.on_rely(self.tools, self.run_anno_tidy)
        self.anno_tidy.on('end',self.end)

    def set_output(self, event):
        """
        prefix = (os.path.basename(self.option("gene_seq").prop['path'])).split("_")[0]
        obj = event['bind_object']
        if event['data'] == 'go':
            go = self.output_dir + "/GO/"
            if not os.path.exists(go):
                os.makedirs(go)
            os.link(obj.output_dir + "/go12level_statistics.xls", go + prefix + "_go12level_statistics.xls")
        if event['data'] == 'kegg':
            kegg = self.output_dir + "/KEGG/"
            if not os.path.exists(kegg):
                os.makedirs(kegg)
            os.link(obj.output_dir + "/kegg_level_stat.xls", kegg + prefix + "_kegg_level_stat.xls")
            self.option("kegg_level", obj.output_dir + "/kegg_level.xls")
            self.option("kegg_xml", self.top_diamond2.option('outxml'))
        elif event['data'] == 'anno_tidy':
            nr = self.output_dir + "/NR/"
            if not os.path.exists(nr):
                os.makedirs(nr)
            os.link(obj.option('tidy_nr').prop['path'], nr + prefix + "_anno_nr.xls")
            swissprot = self.output_dir + "/Swissprot/"
            if not os.path.exists(swissprot):
                os.makedirs(swissprot)
            os.link(obj.option('tidy_swissprot').prop['path'], swissprot + prefix + "_anno_swissprot.xls")
            cog = self.output_dir + "/COG/"
            if not os.path.exists(cog):
                os.makedirs(cog)
            os.link(obj.option('tidy_cog').prop['path'], cog + prefix + "_anno_cog.xls")
            pfam = self.output_dir + "/Pfam/"
            if not os.path.exists(pfam):
                os.makedirs(pfam)
            os.link(obj.option('tidy_pfam').prop['path'], pfam + prefix + "_anno_pfam.xls")
            cazy = self.output_dir + "/CAZy/"
            if not os.path.exists(cazy):
                os.makedirs(cazy)
            os.link(obj.option('tidy_cazy').prop['path'], cazy + prefix + "_anno_cazy.xls")
            kegg = self.output_dir + "/KEGG/"
            if not os.path.exists(kegg):
                os.makedirs(kegg)
            os.link(obj.option('tidy_kegg').prop['path'], kegg + prefix + "_anno_kegg.xls")  # 还差一个KEGG通路图的文件夹
        if event['data'] == 'kegg':
            kegg = self.output_dir + "/KEGG/"
            if not os.path.exists(kegg):
                os.makedirs(kegg)
            os.link(obj.output_dir + "/kegg_level_stat.xls",
                    self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_level_stat.xls")
        """
        obj = event['bind_object']
        if event['data'] == 'anno_tidy':
            if os.path.exists(obj.output_dir + '/NR/'):
                self.linkdir(obj.output_dir + '/NR', 'NR')
                self.option('tidy_nr', self.anno_tidy.option('tidy_nr'))
            if os.path.exists(obj.output_dir + '/Swissprot'):
                self.linkdir(obj.output_dir + '/Swissprot', 'Swissprot')
                self.option('tidy_swissprot', self.anno_tidy.option('tidy_swissprot'))
            if os.path.exists(obj.output_dir + '/GO/'):
                self.linkdir(obj.output_dir + '/GO', 'GO')
                self.option('tidy_go', self.anno_tidy.option('tidy_go'))
            if os.path.exists(obj.output_dir + '/COG/'):
                self.linkdir(obj.output_dir + '/COG', 'COG')
                self.option('tidy_cog', self.anno_tidy.option('tidy_cog'))
            if os.path.exists(obj.output_dir + '/KEGG/'):
                # self.linkdir(obj.output_dir + '/KEGG', 'KEGG')  其中有一个KEGG通路图的文件夹，不允许这样操作
                if os.path.exists(self.output_dir + "/KEGG/"):
                    shutil.rmtree(self.output_dir + "/KEGG/")
                os.rename(obj.output_dir + '/KEGG', self.output_dir + "/KEGG/")
                os.link(self.kegg_anno.output_dir + "/kegg_level_stat.xls",
                        self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_level_stat.xls")
                self.option('tidy_kegg', self.anno_tidy.option('tidy_kegg'))
            if os.path.exists(obj.output_dir + '/Pfam/'):
                self.linkdir(obj.output_dir + '/Pfam', 'Pfam')
                self.option('tidy_pfam', self.anno_tidy.option('tidy_pfam'))
            if os.path.exists(obj.output_dir + '/CAZy/'):
                self.linkdir(obj.output_dir + '/CAZy', 'CAZy')
                self.option('tidy_cazy', self.anno_tidy.option('tidy_cazy'))
            if os.path.exists(obj.output_dir + '/Summary/'):
                self.linkdir(obj.output_dir + '/Summary', 'Summary')
                self.option('tidy_summary', self.anno_tidy.option('tidy_summary'))
        if event['data'] == 'anno_cazy':
            for i in ['gene_cazy_family_stat.xls','gene_cazy_class_stat.xls']:
                if os.path.exists(self.output_dir + '/' + i):
                    os.remove(self.output_dir + '/' + i)
                os.link(self.cazy_anno.output_dir + '/' + i,self.output_dir + '/' + i)
        if event['data'] == 'kegg':
            if os.path.exists(self.output_dir + '/' + 'kegg_pathway_img.tar.gz'):
                os.remove(self.output_dir + '/' + 'kegg_pathway_img.tar.gz')
            os.link(self.kegg_anno.output_dir + '/kegg_pathway_img.tar.gz',self.output_dir + '/' + 'kegg_pathway_img.tar.gz')

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def end(self):
        super(AnnotationModule, self).end()
