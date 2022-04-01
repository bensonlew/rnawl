# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
# last_modifiy =  210325

from biocluster.workflow import Workflow
import os
import gevent
import json
from mainapp.libs.param_pack import group_detail_sort
from mbio.packages.metagenomic.id_convert import id2name
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc

class PersonalAnnoWorkflow(Workflow):
    """
    宏基因组个性化注释
    """

    def __init__(self, wsheet_object):
        super(PersonalAnnoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_type", "type": "string", "default": ""},  # 个性化注释类型
            {"name": "main_col", "type": "string"},  # 注释对应的主表集合名称
            {"name": "main_table_name", "type": "string"},  # 注释对应的主表名称
            {"name": "task_id", "type": "string"},  # 所属项目的task_id
            {"name": "update_info", "type": "string", "default": ""},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.find_one = self.api.api('metagenomic.common_api').find_one

        self.nr_anno_file = ''
        self.database_version = {
            "mvir": "mvirDB_v20110724",
            "cyps": "CYPED_v6.0",
            # "probio" : "probio_v2016.06.20",
            "phi": "phi_v4.9",
            "tcdb": "tcdb_v20200917",
            "sec": "v4.1",
            "t3ss": "v1.0.1",
            "pfam": "pfam_v33.1",
            "go": "go_v1.2",
            "lca": "nr_v20200604",
            "de_unclassified": "nr_v20200604",
            # "qs" : "-" # 没有版本号
        }

    def get_infos(self):
        """
        通过基因集主表获取序列路径
        通过anno_nr 获取项目信息和注释主表需要的其它信息
        """
        gene_set = self.find_one('geneset',
                                 {'task_id': self.option('task_id'),
                                  'name': 'GENESET_Origin'}
                                 )
        self.read_num_file = gene_set['reads_num']
        self.query = self.read_num_file.split('gene_profile')[0] +\
            'uniGeneset/gene.uniGeneset.faa'
        nr_best = self.find_one('anno_nr',
                                {'task_id': self.option('task_id'),
                                 'name': 'NR_Origin'}
                                )
        self.nr_blast_out = '/'.join(nr_best['anno_file'].split('/')[:-1]) +\
            '/nr_align_table.xls'
        params = json.loads(nr_best['params'])
        if 'qs' in self.option('anno_type'):
            self.export_group(params['group_detail'])
        params['database'] = self.option('anno_type')
        self.personal_params = {
            "geneset_id": nr_best['geneset_id'],
            "group_detail": group_detail_sort(params['group_detail']),
            "group_id": params['group_id'],
            "database": "",
            "submit_location": "",
            "task_type": 2
        }
        self.common_param = {
            "specimen": nr_best['specimen'],
            "anno_dir": self._sheet.output
        }

    def export_group(self, group_detail):
        self.group_file = os.path.join(self.work_dir, 'group.txt')
        id_name = id2name(self.option('task_id'), type='task')
        with open(self.group_file, 'w') as w:
            w.write('#sample\tgroup_name\n')
            for k, v in group_detail.items():
                for s in v:
                    w.write('{}\t{}\n'.format(id_name[s], k))

    def insert_main_table(self):
        '''
        运行前导入注释主表，并加入update_info中，基于框架更新此表的状态
        '''
        self.api_dict = {
            'go': "metagenomic.anno_go",  # [col_name, table_name]
            'phi': "metagenomic.mg_anno_phi",
            'mvirdb': "metagenomic.mvirdb",
            'tcdb': "metagenomic.tcdb",
            'qs': "metagenomic.qs_anno",
            'pfam': "metagenomic.mg_anno_pfam",
            'sec': "metagenomic.mg_anno_sec",
            'sec_de': "metagenomic.mg_anno_sec",
            'sec_lca': "metagenomic.mg_anno_sec",
            't3ss': "metagenomic.mg_anno_ttss",
            't3ss_de': "metagenomic.mg_anno_ttss",
            't3ss_lca': "metagenomic.mg_anno_ttss",
            'probio': "metagenomic.probiotics",
            'probio_de': "metagenomic.probiotics",
            'probio_lca': "metagenomic.probiotics",
            'p450': "metagenomic.mg_anno_cyps",
            'nr_de': "metagenomic.mg_anno_nr",
            'nr_lca': "metagenomic.mg_anno_nr"
        }
        self.anno_api = self.api.api(self.api_dict[self.option('anno_type')])
        self.personal_api_anno = self.api.api("metagenomic.personal_anno")
        self.main_table_id = self.personal_api_anno.add_main(
            self.option('anno_type'), self.anno_api,
            self.option('main_table_name'), self.personal_params,
            self.common_param, self.option('task_id')
        )
        update_info = json.loads(self.option('update_info'))
        update_info[str(self.main_table_id)] = self.option('main_col')
        self.sheet.set_option('update_info', json.dumps(update_info))
        self.step.update()

    def run(self):
        self.logger.info("Start Run Workflow")
        self.get_infos()
        self.insert_main_table()
        if self.option('anno_type') in ['nr_de', 'nr_lca']:
            self.diamond_nr = self.add_module('align.meta_diamond')
            self.nr_anno = self.add_module('annotation.mg_common_anno_stat')
            self.diamond_nr.on('end', self.run_nr_anno)
            self.nr_anno.on('end', self.set_db)
            self.run_dimond_nr()
        else:
            self.anno = self.add_module('annotation.personal_anno')
            self.add_event('wait_nr')
            self.on('wait_nr', self.run_anno)
            self.on('start', self.wait_nr)
            self.anno.on('end', self.set_db)
        super(PersonalAnnoWorkflow, self).run()

    def run_dimond_nr(self):
        opts = {
            'query': self.query,
            'query_type': "prot",
            'database': 'nr_v20200604',
            'lines': '50000',
            "target_num": 5,
        }
        self.diamond_nr.set_options(opts)
        self.diamond_nr.run()

    def run_nr_anno(self):
        nr_meth = 'deunclassied' if 'de' in self.option('anno_type') else 'lca'
        opts = {
            'reads_profile_table': self.read_num_file,
            "nr_method": nr_meth,
            "nr_xml_dir": self.diamond_nr.option('outxml_dir'),
            "out_type": 1
        }
        self.nr_anno.set_options(opts)
        self.nr_anno.run()

    def run_anno(self):
        database = self.option('anno_type').split('_')
        set_dic = {
            "database": database[0],
            "query": self.query,
            "reads_profile_table": self.read_num_file,
        }
        if 'sec' in database:
            set_dic['nr_gene_anno'] = self.nr_anno_file
        elif 'de' in database:
            set_dic['nr_gene_anno_de'] = self.nr_anno_file
        elif 'lca' in database:
            set_dic['nr_gene_anno_lca'] = self.nr_anno_file
        elif database[0] in ['t3ss', 'probio']:
            set_dic['nr_gene_anno'] = self.nr_anno_file
        if database[0] == 'qs':
            set_dic['group_table'] = self.group_file
        if database[0] == 'go':
            set_dic['blastout'] = self.nr_blast_out
        self.logger.info("run_anno 参数: {}".format(set_dic))
        self.anno.set_options(set_dic)
        self.anno.run()

    def wait_nr(self):
        if 'de' in self.option('anno_type'):
            table_name = 'NR_Origin_Deunclassified'
        elif 'lca' in self.option('anno_type'):
            table_name = 'NR_Origin_LCA'
        else:
            table_name = 'NR_Origin'
        if 'sec' in self.option('anno_type'):
            table_name = 'NR_Origin'
        wait = True
        nr_info = ''
        if not table_name:
            wait = False
        else:
            self.stop_timeout_check()
        while wait:
            self.logger.info('等待{}完成后，进行后续注释'.format(table_name))
            nr_info = self.find_one('anno_nr',
                                    {'task_id': self.option('task_id'),
                                     'name': table_name}
                                    )
            if not nr_info:
                info = '未进行{}注释，{}注释无法继续'
                self.set_error(info.format(table_name,
                                           self.option('anno_type')))
            elif nr_info['status'] == 'end':
                self.nr_anno_file = nr_info['anno_file']
                break
            elif nr_info['status'] == 'failed':
                info = '{}注释失败，{}注释无法继续'
                self.set_error(info.format(table_name,
                                           self.option('anno_type')))
            gevent.sleep(5)
        self.fire('wait_nr')

    def link(self, source):
        base_dir = os.listdir(source)[0]
        if not os.path.exists(os.path.join(self.output_dir, base_dir)):
            os.mkdir(os.path.join(self.output_dir, base_dir))
        for f in os.listdir(os.path.join(source, base_dir)):
            target = 'output/' + base_dir + '/' + f
            f = os.path.join(source, base_dir, f)
            super(PersonalAnnoWorkflow, self).link(f, target)

    def set_db(self):
        """
        保存结果output，导mongo数据库
        """
        if 'nr' in self.option('anno_type'):
            base_dir = os.listdir(self.nr_anno.output_dir)[0]
            self.link(self.anno_nr.output_dir)
        else:
            base_dir = os.listdir(self.anno.output_dir)[0]
            self.link(self.anno.output_dir)
        if self.option('anno_type') == 'sec_de':
            base_dir = 'Sec_Deunclassifed'
            os.rename(self.output_dir + '/Sec', self.output_dir + '/' + base_dir)
        elif self.option('anno_type') == 'sec_lca':
            base_dir = 'Sec_LCA'
            os.rename(self.output_dir + '/Sec', self.output_dir + '/' + base_dir)
        elif self.option('anno_type') == 'nr_de':
            os.rename(self.output_dir + '/' + base_dir, self.output_dir + '/nr_deunclassied')
            base_dir = 'nr_deunclassied'
        elif self.option('anno_type') == 'nr_lca':
            os.rename(self.output_dir + '/' + base_dir, self.output_dir + '/nr_lca')
            base_dir = 'nr_lca'
        self.personal_api_anno.add_detail(
            self.option('anno_type'), self.anno_api,
            self.main_table_id, os.path.join(self.output_dir, base_dir)
        )
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status and self.option("anno_type") in ["go", "tcdb", "phi", "mvirdb"]:
        # if self.option("save_pdf") and self.option("anno_type") in ["go", "tcdb", "phi", "mvirdb"]:
            # name = get_name(self.option("main_table_id"), self.option("main_col"))
            anno_submit_dic = {"go":"annogo", "tcdb":"annotcdb", "phi":"annophi", "mvirdb":"annomvir"}
            anno_col = {"go": "anno_go", "tcdb": "anno_tcdb", "phi": "anno_phi", "mvirdb": "anno_mvir"}
            main_id = self.find_one(anno_col[self.option("anno_type")], {'task_id': self.option('task_id'), 'name': self.option("main_table_name")})["_id"]
            # submit_loc = get_submit_loc(main_id, self.option("main_col"))
            submit_loc = anno_submit_dic[self.option("anno_type")]
            # self.logger.info(self.option("task_id"),str(main_id),self.option("main_table_name"),submit_loc)
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                # "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "task_id": self.option("task_id"),
                "table_id": str(main_id),
                "table_name": self.option("main_table_name"),
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status and self.option("anno_type") in ["go", "phi", "mvirdb", "tcdb"]:
            pdf_outs = self.figsave.output_dir
            if self.option("anno_type") == "go":
                os.link(os.path.join(pdf_outs, "go_level_bar.pdf.tar.gz"), os.path.join(self.output_dir, "Go", "go_level_bar.pdf.tar.gz"))
            elif self.option("anno_type") == "phi":
                os.link(os.path.join(pdf_outs, "phi_phenotype_bar.pdf.tar.gz"), os.path.join(self.output_dir, "Phi", "phi_phenotype_bar.pdf.tar.gz"))
            elif self.option("anno_type") == "mvirdb":
                os.link(os.path.join(pdf_outs, "Mvirdb_type_bar.pdf.tar.gz"), os.path.join(self.output_dir, "Mvirdb", "Mvirdb_type_bar.pdf.tar.gz"))
            elif self.option("anno_type") == "tcdb":
                os.link(os.path.join(pdf_outs, "TCDB_class_subclass_bar.pdf.tar.gz"), os.path.join(self.output_dir, "Tcdb", "TCDB_class_subclass_bar.pdf.tar.gz"))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["Tcdb", "", "TCDB功能注释结果目录", 0, "120058"],
            ["Tcdb/class_abund_out.xls", "xls", "各样品Class丰度表", 0, "120060"],
            ["Tcdb/subclass_abund_out.xls", "xls", "各样品SubClass丰度表",0,"120297"],
            ["Tcdb/family_abund_out.xls", "xls", "各样品family丰度表",0,"120298"],
            ["Tcdb/tcdb_abund_out.xls", "xls", "各样品tcdb丰度表",0,"120299"],
            ["Tcdb/Class_gene_stat.xls", "xls", "各样品Class基因列表", 0, "120060"],
            ["Tcdb/TCDB/SubClass_gene_stat.xls", "xls", "各样品SubClass基因列表",0,"120300"],
            ["Tcdb/Family_gene_stat.xls", "xls", "各样品family基因列表",0,"120301"],
            ["Tcdb/TCDB_gene_stat.xls", "xls", "各样品tcdb基因列表",0,"120302"],
            ["Tcdb/gene_tcdb_anno.xls", "xls", "每条基因的TCDB功能注释表",0,"120303"],
            ["Tcdb/TCDB_class_subclass_bar.pdf.tar.gz", "pdf", "class_subclass丰度图"],
            ["Mvirdb", "", "MvirDB功能注释结果目录", 0, "120058"],
            ["Mvirdb/Factor_abund_out.xls", "xls", "各样品Factor丰度表", 0, "120060"],
            ["Mvirdb/Type_abund_out.xls", "xls", "各样品Type丰度表",0,"120304"],
            ["Mvirdb/Factor_gene_stat.xls", "xls", "各样品Factor基因列表",0,"120305"],
            ["Mvirdb/Type_gene_stat.xls", "xls", "各样品Type基因列表",0,"120306"],
            ["Mvirdb/gene_mvirdb_anno.xls", "xls", "每条基因的MvirDB功能注释表",0,"120307"],
            ["Mvirdb/Mvirdb_type_bar.pdf.tar.gz", "pdf", "Type丰度图"],
            ["Probio", "", "益生菌功能注释结果目录",0,"120308"],
            ["Probio/gene_probio_anno.xls", "xls", "每条基因的益生菌注释表",0,"120309"],
            ["Go", "", "go功能注释结果目录",0,"120310"],
            ["Go/all.go1.function.xls", "xls", "level1水平各样品丰度表",0,"120311"],
            ["Go/all.go12.function.xls", "xls", "level2水平各样品丰度表",0,"120312"],
            ["Go/all.go123.function.xls", "xls", "level3水平各样品丰度表",0,"120313"],
            ["Go/all.go11234.function.xls", "xls", "level4水平各样品丰度表",0,"120314"],
            ["Go/all.go.annotation.xls", "xls", "各样品go注释gene对应注释表",0,"120315"],
            ["Go/go_level_bar.pdf.tar.gz", "pdf", "GO注释分类统计柱形图"],
            ["P450", "", "P450蛋白功能注释结果目录",0,"120316"],
            ["P450/gene_cyps_anno.xls", "xls", "每条蛋白的CYPS功能注释表",0,"120317"],
            ["P450/cyps_homo_profile.xls", "xls", "各样品CYPS Homologous_family丰度表",0,"120318"],
            ["P450/cyps_super_profile.xls", "xls", "各样品CYPS Superfamily丰度表",0,"120319"],
            ["P450/cyps_sid_profile.xls", "xls", "样品注释的最低层级表",0,"120320"],
            ["P450/cyps_anno_stat.xls", "xls", "pfam注释信息统计表",0,"120321"],
            ["Nr", "", "NR功能注释结果目录", 0, "120035"],
            ["Nr/gene_nr_anno.xls", "xls", "每条基因的物种注释表", 0, "120036"],
            ["Nr/nr_align_table.xls", "", "物种序列比对结果", 0, "120037"],
            ["Nr/tax_d.xls", "xls", "域注释丰度表", 0, "120038"],
            ["Nr/tax_k.xls", "xls", "界注释丰度表", 0, "120039"],
            ["Nr/tax_p.xls", "xls", "门注释丰度表", 0, "120040"],
            ["Nr/tax_c.xls", "xls", "纲注释丰度表", 0, "120041"],
            ["Nr/tax_o.xls", "xls", "目丰注释度表", 0, "120042"],
            ["Nr/tax_f.xls", "xls", "科注释丰度表", 0, "120043"],
            ["Nr/tax_g.xls", "xls", "属注释丰度表", 0, "120044"],
            ["Nr/tax_s.xls", "xls", "种注释丰度表", 0, "120045"],
            ["Pfam", "", "Pfam结构域注释结果目录",0,"120322"],
            ["Pfam/gene_pfam_anno.xls", "xls", "每条基因的Pfam功能注释表",0,"120323"],
            ["Pfam/pfam_type_profile.xls", "xls", "各样品的Pfam Type丰度表",0,"120324"],
            ["Pfam/pfam_acc_profile.xls", "xls", "各样品的Pfam Pfam丰度表/最低层级表",0,"120325"],
            ["Pfam/pfam_clan_profile.xls", "xls", "各样品的Pfam CLAN丰度表",0,"120326"],
            ["Pfam/gene_pfam_anno_stat.xls", "xls", "Pfam 注释基因统计",0,"120327"],
            ["Phi", "", "Phi结果目录",0,"120328"],
            ["Phi/gene_phi_anno.xls", "xls", "PHI功能注释结果信息表", 0, "120329"],
            ["Phi/phi_host_profile.xls", "xls", "各样品Host丰度表", 0, "120330"],
            ["Phi/phi_pathogen_profile.xls", "xls", "各样品Pathogen丰度表", 0, "120331"],
            ["Phi/phi_phenotype_profile.xls", "xls", "各样品Phenotype丰度表", 0, "120332"],
            ["Phi/phi_protein_profile.xls", "xls", "各样品Protein丰度表", 0, "120333"],
            ["Phi/phi_phenotype_bar.pdf.tar.gz", "pdf", "phenotype基因丰度图"],
            ["Qs", "", "QS功能注释结果目录", 0, "120058"],
            ["Qs/qs_class_profile.xls", "xls", "各样品QS的class水平丰度表",0,"120334"],
            ["Qs/gene_qs_anno.xls", "xls", "各样品gene对应注释表",0,"120335"],
            ["Qs/qs_lowest_profile.xls", "xls", "各样品QS的最低水平丰度表",0,"120336"],
            ["Qs/anno_qs_graph.xls", "xls", "各样品QS画图的数据表",0,"120337"],
            ["Sec", "", "分泌蛋白预测结果目录", 0, "120338"],
            ["T3ss", "", "分泌蛋白预测结果目录",0,"120345"]
        ])
        result_dir.add_regexp_rules([
            [r"Sec/.*fisher\.txt", "txt", "分泌蛋白物种注释结果", 0, "120339"],
            [r"Sec/.*summary\.txt", "txt", "分泌蛋白个数统计", 0, "120340"],
            [r"Ttss/.*fisher\.txt", "txt", "Ttss注释结果",0,"120347"],
            [r"Ttss/.*summary\.txt", "txt", "Ttss个数统计",0,"120348"],
            [r"Ttss/.*predict\.txt", "txt", "Ttss基因列表",0,"120349"]
        ])
        super(PersonalAnnoWorkflow, self).end()
