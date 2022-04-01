# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""细菌比较基因组分析导表工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
import os,re
import json
import shutil
from biocluster.config import Config
import time
import datetime
import gevent
import functools
import shutil
import types
from mbio.packages.bac_comp_genome.common_function import link_file,merge_16s,get_sample_from_tree


def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，该阶段导表已进行{}s".format(func_name, end - start))

    return wrapper

class BacComparativeApiWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        比较基因组workflow option参数设置
        """
        self._sheet = wsheet_object
        super(BacComparativeApiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "genome_dir", "type": "infile", "format": "sequence.fasta_dir"},  ###rawdata的文件目录
            {"name": "seq_dir", "type": "string"},  ###样品的分开序列路径
            {"name": "gff_dir", "type": "string"},  ###样品分开的gff序列路径
            {"name": "sample_list", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.remote_dir = self._sheet.output + '/'

    def check_options(self):
        """
        检查参数
        """
        pass

    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        self.samples = self.get_samples()
        self.samples_dict = self.change_samples()
        self.samples_mj, self.samples_up = self.get_majorbiodb_samples()
        gevent.spawn_later(5, self.end)
        super(BacComparativeApiWorkflow, self).run()

    def end(self):
        self.run_api()
        super(BacComparativeApiWorkflow, self).end()

    def run_api(self):
        self.output_dir = "/mnt/ilustre/users/sanger-dev/home/gaohao/compare_genome/mongo_data_20/api_test"
        self.run_data_api()
        self.run_gen_predict_api()
        self.run_core_api()
        self.run_anno_api()
        self.run_prephage_api()
        self.run_antismash_api()
        self.run_island_api()
        self.run_pan_api()
        self.run_tree_api()

    def run_data_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        data_params = {"data": "specimen"}
        self.data_id = api_path.add_main('data', name="data_origin", params=data_params, desc="data数据表!")
        api_path.add_data_specimen("data_specimen", self.samples_dict, self.data_id, "data_id")
        if len(self.samples_mj) > 0:
            api_path.add_data_gene("Majorbiodb", self.samples_mj, "data_gene", {"gene_name": "gene_prefix", "new_gene_name": "gene_prefix", "location": "Genome_accession"}, self.data_id, "data_id")
            api_path.add_data_gene("Majorbiodb", self.samples_mj, "data_detail", {"gene_name": "gene_prefix", "ass_assession": "ass_assession", "location": "Genome_accession","species":"s", "strain ":"strain", "seq_status":"seq_status", "g_size":"total_g_size", "g_location":"g_location", "seq_size":"seq_size", "gc":"gc", "scf_no":"scf_no", "con_no":"con_no", "chr_no":"chr_no", "pla_no":"pla_no", "scf_n50":"n50", "con_n50":"n50"}, self.data_id, "data_id", insert_da={"from_refdb": "T"})
            ##导入样品特殊属性表
            api_path.add_data_gene("Majorbiodb_feature", self.samples_mj, "data_attributes",
                                   { "host": "host", "iso_source": "source_material_id", "geo_loc_name": "geo_loc_name", "lat_lon":"lat_lon"}, self.data_id, "data_id", insert_da={"from_refdb": "T"})
        gene_key = "specimen_id,location,gene_name"
        if len(self.samples_up) > 0:
            for sample in self.samples_up:
                api_path.add_main_gene2(self.output_dir + "/gene_predict/" + sample + "/" + sample + "_stat.xls", "data_gene", self.data_id, gene_key, 'data_id', "gene_name", "new_gene_name", has_head=True)
        ##导基因组基本信息表
        genome_key = "specimen_id,organism"
        if len(self.samples_mj) > 0:
            api_path.add_data_gene("Majorbiodb", self.samples_mj, "genomes_detail", {"rrna": "rrna_no", "species": "s", "cds": "cds_no", "trna":"trna_no", "s16_rna":"s16_num"}, self.data_id, "data_id", convert_dic={"rrna": float, "cds": float, "trna": float, "s16_rna": float})
        if len(self.samples_up) > 0:
            api_path.add_main_detail(self.work_dir + "/data/sample_organism.xls", "genomes_detail", self.data_id, genome_key, main_name='data_id', has_head=True)


    def run_gen_predict_api(self):
        api_stat = self.api.api("bac_comp_genome.gene_stat")
        api_gene = self.api.api("bac_comp_genome.gene_predict")
        params = {
            "soft": {
                "cds": "Glimmer, GeneMark, Prodigal",
                "rrna": "Barrnap",
                "trna": "tRNAscan-SE",
            }
        }
        dir_path = self.output_dir + "/gene_predict"
        cds_path = self.remote_dir + "gene_predict"
        stat_id = api_stat.add_gene_stat(params=params)
        predict_id = api_gene.add_gene_predict(cds_path, params=params)
        seq_type = self.get_seqtype(self.option("sample_list").prop["path"])
        for sample in os.listdir(dir_path):
            type = seq_type[sample]
            sample_dir_path = os.path.join(dir_path, sample)
            sample_path = os.path.join(sample_dir_path, sample + '_sample_stat.xls')
            sample_gff= os.path.join(sample_dir_path, sample + '_CDS.gff')
            api_stat.add_gene_stat_detail(stat_id, sample, sample_path)
            self.logger.info("导入统计表格成功！")
            api_stat.add_genomes_detail(sample, sample_path)
            if type in ["seq", "gff"]:  ##参考库中的数据不再更新，seq->seq或者gff->seq+gff这两种格式的数据需要再重新更新表
                sample_stat_path = os.path.join(sample_dir_path, sample + '_stat.xls')
                api_stat.add_data_detail(sample, sample_stat_path)
            self.logger.info("更新样本信息表和基因组特征表成功！")
            api_gene.add_gene_cds_detail(predict_id, sample, sample_gff)
            self.logger.info("导入cds预测结果成功！")
            rrna_path = os.path.join(sample_dir_path, sample + '_rRNA.gff')
            s16_path = os.path.join(sample_dir_path, sample + '_16S.fna')
            s16_fa = self.remote_dir + "gene_predict/" + sample + "/" + sample + '_16S.fna'
            if os.path.exists(rrna_path):
                if os.path.exists(s16_path):
                    api_gene.add_gene_rrna_detail(predict_id, sample, rrna_path, fnn=rrna_path, fnn_path=s16_fa)
                else:
                    api_gene.add_gene_rrna_detail(predict_id, sample, rrna_path, fnn=rrna_path)
                self.logger.info("导入rrna预测结果成功！")
            trna_path = os.path.join(sample_dir_path, sample + '_tRNA.gff')
            api_gene.add_gene_trna_detail(predict_id, sample, trna_path, type)
            self.logger.info("导入trna预测结果成功！")

    def get_seqtype(self, file):
        seq_type = {}
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                sample = line[0]
                type = line[1]
                seq_type[sample] = type
        return seq_type

    def get_samples(self):
        samples = []
        with open(self.option("sample_list").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                samples.append(lin[0])
        return samples

    def change_samples(self):
        samples_dict = {}
        with open(self.option("sample_list").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                new_name = lin[0].replace(".", "_")
                samples_dict[lin[0]] = new_name
        return samples_dict

    def get_majorbiodb_samples(self):
        samples_mj = []
        samples_up = []
        with open(self.option("sample_list").prop["path"], "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if lin[1] in ["majorbio"]:
                    samples_mj.append(lin[0])
                else:
                    samples_up.append(lin[0])
        return samples_mj,samples_up

    def run_anno_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        cog_params = {"soft": "diamond", "database": "cog"}
        cog_abu_path = self.remote_dir + "/annotation/all_anno/cog"
        cog_id = api_path.add_main('anno_cog', name="cog_origin", params=cog_params, desc="cog注释结果内容!",
                                   others={"cog_path": cog_abu_path})
        cog_key = "gene_id,location,cog,cog_des,type,type_des,category,align_len,identity,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/cog/all.cog_anno.xls', 'anno_cog_detail',
                                 cog_id,
                                 cog_key, has_head=True, main_name='cog_id', convert_dic={'align_len': float, 'identity': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/cog/all.cog_genelist.xls', 'cog_comp_detail',
                                  cog_id, "cog,cog_des,function,function_des,category", has_head=True,
                                  main_name='cog_id', names_convert=self.samples_dict, names_co="False")
        kegg_params = {"soft": "diamond", "database": "kegg"}
        kegg_abu_path = self.remote_dir + "annotation/all_anno/kegg"
        kegg_id = api_path.add_main('anno_kegg', name="kegg_origin", params=kegg_params, desc="kegg注释结果内容!",
                                    others={"kegg_path": kegg_abu_path})
        kegg_key = "gene_id,location,gene_name,ko,ko_des,pathway,enzyme,module,kegg_path,level1,level2,level3,align_len,identity,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/kegg/all.kegg_anno.xls', 'anno_kegg_detail',
                                 kegg_id,
                                 kegg_key, has_head=True, main_name='kegg_id', convert_dic={'align_len': float, 'identity': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/kegg/all.KEGG_genelist.xls',
                                  'kegg_comp_detail',
                                  kegg_id, "ko,gene_name,ko_des,pathway,enzyme,module,kegg_path,level1,level2,level3",
                                  has_head=True,
                                  main_name='kegg_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_main_detail2(self.output_dir + "/annotation/all_anno/kegg/kegg_graph_info.xls", "kegg_comp_graph", kegg_id, "specimen_id,level3,ko,map_id,gene_list,shape,href,title,coords", main_name='kegg_id')
        cazy_params = {"soft": "diamond", "database": "cazy"}
        cazy_id = api_path.add_main('anno_cazy', name="cazy_origin", params=cazy_params, desc="cazy注释结果内容!")
        cazy_key = "gene_id,location,family,evalue,coverage,class,class_des,identity,score,family_des,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/cazy/all.cazy_anno.xls', 'anno_cazy_detail',
                                 cazy_id,
                                 cazy_key, has_head=True, main_name='cazy_id', convert_dic={'score': float, 'identity': float, 'evalue':float, 'coverage':float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/cazy/all.cazy_genelist.xls',
                                  'cazy_comp_detail', cazy_id, "family,class,class_des,family_des",
                                  has_head=True, main_name='cazy_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/cazy/all.cazy_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "cazy")
        card_params = {"soft": "dimond", "database": "card"}
        card_id = api_path.add_main('anno_card', name="card_origin", params=card_params, desc="card注释结果内容!")
        card_key = "gene_id,location,aro_name,accession,des,category,evalue,identity,drug_class,resis_mechanism,score,coverage,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/card/all.card_anno.xls', 'anno_card_detail',
                                 card_id,
                                 card_key, has_head=True, main_name='card_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/card/all.card_genelist.xls',
                                  'card_comp_detail', card_id,
                                  "accession,aro_name,des,category,drug_class,resis_mechanism",
                                  has_head=True, main_name='card_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/card/all.card_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "ant_res")
        phi_params = {"soft": "dimond", "database": "phi"}
        phi_id = api_path.add_main('anno_phi', name="phi_origin", params=phi_params, desc="phi注释结果内容!")
        phi_key = "gene_id,location,phi,protein,gene,tax,pathogen,type,host_des,host_id,host_species,function,identity,evalue,disease,coverage,score,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/phi/all.phi_anno.xls', 'anno_phi_detail',
                                 phi_id,
                                 phi_key, has_head=True, main_name='phi_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/phi/all.phi_genelist.xls',
                                  'phi_comp_detail', phi_id,
                                  "phi,protein,gene,tax,pathogen,type,host_des,host_id,host_species,function,disease",
                                  has_head=True, main_name='phi_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/phi/all.phi_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "phi")
        tcdb_params = {"soft": "dimond", "database": "tcdb"}
        tcdb_id = api_path.add_main('anno_tcdb', name="tcdb_origin", params=tcdb_params, desc="tcdb注释结果内容!")
        tcdb_key = "gene_id,location,tcdb,des,family,subclass,class,identity,evalue,score,align_len,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/tcdb/all.tcdb_anno.xls', 'anno_tcdb_detail',
                                 tcdb_id,
                                 tcdb_key, has_head=True, main_name='tcdb_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'align_len': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/tcdb/all.tcdb_genelist.xls',
                                  'tcdb_comp_detail', tcdb_id, "tcdb,des,family,subclass,class",
                                  has_head=True, main_name='tcdb_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/tcdb/all.tcdb_anno.xls', "Gene ID", "sample",
                                     "genomes_detail", self.data_id, "data_id", self.samples, "tcdb")
        vfdb_params = {"soft": "dimond", "database": "vfdb"}
        vfdb_id = api_path.add_main('anno_vfdb', name="vfdb_origin", params=vfdb_params, desc="vfdb注释结果内容!")
        vfdb_key = "gene_id,location,vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1,identity,evalue,score,coverage,type,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_anno.xls', 'anno_vfdb_detail',
                                 vfdb_id,
                                 vfdb_key, has_head=True, main_name='vfdb_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_predict_genelist.xls',
                                  'vfdb_comp_detail', vfdb_id,
                                  "vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1",
                                  has_head=True, main_name='vfdb_id', insert_d={"type": "predict"},
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_core_genelist.xls',
                                  'vfdb_comp_detail', vfdb_id,
                                  "vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1",
                                  has_head=True, main_name='vfdb_id', insert_d={"type": "core"},
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_all_genelist.xls',
                                  'vfdb_comp_detail', vfdb_id,
                                  "vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1",
                                  has_head=True, main_name='vfdb_id', insert_d={"type": "all"},
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "vir")
        tmhmm_params = {"soft": "dimond", "database": "tmhmm"}
        tmhmm_id = api_path.add_main('anno_tmhmm', name="tmhmm_origin", params=tmhmm_params, desc="tmhmm注释结果内容!")
        tmhmm_key = "gene_id,location,len,num,exp_all,exp_fist,n,topology,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/tmhmm/all.tmhmm_anno.xls', 'anno_tmhmm_detail',
                                 tmhmm_id,
                                 tmhmm_key, has_head=True, main_name='tmhmm_id', convert_dic={'len': float, 'num': float, 'exp_all': float, 'exp_fist': float, 'n': float})
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/tmhmm/all.tmhmm_anno.xls', "Gene ID", "sample",
                                     "genomes_detail", self.data_id, "data_id", self.samples, "tmhmm")
        secretory_params = {"soft": "dimond", "database": "secretory"}
        secretory_id = api_path.add_main('anno_secretory', name="secretory_origin", params=secretory_params,
                                         desc="secretory注释结果内容!")
        secretory_key = "gene_id,location,gene_name,type,ko_id,des,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/secretory/all.secretory_anno.xls',
                                 'anno_secretory_detail',
                                 secretory_id, secretory_key, has_head=True, main_name='secretory_id')
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/secretory/all.secretory_genelist.xls',
                                  'secretory_comp_detail', secretory_id, "type,gene_name,ko_id,des",
                                  has_head=True, main_name='secretory_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/secretory/all.secretory_anno.xls', "Gene ID", 'sample', "genomes_detail", self.data_id, "data_id", self.samples, "sec_sys")
        signalp_params = {"soft": "dimond", "database": "SignalP"}
        signalp_id = api_path.add_main('anno_signalp', name="signalp_origin", params=signalp_params,
                                       desc="signalp注释结果内容!")
        signalp_key = "gene_id,location,specimen_id,cmax,cpos,ymax,ypos,smax,spos,smean,network,d_score"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/signalp/all.Gram+_SignalP.xls',
                                 'anno_signalp_detail',
                                 signalp_id, signalp_key, has_head=True, main_name='signalp_id',
                                 insert_d={"type": "gram+"}, convert_dic={'cpos': float, 'ymax': float, 'smax': float, 'd_score': float, 'ypos': float, 'smean': float, 'cmax': float, 'spos': float})
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/signalp/all.Gram+_SignalP.xls', "Gene ID", 'Sample Name', "genomes_detail", self.data_id, "data_id", self.samples, "signalp_a")
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/signalp/all.Gram-_SignalP.xls',
                                 'anno_signalp_detail', signalp_id, signalp_key, has_head=True, main_name='signalp_id',
                                 insert_d={"type": "gram-"}, convert_dic={'cpos': float, 'ymax': float, 'smax': float, 'd_score': float, 'ypos': float, 'smean': float, 'cmax': float, 'spos': float})
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/signalp/all.Gram-_SignalP.xls', "Gene ID", 'Sample Name', "genomes_detail", self.data_id, "data_id", self.samples, "signalp_b")
        summary_params = {"soft": "diamond", "database": "summary"}
        summary_abu_path = self.remote_dir + "/annotation/all_anno/all/all.anno_summary.xls"
        anno_id = api_path.add_main('anno_summary', name="summary_origin", params=summary_params, desc="summary注释结果内容!",
                                    others={"summary_path": summary_abu_path})
        summary_key = "gene_id,strand,start,end,gene_len,location,cog_id,cog_des,func,func_des,gene_name,ko_id,ko_des,cazy_family,cazy_family_des,cazy_class,cazy_class_des,aro_acc,aro_des,aro_category,phi_id,phi_fun,pathogen_spe,host_spe,vf_id,vfs,vf_des,gram_a,gram_b,tmh_no,tcdb_des,secre_type,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/all/all.anno_summary.xls',
                                 'anno_summary_detail', anno_id, summary_key, has_head=True, main_name='anno_id', convert_dic={'start': int, 'end': int, 'gene_len': float})

    def run_prephage_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "Phage_phinder"}
        pre_id = api_path.add_main("prephage", name="prephage_origin", params=pre_params, desc="prephage的预测！")
        pre_key = "location,ph_id,ph_start,ph_end,strand,len,cds_num,pos_phage,gc,gene_list"
        names = []
        for sample in self.samples:
            if os.path.exists(self.output_dir + "/prephage/" + sample + "_prephage_summary.xls"):
                names.append(sample)
                api_path.add_main_detail3(self.output_dir + "/prephage/" + sample + "_prephage_summary.xls",
                                          self.output_dir + "/prephage/" + sample + "_prephage.fna",
                                          "prephage_gene_detail", pre_id, pre_key, "ph_id", has_head=True,
                                          main_name='prephage_id', insert_d={"specimen_id": sample}, convert_dic={'ph_start': int, 'ph_end': int, 'len': float, 'cds_num': int, 'gc': float})
        api_path.add_main_detail2(self.output_dir + "/prephage/all.prephage_abund.xls", "prephage_abund", pre_id,
                                  "prephage_type", has_head=True, main_name='prephage_id',
                                  names_convert=self.samples_dict)
        api_path.add_main_detail2(self.output_dir + "/prephage/all.prephage_genelist.xls", "prephage_detail", pre_id,
                                  "prephage_type", has_head=True, main_name='prephage_id',
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_num_update(self.output_dir + "/prephage", "_prephage_summary.xls", "genomes_detail", self.data_id, "data_id", "prephage", self.samples)

    def run_island_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        is_params = {"soft": "IslandPath-DIMOB、Islander"}
        is_id = api_path.add_main("island", name="island_origin", params=is_params, desc="island的预测！")
        pre_key = "location,gi_id,island_start,island_end,len,method,cds_num,gene_list"
        names = []
        for sample in self.samples:
            if os.path.exists(self.output_dir + "/island/" + sample + ".GI_summary.xls"):
                names.append(sample)
                api_path.add_main_detail3(self.output_dir + "/island/" + sample + ".GI_summary.xls",
                                          self.output_dir + "/island/" + sample + "_island.fna", "island_gene_detail",
                                          is_id, pre_key, "gi_id", has_head=True, main_name='island_id',
                                          insert_d={"specimen_id": sample}, convert_dic={'island_start': int, 'island_end': int, 'len': float, 'cds_num': int})
        api_path.add_main_detail2(self.output_dir + "/island/all.island_abund.xls", "island_abund", is_id, "Type",
                                  has_head=True, main_name='island_id', names_convert=self.samples_dict)
        api_path.add_main_detail2(self.output_dir + "/island/all.island_genelist.xls", "island_detail", is_id,
                                  "Type", has_head=True, main_name='island_id', names_convert=self.samples_dict, names_co="False")

        api_path.add_num_update(self.output_dir + "/island", ".GI_summary.xls", "genomes_detail",
                                self.data_id, "data_id", "gi", self.samples)

    def run_antismash_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "Antismash"}
        anti_id = api_path.add_main("antismash", name="antismash_origin", params=pre_params, desc="antismash的预测！")
        pre_key = "cluster_id,location,type,start,end,len,cds_num,gene_list,predicted_structure,most_similar_cluster,similarity"
        names = []
        for sample in self.samples:
            if os.path.exists(self.output_dir + "/antismash/" + sample + ".antismash_anno.xls"):
                names.append(sample)
                api_path.add_main_detail3(self.output_dir + "/antismash/" + sample + ".antismash_anno.xls",
                                          self.output_dir + "/antismash/" + sample + "_antismash.fna",
                                          "antismash_gene_detail", anti_id, pre_key, "cluster_id", has_head=True, main_name='anti_id', insert_d={"specimen_id": sample}, convert_dic={'start': int, 'end': int, 'len': float, 'cds_num': int})
        api_path.add_main_detail2(self.output_dir + "/antismash/all.antismash_abund.xls", "antismash_abund", anti_id,
                                  "Type", has_head=True, main_name='anti_id', names_convert=self.samples_dict)
        api_path.add_main_detail2(self.output_dir + "/antismash/all.antismash_genelist.xls", "antismash_detail", anti_id,
                                  "Type", has_head=True, main_name='anti_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_update(self.output_dir + "/antismash", ".antismash_anno.xls", "genomes_detail", self.data_id, "data_id", "snd_metab", self.samples)

    def run_core_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "BLAST", "database": "House-keeping Gene"}
        core_id = api_path.add_main("core_gene", name="core_origin", params=pre_params, desc="core_gene的预测！")
        pre_key = "specimen_id,location,start,end,name,identity,coverage,seq"
        path = self.remote_dir + "core_gene/all.cor_gene.fa"
        api_path.add_main_detail(self.output_dir + "/core_gene/all.coregene.xls", "core_gene_detail", core_id, pre_key, has_head=True, main_table="core_gene", main_name='core_id', update_dic={"seq_path": path}, convert_dic={'start': int, 'end': int, 'identity': float, 'coverage': float})

    def run_tree_api(self):
        self.s16_sampls = ['GCF_000092565.1', 'GCF_000216055.1', 'GCF_000216075.1', 'GCF_000216095.1', 'GCF_000216175.2', 'GCF_000216195.1', 'GCF_000216235.1', 'GCF_000216335.1', 'GCF_000216355.1', 'GCF_000216375.1', 'GCF_000216395.1', 'GCF_000216415.1', 'GCF_000216435.1', 'GCF_000216475.1', 'GCF_000216495.1', 'GCF_000216515.1', 'GCF_000216535.1', 'GCF_000216555.2', 'GCF_000216035.1', 'GCF_000007685.1']
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "MEGA7", "method": "NJ", "Bootstraps": 200, "samples": sorted(self.s16_sampls)}
        tree_id = api_path.add_main("tree", name="core_16s_orign", params=pre_params, desc="16s的NJ进化树！")
        samples = get_sample_from_tree(self.output_dir + "/tree/all.16s_NJ.nwk")
        path = self.remote_dir + "tree/all.16s_NJ.nwk"
        api_path.add_main_tree(self.output_dir + "/tree/all.16s_NJ.nwk", tree_id,
                               update_dic={"type": "s16", "samples_order": samples, "path": path})

    def run_pan_api(self):
        num = len(self.samples)
        samples = ','.join(sorted(self.samples))
        if num >= 60:
            pan_params = {
                "spe_list": samples,
                "cluster_method": 'cdhit+blast+mcl',
                "soft": "Roary",
                "identity": 0.5,
                "coverage": 0.5
            }
        else:
            pan_params = {
                "spe_list": samples,
                "cluster_method": 'blast+mcl',
                "soft": "PGAP",
                "identity": 0.5,
                "coverage": 0.5,
                "pgap_method": "GF",
                "inflation": 1,
                "evalue": "1e-5",
                "score": 40
            }
        pangenomes_dir = os.path.join(self.output_dir, "pangenomes")
        cluster_dir = os.path.join(pangenomes_dir, "cluster")
        category_dir = os.path.join(pangenomes_dir, "category")
        venn_dir = os.path.join(pangenomes_dir, "venn")
        api_pan = self.api.api('bac_comp_genome.pan')
        api_category = self.api.api('bac_comp_genome.pan_category')
        api_venn = self.api.api('bac_comp_genome.pan_venn')
        cluster_data = os.path.join(cluster_dir, 'homologues_cluster.xls')
        variation_data = os.path.join(cluster_dir, 'CDS_variation_analysis.xls')
        anno_file = os.path.join(cluster_dir, 'cluster_anno.xls')
        cluster_path = self.remote_dir + "/pangenomes/cluster/homologues_cluster.xls"

        self.logger.info("正在导泛基因组分析的主表")
        pan_id = api_pan.add_pan(params=pan_params)
        self.logger.info("正在导泛基因组分析的cluster结果表")
        api_pan.add_pan_detail(pan_id, cluster_data, variation_data, anno_file, cluster_path)
        pan_graph = os.path.join(cluster_dir, 'pan_cluster_genome.xls')
        pan_detail = os.path.join(cluster_dir, 'pan_cluster.xls')
        core_graph = os.path.join(cluster_dir, 'core_cluster_genome.xls')
        core_detail = os.path.join(cluster_dir, 'core_cluster.xls')
        new_graph = os.path.join(cluster_dir, 'new_gene_cluster_genome.xls')
        new_detail = os.path.join(cluster_dir, 'new_gene_cluster.xls')
        self.logger.info("正在导公式的结果表")
        api_pan.add_pan_genome(pan_id, 'pan', pan_graph, pan_detail)
        api_pan.add_pan_genome(pan_id, 'core', core_graph, core_detail)
        api_pan.add_pan_genome(pan_id, 'newgene', new_graph, new_detail)
        self.logger.info("正在导变异分析的结果表")
        variation_path = os.path.join(cluster_dir, 'CDS_variation.xls')
        api_pan.add_pan_variation(pan_id, variation_path)
        self.logger.info("正在导分组方案合并的主表")
        category_params = {
            "category_scheme": "number1",
            "pan_id": "{}".format(pan_id)
        }
        pan_category = "core,dispensable,unique"
        category_id = api_category.add_category(pan_id, pan_category, params=category_params)
        self.logger.info("正在导分组方案合并画图的详情表")
        distribution_path = os.path.join(category_dir, "pangenome_distribution.xls")
        api_category.add_category_graph(category_id, distribution_path)
        self.logger.info("正在导分组方案合并的detail表")
        category_path = os.path.join(category_dir, "pangenome_category.xls")
        api_category.add_category_detail(category_id, category_path)
        self.logger.info("正在进行venn图分析")
        venn_path = os.path.join(venn_dir, "venn_table.xls")
        venn_params = {
            "group_detail": {"All": sorted(self.samples)},
            "pan_id": "{}".format(pan_id)
        }
        venn_id = api_venn.add_venn(params=venn_params)
        api_venn.add_venn_detail(venn_id, venn_path)

        category_name1 = "core,dispensable,unique"
        category = {"core": {"max": num, "min": num}, "dispensable": {"max": num, "min": 2},
                    "unique": {"max": 1, "min": 1}}
        value_type = False
        number1 = "number1"
        params1 = {
            "category_scheme": "number1",
        }
        api_category.add_pan_group(category, category_name1, value_type, number1, params=params1)
        category_name2 = "core,softcore,shell,cloud"
        category = {"core": {"max": 100, "min": 100}, "softcore": {"max": 100, "min": 95},
                    "shell": {"max": 95, "min": 5}, "cloud": {"max": 5, "min": 0}}
        value_type = True
        number2 = "number2"
        params2 = {
            "category_scheme": "number2",
        }
        api_category.add_pan_group(category, category_name2, value_type, number2, params=params2)

    def send_files(self):
        """
        结果放置到/upload_results
        """
        pass
