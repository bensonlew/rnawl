# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
from biocluster.config import Config
from Bio import SeqIO
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class TfPredictWorkflow(Workflow):
    """
    TfPredict description
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TfPredictWorkflow, self).__init__(wsheet_object)
        options = list()
        for each in self._sheet.options():
            options.append(dict(name=each, type="string"))
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.tf_predict = self.add_tool("rna.tf_predict")
        self.tf_predict = self.add_module("ref_rna_v2.tf_predict")
        self.dump_tool = self.api.api("whole_transcriptome.tf_predict")
        self.known_tf_trans = dict()
        software_dir = self.config.SOFTWARE_DIR
        self.tfdb = software_dir + "/database/TFDB/"
        self.short_name_dict = dict()
        self.plantname = dict()
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/03 TF_Analysis')
        self.inter_dirs = []

    def run(self):
        self.tf_predict.on("end", self.set_db)
        self.get_run_log()
        self.get_known_tf()
        self.filter_fasta()
        self.run_tf_predict()
        super(TfPredictWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="tf_predict", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(TfPredictWorkflow, self).send_log(data)

    def get_plant_name(self, name):
        species_short_name = self.tfdb + '/plant_tf_pep/species_short_name.list'
        short_name_dict = dict()
        with open(species_short_name) as f:
            for line in f:
                if not line.strip():
                    continue
                short_name, sp = line.strip().split('\t')
                self.plantname[short_name] = sp
        if name in self.plantname:
            return self.plantname[name]
        else:
            return name

    def get_link(self, gene_id):
        species_short_name = self.tfdb + '/plant_tf_pep/species_short_name.list'
        short_name_dict = dict()
        with open(species_short_name) as f:
            for line in f:
                if not line.strip():
                    continue
                short_name, sp = line.strip().split('\t')
                short_name_dict[sp.lower()] = short_name

        if self.option("s") == "plant":
            short_name = self.option("organism")
            ref_link = "http://planttfdb.cbi.pku.edu.cn/tf.php?sp={}&did={}".format(short_name, gene_id)
        else:
            ref_link = "http://bioinfo.life.hust.edu.cn/AnimalTFDB#!/tf_gene_info?tf={}".format(gene_id)
        return ref_link

    def get_species_dict(self):
        species_short_name = self.tfdb + '/plant_tf_pep/species_short_name.list'
        with open(species_short_name) as f:
            for line in f:
                if not line.strip():
                    continue
                short_name, sp = line.strip().split('\t')
                self.short_name_dict[sp.lower()] = short_name

    def get_known_tf(self):
        client_ref = Config().get_mongo_client(mtype="ref_rna", ref=True)
        mongodb = client_ref[Config().get_mongo_dbname("ref_rna", ref=True)]
        known_tfdb_coll = mongodb.known_tf
        s2db = {
            "animal": "AnimalTFDB",
            "plant": "PlantTFDB"
        }
        species_name = self.option('organism')
        if self.option('s') == "plant":
            species_name = self.get_plant_name(species_name)

        self.logger.info("species_name is" + species_name + "db is" + s2db[self.option('s')])
        known_tfs = known_tfdb_coll.find({"specie": species_name, "db": s2db[self.option('s')]})
        if self.option("s") == "plant":
            self.get_species_dict()

        for tf in known_tfs:
            print tf
            self.known_tf_trans[tf["transcript_id"]] = {
                'gene_id': tf['gene_id'],
                'tf_id': tf['tf_id'],
                'family': tf['family']
            }


    def filter_fasta(self):
        seq_records = SeqIO.parse(self.option('seq_db'), 'fasta')
        with open(self.work_dir + "/tf_filter_seq.fa" , 'w') as f, open(self.work_dir + "/known_tf.txt", "w") as known_f:
            known_f.write("query_id\tfamily\tdomain\tdomain_link\tdescription\te_value\tscore\tblast_hit\thit_link\thit_family\thit_pident\thit_evalue\n")
            for seq_record in seq_records:
                seq_seq = seq_record.seq
                seq_name = seq_record.name
                if seq_name not in self.known_tf_trans:
                    f.write('>{}\n{}\n'.format(seq_name, seq_seq))
                else:
                    link = ""
                    if self.option('s') == "plant":
                        link = self.get_link(seq_name)
                    else:
                        link = self.get_link(self.known_tf_trans[seq_name]['gene_id'])

                    known_f.write("\t".join([
                        seq_name,
                        self.known_tf_trans[seq_name]['family'],
                        "",
                        "",
                        "",
                        "",
                        "",
                        self.known_tf_trans[seq_name]['tf_id'],
                        link,
                        self.known_tf_trans[seq_name]['family'],
                        "",
                        ""
                    ]) + "\n")

    def run_tf_predict(self):
        options = dict(
            s=self.option('s'),
            organism=self.option('organism'),
            blast_all=self.option('blast_all'),
            seqfile=self.work_dir + "/tf_filter_seq.fa",
            E=self.option('E'),
            # domE=self.option('domE'),
            evalue=self.option('evalue')
        )
        self.tf_predict.set_options(options)
        self.tf_predict.run()

    def set_db(self):
        """
        dump data to db
        """
        workflow_output = self.get_workflow_output_dir()
        # modify result info
        seq_annot_pd = pd.read_table(self.work_dir + '/seq_annot.xls', header=0, index_col=0)
        known_tf =  pd.read_table(self.work_dir + "/known_tf.txt", index_col=0, header=0)
        # known_tf.set_index("query_id")
        tf_result = pd.read_table(self.tf_predict.output_dir + '/final_tf_predict.xls', header=0, index_col=0)
        tf_result = known_tf.append(tf_result)
        seq_annot_pd = tf_result.join(seq_annot_pd)
        seq_annot_pd.index.name = "transcript_id"
        seq_annot_pd.to_csv(self.tf_predict.output_dir + '/final_tf_predict.xls', header=True, index=True, sep='\t')
        if self.option('s').lower() == 'plant':
            ## modified by shicaiping at 20180809 because of TF_id and meme file are not identical
            TF_id_2_meme = dict()
            TF_id_2_meme_file =  self.config.SOFTWARE_DIR + "/database/TFDB/plant_tf_motif/TF_id_2_MEME.txt"
            with open(TF_id_2_meme_file, "r") as f:
                for line in f:
                    item = line.strip("\n").split("\t")
                    TF_id_2_meme[item[0]] = item[1]
            '''
            known_meme_tf_list_file = self.config.SOFTWARE_DIR + "/database/TFDB/plant_tf_motif/plant.meme.gene.list"
            with open(known_meme_tf_list_file) as f:
                known_meme_tf_list = [x.strip() for x in f]
                known_meme_tf_list = tuple(known_meme_tf_list)
            '''
        else:
            TF_id_2_meme = dict()
        tmp_tuple = zip(seq_annot_pd['gene_id'], seq_annot_pd['blast_hit'], seq_annot_pd['family'])
        tf_select_values = list()
        tf_select_shows = list()
        for x in tmp_tuple:
            if x[1] in TF_id_2_meme.keys():
                show = x[0]+'('+TF_id_2_meme[x[1]]+','+x[2]+')'
                if show not in tf_select_shows:
                    tf_select_shows.append(show)
                    tf_select_values.append(x[0]+"|"+TF_id_2_meme[x[1]])
            '''
            if x[1] not in known_meme_tf_list:
                if '.' in x[1]:
                    tmp_list = x[1].split(".")
                    for i in range(0,len(tmp_list)):
                        tmp_x = '.'.join(tmp_list[0:len(tmp_list)-i-1])
                        if tmp_x in known_meme_tf_list:
                            show = x[0]+'('+x[1]+','+x[2]+')'
                            if show not in tf_select_shows:
                                tf_select_shows.append(show)
                                tf_select_values.append(x[0]+"|"+x[1])
                            break
                        else:
                            continue
                else:
                    continue
            else:
                show = x[0]+'('+x[1]+','+x[2]+')'
                if show not in tf_select_shows:
                    tf_select_shows.append(show)
                    tf_select_values.append(x[0]+"|"+x[1])
            '''
        # dump
        self.dump_tool.add_tf_predict_detail(self.tf_predict.output_dir, self.option("main_id"))
        self.dump_tool.update_db_record('tf_predict', self.option('main_id'),
                                        output_dir=workflow_output, status="end",
                                        tf_select_values=tf_select_values,
                                        tf_select_shows=tf_select_shows,)
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.tf_predict.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.tf_predict.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.tf_predict.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.tf_predict.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/03 TF_Analysis", "", "转录因子分析结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "转录因子预测文件", 0],
            ["./final_tf_predict.xls", "", "转录因子预测详情表", 0],
            ["./predicted_TFs.fa", "", "转录因子序列信息", 0],
            ["./domain_predict.txt", "", "蛋白功能域预测详情表", 0],
            ["./diamond.out.txt", "", "Blast 比对详情表", 0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(TfPredictWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
