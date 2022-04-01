# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'
# 20181130

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
from bson.son import SON
import copy
import shutil
from biocluster.file import download,exists
import os



class CommonApi(Base):
    def __init__(self, bind_object):
        super(CommonApi, self).__init__(bind_object)
        self._project_type = "bacgenome"

    @report_check
    def add_main(self,table, name=None, params=None, others=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Job has been finished",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

        if others:
            for k in others.keys():      #others {mongo字段名: mongo字段值}
                insert_data[k] = others[k]

        collection = self.db[table]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_main_detail(self, infile,detail_table, main_id, mongo_key, has_head =True, main_name='main_id',
                        main_table=None, update_dic=None, other_dic=None):
        data_list = []
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        key_list = mongo_key.split(',')
        num = len(key_list)
        with open(infile, 'rb') as r:
            if has_head:
                r.readline()
            for line in r:
                if other_dic:
                    insert_data = copy.deepcopy(other_dic)
                else:
                    insert_data = {}
                insert_data[main_name]=main_id
                spline = line.strip("\n").split("\t")
                if len(spline) == num:
                    for i in range(num):
                        if key_list[i] !='':
                            insert_data[key_list[i]] = spline[i]
                else:
                    self.bind_object.logger.info(str(len(spline))+' vs '+ str(num))
                    self.bind_object.set_error("data incomplete: %s ", variables=(line), code="52803302")
                data_son = SON(insert_data)
                data_list.append(data_son)
        try:
            collection = self.db[detail_table]
            self.bind_object.logger.info("开始导入%s数据"%infile)
            collection.insert_many(data_list)
            if update_dic:
                main_table = self.db[main_table]
                main_table.update({'_id':main_id},{'$set':update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % e)
            self.bind_object.set_error("import CommonApi data ERROR", code="52803301")
        else:
            self.bind_object.logger.info("导入CommonApi数据成功")

    @report_check
    def update_data(self,main_id,db_table,update_dic):
        """
        功能是跟新数据库数据
        :param download_dir:
        :param task_id:
        :return:
        """
        try:
            collection = self.db[db_table]
            collection.update({'_id': main_id}, {'$set':update_dic})
        except Exception, e:
            self.bind_object.logger.info("导入CommonApi数据出错:%s" % db_table)
            self.bind_object.set_error("import CommonApi update_data ERROR")
        else:
            self.bind_object.logger.info("导入CommonApi数据成功")

    @report_check
    def down_seq_files(self, download_dir, task_id, sample_list):
        """
        功能是根据任务下载组装序列文件
        :param download_dir:
        :param task_id:
        :return:
        """
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)
        else:
            shutil.rmtree(download_dir)
            os.mkdir(download_dir)
        assemble = self.db.assemble
        assemble_seq = self.db.assemble_seq
        assemble_id = assemble.find_one({"task_id" :task_id})["main_id"]
        analysis_type = assemble.find_one({"task_id" :task_id})["analysis_type"]
        if analysis_type == "uncomplete":
            for each in sample_list:
                path = assemble_seq.find_one({"assemble_id" : assemble_id,"specimen_id" : each})["seq_path"]
                s3_path = path + "/" + each + "/assembly_predict/assembly/" + each + "_scaf.fna"
                new_path = os.path.join(download_dir, each + "_scaf.fna")
                download(s3_path, new_path)
        else:
            for sample in sample_list:
                complete_list = []
                complete_list_uniq = []
                if not os.path.exists(download_dir + '/' + sample):
                    os.mkdir(download_dir + '/' + sample)
                path = assemble_seq.find_one({"assemble_id": assemble_id, "specimen_id": sample})["seq_path"]
                sample_data = assemble_seq.find({"assemble_id": assemble_id,"specimen_id":sample})
                for i in sample_data:
                    complete_list.append(i)
                for each in complete_list:
                    if each:
                        if each["type"] not in complete_list_uniq:
                            complete_list_uniq.append(each["type"])
                for a in complete_list_uniq:
                    s3_path = path + "/" + sample + "/assembly_predict/assembly/seq_dir/" + a + ".fasta"
                    new_path = download_dir + '/' + sample + '/' + sample + '_' + a + '_comp.fna'
                    download(s3_path, new_path)
        return (download_dir, analysis_type)


    @report_check
    def down_gene_files(self, gene_dir, task_id, sample_list):
        """
        此功能是下载gff、fnn、faa文件
        :param gene_dir:
        :param task_id:
        :return:
        """
        if not os.path.exists(gene_dir):
            os.mkdir(gene_dir)
        else:
            shutil.rmtree(gene_dir)
            os.mkdir(gene_dir)
        gene_predict = self.db.gene_predict
        file_list = gene_predict.find_one({"task_id" :task_id})["file_path"]
        if file_list:
            for sample in sample_list:
                sample_path = os.path.join(gene_dir, sample)
                if os.path.exists(sample_path):
                    shutil.rmtree(sample_path)
                os.mkdir(sample_path)
                s3_fnn_path = file_list[0] + '/' + sample + file_list[1] + '/' + sample + file_list[2] + 'fnn'
                s3_faa_path = file_list[0] + '/' + sample + file_list[1] + '/' + sample + file_list[2] +'faa'
                s3_gff_path = file_list[0] + '/' + sample + file_list[1] + '/' + sample + file_list[2] +'gff'
                self.bind_object.logger.info(s3_fnn_path)
                self.bind_object.logger.info(s3_faa_path)
                self.bind_object.logger.info(s3_gff_path)
                if not os.path.exists(os.path.join(sample_path, sample + '.fnn')):
                    download(s3_fnn_path, os.path.join(sample_path, sample + '.fnn'))
                if not os.path.exists(os.path.join(sample_path, sample +'.faa')):
                    download(s3_faa_path, os.path.join(sample_path, sample +'.faa'))
                if not os.path.exists(os.path.join(sample_path, sample +'.gff')):
                    download(s3_gff_path, os.path.join(sample_path, sample +'.gff'))
        return gene_dir

    @report_check
    def get_seq_files(self, download_dir, task_id):
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)
        else:
            shutil.rmtree(download_dir)
            os.mkdir(download_dir)
        self.faa_dir = download_dir
        all_data_detail = []
        sample_id = []
        assemble = self.db.assemble
        assemble_seq = self.db.assemble_seq
        assemble_id = assemble.find_one({"task_id": task_id})["main_id"]
        analysis_type = assemble.find_one({"task_id": task_id})["analysis_type"]
        all_data = assemble_seq.find({"assemble_id": assemble_id})
        if analysis_type == "uncomplete":
            for i in all_data:
                all_data_detail.append(i)
            for i in range(len(all_data_detail)):
                if all_data_detail[i]["specimen_id"] not in sample_id:
                    sample_id.append(all_data_detail[i]["specimen_id"])
            for each in sample_id:
                path = assemble_seq.find_one({"assemble_id": assemble_id, "specimen_id": each})["seq_path"]
                s3_path = path + "/" + each + "/assembly_predict/assembly/" + each + "_scaf.fna"
                new_path = os.path.join(self.faa_dir, each + "_scaf.fna")
                # print s3_path, new_path
                download(s3_path, new_path)
            # os.link("/mnt/ilustre/users/sanger-dev/home/zhaozhigang/BacgenomeV3/repeatmasker/workflow-test/FJXPY27_L6_scaf.fna",new_path)
        else:
            complete_list = []
            complete_list_uniq = []
            for i in all_data:
                all_data_detail.append(i)
            for i in range(len(all_data_detail)):
                if all_data_detail[i]["specimen_id"] not in sample_id:
                    sample_id.append(all_data_detail[i]["specimen_id"])
            for sample in sample_id:
                if not os.path.exists(self.faa_dir + '/' + sample):
                    os.mkdir(self.faa_dir + '/' + sample)
                path = assemble_seq.find_one({"assemble_id": assemble_id, "specimen_id": sample})["seq_path"]
                sample_data = assemble_seq.find({"assemble_id": assemble_id, "specimen_id": sample})
                for i in sample_data:
                    complete_list.append(i)
                for each in complete_list:
                    if each:
                        if each["type"] not in complete_list_uniq:
                            complete_list_uniq.append(each["type"])
                for a in complete_list_uniq:
                    s3_path = path + "/" + sample + "/assembly_predict/assembly/seq_dir/" + a + ".fasta"
                    new_path = self.faa_dir + '/' + sample + '/' + sample + '_' + a + '_comp.fna'
                    # print s3_path,new_path
                    download(s3_path, new_path)
                complete_list = []
                complete_list_uniq = []
                # os.link("/mnt/ilustre/users/sanger-dev/home/zhaozhigang/BacgenomeV3/repeatmasker/workflow-test/Chromosome.fasta",self.faa_dir + '/' + sample + '/' + sample + '_' + "Chromosome" + '_comp.fna')
                # os.link("/mnt/ilustre/users/sanger-dev/home/zhaozhigang/BacgenomeV3/repeatmasker/workflow-test/Plasmid.fasta",self.faa_dir + '/' + sample + '/' + sample + '_' + "Plasmid" + '_comp.fna')

        return analysis_type, sample_id

    @report_check
    def get_gene_gff(self, gff_dir, task_id, sample):
        """
        此功能是下载gff文件
        :param gff_dir:
        :param task_id:
        :return:
        """
        if not os.path.exists(gff_dir):
            os.mkdir(gff_dir)
        else:
            shutil.rmtree(gff_dir)
            os.mkdir(gff_dir)
        gene_predict = self.db.gene_predict
        file_list = gene_predict.find_one({"task_id": task_id})["file_path"]
        s3_gff_path = file_list[0] + '/' + sample + file_list[1] + '/' + sample + file_list[2] + 'gff'
        self.bind_object.logger.info(s3_gff_path)
        if not os.path.exists(os.path.join(gff_dir, sample + '.gff')):
            download(s3_gff_path, os.path.join(gff_dir, sample + '.gff'))
        return os.path.join(gff_dir, sample + '.gff')

    @report_check
    def get_gene_ngloc(self, dir, task_id, sample):
        """
        此功能是下载faa文件和pfam文件
        :param gff_dir:
        :param task_id:
        :return:
        """
        gene_predict = self.db.gene_predict
        file_list = gene_predict.find_one({"task_id": task_id})["file_path"]
        s3_faa_path = file_list[0] + '/' + sample + file_list[1] + '/' + sample + file_list[2] + 'faa'
        self.bind_object.logger.info(s3_faa_path)
        if not os.path.exists(os.path.join(dir, sample + '.faa')):
            download(s3_faa_path, os.path.join(dir, sample + '.faa'))
        assemble = self.db.assemble
        analysis_type = assemble.find_one({"task_id": task_id})["analysis_type"]
        if analysis_type == "uncomplete":
            s3_pfam_path = file_list[0] + '/' + sample + "/annotation/Pfam/" +sample+ "_anno_pfam.xls"
        else:
            s3_pfam_path = file_list[0] + '/' + sample + "/annotation/Pfam/" + sample + "_whole_genome_anno_pfam.xls"
        if not os.path.exists(os.path.join(dir, sample + '.pfam.xls')):
            download(s3_pfam_path, os.path.join(dir, sample + '.pfam.xls'))
        return os.path.join(dir, sample + '.faa'),os.path.join(dir, sample + '.pfam.xls')

    @report_check
    def get_table_file(self, dir, task_id, sample, type):
        """
        此功能是下载type类型的文件
        :param gff_dir:
        :param task_id:
        :return:
        """
        gene_predict = self.db.gene_predict
        file_list = gene_predict.find_one({"task_id": task_id})["file_path"]
        if type == "prephage":
            s3_pfam_path = file_list[0] + sample + "/mobile_elements/prephage/" + sample + "_prephage_summary.xls"
            self.bind_object.logger.info(s3_pfam_path)
            if not os.path.exists(os.path.join(dir, sample + '.prephage_summary.xls')):
                download(s3_pfam_path, os.path.join(dir, sample + '.prephage_summary.xls'))
            return  os.path.join(dir, sample + '.prephage_summary.xls')
        elif type == "island":
            s3_pfam_path = file_list[0]  + sample + "/mobile_elements/Genomic_Islands/" + sample + "_GI_summary.xls"
            if not os.path.exists(os.path.join(dir, sample + '.island_summary.xls')):
                download(s3_pfam_path, os.path.join(dir, sample + '.island_summary.xls'))
            return os.path.join(dir, sample + '.island_summary.xls')

    @report_check
    def get_anno_file(self, dir, task_id, samples, type):
        """
        此功能是下载type类型的文件
        :param gff_dir:
        :param task_id:
        :return:
        """
        gene_predict = self.db.gene_predict
        file_list = gene_predict.find_one({"task_id": task_id})["file_path"]
        for sample in list(samples.keys()):
            if type == "COG":
                s3_pfam_path = file_list[0] + '/' + sample + "/annotation/COG/" + sample + "_cog_summary.xls"
                if not os.path.exists(os.path.join(dir, samples[sample] + '.cog.xls')):
                    download(s3_pfam_path, os.path.join(dir, samples[sample] + '.cog.xls'))
            elif type == "KEGG":
                s3_pfam_path = file_list[0] + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_level_stat.xls"
                if not os.path.exists(os.path.join(dir, samples[sample] + '.kegg.xls')):
                    download(s3_pfam_path, os.path.join(dir, samples[sample] + '.kegg.xls'))
            elif type == "GO":
                s3_pfam_path = file_list[0] + '/' + sample + "/annotation/GO/" + sample + "_go_statistics.xls"
                if not os.path.exists(os.path.join(dir, samples[sample] + '.go.xls')):
                    download(s3_pfam_path, os.path.join(dir, samples[sample] + '.go.xls'))
        return dir

    def add_geneset_ppi(self, main_id):
        try:
            main_id = ObjectId(main_id)
        except:
            pass
        main_table = self.db['sg_tool_lab_geneset_ppi']
        main_table.update({'_id': main_id}, {'$set': {"status":"end", "main_id":main_id}})