# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 20180409

from api_base import ApiBase
import re


class CnvCompare(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目导表
        """
        super(CnvCompare, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_sg_cnv_compare_detail(self, compare_id, file_path, analysis_name=None, analysis_mode='multiple'):
        """
        AH03_AH19_same.xls
        :param compare_id:
        :param file_path:
        :param analysis_name:
        :param analysis_mode:
        :return:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        header = []
        header_dict = {}
        title = open(file_path, 'r').readlines()[0].strip().split("\t")
        for i in range(5, len(title) - 2):
            sam = '_'.join(title[i].split('_')[:-1])
            if re.match('.*_genotype$', title[i]):
                header.append("{} Genotype".format(sam))
                header_dict[title[i]] = "{} Genotype".format(sam)
            elif re.match('.*_pvalue$', title[i]):
                header.append("{} Pvalue".format(sam))
                header_dict[title[i]] = "{} Pvalue".format(sam)
        with open(file_path, 'r') as r:
            lines = r.readlines()[1:]
            for line in lines:
                tmp = line.strip().split('\t')
                insert_data = {
                    "cnv_id": tmp[0] + ':' + tmp[1] + '-' + tmp[2],
                    "compare_id": compare_id,
                    "chr": tmp[0],
                    "start": int(tmp[1]),
                    "end": int(tmp[2]),
                    "len": int(float(tmp[3])),
                    "type": tmp[4],
                    "gene_num": tmp[-2],
                    "gene": tmp[-1],
                    "name": analysis_name,
                    "search_id": tmp[0] + ':' + tmp[1] + '-' + tmp[2]
                }
                for i in range(5, len(title) - 2):
                    insert_data.update({title[i]: tmp[i]})
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_cnv_compare_detail", data_list)
            if analysis_mode == 'multiple':
                self.update_db_record('sg_cnv_compare', {"_id": compare_id}, {"header": header,
                                                                              "header_dict": header_dict})

    def add_sg_cnv_compare_stat(self, compare_id, file_path, analysis_name=None):
        """
        AH03_AH19_same.stat.xls
        :param compare_id:
        :param file_path:
        :param analysis_name:
        :return:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            r.next()
            for line in r:
                tmp = line.strip().split('\t')
                insert_data = {
                    "chr": tmp[0],
                    "compare_id": compare_id,
                    "del": tmp[1],
                    "dup": tmp[2],
                    "gene": tmp[3],
                    "name": analysis_name
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_cnv_compare_stat", data_list)

    def update_subname(self, compare_id, subname, vcf_path):
        compare_id = self.check_objectid(compare_id)
        self.update_db_record('sg_cnv_compare', {"_id": compare_id}, {"subname": subname, 'vcf_path': vcf_path})

if __name__ == "__main__":
    a = CnvCompare(None)
    a.add_sg_cnv_compare_detail("1q11", "/mnt/ilustre/users/sanger-dev/workspace/20190415/CnvCompare_sanger_85433_0415175703623375_2781/CnvCompare/output/multiple_same.xls")
