# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import os
import pandas as pd
import unittest
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob
import shutil
import traceback

class SangerChartWorkflow(Workflow):
    """
    下载sanger在线交互页面的pdf图片
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SangerChartWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='url', type='string'),
            dict(name="userpass", type="string"),
            dict(name="mode", type="string", default="all"),
            dict(name="imageversion", type="string", default="latest"),
            dict(name="chart_main_id", type='string'),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("ref_rna_v2.download_sangerpdf")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/00 Charts')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]).encode('utf-8'),
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
        super(SangerChartWorkflow, self).send_log(data)

    def run(self):
        self.tool.on("end", self.set_db)
        self.run_tool()
        super(SangerChartWorkflow, self).run()

    def set_db(self):
        """
        更新主表状态
        """
        api_base = self.api.api("ref_rna_v2.api_base")
        api_base.update_db_record('sanger_chart', self.option('chart_main_id'), status="end")
        self.end()



    def move_file(self, src_path, dst_path, file):
    #     print 'from : ',src_path
    #     print 'to : ',dst_path
        try:
            # cmd = 'chmod -R +x ' + src_path
            # os.popen(cmd)
            f_src = os.path.join(src_path, file)
            if not os.path.exists(dst_path):
                os.mkdir(dst_path)
            f_dst = os.path.join(dst_path, file)
            shutil.move(f_src, f_dst)
        except Exception as e:
            print 'move_file ERROR: ',e
            traceback.print_exc()


    def split_dirs_files_to_seperate_child_dirs(self, t_dir):


        # t_dir = "./exp_graph/"
        name_detail_dic = {}
        with open(os.path.join(t_dir,'name_detail.txt')) as f:
            for line in f:
                pdf_info = re.split('\s',line.strip())
                pdf_name = pdf_info[0].strip()
                name_detail_dic[pdf_name] = "no_table_id"
                if pdf_name.endswith('.pdf'):
                    for ele in pdf_info:
                        if ele.strip().startswith('table_id-'):
                            name_detail_dic[pdf_name] = ele.lstrip('table_id-')

        if len(list(set(name_detail_dic.values()))) > 1:
            for file,path in name_detail_dic.items():
                if path is not "no_table_id":
                    self.move_file(t_dir, os.path.join(t_dir, path), file)

            origin_name_detail = open(os.path.join(t_dir,'name_detail.txt')).read().split('\n')
            origin_name_detail_cp = open(os.path.join(t_dir,'name_detail.txt')).read().split('\n')
            for p in list(set(name_detail_dic.values())):
                if p is not "no_table_id":
                    child_name_detail = os.path.join(t_dir, p, 'name_detail.txt')
                    with open(child_name_detail,'w') as w:
                        for line in origin_name_detail:
                            if p in line:
                                w.write(line+'\n')
                                origin_name_detail_cp.remove(line)
            origin_name_detail_cp = [x for x in origin_name_detail_cp if x != '']
            if origin_name_detail_cp:
                with open(os.path.join(t_dir,'name_detail.txt'),'w') as w:
                    w.write('\n'.join(origin_name_detail_cp))
            else:
                os.remove(os.path.join(t_dir,'name_detail.txt'))


    def filePath_to_fileInfo(self, pdf_file_path):
        pdf_file_path_dir = os.path.split(pdf_file_path)[0]
        pdf_file_path_file = os.path.split(pdf_file_path)[1]
        name_detail_list = open(os.path.join(pdf_file_path_dir, "name_detail.txt")).read().split('\n')
        pdf_detail = ""
        for line in name_detail_list:
            if line.startswith(pdf_file_path_file):
                pdf_detail = line
        return pdf_detail


    def end(self):
        if self.option('mode') == "single":
            #-------------------------- 目录重新设置层级结构 -------------------------------
            #-------------------------- 目录重新设置层级结构 -------------------------------
            target_dir = os.path.join(self.tool.output_dir,'pdf')
            # target_dir = './pdf_5/'

            for d in os.listdir(target_dir):
                #print os.path.join(target_dir, d)
                self.split_dirs_files_to_seperate_child_dirs(os.path.join(target_dir, d))

            dirr_first_child_dir = os.path.join(target_dir, os.listdir(target_dir)[0])
            #print dirr_first_child_dir
            for dirr_second_child_dir in os.listdir(dirr_first_child_dir):
                dirr_second_child_dir_path = os.path.join(dirr_first_child_dir, dirr_second_child_dir)
                #print dirr_second_child_dir_path,target_dir
                os.system("mv {} {}".format(dirr_second_child_dir_path.replace(' ','\ '), target_dir.replace(' ','\ ')))
            shutil.rmtree(dirr_first_child_dir)

            #-------------------------- 目录重新设置层级结构 -------------------------------
            #-------------------------- 目录重新设置层级结构 -------------------------------


            #-------------------------- 生成 add_relpath_rules_list -------------------------------
            all_files = []
            all_name_detail = {}
            all_pdf_full_path = []
            all_pdf_full_path_and_info = {}
            for root, dirs, files in os.walk(target_dir):
                if "name_detail.txt" in files:
                    all_pdf_full_path_and_info[os.path.join(root, "name_detail.txt").split(target_dir)[1].strip('/')] = "项目结果图片说明文档"
                    for file in os.listdir(root):
                        if file.endswith('.pdf'):
                            #print os.path.join(root, file)
                            all_pdf_full_path.append(os.path.join(root, file))


            for pdf in all_pdf_full_path:
                all_pdf_full_path_and_info[pdf.split(target_dir)[1].strip('/')] = self.filePath_to_fileInfo(pdf)
                
            add_relpath_rules_list_2 = []
            for pdf_path,pdf_info in all_pdf_full_path_and_info.items():
                pdf_info_dropTableid = '\t'.join(filter(None, [ele for ele in re.split('\s',pdf_info.strip()) if not ele.strip().startswith('table_id-') and not ele.strip().endswith('.pdf')]))
                add_relpath_rules_list_2.append([pdf_path, 'pdf', pdf_info_dropTableid.replace('\t','____').replace(' ','____').strip('-').strip('____'), 0])
                
            #-------------------------- 生成 add_relpath_rules_list -------------------------------

            result_dir = self.add_upload_dir(os.path.join(self.tool.output_dir,'pdf'))
            self.inter_dirs = [
                ["00 Charts", "", "项目结果图片",0],
            ]
            add_relpath_rules_list = [
                ['.', '', '项目结果图片', 0],
            ]
            # self.logger.info("xxxxxxxx____add_relpath_rules_list_2:\n"+str(add_relpath_rules_list_2))
            add_relpath_rules_list.extend(add_relpath_rules_list_2)
            result_dir.add_relpath_rules(add_relpath_rules_list)
            super(SangerChartWorkflow, self).end()

            # all_files = []
            # all_name_detail = {}
            # for root, dirs, files in os.walk(os.path.join(self.tool.output_dir,'pdf')):
            #     if dirs == []:
            #         for file_name in files:
            #             if file_name[-3:] == 'pdf':
            #                 all_files.append([root,file_name])
            #             if file_name == 'name_detail.txt':
            #                 name_detail = open(os.path.join(root,file_name)).read().strip().split('\n')
            #                 for i in name_detail:
            #                     if i.split('\t')[0]:
            #                         all_name_detail[root.split('/')[-1]+'/'+i.split('\t')[0]] = ''.join(i.split('\t')[1:]).replace('__','  ').replace('-  ','')
            # add_relpath_rules_list = []
            # add_relpath_rules_list.append(["name_detail.txt", 'txt', "项目结果图片说明文档", 0])
            # for i in all_files:
            #     pdf_name = str(i[0].split('/')[-1]+'/'+i[1])
            #     if pdf_name in all_name_detail.keys():
            #         this_pdf_info = all_name_detail[pdf_name]
            #         this_pdf_info_drop_tableid = '\t'.join(filter(None, [ele for ele in re.split('\s',this_pdf_info.strip()) if not ele.strip().startswith('table_id-')]))
            #         add_relpath_rules_list.append([pdf_name, 'pdf', this_pdf_info_drop_tableid, 0])


        if self.option('mode') == "all":
            #-------------------------- 目录重新设置层级结构 -------------------------------
            #-------------------------- 目录重新设置层级结构 -------------------------------
            target_dir = os.path.join(self.tool.output_dir,'pdf')

            dirs_dict = {
            "specimen_before" : "02_QC",
            "specimen_after" : "02_QC",
            "assessmentsaturation" : "03_Align",
            "assessmentcoverage" : "03_Align",
            "assessmentdistribution" : "03_Align",
            "assessmentchrom" : "03_Align",
            "transcripts" : "04_Assemble",
            "annotationstat" : "05_Annotation",
            "expgraph" : "06_Express",
            "expvenn" : "06_Express",
            "expcorr" : "06_Express",
            "exppca" : "06_Express",
            "expbatch" : "06_Express",
            "diff_detail" : "07_DiffExpress",
            "genesetvenn" : "08_GeneSet",
            "genesetcluster" : "08_GeneSet",
            "genesetcog" : "08_GeneSet",
            "genesetgo" : "08_GeneSet",
            "genesetkegg" : "08_GeneSet",
            "genesetgo_rich" : "08_GeneSet",
            "genesetkegg_rich" : "08_GeneSet",
            "genesetcirc" : "08_GeneSet",
            "genesetgo_acyclic" : "08_GeneSet",
            "genesetcorrsf" : "08_GeneSet",
            "splicingrmats_count" : "09_AS",
            "splicingrmats_stat" : "09_AS",
            "snp_collect" : "10_SNP",
            "wgcnaprepare" : "11_WGCNA",
            "wgcnamodule" : "11_WGCNA",
            "wgcnarelate" : "11_WGCNA",
            "stem" : "12_Time_series_Analysis",
            "tfstat" : "13_TF_Analysis",
            "ppinetwork" : "14_PPI_Analysis"
            }

            # target_dir = './pdf_5/'

            for d in os.listdir(target_dir):
                if d not in dirs_dict.keys():
                    shutil.rmtree(os.path.join(target_dir, d))

            for d in os.listdir(target_dir):
            #     print os.path.join(target_dir, d)
                self.split_dirs_files_to_seperate_child_dirs(os.path.join(target_dir, d))
                
            for old_dir,new_dir in dirs_dict.items():
                old_path = os.path.join(target_dir, old_dir)
                new_path = os.path.join(target_dir, new_dir)
                if not os.path.exists(new_path):
                    os.mkdir(new_path)
                #shutil.move(old_path, new_path)
                os.system("mv {} {}".format(old_path, new_path.replace(' ','\ ')))

            need_minus_one_layer_dirs = ["04_Assemble", "05_Annotation", "07_DiffExpress", "10_SNP", "12_Time_series_Analysis", "13_TF_Analysis", "14_PPI_Analysis"]
            for directory in need_minus_one_layer_dirs:
                dirr = os.path.join(target_dir, directory)
                if len(os.listdir(dirr)) > 0:
                    dirr_first_child_dir = os.path.join(dirr, os.listdir(dirr)[0])
                    #print dirr_first_child_dir
                    for dirr_second_child_dir in os.listdir(dirr_first_child_dir):
                        dirr_second_child_dir_path = os.path.join(dirr_first_child_dir, dirr_second_child_dir)
                        #print dirr_second_child_dir_path,dirr
                        os.system("mv {} {}".format(dirr_second_child_dir_path.replace(' ','\ '), dirr.replace(' ','\ ')))
                    shutil.rmtree(dirr_first_child_dir)

            dirs_need_rename = {
                "02_QC/specimen_before" : "02_QC/Rawdata",
                "02_QC/specimen_after" : "02_QC/Cleandata",
                "06_Express/expgraph" : "06_Express/expdistribution",
                "08_GeneSet/genesetvenn" : "08_GeneSet/01_GenesetVenn",
                "08_GeneSet/genesetcluster" : "08_GeneSet/02_Cluster_Analysis",
                "08_GeneSet/genesetcog" : "08_GeneSet/03_COG_Annotation",
                "08_GeneSet/genesetgo" : "08_GeneSet/04_GO_Annotation",
                "08_GeneSet/genesetkegg" : "08_GeneSet/05_KEGG_Annotation",
                "08_GeneSet/genesetgo_rich" : "08_GeneSet/06_GO_Enrich",
                "08_GeneSet/genesetkegg_rich" : "08_GeneSet/07_KEGG_Enrich",
                "08_GeneSet/genesetcirc" : "08_GeneSet/08_Enrich_Circ",
                "08_GeneSet/genesetgo_acyclic" : "08_GeneSet/09_GO_Dag",
                "08_GeneSet/genesetcorrsf" : "08_GeneSet/10_Exp_Network",
                "09_AS/splicingrmats_count" : "09_AS/AS_count",
                "09_AS/splicingrmats_stat" : "09_AS/AS_stat",
            }
            for src,dst in dirs_need_rename.items():
                src_ = os.path.join(target_dir, src)
                dst_ = os.path.join(target_dir, dst)
                # print src_,dst_
                try:
                    os.rename(src_, dst_)
                except Exception,err:
                    print err
            #-------------------------- 目录重新设置层级结构 -------------------------------
            #-------------------------- 目录重新设置层级结构 -------------------------------


            #-------------------------- 生成 add_relpath_rules_list -------------------------------
            all_files = []
            all_name_detail = {}
            all_pdf_full_path = []
            all_pdf_full_path_and_info = {}
            for root, dirs, files in os.walk(target_dir):
                if "name_detail.txt" in files:
                    all_pdf_full_path_and_info[os.path.join(root, "name_detail.txt").split(target_dir)[1].strip('/')] = "项目结果图片说明文档"
                    for file in os.listdir(root):
                        if file.endswith('.pdf'):
                            #print os.path.join(root, file)
                            all_pdf_full_path.append(os.path.join(root, file))


            for pdf in all_pdf_full_path:
                all_pdf_full_path_and_info[pdf.split(target_dir)[1].strip('/')] = self.filePath_to_fileInfo(pdf)



            add_relpath_rules_list_2 = []
            for pdf_path,pdf_info in all_pdf_full_path_and_info.items():
                pdf_info_dropTableid = '\t'.join(filter(None, [ele for ele in re.split('\s',pdf_info.strip()) if not ele.strip().startswith('table_id-') and not ele.strip().endswith('.pdf')]))
                add_relpath_rules_list_2.append([pdf_path.encode('utf-8'), 'pdf', pdf_info_dropTableid.replace('\t','____').replace(' ','____').strip('-').strip('____').encode('utf-8'), 0])
            
            #-------------------------- 生成 add_relpath_rules_list -------------------------------

            result_dir = self.add_upload_dir(os.path.join(self.tool.output_dir,'pdf'))
            self.inter_dirs = [
                ["00 Charts", "", "项目结果图片",0],
            ]
            add_relpath_rules_list = [
                ['.', '', '项目结果图片', 0],
                ['02_QC', '', '测序数据质控', 0],
                ['02_QC/Rawdata', '', '原始测序数据', 0],
                ['02_QC/Cleandata', '', '质控测序数据', 0],
                ['03_Align', '', '序列比对', 0],
                ['03_Align/assessmentsaturation', '', '测序饱和度分析', 0],
                ['03_Align/assessmentcoverage', '', '测序覆盖度分析', 0],
                ['03_Align/assessmentdistribution', '', '不同区域Reads分布', 0],
                ['03_Align/assessmentchrom', '', '不同染色体Reads分布', 0],
                ['04_Assemble', '', '转录本组装', 0],
                ['05_Annotation', '', '功能注释与查询', 0],
                ['06_Express', '', '表达量分析', 0],
                ['06_Express/expdistribution', '', '表达量分布', 0],
                ['06_Express/expvenn', '', '样本间Venn分析', 0],
                ['06_Express/expcorr', '', '样本间相关性分析', 0],
                ['06_Express/exppca', '', '样本间PCA分析', 0],
                ['06_Express/expbatch', '', '批次效应评估', 0],
                ['07_DiffExpress', '', '表达量差异分析', 0],
                ['08_GeneSet', '', '基因集分析', 0],
                ['08_GeneSet/01_GenesetVenn', '', '基因集Venn分析', 0],
                ['08_GeneSet/02_Cluster_Analysis', '', '基因集聚类分析', 0],
                ['08_GeneSet/03_COG_Annotation', '', '基因集COG注释', 0],
                ['08_GeneSet/04_GO_Annotation', '', '基因集GO注释', 0],
                ['08_GeneSet/05_KEGG_Annotation', '', '基因集KEGG注释', 0],
                ['08_GeneSet/06_GO_Enrich', '', '基因集GO富集分析', 0],
                ['08_GeneSet/07_KEGG_Enrich', '', '基因集KEGG富集分析', 0],
                ['08_GeneSet/08_Enrich_Circ', '', '基因集富集弦图', 0],
                ['08_GeneSet/09_GO_Dag', '', '基因集GO有向无环图', 0],
                ['08_GeneSet/10_Exp_Network', '', '基因集表达相关性分析', 0],
                ['09_AS', '', '可变剪切分析', 0],
                ['09_AS/AS_count', '', '可变剪切类型统计', 0],
                ['09_AS/AS_stat', '', '差异可变剪切类型统计', 0],
                ['10_SNP', '', 'SNP/InDel分析', 0],
                ['11_WGCNA', '', 'WGCNA分析', 0],
                ['11_WGCNA/wgcnaprepare', '', 'WGCNA预处理分析', 0],
                ['11_WGCNA/wgcnamodule', '', 'WGCNA模块识别分析', 0],
                ['11_WGCNA/wgcnarelate', '', 'WGCNA模块与表型分析', 0],
                ['12_Time_series_Analysis', '', '时序分析', 0],
                ['13_TF_Analysis', '', '转录因子分析', 0],
                ['14_PPI_Analysis', '', '蛋白互作网络分析', 0],
            ]
            # self.logger.info("xxxxxxxx____add_relpath_rules_list_2:\n"+str(add_relpath_rules_list_2))
            add_relpath_rules_list.extend(add_relpath_rules_list_2)
            result_dir.add_relpath_rules(add_relpath_rules_list)
            super(SangerChartWorkflow, self).end()

    def run_tool(self):
        options = dict(
            url=self.option('url'),
            userpass=self.option('userpass'),
            mode=self.option('mode'),
            imageversion=self.option('imageversion'),
        )
        self.tool.set_options(options)
        self.tool.run()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.ref_rna_v2.report.sanger_chart import SangerChartWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'sanger_chart_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'ref_rna_v2.report.sanger_chart',
            'options': {
                'url' : 'http://report.sanger.com/wholerna/expcorr/task_id/sanger_316281.html',
                'userpass': 'sgtest@majorbio.com____username_password____test01',
                'mode': 'single',
                'imageversion': 'latest',
                'chart_main_id': '000000000000000',
            }
        }
        wsheet = Sheet(data=data)
        wf =SangerChartWorkflow(wsheet)
        wf.sheet.id = 'sanger_chart'
        wf.sheet.project_sn = 'sanger_chart'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
