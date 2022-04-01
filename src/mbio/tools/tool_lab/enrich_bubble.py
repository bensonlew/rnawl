# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import os
import re
import xlrd
import shutil


class EnrichBubbleAgent(Agent):
    """
    Go富集气泡图（单组）
    """

    def __init__(self, parent):
        super(EnrichBubbleAgent, self).__init__(parent)
        options = [
            {"name": "table_file1", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file2", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file3", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file4", "type": "infile", "format": "small_rna.common"},
            {"name": "table_file5", "type": "infile", "format": "small_rna.common"},
            {"name": "group_num", "type":"int"},
            {"name": "p_value", "type": "float"},
            {"name": "result_num", "type": "int"},
            {"name": "output_form", "type":"string"},
            {"name": "ids","type":"string"},
            {"name": "show_order","type":"string"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("group_num"):
            if not self.option("table_file1"):
                raise OptionError("必须设置作图数据")
            if not self.option("p_value"):
                raise OptionError("必须设置P_value值")
            if not self.option("result_num"):
                raise OptionError("result_num")
        else:
            if not self.option("group_num"):
                raise OptionError("必须提供表格数")
            for i in range(1,self.option("group_num")):
                if not self.option("table_file{}".format(i)):
                    raise OptionError("缺少输入表格{}".format(i))
            if self.option("output_form") not in ["union","intersection","ids"]:
                raise OptionError("需要提供输出形式")
            if self.option("output_form") == "interation":
                if not self.option("show_order"):
                    raise OptionError("当结果以交集展示时，需要指定显示方式")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "筛选气泡图数据"],
            ["./bubble.csv", "csv", "气泡图数据"]
        ])
        super(EnrichBubbleAgent, self).end()


class EnrichBubbleTool(Tool):
    """
    富集气泡图
    """

    def __init__(self, config):
        super(EnrichBubbleTool, self).__init__(config)
        self._version = 'v1.0'
        self.out_table = {}
        self.p_value = self.option("p_value")
        # self.select_num = self.option("result_num")


    def run(self):
        """
        运行小工具
        """ 
        super(EnrichBubbleTool, self).run()
        if not self.option("group_num"):
            table_file = self.option("table_file1").prop['path']
            term_p, term_all = self.run_get_data(table_file)
            sort_pvalue = sorted(zip(term_p.itervalues(),term_p.iterkeys()),key=lambda x:(x[0]["pvalue"],x[0]["id_name"]))
            sort_term= []
            num = 0
            for i in sort_pvalue:
                if num < self.option("result_num"):
                    sort_term.append(i[1])
                    num += 1
            self.run_set_table(term_all,"table_file1",sort_term)
        else:
            output_form = self.option("output_form")
            if output_form == "intersection":
                self.group_intersection()
            elif output_form == "union":
                self.group_union()
            else:
                self.group_target(self.option("ids"))

        self.set_output()
        self.end()

    def group_intersection(self):
        term_p = {}
        table_info = {}
        # self.table_info = {}
        sort_term = []
        a = []
        for i in range(1,self.option("group_num") + 1):
            table_name = "table_file{}".format(i)
            table_file = self.option(table_name).prop['path']
            term_p[table_name],table_info[table_name] = self.run_get_data(table_file)
            if i == 1 :
                a = term_p[table_name].keys()
            else:
                a = set(a).intersection(term_p[table_name].keys())
        temp_inter = {}
        for key in a:
            temp_inter[key] = term_p[self.option("show_order")][key]
        sort_pvalue = sorted(zip(temp_inter.itervalues(),temp_inter.iterkeys()),key=lambda x:(x[0]["pvalue"],x[0]["id_name"]))# sort_term = []
        n = 0
        for i in sort_pvalue:
            if n >= self.option("result_num"):
                break
            sort_term.append(i[1])
            n += 1
            if i[1] not in a:
                continue
            # for i in range(1,self.option("group_num") + 1):
            #     table_name = "table_file{}".format(i)
            #     self.table_info[table_name][i[1]] = table_info[table_name][i[1]] if table_info[table_name].has_key(i[i]) else 0
        for i in range(1,self.option("group_num") + 1):
            table_name = "table_file{}".format(i)
            self.run_set_table(table_info[table_name],table_name,sort_term)

    def group_union(self):
        term_p = {}
        table_info = {}
        table_all_info = {}
        self.table_info = {}
        a = []
        term_id = {}
        for i in range(1,self.option("group_num") + 1):
            table_name = "table_file{}".format(i)
            self.table_info[table_name] = {}
            table_file = self.option(table_name).prop['path']
            term_p[table_name],table_all_info[table_name] = self.run_get_data(table_file)
            for i in table_all_info[table_name].keys():
                term_id[i] = table_all_info[table_name][i]["id_name"]
            # _, = self.run_get_data(table_file,num=-1)
            if i == 1 :
                a =  term_p[table_name].keys()
            else:
                a = list(set(a).union(term_p[table_name].keys()))
        term_p1 = {}
        for term in a:
            for i in range(1,self.option("group_num") + 1):
                table_name = "table_file{}".format(i)
                table_null = {"term_description": term,"number":"0","pvalue":float(0),"enrich_factor":"0","id_name":term_id[term]}
                self.table_info[table_name][term] = table_all_info[table_name][term] if table_all_info[table_name].has_key(term) else table_null
            term_p1[term] = self.table_info["table_file1"][term]
        sort_pvalue = sorted(zip(term_p1.itervalues(),term_p1.iterkeys()),key=lambda x:(x[0]["pvalue"],x[0]["id_name"]))
        sort_term = []
        for i in sort_pvalue:
            sort_term.append(i[1])
            if i[1] not in a:
                continue
        for i in range(1,self.option("group_num") + 1):
            table_name = "table_file{}".format(i)
            self.run_set_table(self.table_info[table_name],table_name,sort_term)

    def group_target(self,target):
        table_info = {}
        table_all_info = {}
        ids = target.split(";")
        term_p = {}
        for i in range(1,self.option("group_num") + 1):
            table_name = "table_file{}".format(i)
            table_info[table_name] = {}
            table_file = self.option(table_name).prop['path']
            _,table_all_info[table_name] = self.run_get_data(table_file)
        for term_id in ids:
            term = self.id_term[term_id]
            for i in range(1,self.option("group_num") + 1):
                table_name = "table_file{}".format(i)
                table_null = {"term_description": term,"number":"0","pvalue":float(0),"enrich_factor":"0"}
                table_info[table_name][term] = table_all_info[table_name][term] if table_all_info[table_name].has_key(term) else table_null
            term_p[term] = table_info["table_file1"][term]["pvalue"] 
        sort_pvalue = sorted(zip(term_p.itervalues(),term_p.iterkeys()),key=lambda x:(x[0]["pvalue"],x[0]["id_name"]))
        sort_term = []
        for i in sort_pvalue:
            sort_term.append(i[1])
        for i in range(1,self.option("group_num") + 1):
            table_name = "table_file{}".format(i)
            self.run_set_table(table_info[table_name],table_name,sort_term)

    def run_get_data(self,table_file):
        """
        读取表格，重组分析数据
        """
        self.id_term = {}
        term_p = {}
        term_all = {}
        title_name = {}
        txt_file = os.path.join(self.work_dir,"{}.txt".format(table_file))
        self.xls2txt(table_file,txt_file)
        self.p_term = ""
        self.description_term = ""
        p_name = ""
        term_name = ""
        id_name = ""
        with open(txt_file, "r") as tbf:
            line = tbf.readline()
            titles = line.rstrip().split("\t")
            for title in titles:
                title_name[title.lower()] = titles.index(title)
                if re.match("pvalue", title.lower()) or re.match("padjust", title.lower()):
                    p_name = title.lower()
                    self.p_term = title.lower()
                if re.match("term description", title.lower()) or re.match("description", title.lower()) or re.match("pathway description", title.lower()):
                    term_name = title.lower()
                    self.description_term = title.lower()
                if re.match("go id", title.lower()) or re.match("pathway id", title.lower()):
                    id_name = title.lower()
                    # self.id = title.lower()
            while 1:
                term_info = {}
                line = tbf.readline()
                if not line:
                    break
                # if num == 0:
                #     break
                terms = line.rstrip().split("\t")
                term_description = terms[title_name[term_name]]
                number = terms[title_name["number"]] if title_name.has_key("number") else terms[title_name["num"]]
                ratio_in_study = float(
                    terms[title_name["ratio_in_study"]].split('/')[0])
                ratio_in_pop = float(
                    terms[title_name["ratio_in_pop"]].split('/')[0])
                pvalue = float(terms[title_name[p_name]])
                enrich_factor = ratio_in_study/ratio_in_pop
                term_info["term_description"] = term_description
                term_info["number"] = number
                term_info["pvalue"] = pvalue
                term_info["enrich_factor"] = enrich_factor
                term_info["id_name"] = terms[title_name[id_name]]
                self.id_term[terms[title_name[id_name]]] = term_description
                if pvalue < self.option("p_value"):
                    # num = num - 1
                    # term_p[term_description] = pvalue
                    term_p[term_description] = term_info
                    term_all[term_description] = term_info
        if self.option("group_num") and self.option("output_form") == "union":
            sort_pvalue = sorted(zip(term_p.itervalues(),term_p.iterkeys()),key=lambda x:(x[0]["pvalue"],x[0]["id_name"]))
            sort_term = []
            for i in sort_pvalue:
                sort_term.append(i[1])
            num = 1
            temp_p={}
            for i in sort_term:
                if num > self.option("result_num"):
                    break
                temp_p[i]=term_p[i]
                num += 1
            term_p=temp_p
        return term_p, term_all


    def xls2txt(self,xls_file, txt_file):
        try:
            file = xlrd.open_workbook(xls_file.decode('utf-8'))
            sheet = file.sheet_names()[0]
            df = pd.read_excel(xls_file.decode('utf-8'), sheet_name=sheet, header=None)		# 使用pandas模块读取数据
            # print('开始写入txt文件...')
            df.to_csv(txt_file, header=None, sep='\t', index=False)		# 写入，逗号分隔
            # print('文件写入成功!')
        except:
            shutil.copyfile(xls_file,txt_file)

    def run_set_table(self,table_info,table_file,sort_term):
        """
        重新生成气泡图所需要数据的表格
        """
        num = 0
        with open("{}_bubble.xls".format(table_file), "w") as bb:
            bb.write("X\tY\tsize\tcolor\tnum\n")
            text = ""
            for term in sort_term:
                if int(table_info[term]["number"]) == 0:
                    continue
                else:
                    num += 1
                    info = table_info[term]
                    tmp_text = "{}\t{}\t{}\t{}\t{}\n".format(
                        info["enrich_factor"], info["term_description"], info["number"], info["pvalue"], num)
                    text = text + tmp_text
            bb.write(text)

    def set_output(self):
        """
        设置输出文件
        """
        with open(os.path.join(self.output_dir,"term.txt"),"w") as tt:
            tt.write("{}\t{}".format(self.p_term,self.description_term))
        if not self.option("group_num"):
            input_name = os.path.basename(self.option("table_file1").prop["path"]).split(".")[0]
            path1 = self.output_dir + "/{}_bubble.xls".format(input_name)
            if os.path.exists(path1):
                os.remove(path1)
            os.link(self.work_dir + "/table_file1_bubble.xls", self.output_dir + "/{}_bubble.xls".format(input_name))
        else:
            for i in range(1,self.option("group_num") + 1):
                table_name = "table_file{}".format(i)
                input_name = os.path.basename(self.option(table_name).prop["path"]).split(".")[0]
                path1 = self.output_dir + "/{}_bubble.xls".format(input_name)
                if os.path.exists(path1):
                    os.remove(path1)
                os.link(self.work_dir + "/{}_bubble.xls".format(table_name), self.output_dir + "/{}_bubble.xls".format(input_name))
        
class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        # data = {
        #     "id": "gobubble_" + str(random.randint(1, 10000)),
        #     "type": "tool",
        #     "name": "tool_lab.go_enrich_bubble_mult",
        #     "options": {
        #         "table_file1": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue.xls",
        #         "table_file2": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue1.xls",
        #         "group_num":2,
        #         "output_form":"intersection",
        #         "p_value" : 0.05, 
        #         "result_num":20,
        #         "show_order":"table_file1"
        #         }
        #     }
        data = {
            "id": "enrichbubble_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.enrich_bubble",
            "options": {
                "table_file1": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/KEGG-1.xls",
                "table_file2": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/KEGG-2.xls",
                "table_file3": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/KEGG-3.xls",
                "group_num":3,
                "output_form":"intersection",
                "p_value" : 0.05, 
                "result_num":20,
                "show_order":"table_file1"
                }
        }
        # data = {
        #     "id": "gobubble_" + str(random.randint(1, 10000)),
        #     "type": "tool",
        #     "name": "tool_lab.go_enrich_bubble_mult",
        #     "options": {
        #         "table_file1": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue.xls",
        #         "table_file2": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/enrich_bubble/test_data/go_pvalue1.xls",
        #         "group_num":2,
        #         "output_form":"terms",
        #         "p_value" : 0.05, 
        #         "result_num":20,
        #         "show_order":"table_file1",
        #         "terms":"negative regulation of carbohydrate metabolic process1;cellular hormone metabolic process;response to testosterone1"
        #         }
        #     }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

   

