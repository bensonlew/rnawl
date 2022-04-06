# -*- coding: utf-8 -*-

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.itraq_and_tmt.go_graph import draw_GO
import glob
import unittest
import pandas as pd
import json

class GoDagAgent(Agent):
    def __init__(self, parent):
        super(GoDagAgent, self).__init__(parent)
        options = [
            {"name": "go_enrich_detail", "type": "string", "default": None},
            {"name": "go_list", "type": "string", "default": None},
            {"name": "top_num", "type": "int", "default": 20},
            {"name": "significant_diff", "type": "string", "default": None},
            {"name": "significant_value", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.step.add_steps("go_dag")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.go_dag.start()
        self.step.update()

    def stepfinish(self):
        self.step.go_dag.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option("go_enrich_detail")):
            raise OptionError("%s not exist", variables = (self.option("go_enrich_detail")), code = "33706001")
        if self.option("go_list") is None and self.option("significant_diff") is None:
            raise OptionError("go_list, significant_diff  one of them should not be None at least", code = "33706002")

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["go_dag.png", "png", "go有向无环图"]
        ])
        super(GoDagAgent,self).end()

class GoDagTool(Tool):
    def __init__(self,config):
        super(GoDagTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.go_graph = self.config.PACKAGE_DIR + "/ref_rna_v2/go_graph.py"
        # self.obo = self.config.SOFTWARE_DIR + '/database/GO/go-basic.obo'
        self.obo = self.config.SOFTWARE_DIR + '/database/Annotation/other2019/go.obo'

    def generate_go_dict(self):
        #self.logger.info(self.option("go_list"))
        #self.option(8888888888)
        go_p_dict = {}
        with open(self.option("go_enrich_detail")) as f:
            _ = f.readline().strip().split('\t')
            # pvidx = _.index('p_uncorrected')
            # paidx = _.index('p_corrected')
            if self.option("significant_value") is None:
                for line in f:
#                    if self.option("significant_diff").lower() == "pvalue":
#                        significant_diff = line.split("\t")[6]
#                    else:
#                        significant_diff = line.split("\t")[9]
                    for items in self.option("go_list").split(";"):
                        if line.startswith(items):
                            go_p_dict[line.split("\t")[0]] = float(line.split("\t")[9])
            elif self.option("go_list") is None:
                df = pd.read_table(self.option('go_enrich_detail'))
                if self.option('significant_diff').lower() == 'pvalue':
                    df = df.sort_values('p_uncorrected')
                    if "GO" in df.columns:
                        subdf = df.head(self.option('top_num')).reindex(['GO', 'p_uncorrected'], axis=1)
                    else:
                        subdf = df.head(self.option('top_num')).reindex(['go_id', 'p_uncorrected'], axis=1)
                else:
                    # padjust = df.columns[9]
                    padjust = "p_corrected"
                    df["r_id"] = df.index
                    df = df.sort_values(by=["p_corrected","r_id"],ascending=[True,True])
                    df = df.drop("r_id", 1)
                    if "GO" in df.columns:
                        subdf = df.head(self.option('top_num')).reindex(['go_id', padjust], axis=1)
                    else:
                        subdf = df.head(self.option('top_num')).reindex(['go_id', padjust], axis=1)
                go_p_dict = {row[0]: row[1] for i, row in subdf.iterrows()}
                # count = 0
                # for line in f:
                #     if self.option("significant_diff").lower() == "pvalue":
                #         significant_diff = line.split("\t")[pvidx]
                #     else:
                #         significant_diff = line.split("\t")[paidx]
                #     if count < self.option("top_num"):
                #         if float(significant_diff) < float(self.option("significant_value")):
                #             #self.logger.info(significant_diff)
                #             go_p_dict[line.split("\t")[0]] = float(significant_diff)
                #     else:
                #         break
                #     count = count + 1
        #self.logger.info("666666")
        #self.logger.info(go_p_dict.keys())
        return go_p_dict

    def run_dag(self):
        go_p_dict = self.generate_go_dict()
        self.logger.info("go_p_dict is {}".format(go_p_dict))
        node_dict, edge_dict = draw_GO(go_p_dict, out="go_dag", obo=self.obo, parse_altid = True)

        # 获取相关关系，动态图
        # 因为图形插件和go 数据库方向相反，所以导致child parent的关系

        enrich_df = pd.read_table(self.option("go_enrich_detail"), header=0)
        enrich_df.set_index("go_id", inplace=True)


        with open("relation.tsv", 'w') as fo:
            fo.write("child\tparent\trelation\tcategory\tdetail\tp_value\tnumber_ratio\tchoosed\n")
            par_list = [edge_dict[node][0] for node in edge_dict]
            edge_dict2 = dict()
            for node1, value in edge_dict.items():
                edge_dict2[value[0]] = [node1, value[1], ]
            print "par_list is {}".format(par_list)
            node_set2 = set()
            node_set1 = set()
            dag_list = list()
            print "edge_dict is {}".format(edge_dict)

            for node, node_values in edge_dict.items():
                [node1, node2] = node.split("|")

                node_set2.add(node2)
                node_set1.add(node1)

                eles = [
                    node2,
                    node1,
                    node_values,
                    node_dict[node2][2],
                    node_dict[node2][1]]
                dag_list.append(eles)

            for node in node_set1 - node_set2:
                eles = [node,
                        "",
                        "is_a",
                        node_dict[node][2],
                        node_dict[node][1]]
                dag_list.append(eles)

            for eles in dag_list:
                try:
                    dic = dict(enrich_df.loc[eles[0]])
                    eles.append(str(dic.get("p_uncorrected", "")))
                    eles.append(dic.get("ratio_in_study", ""))
                except:
                    eles.append("")
                    eles.append("")

                if eles[0] in go_p_dict:
                    eles.append("yes")
                else:
                    eles.append("no")
                fo.write("\t".join(eles) + "\n")
        relations = set([e for e in edge_dict.values()])
        relationship_dict = {'is_a': ('#130c0e', 'solid'),
                             'part_of': ('#2a5caa', 'solid'),
                             'negatively_regulates': ('#d71345', 'solid'),
                             'positively_regulates': ('#1d953f', 'solid'),
                             'regulates': ('#ffc20e', 'solid'),
                             'occurs_in': ('#008792', 'solid'),
                             'capable_of': ('#33a3dc', 'dashed'),
                             'capable_of_part_of': ('#f36c21', 'dashed')}
        colors = [relationship_dict[r][0] for r in relations]
        visualMap = [
            {
                "dataType": "categories",
                "visualType": "color",
                "visualValue": colors,
                "pieces": list(relations),
                "data": "relation",
                "legendIndex": 0,
            },
            {
                "dataType": "categories",
                "visualType": "color",
                "visualValue": ["none", "#F08080"],
                "pieces": ["no", "yes"],
                "data": "choosed",
                "legendIndex": 0,
            }
        ]
        jsonout = open("relation_visual.json",'w')
        jsonout.write(json.dumps(visualMap, indent=4))




    def set_output(self):
        go_dag = glob.glob(self.work_dir + '/go_dag.*')
        for each in go_dag:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(GoDagTool, self).run()
        self.run_dag()
        self.set_output()
        self.end()
