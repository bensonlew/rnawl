# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180508

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os


class SummarySetAgent(Agent):
    """
     用于对组装后的结果文件Lands.anno.summary，与blast结果文件进行处理，获取拼接后的每个基因在ref上的chr end start
     注意在blast结果中，每个基因取哪条记录，首先是判断长度最长的，如果长度一样，去identity值最大，最后去evalue最小的
    """
    def __init__(self, parent):
        super(SummarySetAgent, self).__init__(parent)
        options = [
            {"name": "anno_summary", "type": "string"},  # 每个样本组装后序列的注释文件
            {"name": "blast_result", "type": "string"},   # 每个样本组装后序列的blast结果文件
            {"name": "sample_id", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('gatk')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gatk.start()
        self.step.update()

    def step_end(self):
        self.step.gatk.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("anno_summary"):
            raise OptionError("缺少anno_summary参数", code="34506701")
        if not self.option("blast_result"):
            raise OptionError("缺少blast_result参数", code="34506702")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '50G'
        
    def end(self):
        super(SummarySetAgent, self).end()


class SummarySetTool(Tool):
    def __init__(self, config):
        super(SummarySetTool, self).__init__(config)
        
    def get_blast_info(self):
        """
        获取每个gene对应的信息，格式如下：[{C7: {"chr": "chr1", "start": "1111", "end": "2222", "length": "30", "identity": "100"}}]
        :return:
        """
        # genes = []
        # blast_info = []
        # best_genes = []
        item = {}
        with open(self.option("blast_result"), 'r') as r:
            # data = r.readlines()
            for line in r:
                temp = line.strip().split("\t")
                # if temp[0] not in genes:
                #     genes.append(temp[0])
                if temp[0] not in item.keys():
                    item[temp[0]] = {'chr': temp[1], "identity": float(temp[2]), "length": int(temp[3]),
                                     "start": temp[8], "end": temp[9], "score": float(temp[11])}
                else:
                    if int(temp[3]) > item[temp[0]]['length']:
                        item[temp[0]] = {'chr': temp[1], "identity": float(temp[2]), "length": int(temp[3]),
                                         "start": temp[8], "end": temp[9], "score": float(temp[11])}
                    else:
                        if float(temp[2]) > item[temp[0]]['identity'] or float(temp[11]) > item[temp[0]]['score']:
                            item[temp[0]] = {'chr': temp[1], "identity": float(temp[2]), "length": int(temp[3]),
                                             "start": temp[8], "end": temp[9], "score": float(temp[11])}
        # self.logger.info("item:{}".format(item))
        return item

        #         blast_info.append(item)
        # # 获取最佳序列
        # for gene in genes:
        #     length = 0
        #     identity = 0
        #     score = 0
        #     best_item = {}
        #     for info in blast_info:
        #         if info.keys()[0] == gene:
        #             if info.values()[0]["length"] > length:
        #                 length = info.values()[0]["length"]
        #                 identity = info.values()[0]["identity"]
        #                 score = info.values()[0]["score"]
        #                 best_item['chr'] = info.values()[0]['chr']
        #                 best_item['start'] = info.values()[0]['start']
        #                 best_item['end'] = info.values()[0]['end']
        #             elif info.values()[0]["length"] == length:
        #                 if info.values()[0]["identity"] > identity:
        #                     length = info.values()[0]["length"]
        #                     identity = info.values()[0]["identity"]
        #                     score = info.values()[0]["score"]
        #                     best_item['chr'] = info.values()[0]['chr']
        #                     best_item['start'] = info.values()[0]['start']
        #                     best_item['end'] = info.values()[0]['end']
        #                 elif info.values()[0]["identity"] == identity and info.values()[0]["score"] > score:
        #                     length = info.values()[0]["length"]
        #                     identity = info.values()[0]["identity"]
        #                     score = info.values()[0]["score"]
        #                     best_item['chr'] = info.values()[0]['chr']
        #                     best_item['start'] = info.values()[0]['start']
        #                     best_item['end'] = info.values()[0]['end']
        #     result = {gene: best_item}
        #     best_genes.append(result)
        # return best_genes, genes

    def set_summary_start_end(self):
        if self.option("sample_id"):
            file_path = self.output_dir + "/{}.final.summary".format(self.option("sample_id"))
        else:
            file_path = self.output_dir + "/{}.final.summary"\
                .format(os.path.basename(self.option("anno_summary")).split('.')[0])
        # best_genes, genes = self.get_blast_info()
        best_genes = self.get_blast_info()
        with open(self.option("anno_summary"), 'r') as r, open(file_path, 'w') as w:
            # data = r.readlines()
            for line in r:
                if re.match(r'^#', line) or re.match(r'^##', line):
                    w.write(line.rstrip('\n') + "\tChr\tStart\tEnd\n")
                else:
                    temp = line.strip().split("\t")
                    if temp[0] in best_genes.keys():
                        w.write(line.rstrip('\n') + "\t{}\t{}\t{}\n".format(best_genes[temp[0]]['chr'],
                                                                            best_genes[temp[0]]['start'],
                                                                            best_genes[temp[0]]['end']))
                        # for m in best_genes:
                        #     if m.keys()[0] == temp[0]:
                        #         w.write(line.rstrip('\n') + "\t{}\t{}\t{}\n".format(m.values()[0]['chr'],
                        #                                                             m.values()[0]['start'],
                        #                                                             m.values()[0]['end']))
                    else:
                        w.write(line.rstrip('\n') + "\t--\t--\t--\n")

    def run(self):
        super(SummarySetTool, self).run()
        self.logger.info("开始进行summary文件处理！")
        self.set_summary_start_end()
        self.logger.info("进行summary文件处理完成！")
        self.end()
