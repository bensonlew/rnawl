# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import regex
import re
import os
import xml.etree.ElementTree as ET


class gene2trans(object):
    def __init__(self):
        """
        根据转录本基因的对应关系提取信息
        """
        self.tran_gene = {}  # 转录本ID和对应的基因ID
        self.tran_list = []  # 所有的转录本ID
        self.gene_dict = {}

    def get_gene_transcript(self, gtf_path):
        """
        得到转录本ID和对应的基因ID
        """
        tran_list = []
        gene_list = []
        tran_gene = {}
        for line in open(gtf_path):
            txpt_id = ''
            gene_id = ''
            content_m = regex.match(
                r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                line.strip())
            if content_m:
                if 'transcript_id' in content_m.captures(6):
                    txpt_id = content_m.captures(7)[0]
                    gene_id = content_m.captures(10)[0]
                else:
                    txpt_id = content_m.captures(10)[0]
                    gene_id = content_m.captures(7)[0]
            if txpt_id not in tran_list:
                tran_list.append(txpt_id)
                gene_list.append(gene_id)
                tran_gene[txpt_id] = gene_id
        return tran_gene, tran_list, gene_list

    def get_gene_transcript_denovo(self, trans2gene):
        """
        得到转录本ID和对应的基因ID
        """
        tran_list = []
        gene_list = []
        tran_gene = {}
        for line in open(trans2gene):
            line = line.strip().split("\t")
            txpt_id = line[0]
            gene_id = line[1]

            if not tran_gene.has_key(txpt_id):
                tran_list.append(txpt_id)
                if line[2] == "yes":
                    gene_list.append(txpt_id)
                tran_gene[txpt_id] = gene_id
        return tran_gene, tran_list, gene_list

    def get_gene_blast_xml(self, tran_list, tran_gene, xml_path, gene_xml_path):
        """
        根据提供的基因和转录本对应关系，查找xml中的查询序列，将转录本ID替换成基因ID,生成新的xml
        """
        xml = ET.parse(xml_path)
        root = xml.getroot()
        BlastOutput_iterations = root.find('BlastOutput_iterations')
        for one_query in BlastOutput_iterations.findall('Iteration'):
            query_def = one_query.find('Iteration_query-def')
            query_def_split = re.split(r'\s', query_def.text, maxsplit=1)
            query_ID = query_def_split[0]
            if query_ID in tran_list:
                query_ID = re.sub(r'{}'.format(query_ID), tran_gene[query_ID], query_def_split[0])
                if len(query_def_split) == 2:
                    query_def.text = query_ID + ' ' + query_def_split[1]
                else:
                    query_def.text = query_ID
            else:
                BlastOutput_iterations.remove(one_query)
        xml.write('tmp.txt')
        with open('tmp.txt', 'rb') as f, open(gene_xml_path, 'wb') as w:
            lines = f.readlines()
            a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
            w.write(a)
            w.writelines(lines)
        # os.remove('tmp.txt')

    def select_xml_buffer(self, xml_iter_buffer):
        '''
        判断xml中的iter 是否为gene 并修改id
        '''
        iter_node = ET.fromstring(xml_iter_buffer)
        tran_id = iter_node.find('Iteration_query-def').text
        if tran_id in self.gene_dict:
            gene_id = self.tran_gene[tran_id]
            return xml_iter_buffer.replace(tran_id, gene_id)
        else:
            return ""


    def get_gene_blast_xml_from_trans(self, tran_list, tran_gene, xml_path, gene_xml_path):
        """
        根据提供的基因(代表转录本id列表)和转录本基因对应关系，查找xml中的查询序列，将转录本ID替换成基因ID,生成新的xml xml文件大时ET.parse效率低
        """
        self.tran_gene = tran_gene
        self.gene_dict = dict(zip(tran_list, xrange(len(tran_list))))
        iter_buffer = ''
        with open(xml_path, 'rb') as f, open(gene_xml_path, 'wb') as w:
            a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
            w.write(a)
            f.readline()
            f.readline()
            append = False
            for line in f:
                if line.startswith('<Iteration>'):
                    iter_buffer = line
                    append = True
                elif line.startswith('</Iteration>'):
                    iter_buffer += line
                    w.write(self.select_xml_buffer(iter_buffer))
                    append = False
                elif append == True:
                    iter_buffer += line
                else:
                    w.write(line)

        # os.remove('tmp.txt')

    def get_gene_blast_table(self, tran_list, tran_gene, table_path, gene_table_path):
        """
        根据提供的基因和转录本对应关系，查找table中的查询序列，将转录本ID替换成基因ID,生成新的table
        """
        with open(table_path, "rb") as f, open("tmp.xls", "wb") as w:
            lines = f.readlines()
            w.write(lines[0])
            for line in lines[1:]:
                item = line.strip().split("\t")
                query_id = item[5]
                if query_id in tran_list:
                    query = re.sub(r'{}'.format(query_id), tran_gene[query_id], line)
                    w.write(query)
        with open("tmp.xls", "rb") as f, open(gene_table_path, "wb") as w:
            lines = f.readlines()
            w.writelines(lines)

    def get_gene_go_list(self, tran_list, tran_gene, go_list, gene_go_list):
        """
        根据提供的基因和转录本对应关系，将go注释的go.list转录本ID替换成基因ID
        """
        with open(go_list, "rb") as f, open("tmp.list", "wb") as w:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                if item[0] in tran_list:
                    gene = re.sub(r"{}".format(item[0]), tran_gene[item[0]], line)
                    w.write(gene)
        with open("tmp.list", "rb") as f, open(gene_go_list, "wb") as w:
            lines = f.readlines()
            w.writelines(lines)
        # os.remove("tmp.list")


if __name__ == "__main__":
    a = gene2trans()
    a.get_gene_transcript()
    a.get_gene_blast_xml()
    a.get_gene_go_list()
    a.get_gene_blast_table()
