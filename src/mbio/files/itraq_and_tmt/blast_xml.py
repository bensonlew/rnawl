# -*- coding: utf-8 -*-
# __author__ = '刘彬旭'
from biocluster.iofile import File
import os
from biocluster.core.exceptions import FileError
from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET
import re


class BlastXmlFile(File):
    """
    定义blast+比对输出类型为5结果文件的xml格式， 测试blast+为2.3.0版本
    """

    def __init__(self):
        super(BlastXmlFile, self).__init__()

    def get_info(self):
        """
        获取文件属性
        """
        super(BlastXmlFile, self).get_info()
        blast_info = self.get_xml_info()
        self.set_property('query_num', blast_info[0])
        self.set_property('hit_num', blast_info[1])
        self.set_property('hsp_num', blast_info[2])
        self.set_property('hit_query_list', blast_info[3])

    def get_xml_info(self):
        """
        获取blast结果xml的信息

        :return
        """
        blastxml = NCBIXML.parse(open(self.path))
        try:
            query_count = 0
            align_count = 0
            hsp_count = 0
            query_list = []
            for query in blastxml:
                query_count += 1
                if query.alignments:
                    query_list.append(re.split(' ', query.query, maxsplit=1)[0])
                for align in query.alignments:
                    align_count += 1
                    for hsp in align.hsps:
                        hsp_count += 1
        except ValueError:
            raise FileError('传入文件类型不正确，无法解析，请检查文件是否正确，或者生成文件的blast版本不正确，本程序测试版本blast+2.3.0', code="42500301")
        except Exception as e:
            raise FileError('未知原因导致blastxml文件解析失败:%s', variables=(e), code="42500302")
        return query_count, align_count, hsp_count, query_list

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        if super(BlastXmlFile, self).check():
            # 父类check方法检查文件路径是否设置，文件是否存在，文件是否为空
            blastxml = NCBIXML.parse(open(self.path))
            try:
                blastxml.next()
            except ValueError:
                raise FileError('传入文件类型不正确，无法解析，请检查文件是否正确，或者生成文件的blast版本不正确，本程序测试版本blast+2.3.0', code="42500303")
            except Exception as e:
                raise FileError('未知原因导致blastxml文件解析失败:%s', variables=(e), code="42500304")
            return True

    def convert2table(self, outfile):
        """调用packages中的xml2table方法来转换到table格式，结果文件有表头"""
        from mbio.packages.align.blast.xml2table import xml2table
        xml2table(self.path, outfile)

    def convert2blast6default(self, outfile):
        """调用packages中的xml2blast6方法将xml转换成blast默认的outfmt 6输出格式，结果没有表头"""
        from mbio.packages.align.blast.xml2table import xml2blast6
        xml2blast6(self.path, outfile)

    def sub_blast_xml(self, genes, new_fp, trinity_mode=False):
        """
        根据提供的基因列表，查找xml中的查询序列，生成并集新的xml
        trinity_mode用于在新生成的xml的queryID是去除结尾的_i(数字) 的
        """
        genes = dict(zip(genes, xrange(len(genes))))
        xml = ET.parse(self.path)
        root = xml.getroot()
        BlastOutput_iterations = root.find('BlastOutput_iterations')
        for one_query in BlastOutput_iterations.findall('Iteration'):
            query_def = one_query.find('Iteration_query-def')
            query_def_split = re.split(r'\s', query_def.text, maxsplit=1)
            query_ID = query_def_split[0]
            if query_ID in genes:
                if trinity_mode:
                    query_ID = re.sub(r'_i[0-9]+$', '', query_def_split[0])
                    if len(query_def_split) == 2:
                        query_def.text = query_ID + ' ' + query_def_split[1]
                    else:
                        query_def.text = query_ID
            else:
                BlastOutput_iterations.remove(one_query)
        xml.write('tmp.txt')
        with open('tmp.txt', 'rb') as f, open(new_fp, 'wb') as w:
            lines = f.readlines()
            a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
            w.write(a)
            w.writelines(lines)
        os.remove('tmp.txt')

    def filter_blast_xml(self, filter_xml, evalue=0.001, identity=0, similarity=0):
        """
        根据evalue, identity, similarty三个阈值对xml文件进行过滤
        """
        xml = ET.parse(self.path)
        root = xml.getroot()
        BlastOutput_iterations = root.find('BlastOutput_iterations')
        for one_query in BlastOutput_iterations.findall('Iteration'):
            for one_hit in one_query.find('Iteration_hits').findall('Hit'):
                hsp_evalue = float(one_hit.find('Hit_hsps').find('Hsp').find('Hsp_evalue').text)
                hsp_identity = float(one_hit.find('Hit_hsps').find('Hsp').find('Hsp_identity').text)
                hsp_similarity = float(one_hit.find('Hit_hsps').find('Hsp').find('Hsp_positive').text)
                hsp_len = float(one_hit.find('Hit_hsps').find('Hsp').find('Hsp_align-len').text)
                hsp_identity_rate = hsp_identity/hsp_len * 100
                hsp_similarity_rate = hsp_similarity/hsp_len * 100

                if hsp_evalue > evalue or hsp_identity_rate < identity or hsp_similarity_rate < similarity:
                    one_query.find('Iteration_hits').remove(one_hit)
                else:
                    pass
        xml.write('tmp.txt')
        with open('tmp.txt', 'rb') as f, open(filter_xml, 'wb') as w:
            lines = f.readlines()
            a = '<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n'
            w.write(a)
            w.writelines(lines)
        os.remove('tmp.txt')


    def change_blast_version(self, fp, version='2.2.25+'):
        """
        特殊为blast2go准备，程序特定不能做当前2.3.0+版本的blast+
        """
        with open(self.path) as f, open(fp, 'w') as w:
            for line in f:
                l_str = line.strip()
                if l_str.startswith('<BlastOutput_version>'):
                    if '2.3.0+' not in line:
                        self.set_error('blast 程序版本不是2.3.0版本，此处报错不是必须用2.3.0版本，而是b2gPipe程序需要版本为2.2.25，此处特殊改为2.2.25，依然可以进行blast2go，\
                            但是不代表后续版本仍然可以使用，既然blast版本修改，请重新检查')
                    line = line.replace('2.3.0+', version)
                    w.write(line)
                    break
                else:
                    w.write(line)
            for line in f:
                w.write(line)
        return fp


if __name__ == '__main__':  # for test
    a = BlastXmlFile()
    # a.set_path('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\annotation\\annotation\\NR\\transcript.fa_vs_nr.blasttable.xls')
    a.set_path("/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/annotation/output/blast_xml/Trinity_vs_animal.xml")
    a.check()
    a.get_info()
    # a.convert2table('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\test.xls')
    a.filter_blast_xml("new_filter.xml", 0.00001, 90, 90)
