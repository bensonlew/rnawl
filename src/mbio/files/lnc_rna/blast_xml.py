# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import File
import os
from biocluster.core.exceptions import FileError
from Bio.Blast import NCBIXML
from mbio.files.ref_rna_v2.common import CommonFile
import xml.etree.ElementTree as ET
import re
from biocluster.config import Config
import sqlite3

class BlastXmlFile(File):
    '''
    Define XML format for blast+ alignment output of type 5
    Test version: blast+ 2.3.0
    '''

    def __init__(self):
        self.db_path = os.path.join(Config().SOFTWARE_DIR, "database/Annotation/latest")
        super(BlastXmlFile, self).__init__()

    def get_info(self):
        '''
        Get file properties
        '''
        super(BlastXmlFile, self).get_info()
        blast_info = self.get_xml_info()
        self.set_property('query_num', blast_info[0])
        self.set_property('hit_num', blast_info[1])
        self.set_property('hsp_num', blast_info[2])
        self.set_property('hit_query_list', blast_info[3])

    def convert_xml2go(self, new_xml):
        '''
        Modify the annotation format for blast2go input, blast attribute information was labeled by gene ID
        '''
        nracc2gi = self.db_path + "/nr_acc2gi.db"
        conn = sqlite3.connect(nracc2gi)
        cursor = conn.cursor()

        xml = ET.parse(self.path)
        root = xml.getroot()
        blast_version=root.find('BlastOutput_version')
        blast_version.text = 'BLAST 2.2.25+'
        BlastOutput_iterations = root.find('BlastOutput_iterations')
        BlastOutput_iterations_changed = BlastOutput_iterations.copy()
        BlastOutput_iterations_changed.clear()
        iter_id = 1
        for one_query in BlastOutput_iterations.findall('Iteration'):
            all_hit = one_query.find('Iteration_hits').findall('Hit')
            for one_hit in all_hit:
                # 由每个基因比对一次，改为每个基因和结果存在一次记录，为了标记结果中的数据库序列来源
                # one_hit.find('Hit_id').text
                one_query_new = one_query.copy()
                one_query_new.clear()

                iter_new = one_query.find('Iteration_iter-num').copy()
                iter_new.text = str(iter_id)
                one_query_new.append(iter_new)

                query_id_old = one_query.find('Iteration_query-ID')
                query_id_new = query_id_old.copy()
                query_id_new.text = "Query_" + str(iter_id)
                one_query_new.append(query_id_new)

                query_def_old = one_query.find('Iteration_query-def')
                query_def_new = query_def_old.copy()
                acc_id = one_hit.find('Hit_id')

                acc_id_new = acc_id
                # 根据数据库获得acc对应的gi编号
                if len(acc_id.text.split("|")) == 1:
                    acc = acc_id.text
                    gi = 0
                    try:
                        cursor.execute('select * from acc2gi where acc="{}"'.format(acc))
                        gi = cursor.fetchall()[0][1]
                    except:
                        gi = 0
                    db = "gb"
                    if acc.startswith("XP_") or acc.startswith("WP_"):
                        db = "ref"
                    acc_id_new  =  "gi|{}|{}|{}|".format(str(gi), db, acc)
                '''
                if one_hit.find('Hit_id').text.split("|")[-1]:
                    acc_id = one_hit.find('Hit_id').text.split("|")[-1]
                else:
                    acc_id = one_hit.find('Hit_id').text.split("|")[-2]
                '''
                one_hit.find('Hit_num').text = "1"
                one_hit.find('Hit_id').text = acc_id_new

                query_def_new.text = query_def_old.text.split()[0] + '__id__' + str(acc_id.text)
                one_query_new.append(query_def_new)

                one_query_new.append(one_query.find('Iteration_query-len'))
                hit_iter_new=one_query.find("Iteration_hits").copy()
                hit_iter_new.clear()
                hit_iter_new.append(one_hit)
                hit_iter_new.text = "\n"
                hit_iter_new.tail = "\n"
                one_query_new.append(hit_iter_new)

                one_query_new.append(one_query.find('Iteration_stat'))
                # 只保留一次比对结果
                iter_id = iter_id + 1
                one_query_new.text = "\n"
                one_query_new.tail = "\n"
                BlastOutput_iterations_changed.append(one_query_new)
                # print BlastOutput_iterations_changed.getchildren()

        root.remove(BlastOutput_iterations)
        BlastOutput_iterations_changed.text = "\n"
        BlastOutput_iterations_changed.tail = "\n"
        root.append(BlastOutput_iterations_changed)

        xml.write(new_xml, encoding='utf-8', xml_declaration=True)


    def get_xml_info(self):
        '''
        Get information about the blast XML result
        '''
        blastxml = NCBIXML.parse(open(self.path))
        try:
            query_count = 0
            align_count = 0
            hsp_count = 0
            query_list = []
            if os.path.getsize(self.path) > 0:
                for query in blastxml:
                    query_count += 1
                    if query.alignments:
                        query_list.append(re.split(' ', query.query, maxsplit=1)[0])
                    for align in query.alignments:
                        align_count += 1
                        for hsp in align.hsps:
                            hsp_count += 1
        except ValueError:
            raise FileError('传入文件类型不正确，无法解析，请检查文件是否正确，或者生成文件的blast版本不正确，本程序测试版本blast+2.3.0', code = "43700801")
        except Exception as e:
            raise FileError('未知原因导致blastxml文件解析失败：%s', variables = (e), code = "43700802")
        return query_count, align_count, hsp_count, query_list

    def check(self):
        '''
        Check whether the documents meet the requirements
        '''
        if super(BlastXmlFile, self).check():
            # The parent class check method checks if the file path is set
            blastxml = NCBIXML.parse(open(self.path))
            if os.path.getsize(self.path) > 0:
                try:
                    blastxml.next()
                except ValueError:
                    raise FileError('传入文件类型不正确，无法解析，请检查文件是否正确，或者生成文件的blast版本不正确，本程序测试版本blast+2.3.0', code = "43700803")
                except Exception as e:
                    raise FileError('未知原因导致blastxml文件解析失败：%s', variables = (e), code = "43700804")
            return True

    def convert2table(self, outfile):
        '''
        Call the xml2table method in packages to convert to the table format, and the resulting file has a table header
        '''
        from mbio.packages.align.blast.xml2table import xml2table
        xml2table(self.path, outfile)

    def convert2blast6default(self, outfile):
        """
        Call the xml2blast6 method in packages to transform the XML into blast outfmt 6 with no header
        """
        from mbio.packages.align.blast.xml2table import xml2blast6
        xml2blast6(self.path, outfile)

    def sub_blast_xml(self, genes, new_fp, trinity_mode=False):
        '''
        According to the gene list provided, search the query sequence in XML and generate new XML
        '''
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
        Filter XML files according to evalue, identity and similarty
        """
        xml = ET.parse(self.path)
        root = xml.getroot()
        BlastOutput_iterations = root.find('BlastOutput_iterations')
        for one_query in BlastOutput_iterations.findall('Iteration'):
            for one_hit in one_query.find('Iteration_hits').findall('Hit'):
                hsp_evalue = float(one_hit.find('Hit_hsps').find('Hsp').find('Hsp_evalue').text) #
                hsp_identity = float(one_hit.find('Hit_hsps').find('Hsp').find('Hsp_identity').text) #
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
        '''
        Blast + is specific for blast2go, and specific programs cannot blast+ the current 2.3.0+ version
        '''
        with open(self.path) as f, open(fp, 'w') as w:
            for line in f:
                l_str = line.strip()
                if l_str.startswith('<BlastOutput_version>'):
                    if '2.3.0+' not in line:
                        raise FileError('blast 程序版本不是2.3.0版本，此处报错不是必须用2.3.0版本，而是b2gPipe程序需要版本为2.2.25，此处特殊改为2.2.25，依然可以进行blast2go，但是不代表后续版本仍然可以使用，既然blast版本修改，请重新检查', code = "43700805")
                    line = line.replace('2.3.0+', version)
                    w.write(line)
                    break
                else:
                    w.write(line)
            for line in f:
                w.write(line)
        return fp

    def change_blast_version2(self, fp, version='2.2.25+', sub_num=100000):
        '''
        Blast + is specific for blast2go, and specific programs cannot blast+ the current 2.3.0+ version
        Add the file-splitting feature, blast2go only has as many sub_num per input
        '''
        iter_num = 0
        headers = []
        count = 0
        result = []
        split_file = fp + 'split001' + '.xml'
        w=open(split_file, 'w')
        with open(self.path) as f:
            for line in f:
                l_str = line.strip()
                if l_str.startswith('<BlastOutput_version>'):
                    if '2.3.0+' in line:
                        line = line.replace('2.3.0+', version)
                    else:
                        pass
                    #w.write(line)
                    #break
                else:
                    pass
                    #w.write(line)
                if  l_str.startswith('<BlastOutput_iterations>'):
                    headers.append(line)
                    iter_num = 1
                    continue
                else:
                    pass
                if iter_num == 0:
                    headers.append(line)
                else:
                    if(l_str.startswith('<Iteration>')):
                        if count % sub_num == 0:
                            if count != 0:
                                w.write("  </BlastOutput_iterations>\n" + "</BlastOutput>\n")
                                w.close()
                            else:
                                pass
                            iter_num += 1
                            split_file = fp + 'split00' + str(iter_num) + '.xml'
                            w=open(split_file, 'w')
                            result.append(split_file)
                            w.write("".join(headers))
                            w.write(line)
                        else:
                            w.write(line)
                        count += 1
                    else:
                        w.write(line)
            w.close()
        return result
