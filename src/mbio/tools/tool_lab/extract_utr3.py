# coding=utf-8
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import glob
import shutil
import regex
from Bio import SeqIO
from Bio.Alphabet import IUPAC
__author__ = 'liubinxu'


class ExtractUtr3Agent(Agent):
    """
    predict utr3 sequence
    1 extract utr3 from gtf and reference by primer4_utr tag
    2 extract utr3 from gtf and reference by exon - cds
    3 extract utr3 from gtf and reference by cds down 500 bg
    4 extract utr3 from gtf and reference by transcript orf predict
    5 extract utr3 from biomart
    species: vertebrates,plants,insects,fungi
    pro_dis: promoter distance before translation start site
    """
    def __init__(self, parent):
        super(ExtractUtr3Agent, self).__init__(parent)
        options = [
            {"name": "ref", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "transcript", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "gtf", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "utr_length", "type": "int", "default": 500},
            {"name": "Markov_length", "type": "int", "default": 2000},  # 马尔科夫训练长度
            {"name": "p_length", "type": "int", "default": 100},  # 最小蛋白长度
            {'name': "type", "type": 'string', 'default': "genome"}, # genome or biomart
            {'name': 'species_name', 'type': 'string', 'default': ''},
            {'name': 'species_class', 'type': 'string', 'default': ''},
            {'name': 'gene_list', 'type': 'string', 'default': ""}
        ]
        self.add_option(options)

    def check_options(self):

        # if self.option("transcript").is_set:
        #     pass
        # elif self.option("ref").is_set and self.option("gtf").is_set:
        #     pass
        # else:
        #     raise OptionError("必须设置参数 ref 和 gtf  或 transcript")
        pass

    def set_resource(self):
        self._cpu = 10
        self._memory = "{}G".format('30')

    def end(self):
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(ExtractUtr3Agent, self).end()

class ExtractUtr3Tool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(ExtractUtr3Tool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        self.perl_path = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin'
        self.python = '/miniconda2/bin/python'
        self.bedtool_path = self.config.SOFTWARE_DIR + ''
        self.biomart_utr = self.config.PACKAGE_DIR + '/tool_lab/biomart_utr.sh'
        self.MOODS = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/MOODS/python/build/lib.linux-x86_64-2.7'
        self.mood_script = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/MOODS/python/scripts/moods_dna.py'
        self.jaspar_db = self.config.SOFTWARE_DIR + '/database/JASPAR_CORE'
        self.bedtool_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
        self.transmir_db = self.config.SOFTWARE_DIR + '/database/transmir'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PYTHONPATH=self.MOODS)
        self.set_environ(LD_LIBRARY_PATH=self.gcc_lib)
        self.set_environ(PATH=self.perl_path)
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.parafly = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'
        self.transdecoder_path = "bioinfo/gene-structure/TransDecoder-3.0.0/"
        self.pre = dict()
        self.trans_mir = dict()

    def get_gtf_utr(self):
        with open(self.option("gtf").prop['path'], 'r') as gtf_f, open("gtf.utr3.gtf", 'w') as utr_gtf_f:
            utr_num, cds_num, exon_num = 0, 0, 0


            if len(self.choosed_list) == 0:
                for line in gtf_f.readlines():
                    cols = line.split("\t")
                    if len(cols) > 8 and cols[2].lower() in ["three_prime_utr", "UTR3", "UTR_3"]:
                        utr_gtf_f.write(line)
                        utr_num += 1
                    if len(cols) > 8 and cols[2].lower() == "cds":
                        cds_num += 1
                    if len(cols) > 8 and cols[2].lower() == "exon":
                        exon_num += 1
            else:
                content_prog = regex.compile(
                    r'^([^#][^\t]*?)\t+(([^\t]+)\t+){7}(.*;)*\s*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$'
                )


                for line in gtf_f.readlines():
                    cols = line.split("\t")
                    content_m = content_prog.match(line)
                    if content_m:
                        if 'transcript_id' in content_m.captures(6):
                            txpt_id = content_m.captures(7)[0]
                            gene_id = content_m.captures(10)[0]
                        else:
                            txpt_id = content_m.captures(10)[0]
                            gene_id = content_m.captures(7)[0]
                        if gene_id in self.choosed_list or txpt_id in self.choosed_list:
                            if len(cols) > 8 and cols[2].lower() in ["three_prime_utr", "UTR3", "UTR_3"]:
                                utr_gtf_f.write(line)
                                utr_num += 1
                            if len(cols) > 8 and cols[2].lower() == "cds":
                                cds_num += 1
                            if len(cols) > 8 and cols[2].lower() == "exon":
                                exon_num += 1

        return utr_num, cds_num, exon_num

    def gtf_to_bed(self):
        self.logger.info("转换gtf文件为bed文件")
        new_gtf = self.work_dir + "/" + os.path.basename(self.option("gtf").prop["path"])
        if os.path.exists(new_gtf):
            os.remove(new_gtf)
        os.link(self.option("gtf").prop["path"], new_gtf)
        self.option("gtf").set_path(new_gtf)
        self.option("gtf").to_bed()

    def filter_bed(self):
        with open(self.option("gtf").prop["path"] + ".bed", "r") as f, open(self.option("gtf").prop["path"] + ".filter.bed", "w") as w:
            for line in f:
                if int(line.split("\t")[2]) < 500000000:
                    w.write(line)
        self.option("bed").set_path(self.option("gtf").prop["path"] + ".filter.bed")

    def get_bed_utr(self):
        utr_num = 0
        with open(self.option("gtf").prop["path"] + ".bed", "r") as f, open("temp.utr3.bed", "w") as utr_w:
            for line in f:
                cols = line.strip("\n").split("\t")
                if cols[5] == "+":
                    if int(cols[2]) - int(cols[7]) > 10:
                        # extract utr exon sequence
                        exon_len = cols[10].strip(",").split(",")
                        exon_start = cols[11].strip(",").split(",")
                        utr_len = []
                        utr_start = []
                        for n,start in enumerate(exon_start):
                            g_start = int(cols[1]) + int(start)
                            g_end = int(cols[1]) + int(start) + int(exon_len[n])
                            if g_start <= int(cols[7]) and g_end > int(cols[7]):
                                utr_start.append(str(0))
                                utr_len.append(str(int(g_end) - int(cols[7])))
                            elif g_start > int(cols[7]) :
                                utr_start.append(str(g_start - int(cols[7])))
                                utr_len.append(exon_len[n])

                        num = str(len(utr_len))
                        utr_w.write("\t".join([
                            cols[0], cols[7], cols[2], cols[3], cols[4], cols[5], \
                            cols[7], cols[2], "0", num, ",".join(utr_len) + ",", ",".join(utr_start) + ","
                        ]) + "\n")
                        utr_num += 1
                else:
                    if int(cols[6]) - int(cols[1]) > 10:
                        # extract utr exon sequence
                        exon_len = cols[10].strip(",").split(",")
                        exon_start = cols[11].strip(",").split(",")
                        utr_len = []
                        utr_start = []
                        for n,start in enumerate(exon_start):
                            g_start = int(cols[1]) + int(start)
                            g_end = int(cols[1]) + int(start) + int(exon_len[n])
                            end = int(start) + int(exon_len[n])
                            if g_start < int(cols[6]) and g_end > int(cols[6]):
                                utr_start.append(str(g_start - int(cols[1])))
                                utr_len.append(str(int(cols[6]) - g_start))
                            elif g_end < int(cols[6]):
                                utr_start.append(str(g_start - int(cols[1])))
                                utr_len.append(exon_len[n])

                        num = str(len(utr_len))
                        utr_w.write("\t".join([
                            cols[0], cols[1], cols[6], cols[3], cols[4], cols[5], \
                            cols[6], cols[1], "0", num, ",".join(utr_len) + ",", ",".join(utr_start) + ","
                        ]) + "\n")
                        utr_num += 1
        return utr_num

    def get_bed_cds_utr(self):
        utr_num = 0
        with open(self.option("gtf").prop["path"] + ".bed", "r") as f, open("temp.utr3.bed", "w") as utr_w:
            for line in f:
                cols = line.strip("\n").split("\t")
                if cols[5] == "+":
                    utr_w.write("\t".join([
                        cols[0], cols[7], str(int(cols[7]) + self.option("utr_length")), cols[3], cols[4], cols[5]
                    ]) + "\n")
                    utr_num += 1
                else:
                    cols_start = int(cols[6]) - 500
                    if cols_start < 0:
                        cols_start = 0
                    if int(cols[6]) < 20:
                        continue
                    utr_w.write("\t".join([
                        cols[0], str(cols_start), cols[6], cols[3], cols[4], cols[5]
                    ]) + "\n")
                    utr_num += 1
        return utr_num


    def run_gtf_to_fa(self, gtf, ref, fasta):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """

        cmd = self.cufflinks_path + "gffread %s -g %s -w %s" % (
            gtf, ref, fasta)
        self.logger.info('运行gtf_to_fasta，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

    def td_longorfs(self, transcript):
        self.logger.info(self.option("p_length"))
        cmd = "{}TransDecoder.LongOrfs -t {} -m {}".format(self.transdecoder_path, transcript, self.option("p_length"))
        print(cmd)
        self.logger.info("开始提取长orf")
        command = self.add_command("transdecoder_longorfs", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("提取结束！")
        else:
            self.set_error("提取orf过程出错")

    def td_predict(self, transcript, hmm_out=None):
        if hmm_out is None:
            hmm_out = ""
        else:
            hmm_out = "--retain_pfam_hits {}".format(hmm_out)
        cmd = "{}TransDecoder.Predict -t {} -T {}".format(
            self.transdecoder_path, transcript, self.option("Markov_length"))
        print(cmd)
        self.logger.info("开始预测编码区域")
        command = self.add_command("transdecoder_predict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("预测编码区域完成！")
        else:
            self.set_error("预测过程过程出错")


    def predict_orf(self, transcript):
        """
        运行
        """
        self.td_longorfs(transcript)
        self.td_predict(transcript)


    def get_orf_utr(self, orf_bed):
        with open(orf_bed, "r") as f, open("temp.utr3.bed", "w") as utr_w:
            trans = []
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) > 11:
                    if cols[0] in trans:
                        pass
                    else:
                        trans.append(cols[0])
                        if "complete_len" in cols[3] or "5prime_partial_len" in cols[3]:
                            if cols[5] == "+":
                                utr_w.write("\t".join([
                                    cols[0], cols[7], cols[2], cols[0], cols[4], cols[5]
                                ]) + "\n")
                            else:
                                utr_w.write("\t".join([
                                    cols[0], cols[1], cols[6], cols[0], cols[4], cols[5]
                                ]) + "\n")
                else:
                    pass

    def get_bed_exon_utr(self):
        utr_num = 0

        self.run_gtf_to_fa(self.option("gtf").prop["path"],
                           self.option("ref").prop["path"],
                           self.work_dir + "/transcript.fa"
        )

        self.predict_orf(self.work_dir + "/transcript.fa")
        self.get_orf_utr("transcript.fa.transdecoder.bed")


    def get_fasta_from_bed(self, bed, fasta, utr):
        cmd = "{bedtool_path} getfasta -fi {fna} -bed {bed} -s -name -fo {utr3}".format(
            bedtool_path=self.bedtool_path, fna=fasta, bed=bed, utr3=utr)
        os.system(cmd)

    def set_output(self):
        if os.path.exists(self.output_dir + '/utr3.fa'):
            os.remove(self.output_dir + '/utr3.fa')
        os.link(self.work_dir + '/utr3.fa', self.output_dir + '/utr3.fa')

    def run_utr3_biomart(self):
        species_name = self.option("species_name")
        db_class = self.option("species_class")

        spes = species_name.split(" ")
        a = spes[0][0]
        b = "".join(spes[1:])
        abr = (a+b).lower()

        mart = db_class + '_mart'
        bio_class = db_class

        if db_class == "vertebrates" or db_class == "":
            mart = "default"
            abr = abr
            bio_class = "asia"

        get1 = "/sh {} {} {} {}".format(self.biomart_utr, abr, mart, bio_class)
        outfile = "utr3.txt"


        self.logger.info('下载，fasta文件')
        command = self.add_command("biomart_download", get1, ignore_error=True)
        command.software_dir = "/bin"
        if os.path.exists(self.work_dir + "/" + abr+"_gene.txt"):
            pass
        else:
            command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("biomar完成")
        else:
            if command.return_code == 4:
                self.logger.info("useast 连接")
                get1 = "/sh {} {} {} {}".format(self.biomart_utr, abr, mart, "useast")
                command = self.add_command("biomart_download_retry", get1)
                command.software_dir = "/bin"
                command.run()
                self.wait(command)
                if command.return_code == 0:
                    pass
                else:
                    self.set_error("biomart useast 无法连接!")
            elif os.path.getsize(self.work_dir + "/" + abr+"_gene.txt") > 0:
                pass
            else:
                self.set_error("biomart无法连接!")
        # if os.path.exists(self.rel_dir + 'get_utr' + '.finished'):
        #     print "go 下载成功， 跳过此步"
        # else:
        #     while 1:
        #         f = self.run_cmd(get1)
        #         if f != 130:
        #             break

        with open("utr3.fa", 'w') as f:
            for seq in SeqIO.parse(abr+"_gene.txt", "fasta"):
                ids = seq.id.split("|")
                if ids[0] in self.choosed_list or ids[2] in self.choosed_list or len(self.choosed_list) == 0:
                    f.write(">{} {}\n{}\n".format(ids[2], ids[0], seq.seq))

    def get_choose_list(self):
        self.choosed_list = []
        if  os.path.exists(self.option("gene_list")):
            with open(self.option("gene_list"), 'r') as f:
                for line in f:
                    self.choosed_list.append(line.strip())

    def run(self):
        super(ExtractUtr3Tool, self).run()
        self.get_choose_list()
        if self.option("transcript").is_set:
            if os.path.exists(os.path.join(self.work_dir,"transcript.fa")):
                os.remove(os.path.join(self.work_dir,"transcript.fa"))
            os.link(self.option("transcript").prop['path'], os.path.join(self.work_dir,"transcript.fa"))
            self.predict_orf(os.path.join(self.work_dir,"transcript.fa"))
            self.get_orf_utr(os.path.join(self.work_dir,"transcript.fa.transdecoder.bed"))
            self.get_fasta_from_bed("temp.utr3.bed", "transcript.fa", "utr3.fa")
        elif self.option("type") == "genome":
            self.gtf_to_bed()
            utr_num, cds_num, exon_num = self.get_gtf_utr()
            self.run_gtf_to_fa(self.work_dir + "/gtf.utr3.gtf", self.option("ref").prop['path'], "utr3.fa")
        elif self.option("type") == "biomart":
            self.run_utr3_biomart()



            '''
            if utr_num > 1000:
                self.run_gtf_to_fa(self.work_dir + "/gtf.utr3.gtf", self.option("ref").prop['path'], "utr3.fa")


            elif cds_num > 1000:
                utr_num = self.get_bed_utr()
                if utr_num > 1000:
                    self.get_fasta_from_bed("temp.utr3.bed", self.option("ref").prop['path'], "utr3.fa")
                    pass
                else:
                    utr_num = self.get_bed_cds_utr()
                    if utr_num > 1000:
                        self.get_fasta_from_bed("temp.utr3.bed", self.option("ref").prop['path'], "utr3.fa")
                        pass
                    else:
                        self.get_bed_exon_utr()
                        self.get_fasta_from_bed("temp.utr3.bed", "transcript.fa", "utr3.fa")

            else:
                self.get_bed_exon_utr()
                self.get_fasta_from_bed("temp.utr3.bed", "transcript.fa", "utr3.fa")
            '''

        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data5'
        data = {
            "id": "extract_utr" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.extract_utr3",
            "instant": False,
            "options": dict(
                ref = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/Ensemble_release_36/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa",
                gtf = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/Ensemble_release_36/gtf/Arabidopsis_thaliana.TAIR10.36.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
