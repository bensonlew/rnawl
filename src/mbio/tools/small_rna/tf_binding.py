# coding=utf-8
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import glob
import shutil
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
__author__ = 'liubinxu'


class TfBindingAgent(Agent):
    """
    predict small_rna TF binding site
    db: vertebrates,plants,insects,fungi
    species: species for trans mir
    pro_dis: promoter distance before translation start site
    """
    def __init__(self, parent):
        super(TfBindingAgent, self).__init__(parent)
        options = [
            {'name': 'species', 'type': 'string', 'default': 'vertebrates'},
            {'name': 'sub_species', 'type': 'string', 'default': None},
            {"name": "known_pre", "type": "infile", "format": "small_rna.common"},
            {"name": "novol_pre", "type": "infile", "format": "small_rna.common"},
            {'name': 'tss_up', 'type': 'string', 'default': 5000},
            {'name': 'tss_down', 'type': 'string', 'default': 1000},
            {"name": "ref", "type": "infile", "format": "small_rna.fasta"},
            {"name": "method", "type": "string", "format": "meme"},
            {"name": "pro_dis", "type": "int", "default": 500},
            {"name": "mood_pvalue", "type": "float", "default": 0.05},
            {"name": "thresh", "type": "float", "default": 0.0001},
            {"name": "known_species", "type": "string", "default": None},
            {"name": "db", "type": "string", "default": None},
            {'name': 'gene_link', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'gene_detail', 'type': 'infile', 'format': 'small_rna.common'},
            # {"name": "mood_score", "type": "float", "default": 10}
        ]
        self.add_option(options)



    def check_options(self):
        if (not self.option("known_pre").is_set) and (not self.option("novol_pre").is_set):
            raise OptionError("必须设置参数known_pre或者novol_pre]")
        if not self.option("ref").is_set:
            raise OptionError("必须设置参数 ref")
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
        super(TfBindingAgent, self).end()

class TfBindingTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(TfBindingTool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.python = '/program/Python/bin/python'
        self.bedtool_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
        self.MOODS = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/MOODS/python/build/lib.linux-x86_64-2.7'
        self.mood_script = self.config.SOFTWARE_DIR + '/bioinfo/miRNA/MOODS/python/scripts/moods_dna.py'
        self.jaspar_db = self.config.SOFTWARE_DIR + '/database/JASPAR_CORE'
        self.fimo = '/bioinfo/rna/meme_4.12.0/bin/fimo'
        self.transmir_db = self.config.SOFTWARE_DIR + '/database/transmir'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PYTHONPATH=self.MOODS)
        self.set_environ(LD_LIBRARY_PATH=self.gcc_lib)
        self.parafly = '/bioinfo/denovo_rna_v2/trinityrnaseq-Trinity-v2.5.0/trinity-plugins/ParaFly-0.1.0/bin/ParaFly'

        self.pre = dict()
        self.trans_mir = dict()
        self.tss = dict()
        self.species2trans_mir = {
            "Homo_sapiens": "hsa",
            "Rattus_norvegicus": "rno",
            "Mus_musculus": "mmu",
            "Arabidopsis_thaliana": "ath",
            "Oryza_sativa": "oth",
            "Gallus_gallus": "gga",
            "Danio_rerio": "dre",
            "Drosophila_melanogaster": "dme",
            "Caenorhabditis_elegans": "cel"
        }

    def get_literature(self):
        """
        获取数据库关系对
        """
        db = self.transmir_db + "/literature/" + self.species2trans_mir[self.option("species")] + ".tsv"
        with open(db, 'r') as db_f:
            for line in db_f.readlines():
                cols = line.strip().split("\t")
                bind = cols[1] + " " + cols[0]
                self.trans_mir[bind] = {
                    "regulate": cols[4],
                    "pubmed_id": cols[5]
                }
    def get_name2link(self):
        name2link = dict()
        with open(self.option("gene_link").prop['path'], 'r') as gene_link_f:
            gene_link_f.readline()
            gene2link = {line.strip("\n").split("\t")[0]: line.strip("\n").split("\t")[1] for line in gene_link_f}
        with open(self.option("gene_detail").prop['path'], 'r') as gene_name_f:
            gene_name_f.readline()
            for line in gene_name_f:
                name = line.strip("\n").split("\t")[2].lower()
                gene = line.strip("\n").split("\t")[1]
                if name in name2link:
                    if gene in gene2link:
                        name2link[name].update({gene: gene2link[gene]})
                else:
                    if gene in gene2link:
                        name2link[name] = {gene: gene2link[gene]}
        return name2link

    def get_bed(self):
        """
        获取启动子区域bed文件
        """

        distance = self.option("pro_dis")
        with open("small_rna.promoter.bed", 'w') as pro_w:
            if self.option("known_pre").is_set:
                with open(self.option("known_pre").prop['path'], 'r') as known_f:
                    known_f.readline()
                    for line in known_f.readlines():
                        cols = line.strip().split("\t")
                        mirna_id = cols[0]
                        pre_id = cols[3]
                        if pre_id in self.pre:
                            self.pre[pre_id].append(mirna_id)
                            continue
                        else:
                            self.pre[pre_id] = [mirna_id]
                        chr_id = cols[-1].split("(")[0]
                        if chr_id.startswith("chr"):
                            chr_id = chr_id[3:]
                        try:
                            strand = cols[-1].split(")")[0].split("(")[1]
                        except:
                            continue
                        strand = cols[-1].split(")")[0].split("(")[1]
                        start, end = cols[-1].split(":")[1].split("-")
                        if strand == "+":
                            pro_start = str(max(int(start) - int(self.option("tss_up")), 0))
                            pro_end = str(int(start) + int(self.option("tss_down")))
                            self.tss[pre_id] = {"strand": strand, "pos": int(start), "chr": chr_id, "begin": pro_start}
                        else:
                            pro_start = str(max(int(end) - int(self.option("tss_down")), 0))
                            pro_end = str(int(end) + int(self.option("tss_up")))
                            self.tss[pre_id] = {"strand": strand, "pos": int(end), "chr": chr_id, "begin": pro_start}
                        pro_w.write("\t".join([chr_id, pro_start, pro_end, pre_id, "0", strand]) + "\n")

            if self.option("novol_pre").is_set:
                with open(self.option("novol_pre").prop['path'], 'r') as known_f:
                    known_f.readline()
                    for line in known_f.readlines():
                        cols = line.strip().split("\t")
                        mirna_id = cols[0]
                        pre_id = cols[3]
                        if pre_id in self.pre:
                            self.pre[pre_id].append(mirna_id)
                            continue
                        else:
                            self.pre[pre_id] = [mirna_id]
                        chr_id = cols[6].split("(")[0]
                        if chr_id.startswith("chr"):
                            chr_id = chr_id[3:]
                        strand = cols[6].split(")")[0].split("(")[1]
                        start, end = cols[6].split(":")[1].split("-")
                        if strand == "+":
                            pro_start = str(max(int(start) - int(self.option("tss_up")), 0))
                            pro_end = str(int(start) + int(self.option("tss_down")))
                            self.tss[pre_id] = {"strand": strand, "pos": int(start), "chr": chr_id, "begin": pro_start}
                        else:
                            pro_start = str(max(int(end) - int(self.option("tss_down")), 0))
                            pro_end = str(int(end) + int(self.option("tss_up")))
                            self.tss[pre_id] = {"strand": strand, "pos": int(end), "chr": chr_id, "begin": pro_start}


                        pro_w.write("\t".join([chr_id, pro_start, pro_end, pre_id, "0", strand]) + "\n")


    def get_fasta_from_bed(self):
        cmd = "{bedtool_path} getfasta -fi {fna} -bed {work_dir}/small_rna.promoter.bed -s -name -fo {work_dir}/smallrna.promoter.fa".format(
            bedtool_path=self.bedtool_path, fna=self.option("ref").prop['path'], work_dir = self.work_dir)
        os.system(cmd)
        '''
        with open(self.work_dir + '/Rockhopper_Results/cds.fa', 'r') as fa_r, \
                open(self.work_dir + '/Rockhopper_Results/cds.faa', 'w') as faa_w:
            for block in fa_r.read().split('\n>'):
                block = block.lstrip('>').split('\n')
                coding_dna = Seq(''.join(block[1:]), IUPAC.ambiguous_dna)
                protein = coding_dna.translate()
                faa_w.write('>' + block[0].strip() + '\n' + str(protein) + '\n')

        '''
    def binding_predict(self):
        if self.option("db") in ['vertebrate', 'plant', 'insect']:
            db_dir = self.option("db") + "s"
        else:
            db_dir = self.option("db")
        pfms =  glob.glob("{}/{}/*.pfm".format(self.jaspar_db, db_dir))
        with open("mood.sh", 'w') as para_f:
            for pfm in pfms:
                cmd = "{} {}  -m {} -s {} -p {} -o {}.binding.csv\n".format(self.python_path, self.mood_script, pfm, "smallrna.promoter.fa", self.option("mood_pvalue"), os.path.basename(pfm))
                para_f.write(cmd)

        cmd = '{} '.format(self.parafly)
        cmd += '-{} {} '.format("c", "mood.sh")
        cmd += '-{} {} '.format("CPU", '10')
        cmd += '-v -shuffle'
        cmd_name = 'para_mood'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()

        os.system("cat *.binding.csv > binding.csv")
        self.logger.info("开始运行mood")
        mood = self.add_command("mood", cmd).run()
        self.wait(mood)
        if mood.return_code == 0:
            self.logger.info("mood")
        else:
            self.set_error("mood 运行错误")

        with open("binding.csv", 'r') as mood_f, open("binding.xls", 'w') as mood_w:
            if self.option("species") in self.species2trans_mir:
                mood_w.write("small_rna\tsmall_rna_pre\tTF\tbinding_site\tstrand\tscore\ttransmir_reg\ttransmir_pubmed\n")
            else:
                mood_w.write("small_rna\tsmall_rna_pre\tTF\tbinding_site\tstrand\tscore\n")
            tf_binding = dict()
            for line in mood_f.readlines():
                cols = line.strip().split(",")
                cols[1] = cols[1].split(".")[0]
                bind = cols[0] + " " + cols[1]
                if bind in tf_binding:
                    if tf_binding[bind]['score'] < float(cols[4]):
                        tf_binding[bind]['score'] = float(cols[4])
                        tf_binding[bind]['pos'] = cols[2]
                        tf_binding[bind]['strand'] = cols[3]
                else:
                    tf_binding[bind] = {
                        'score': float(cols[4]),
                        'pos': cols[2],
                        'strand': cols[3]
                    }
            # print tf_binding
            if self.option("species") in self.species2trans_mir:
                # print tf_binding
                print self.trans_mir
                print self.tss.keys()
                for bind in tf_binding:
                    if bind in self.trans_mir:
                        transmir_reg = self.trans_mir[bind]["regulate"]
                        pubmed = self.trans_mir[bind]["pubmed_id"]
                    else:
                        transmir_reg = ""
                        pubmed = ""

                    mood_w.write("\t".join([
                        ";".join(self.pre[bind.split(" ")[0]]),
                        bind.split(" ")[0],
                        bind.split(" ")[1],
                        tf_binding[bind]['pos'],
                        tf_binding[bind]['strand'],
                        str(tf_binding[bind]['score']),
                        transmir_reg,
                        pubmed
                    ]) + "\n")

                known_binding = set(self.trans_mir.keys())
                predict_binding = set(tf_binding.keys())
                for bind in known_binding - predict_binding:
                    if bind.split(" ")[0] in self.pre:
                        if bind in self.trans_mir:
                            transmir_reg = self.trans_mir[bind]["regulate"]
                            pubmed = self.trans_mir[bind]["pubmed_id"]
                        else:
                            transmir_reg = ""
                            pubmed = ""
                        mood_w.write("\t".join([
                            ";".join(self.pre[bind.split(" ")[0]]),
                            bind.split(" ")[0],
                            bind.split(" ")[1],
                            "",
                            "",
                            "",
                            transmir_reg,
                            pubmed
                        ]) + "\n")
            else:
                for bind in tf_binding:
                    mood_w.write("\t".join([
                        ";".join(self.pre[bind.split(" ")[0]]),
                        bind.split(" ")[0],
                        bind.split(" ")[1],
                        tf_binding[bind]['pos'],
                        tf_binding[bind]['strand'],
                        str(tf_binding[bind]['score']),
                    ]) + "\n")

    def get_sub_db(self):
        '''
        获取单物种的文件
        '''
        spe = self.jaspar_db + '/{}.species.tsv'.format(self.option("db"))
        spe_dir = self.jaspar_db + '/{}/'.format(self.option("db"))

        '''
        if self.option("db") == 'vertebrates':
            spe = self.jaspar_db + '/vertebrates.species.tsv'
            spe_dir = self.jaspar_db + '/vertebrates/'
        elif self.option("db") == 'plants':
            spe = self.jaspar_db + '/plants.species.tsv'
            spe_dir = self.jaspar_db + '/plants/'
        '''
        tf_ids = []
        with open(spe, 'r') as f:
            for line in f:
                cols = line.split("\t")
                if self.option("sub_species") in cols[2] or self.option("sub_species") in ["all", ""]:
                    tf_ids.append(cols[0])
        for tf_id in tf_ids:
            os.system("cat {} >> sub_spe.meme".format(spe_dir + tf_id + '.meme'))
        return "sub_spe.meme"



    def binding_predict_meme(self):
        name2link = self.get_name2link()

        # self.option("sub_species"):
        db_dir = self.get_sub_db()
        '''
            else:
            if self.option("db") in ['vertebrate', 'plant', 'insect']:
                db_dir = self.jaspar_db + "/" + self.option("db") + "s" + '/{}.meme'.format(self.option("db") + "s")
            else:
                db_dir = self.option("db")
        '''

        cmd = "{} ".format(self.fimo)
        cmd += "--bgfile --motif-- "
        cmd += "--max-strand "
        cmd += "--thresh {} ".format(self.option("thresh"))
        cmd += "{} {}".format(db_dir, "smallrna.promoter.fa")

        cmd_name = "fimo"
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()

        with open("fimo_out/fimo.txt", 'r') as fimo_f, open("binding.xls", 'w') as fimo_w:
            if self.option("species") in self.species2trans_mir:
                fimo_w.write("small_rna\tsmall_rna_pre\tTF\ttf_gene_id\tpvalue\tbinding_site\tstrand\ttss_postion\ttf_link\tgene_link\ttransmir_reg\ttransmir_pubmed\n")
            else:
                fimo_w.write("small_rna\tsmall_rna_pre\tTF\ttf_gene_id\tpvalue\tbinding_site\tstrand\ttss_postion\ttf_link\tgene_link\n")
            tf_binding = dict()
            fimo_f.readline()
            for line in fimo_f.readlines():
                cols = line.strip().split("\t")
                # cols[1] = cols[1].split(".")[0]
                bind = cols[2] + " " + cols[1]
                if bind in tf_binding:
                    pass
                else:
                    tf_binding[bind] = {
                        'motif_id': cols[0],
                        'score': float(cols[6]),
                        'pvalue': cols[7],
                        'pos': cols[3],
                        'strand': cols[5],
                        'start': int(cols[3]),
                        'end': int(cols[4]),
                    }
            # print tf_binding.keys()
            if self.option("species") in self.species2trans_mir:
                # print tf_binding
                # print self.trans_mir
                print self.tss.keys()
                for bind in tf_binding:
                    tss_info = self.tss[bind.split(" ")[0]]
                    tss_start = "{}:{}".format(tss_info["chr"], tss_info["pos"])
                    tf_link = "http://jaspar.genereg.net/matrix/{}/".format(tf_binding[bind]['motif_id'])
                    if bind in self.trans_mir:
                        transmir_reg = self.trans_mir[bind]["regulate"]
                        pubmed = self.trans_mir[bind]["pubmed_id"]
                    else:
                        transmir_reg = ""
                        pubmed = ""

                    name_str = bind.split(" ")[1]
                    names = name_str.split("(")[0].split("::")
                    gene_list = list()
                    link_list = list()
                    for name in names:
                        if name.lower() in name2link:
                            gene_list.extend(name2link[name.lower()].keys())
                            link_list.extend(name2link[name.lower()].values())

                    fimo_w.write("\t".join([
                        ";".join(self.pre[bind.split(" ")[0]]),
                        bind.split(" ")[0],
                        bind.split(" ")[1],
                        ";".join(gene_list),
                        tf_binding[bind]['pvalue'],
                        "{}-{}".format(int(tss_info["begin"]) + int(tf_binding[bind]['start']), int(tss_info["begin"]) + int(tf_binding[bind]['end'])),
                        tf_binding[bind]['strand'],
                        tss_start,
                        tf_link,
                        ";".join(link_list),
                        transmir_reg,
                        pubmed
                    ]) + "\n")

                known_binding = set(self.trans_mir.keys())
                predict_binding = set(tf_binding.keys())
                for bind in known_binding - predict_binding:
                    name_str = bind.split(" ")[1]
                    names = name_str.split("(")[0].split("::")
                    gene_list = list()
                    link_list = list()
                    for name in names:
                        if name.lower() in name2link:
                            gene_list.extend(name2link[name.lower()].keys())
                            link_list.extend(name2link[name.lower()].values())
                    if bind.split(" ")[0] in self.pre:
                        if bind in self.trans_mir:
                            transmir_reg = self.trans_mir[bind]["regulate"]
                            pubmed = self.trans_mir[bind]["pubmed_id"]
                        else:
                            transmir_reg = ""
                            pubmed = ""
                        fimo_w.write("\t".join([
                            ";".join(self.pre[bind.split(" ")[0]]),
                            bind.split(" ")[0],
                            bind.split(" ")[1],
                            ";".join(gene_list),
                            "",
                            "",
                            "",
                            "",
                            "",
                            ";".join(link_list),
                            transmir_reg,
                            pubmed
                        ]) + "\n")
            else:
                for bind in tf_binding:
                    tss_info = self.tss[bind.split(" ")[0]]
                    tss_start = "{}:{}".format(tss_info["chr"], tss_info["pos"])
                    tf_link = "http://jaspar.genereg.net/matrix/{}/".format(tf_binding[bind]['motif_id'])

                    name_str = bind.split(" ")[1]
                    names = name_str.split("(")[0].split("::")
                    gene_list = list()
                    link_list = list()
                    for name in names:
                        if name.lower() in name2link:
                            gene_list.extend(name2link[name.lower()].keys())
                            link_list.extend(name2link[name.lower()].values())


                    fimo_w.write("\t".join([
                        ";".join(self.pre[bind.split(" ")[0]]),
                        bind.split(" ")[0],
                        bind.split(" ")[1],
                        ";".join(gene_list),
                        tf_binding[bind]['pvalue'],
                        "{}-{}".format(int(tss_info["begin"]) + int(tf_binding[bind]['start']), int(tss_info["begin"]) + int(tf_binding[bind]['end'])),
                        tf_binding[bind]['strand'],
                        tf_binding[bind]['pos'],
                        tss_start,
                        tf_link,
                        ";".join(link_list),
                    ]) + "\n")


    def set_output(self):
        if os.path.exists(self.output_dir + '/binding.xls'):
            os.remove(self.output_dir + '/binding.xls')
        os.link(self.work_dir + '/binding.xls', self.output_dir + '/binding.xls')

    def run(self):
        super(TfBindingTool, self).run()
        if self.option("species") in self.species2trans_mir:
            self.get_literature()
        self.get_bed()
        self.get_fasta_from_bed()
        self.binding_predict_meme()
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
            "id": "tf_binding" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.tf_binding",
            "instant": False,
            "options": dict(
                ref = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Arabidopsis_thaliana/Ensemble_release_36/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa",
                known_pre = test_dir + "/" + "known_mirna_detail",
                novol_pre = test_dir + "/" + "novel_miR_mature_infor.xls",
                species = "ath",
                db = "plant"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
