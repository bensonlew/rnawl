# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os, glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import requests, sys
import unittest


class GsekeggAgent(Agent):
    """
    GSEA using KEGG pathway
    """

    def __init__(self, parent):
        super(GsekeggAgent, self).__init__(parent)
        options = [
            {"name": "geneList", "type": "infile", "format": "ref_rna_v2.common"},
            # A csv file contains two columns, one for gene ID (no duplicated allowed) and another one for fold change.
            {"name": "organism", "type": "string", "default": ""},
            # such as 'Arabidopsis thaliana', 'arabidopsis_thaliana', 'Arabidopsis_thaliana'
            {"name": "id_type", "type": "string", "default": "EntrezGene"},
            # EntrezGene or EnsemblGene
            {"name": "nPerm", "type": "int", "default": 1000},
            # The number of permutations, default 1000.
            {"name": "minGSSize", "type": "int", "default": 10},
            # Gene sets smaller than this number are EXLCUDED from the analysis.
            {"name": "maxGSSize", "type": "int", "default": 500},
            # Gene sets bigger than this number are EXLCUDED from the analysis.
            {"name": "pvalueCutoff", "type": "float", "default": 0.05},
            {"name": "nGenesets", "type": "int", "default": 5},
            # The number of gene sets to be displayed on the same figure.
        ]
        self.add_option(options)
        self.step.add_steps("gsekegg")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gsekegg.start()
        self.step.update()

    def stepfinish(self):
        self.step.gsekegg.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('geneList').is_set:
            raise OptionError('基因集文件必须输入')
        if not self.option('organism'):
            raise OptionError('物种信息必须输入')
        if self.option('id_type').lower() not in ['entrezgene', 'ensemblgene']:
            raise OptionError('目前只支持EntrezGene或EnsemblGene两种类型的Gene ID！')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(GsekeggAgent, self).end()


class GsekeggTool(Tool):
    def __init__(self, config):
        super(GsekeggTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self._LD_LIBRARY_PATH = software_dir + "/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/lib64/:$LD_LIBRARY_PATH"
        self._PATH = software_dir + "/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/:$PATH"
        self._C_INCLUDE_PATH = software_dir + "/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/include/:$C_INCLUDE_PATH"
        self.set_environ(PATH=self._PATH, C_INCLUDE_PATH=self._C_INCLUDE_PATH, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/Rscript',
        }
        self.script = {
            'gsekegg': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/gsekegg/gsekegg.r'),
            'get_go_info': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/gsekegg/get_enterz_info.sh'),
            'kegg_organism': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/gsekegg/kegg_organism.txt')
        }
        self.division = {
            'EnsemblPlants': 'plants',
            'EnsemblVertebrates': 'asia',
            'EnsemblFungi': 'fungi',
            'EnsemblProtists': 'protists'
        }
        self.virtualSchemaName = {
            'EnsemblPlants': 'plants_mart',
            'EnsemblVertebrates': 'default',
            'EnsemblFungi': 'fungi_mart',
            'EnsemblProtists': 'protists_mart'
        }
        self.datasetName = {
            'EnsemblPlants': '_eg_gene',
            'EnsemblVertebrates': '_gene_ensembl',
            'EnsemblFungi': '_eg_gene',
            'EnsemblProtists': '_eg_gene'
        }

    def run(self):
        """
        运行
        :return:
        """
        super(GsekeggTool, self).run()
        if self.option('id_type').lower() == 'ensemblgene':
            download = self.get_organism_xrefs()
            if not download:
                self.set_error("请检查物种名是否填写正确")
            self.get_xrefs(download)
        self.run_gsekegg()
        #self.convert_pdf_to_png()
        self.set_output()
        self.end()

    def convert_pdf_to_png(self):
        self.image_magick = '/program/ImageMagick/bin/convert'
        pdfs = glob.glob(self.work_dir + "/*.pdf")
        num = 0
        for pdf in pdfs:
            num += 1
            png = os.path.basename(pdf).replace("pdf", "png")
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + pdf + ' ' + png
            self.logger.info(cmd)
            command = self.add_command('convert_pdf_to_png_{}'.format(num), cmd)
            command.run()
        self.wait()
        if command.return_code == 0:
            pass
        else:
            self.set_error("PDF转换PNG出错!")

    def get_organism_xrefs(self):
        server = "https://rest.ensembl.org"
        if " " in self.option('organism'):
            organism = self.option('organism').replace(" ", "_")
        else:
            organism = self.option('organism')
        ext = "/info/genomes/{}?".format(organism)
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        if decoded:
            self.logger.info(decoded)
            web = self.division[decoded['division']]
            database = self.virtualSchemaName[decoded['division']]
            species = "".join([decoded['name'].split("_")[0][0], decoded['name'].split("_")[1], self.datasetName[decoded['division']]])
            cmd = "wget -O {}_entrez.txt \'http://{}.ensembl.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"{}\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" completionStamp = \"1\"><Dataset name = \"{}\" interface = \"default\" ><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"entrezgene_id\" /></Dataset></Query>\'".format(species, web, database, species)
            outfile = species + "_entrez.txt"
            if os.path.exists('get_entrez' + '.finished'):
                print("enterz下载成功跳过此步")
            else:
                while 1:
                    self.logger.info("命令内容为：{}.".format(cmd))
                    f = os.system(cmd)
                    if f != 130:
                        break
                with open(outfile, 'r') as out:
                    lines = out.readlines()
                    if lines[-1].strip() == "[success]":
                        pass
                    else:
                        print(lines)
            return outfile
        else:
            return False

    def get_xrefs(self, file):
        self.ensembl2ncbi = dict()
        with open(file, "r") as f:
            for line in f:
                if line.strip() == '[success]':
                    pass
                else:
                    items = line.strip().split("\t")
                    if len(items) == 2 and items[0] and items[1]:
                        self.ensembl2ncbi[items[0]] = items[1]
        with open(self.option("geneList").path, "r") as f, open(os.path.join(self.work_dir, 'geneList.txt'), "w") as w:
            for line in f:
                items = line.strip().split()
                gene_id = items[0]
                if gene_id in self.ensembl2ncbi:
                    w.write(self.ensembl2ncbi[gene_id] + "\t" + items[1] + "\n")
        self.option('geneList', os.path.join(self.work_dir, 'geneList.txt'))

    def run_gsekegg(self):
        organism = self.option('organism')
        with open(self.script['kegg_organism'], "r") as f:
            for line in f:
                items = line.strip().split("\t")
                if self.option('organism') == items[1]:
                    break
                elif self.option('organism').lower() == items[2].split(" (")[0].lower():
                    organism = items[1]
                elif "_" in self.option('organism') and self.option('organism').replace("_", " ").lower() == items[2].split(" (")[0].lower():
                    organism = items[1]
        cmd = '{} {}'.format(self.program['rscript'], self.script['gsekegg'])
        cmd += ' -g {}'.format(self.option('geneList').prop["path"])
        cmd += ' -o {}'.format(organism)
        cmd += ' -p {}'.format(self.option('pvalueCutoff'))
        cmd += ' -n {}'.format(self.option('nPerm'))
        cmd += ' -m {}'.format(self.option('minGSSize'))
        cmd += ' -s {}'.format(self.option('nGenesets'))
        cmd_name = 'run_gsekegg'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        file_list = ['gseaplot.pdf', 'gseaupset.pdf', 'gsea.txt']
        for file in file_list:
            if os.path.exists(os.path.join(self.work_dir, file)):
                if os.path.exists(os.path.join(self.output_dir, file)):
                    os.remove(os.path.join(self.output_dir, file))
                os.link(os.path.join(self.work_dir, file), os.path.join(self.output_dir, file))