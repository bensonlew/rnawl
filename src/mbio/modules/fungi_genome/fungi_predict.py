# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modify:2018.5.26

from biocluster.module import Module
import os
from biocluster.core.exceptions import OptionError


class FungiPredictModule(Module):
    def __init__(self, work_id):
        super(FungiPredictModule, self).__init__(work_id)
        options = [
            {"name": "masked", "type": "infile", "format": "sequence.fasta"},
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "cegma_gff", "type":"string", "default": "" },
            {"name": "sample", "type": "string","default":""}, # 样品名称
            {"name": "ref_protein","type": "infile", "format":"sequence.fasta"} #参考蛋白

        ]
        self.add_option(options)

        self.fungi_cds = self.add_module('fungi_genome.fungi_cds')
        self.rna_predict = self.add_module('fungi_genome.dnafungi_predict')
        #self.updata_fasta_info = self.add_tool('fungi_genome.')




    def check_options(self):
        if not self.option("masked").is_set:
            raise OptionError("必须设置参数masked", code="22101001")
        if not self.option("ref_protein").is_set:
            raise OptionError("必须设置参数ref_protein", code="22101002")
        if not self.option("cegma_gff"):
            raise OptionError("必须设置参数cegma_gff", code="22101003")
        if not self.option("sample"):
            raise OptionError("必须设置参数sample", code="22101004")
        return True


    def run_fungi_cds(self):
        self.fungi_cds.set_options({
            "masked": self.option("masked"),
            "fasta" : self.option("fasta"),
            "cegma_gff" : self.option("cegma_gff"),
            "sample" : self.option("sample"),
            "ref_protein" : self.option("ref_protein")
        })
        self.fungi_cds.on("end",self.run_rna_predict)
        self.fungi_cds.run()

    def run_rna_predict(self):
        self.rna_predict.set_options({
            "gene_gff" : os.path.join(self.fungi_cds.output_dir,"{}_CDS.gff".format(self.option('sample'))),
            "genome" : self.option("fasta"),
            "sample" : self.option("sample"),
            "gene_faa" : os.path.join(self.fungi_cds.output_dir,"{}_CDS.faa".format(self.option('sample'))),
            "gene_ffn" : os.path.join(self.fungi_cds.output_dir,"{}_CDS.fnn".format(self.option('sample')))
        })
        self.rna_predict.on('end',self.set_output)
        self.rna_predict.run()





    def set_output(self):
        files_list = ['CDS_predict','rRNA','tRNA']
        for f in files_list:
            f_path = os.path.join(self.output_dir,f)
            if not os.path.exists(f_path):
                os.mkdir(f_path)

        link_file_str = "{1}/CDS_predict/{0}_CDS.gff {1}/CDS_predict/{0}_CDS.faa {1}/CDS_predict/{0}_CDS.fnn {1}/rRNA/{0}_rRNA.fnn {1}/rRNA/{0}_rRNA.gff {1}/tRNA/{0}_tRNA.fnn {1}/tRNA/{0}_tRNA.gff {1}/tRNA/{0}_tRNA.struc".format(self.option('sample'),self.rna_predict.output_dir)
        link_files = link_file_str.split(' ')
        new_files_str = "{0}/CDS_predict/{1}_CDS.gff {0}/CDS_predict/{1}_CDS.faa {0}/CDS_predict/{1}_CDS.fnn {0}/rRNA/{1}_rRNA.fnn {0}/rRNA/{1}_rRNA.gff {0}/tRNA/{1}_tRNA.fnn {0}/tRNA/{1}_tRNA.gff {0}/tRNA/{1}_tRNA.struc".format(self.output_dir,self.option('sample'))
        new_files = new_files_str.split(' ')
        for i in link_files :
            if os.path.exists(new_files[link_files.index(i)]):
                os.remove(new_files[link_files.index(i)])
            os.link(i,new_files[link_files.index(i)])

        new_files2 = "{1}/CDS_predict/{0}_CDS_length.xls {1}/CDS_predict/{0}_CDS_statistics.xls".format(self.option("sample"),self.output_dir).split(" ")
        old_files2 = "{1}/{0}_CDS_length.xls {1}/{0}_CDS_statistics.xls".format(self.option("sample"),self.fungi_cds.output_dir).split(" ")
        for file in new_files2:
            if os.path.exists(file) :
                os.remove(file)
            os.link(old_files2[new_files2.index(file)],file)


        self.end()




    def run(self):
        super(FungiPredictModule, self).run()
        self.run_fungi_cds()


    def end(self):
        super(FungiPredictModule, self).end()
