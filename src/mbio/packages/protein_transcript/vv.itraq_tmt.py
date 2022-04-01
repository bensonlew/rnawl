# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import ConfigParser
import time
import os


class ItraqTmt(Controller):
    def __init__(self, config):
        super(ItraqTmt, self).__init__()
        self.bin = "/mnt/ilustre/users/ting.kuang/ITRAQ/bin"
        self.db = "/mnt/ilustre/users/ting.kuang/ITRAQ/db"
        self.config = config
        if not os.path.exists(self.config):
            cmd = 'python %s/creat_ITRAQ_seq_config.py %s %s' % (bin, self.config, os.getcwd())
            exit('您传入的配置文件不存在，可以用以下命令生成:\n%s'%cmd)
        self.get_global()

    def get_global(self):
        cfg = ConfigParser.ConfigParser()
        cfg.read(self.config)
        if cfg.get("Global", "Confirm") == 'F':
            exit('请确认配置文件后,将Globle的Confirm参数值改为T,再提交运行！')
        # 各类参数的读入
        self.output = cfg.get("Global", "output")
        if not os.path.exists(self.output):
            os.mkdir(self.output)
        self.raw_data = cfg.get("Global", "raw_data")
        if not os.path.exists(self.raw_data) or not os.path.isdir(self.raw_data):
            exit('请传入正确的原始数据路径，%s不存在或者不是目录'%self.raw_data)
        self.go_obo = cfg.get("GO", "go_obo")
        self.blast2go = cfg.get("GO", "blast")
        self.go_identity = cfg.get("GO", "identity")
        self.blast2kegg = cfg.get("KEGG", "blast")
        self.uniprot2go = cfg.get("GO", "uniprot2go")
        self.KEGG_dir = cfg.get("KEGG", "KEGG_dir")
        self.uniprot2kegg = cfg.get("KEGG", "uniprot2kegg")
        # time_ = time.strftime("%Y-%m-%d-%H-%M", time.localtime())
        self.dup = cfg.get("Sample_Info", "dup")
        self.cog = cfg.get("COG", "COG")
        self.COGclassification = cfg.get("COG", "COGclassification")
        self.KOGclassification = cfg.get("COG", "KOGclassification")
        self.uniprot2cog = cfg.get("COG", "uniprot2cog")
        self.uniprot2kog = cfg.get("COG", "uniprot2kog")
        self.up = cfg.get("Sample_Info", "up")
        self.down = cfg.get("Sample_Info", "down")
        self.pvalue = cfg.get("Sample_Info", "pvalue")
        self.mb = cfg.get("Heatmap", "mb")
        self.ml = cfg.get("Heatmap", "ml")
        self.clt = cfg.get("Heatmap", "clt")
        self.orgC = cfg.get("KEGG", "orgC")
        self.speics = cfg.get("STRING", "speics")
        self.go_blast_db = cfg.get("GO", "blast_db")
        self.KEGG_dir = self.KEGG_dir + '.' + self.orgC

    def add_commands(self):
        self.diff = self.add_command(name='diff')
        self.venn = self.add_command(name='venn')
        self.qc = self.add_command(name='qc')
        self.cluster = self.add_command(name='cluster')

    def how_run(self):
        self.after_one(self.diff, self.run_venn)
        self.after_one(self.diff, self.run_cluster)
        self.after_one(self.venn, self.run_qc)
        self.after_some([self.diff, self.venn, self.qc], self.end)

    def end(self):
        super(ItraqTmt, self).end(out='your program finished')

    def run_diff(self):
        cmd = r'''
${bin}/run_diff_analysis.pl -i ${cwd}/exp.txt -g ${cwd}/sample.config -u ${up} -d ${down}
for i in `ls *.diff.exp.xls | awk -F ".diff.exp.xls" '{print $1}'`
do
    less $i.diff.exp.xls | awk '{if($6<${pvalue} && ($4>${up} || $4<${down})){print $1}}' > $i.DE.list
    less $i.diff.exp.xls | awk '{if($6<${pvalue} && $4>${up}){print $1}}' > $i.up.list
    less $i.diff.exp.xls | awk '{if($6<${pvalue} && $4<${down}){print $1}}' > $i.down.list
    less $i.diff.exp.xls | awk -F "\t" '{printf $1"\t"$5"\t"$6"\t"; if(NR==1){print "sig"}else if($6<0.05){if($4>${up}){printf "up"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else if($4<${down}){printf "down"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else{print "nosig"}}else{print "nosig"}}' > $i.volcano
    less $i.diff.exp.xls | awk -F "\t" '{printf $1"\t"$2"\t"$3"\t"; if(NR==1){print "sig"}else if($6<0.05){if($4>${up}){printf "up"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else if($4<${down}){printf "down"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else{print "nosig"}}else{print "nosig"}}' > $i.scatter
    ${bin}/Highchart_for_ITRAQ.pl -yAxis_log -yAxis_min 0.0001 -yAxis_max 1 -type scatter -t $i.volcano -yAxis_reversed -scatter_series down-p-0.01,down-p-0.05,nosig,up-p-0.05,up-p-0.01 -width 700 -height 500 -color_the "'#2222FF','#22CCFF','#222222','#FFCC22','#FF2222'" -scatter_size 2,2,1,2,2 -scatter_symbol "'triangle-down','triangle-down','diamond','triangle','triangle'"  
    ${bin}/Highchart_for_ITRAQ.pl -type scatter -t $i.scatter -xAxis_log -yAxis_log -scatter_series down-p-0.01,down-p-0.05,nosig,up-p-0.05,up-p-0.01 -width 700 -height 500 -color_the "'#2222FF','#22CCFF','#222222','#FFCC22','#FF2222'" -scatter_size 2,2,2,2,2 -scatter_symbol "'triangle-down','triangle-down','diamond','triangle','triangle'"  
done
${bin}/get_diff_up_down.py
'''
        cmd = Template(cmd)
        cmd = cmd.render(bin = self.bin,
                         cwd = self.raw_data,
                         up = self.up,
                         down = self.down,
                         pvalue = self.pvalue)
        params = dict(
            cmd = cmd,
            node = 3,
            memory = 6
        )
        self.diff.set_params(params)
        self.diff.run()

    def run_cluster(self):
        cmd = '''
cat ${diff}/*.DE.list | sort | uniq | awk 'BEGIN {print "Accession"} {print $0}' > ${cluster}/All.DE.list
${bin}/get_exp_from_list.pl All.DE.list ${cwd}/exp.txt all.diffexp.txt
for i in 5 10 20 40 80
do    
    mkdir tmp_$i
    cd tmp_$i
    cp ${cluster}/all.diffexp.txt ./
    ${bin}/plot_heatmap_trendline.pl -i all.diffexp.txt -clt ${clt} -n $i -o $i.Heatmap -mr ${mb} -ml ${ml} 
    cp -r subclusters* ${cluster}/
    cp *pdf ${cluster}/
done
'''
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         cwd=self.raw_data,
                         clt=self.clt,
                         mb=self.mb,
                         cluster=self.cluster.work_dir,
                         diff=self.diff.work_dir,
                         ml=self.ml)
        params = dict(
            cmd=cmd,
            node=3,
            memory=6
        )
        self.cluster.set_params(params)
        self.cluster.run()

    def run_venn(self):
        cmd = '''
cp ${diff}/*.DE.list ./
${bin}/combine_venn.py 
'''
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         diff=self.diff.work_dir
                         )
        params = dict(
            cmd=cmd,
            node=3,
            memory=6
        )
        self.venn.set_params(params)
        self.venn.run()

    def run_qc(self):
        cmd = '''
cp ${cwd}/psm.xlsx ./
cp ${cwd}/peptide.xlsx ./
cp ${cwd}/protein.xlsx ./
cp ${cwd}/Protein_information.xls ./
${bin}/extract_special_col_from_xlsx.py -i psm.xlsx -o psm.xls
${bin}/extract_special_col_from_xlsx.py -i peptide.xlsx -o peptide.xls
${bin}/extract_special_col_from_xlsx.py -i protein.xlsx -o protein.xls

${bin}/dMass_plot.pl -i psm.xls -s "m/z [Da],DeltaM [ppm]"
${bin}/Peptide_length_distribution.pl -i peptide.xls -s "Annotated Sequence"
${bin}/Protein_molecular_weight_distribution.pl -i protein.xls -s "Accession,MW [kDa]"
${bin}/Peptide_number_distribution.pl -i protein.xls -s "Accession,# Peptides"
${bin}/Protein_Coverage_distribution.pl -i protein.xls -s "Accession,Coverage"
${bin}/Protein_information.pl -i Protein_information.xls
'''
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         cwd=self.raw_data,
                         )
        params = dict(
            cmd=cmd,
            node=3,
            memory=6
        )
        self.qc.set_params(params)
        self.qc.run()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_diff()
        super(ItraqTmt, self).fire()


if __name__ == '__main__':
    config = sys.argv[1]
    itraq = ItraqTmt(config)
    itraq.run()