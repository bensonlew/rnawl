# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import ConfigParser
import os
import glob


class V10Labelfree(Controller):
    def __init__(self, config):
        super(V10Labelfree, self).__init__()
        self.bin = "/mnt/ilustre/users/ting.kuang/LABELFREE/bin"
        self.bin_v = "/mnt/ilustre/users/yitong.feng/scripts/itraq_v10"
        self.bin_l = "/mnt/ilustre/users/yitong.feng/scripts/labelfree_v10"
        self.db = "/mnt/ilustre/users/ting.kuang/ITRAQ/db"
        self.config = config
        if not os.path.exists(self.config):
            cmd = 'python %s/creat_labelfree_v10_config.py %s %s' % (self.bin_l, self.config, os.getcwd())
            os.system(cmd)
            super(V10Labelfree, self).end(normal=1, out='您传入的配置文件不存在，可以用以下命令生成:%s\n而我已经贴心的为你生成了一个\n'%cmd)
        self.get_global()
        self.check_dbs()

    def get_global(self):
        cfg = ConfigParser.ConfigParser()
        cfg.read(self.config)
        if cfg.get("Global", "Confirm") == 'F':
            super(V10Labelfree, self).end(normal=1, out='请确认配置文件后,将Globle的Confirm参数值改为T,再提交运行！')
        # 各类参数的读入
        self.output = cfg.get("Global", "output")
        if not os.path.exists(self.output):
            try:
                os.mkdir(self.output)
            except:
                super(V10Labelfree, self).end(normal=1, out='输出文件夹路径不对')
        self.raw_data = cfg.get("Global", "raw_data")
        self.run_ppi_ = False
        if cfg.get("Global", "run_ppi") == 'T':
            self.run_ppi_ = True
        if not os.path.exists(self.raw_data) or not os.path.isdir(self.raw_data):
            super(V10Labelfree, self).end(normal=1, out='请传入正确的原始数据路径，%s不存在或者不是目录'%self.raw_data)
        self.go_obo = cfg.get("GO", "go_obo")
        self.blast2go = cfg.get("GO", "blast")
        self.go_identity = cfg.get("GO", "identity")
        self.blast2kegg = cfg.get("KEGG", "blast")
        self.uniprot2go = cfg.get("GO", "uniprot2go")
        self.KEGG_dir = cfg.get("KEGG", "KEGG_dir")
        self.kegg_identity = cfg.get("KEGG", "identity")
        self.uniprot2kegg = cfg.get("KEGG", "uniprot2kegg")
        # time_ = time.strftime("%Y-%m-%d-%H-%M", time.localtime())
        self.dup = cfg.get("Sample_Info", "dup")
        self.diff_method = cfg.get("Sample_Info", "diff_method")
        self.diff_pcorrect_method = cfg.get("Sample_Info", "diff_pcorrect_method")
        self.alternative = cfg.get("Sample_Info", "alternative")
        self.cutoffs = cfg.get("Sample_Info", "cutoffs")
        self.cog = cfg.get("COG", "COG")
        self.cog_db = cfg.get("COG", "COG_db")
        self.COGclassification = cfg.get("COG", "COGclassification")
        self.KOGclassification = cfg.get("COG", "KOGclassification")
        self.uniprot2cog = cfg.get("COG", "uniprot2cog")
        self.uniprot2kog = cfg.get("COG", "uniprot2kog")
        self.ppi_identity = cfg.get("PPI", "identity")
        self.ppi_db = cfg.get("PPI", "ppi_db")
        self.up = cfg.get("Sample_Info", "up")
        self.down = cfg.get("Sample_Info", "down")
        self.pvalue = cfg.get("Sample_Info", "pvalue")
        try:
            self.average = cfg.get("Sample_Info", "average")
        except:
            self.average = 'false'
        self.mb = cfg.get("Heatmap", "mb")
        self.ml = cfg.get("Heatmap", "ml")
        self.clt = cfg.get("Heatmap", "clt")
        self.orgC = cfg.get("KEGG", "orgC")
        self.speics = cfg.get("STRING", "speics")
        self.go_blast_db = cfg.get("GO", "blast_db")
        self.KEGG_dir = self.KEGG_dir + '.' + self.orgC
        self.single = cfg.get("KEGG", "single")
        if self.single == 'F':
            self.KEGG_dir = cfg.get("KEGG", "map")

    def check_dbs(self):
        if not os.path.isfile(self.go_obo):
            super(V10Labelfree, self).end(normal=1, out='NEED TO DO FIRST:\n{bin}/go_obo_downloader.py {go_obo}'.format(bin=self.bin, go_obo=self.go_obo))
        if self.single == 'T' and not os.path.isfile('%s/%s.pathway.xls' % (self.KEGG_dir, self.orgC)):
            super(V10Labelfree, self).end(normal=1, out='NEED TO DO FIRST:\n{bin}/download_kegg_singlespecies.py {orgC} {KEGG_dir}'.format(bin=self.bin, orgC=self.orgC,
                                                                                                    KEGG_dir=self.KEGG_dir))
        if self.single == 'T' and not os.path.isfile('%s/%s.aa.KEGG.fasta' % (self.KEGG_dir, self.orgC)):
            super(V10Labelfree, self).end(normal=1, out=
                'NEED TO DO FIRST:\ncd {KEGG_dir} && {bin}/get_aaseq_from_KEGG_by_org.py {orgC} && cd -'.format(bin=self.bin,
                                                                                                                orgC=self.orgC,
                                                                                                                KEGG_dir=self.KEGG_dir))

    def add_commands(self):
        params = dict(run_wd=self.output)
        params_ = dict(run_wd=self.output, qsub=0)
        if u'/centos7users/' in self.output:
            self.annot_diomand = self.add_command(name='annot_diomand', params=params_)
            self.diff = self.add_command(name='diff', params=params_)
            self.cluster = self.add_command(name='cluster', params=params_)
            self.annot_enrich = self.add_command(name='diff_annot_enrich', params=params_)
            self.go_result_nr = self.add_command(name='go_result_nr', params=params_)
            self.ipath_picture = self.add_command(name='ipath_picture', params=params_)
            if self.run_ppi_:
                self.ppi = self.add_command(name='ppi', params=params_)
        else:
            self.annot_diomand = self.add_command(name='annot_diomand', params=params)
            self.diff = self.add_command(name='diff', params=params)
            self.cluster = self.add_command(name='cluster', params=params)
            self.annot_enrich = self.add_command(name='diff_annot_enrich', params=params)
            self.go_result_nr = self.add_command(name='go_result_nr', params=params)
            self.ipath_picture = self.add_command(name='ipath_picture', params=params)
            if self.run_ppi_:
                self.ppi = self.add_command(name='ppi', params=params)
        self.go_result = self.add_command(name='go_result', params=params)
        self.ipath = self.add_command(name='ipath', params=params)
        self.kegg_result = self.add_command(name='kegg_result', params=params)
        self.cog_result = self.add_command(name='cog_result', params=params)
        self.venn = self.add_command(name='venn', params=params)
        self.qc = self.add_command(name='qc', params=params)
        self.results_collect = self.add_command(name='results_collect', params=params)

    def how_run(self):
        self.after_one(self.annot_diomand, self.run_go_result_nr)
        self.after_one(self.go_result_nr, self.run_go_result)
        self.after_one(self.annot_diomand, self.run_kegg_result)
        self.after_one(self.annot_diomand, self.run_cog_result)
        if self.run_ppi_:
            self.after_some([self.annot_diomand, self.diff], self.run_ppi)
        self.after_one(self.diff, self.run_cluster)
        self.after_one(self.cluster, self.run_venn)
        self.after_some([self.diff, self.annot_enrich], self.run_ipath)
        self.after_some([self.diff, self.kegg_result], self.run_ipath_picture)
        self.after_some([self.go_result, self.kegg_result, self.cog_result, self.diff], self.run_annot_enrich)
        col_list = [self.annot_enrich, self.qc, self.venn, self.ipath, self.ipath_picture]
        if self.run_ppi_:
            col_list.append(self.ppi)
        self.after_some(col_list, self.run_results_collect)
        self.after_one(self.results_collect, self.end)

    def end(self):
        super(V10Labelfree, self).end(out='your program finished')

    def run_annot_diomand(self):
        if self.single == 'F':
            kegg_db = '/mnt/ilustre/users/yitong.feng/scripts/annotion/ko.pep.dmnd'
            if u'/centos7users/' in os.getcwd():
                kegg_db = '/mnt/ilustre/centos7users/jun.yan/fengyitong/database/kegg.pep.dmnd'
        else:
            kegg_db = os.path.join(self.annot_diomand.work_dir, 'kegg_db.fasta')
            try:
                # os.link(os.path.join(self.KEGG_dir, self.orgC + '.aa.KEGG.fasta'), kegg_db)
                #改成cp主要是因为新服务器link不过来。。。
                cp_c = "cp " + os.path.join(self.KEGG_dir, self.orgC + '.aa.KEGG.fasta') + " " + kegg_db
                os.system(cp_c)
            except:
                if not os.path.exists(kegg_db):
                    super(V10Labelfree, self).end(normal=1, out='kegg的blast用的db文件不存在')
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/annotion/run_annotation.py -fasta ' + os.path.join(self.raw_data, 'exp.fasta ')
        cmd += '-go_db ' + self.go_blast_db
        cmd += ' -kegg_db ' + kegg_db
        cmd += ' -string_db ' + self.cog_db
        cmd += ' -out ' + os.path.join(self.annot_diomand.work_dir, 'diomand_results')
        params = dict(
            cmd=cmd,
            node=2,
            memory=6
        )
        self.annot_diomand.set_params(params)
        self.annot_diomand.run()

    def run_go_result(self):
        nrgo = os.path.join(self.go_result.work_dir, 'nr.GO.list')
        try:
            os.link(os.path.join(self.go_result_nr.work_dir, 'nr.GO.list'), nrgo)
            os.link(os.path.join(self.raw_data, 'exp.list'), os.path.join(self.go_result.work_dir, 'exp.list'))
        except:
            if not os.path.exists(nrgo):
                super(V10Labelfree, self).end(normal=1, out='GO报错，nr数据库mapping没有正常完成')

        cmd = '{bin}/get_annot_list.pl exp.list {uniprot2go} pir.GO.list\n'.format(bin=self.bin,
                                                                                   uniprot2go=self.uniprot2go)
        cmd += '{bin}/merge_2tab_file.py nr.GO.list,pir.GO.list GO.list\n'.format(bin=self.bin)
        cmd_ = """${bin}/go_annotation.py -go_obo ${go_obo} -i GO.list -o GO.level.xls 

${bin}/go_level_bar.pl -i GO.level.xls -l 2 -o GO.level2.bar.pdf -w 12 -h 8
${bin}/go_level_pie.py -i GO.level.xls -lv 2 -o GO.level2.pie.pdf
${bin}/go_level_pie_percent.py -i GO.level.xls -lv 2 -o GO.level2.pie_percent.pdf

${bin}/go_level_bar.pl -i GO.level.xls -l 3 -o GO.level3.bar.pdf -w 18 -h 12
${bin}/go_level_pie.py -i GO.level.xls -lv 3 -o GO.level3.pie.pdf
${bin}/go_level_pie_percent.py -i GO.level.xls -lv 3 -o GO.level3.pie_percent.pdf

${bin}/go_level_bar.pl -i GO.level.xls -l 4 -o GO.level4.bar.pdf -w 24 -h 16
${bin}/go_level_pie.py -i GO.level.xls -lv 4 -o GO.level4.pie.pdf
${bin}/go_level_pie_percent.py -i GO.level.xls -lv 4 -o GO.level4.pie_percent.pdf
"""
        cmd_ = Template(cmd_)
        cmd_ = cmd_.render(bin=self.bin,
                           go_obo=self.go_obo,
                           )
        cmd += cmd_
        params = dict(
            cmd=cmd,
            node=4,
            memory=24
        )
        self.go_result.set_params(params)
        self.go_result.run()

    def run_go_result_nr(self):
        vsnr = os.path.join(self.go_result_nr.work_dir, 'exp.fasta_vs_nr.fast.xls')
        try:
            os.link(os.path.join(self.annot_diomand.work_dir, 'diomand_results', 'exp.fasta_vs_nr.fast.xls'), vsnr)
            os.link(os.path.join(self.raw_data, 'exp.list'), os.path.join(self.go_result.work_dir, 'exp.list'))
        except:
            if not os.path.exists(vsnr):
                super(V10Labelfree, self).end(normal=1, out='GO报错，diomand没有正常完成')
        # cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/annotion/get_GO_from_blast_by_nr.py ' + vsnr + ' /mnt/ilustre/users//bing.yang/DB/GO/idmapping.tb ' + self.go_identity + ' nr.GO.list\n'
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/annotion/get_go_from_nr.py'
        cmd += ' -idmapping ' + '/mnt/ilustre/users//bing.yang/DB/GO/idmapping.tb'
        cmd += ' -fasta_vs_nr ' + 'exp.fasta_vs_nr.fast.xls'
        cmd += ' -go_identity ' + self.go_identity
        cmd += ' -out ' + os.path.join(self.go_result.work_dir, 'nr.GO.list') + '\n'
        params = dict(
            cmd=cmd,
            node=3,
            memory=9
        )
        self.go_result_nr.set_params(params)
        self.go_result_nr.run()

    def run_kegg_result(self):
        vskegg = os.path.join(self.kegg_result.work_dir, 'exp.fasta_vs_kegg.xls')
        try:
            os.link(os.path.join(self.annot_diomand.work_dir, 'diomand_results', 'exp.fasta_vs_KEGG.xls'), vskegg)
            os.link(os.path.join(self.raw_data, 'exp.list'), os.path.join(self.kegg_result.work_dir, 'exp.list'))
        except:
            if not os.path.exists(vskegg):
                super(V10Labelfree, self).end(normal=1, out='kegg报错，diomand没有正常完成')
        if self.single == 'F':
            kegg_p = os.path.join(self.KEGG_dir, 'ko.pathway.xls')
        else:
            kegg_p = os.path.join(self.KEGG_dir, self.orgC + '.pathway.xls')
        cmd = r"""${bin}/get_annot_list.pl exp.list ${uniprot2kegg} KEGG.txt
${bin_v}/combine_keggblast_uniprot.py exp.fasta_vs_kegg.xls KEGG.txt ${orgC} pathway.txt ${kegg_identity}
${bin}/merge_2tab_file.py pathway.txt pathway.txt1
mv pathway.txt1 pathway.txt
${bin}/KEGG_all_annotation_col.py -paths ${KEGG_dir} -pathwaytxt pathway.txt -out pathways
awk 'BEGIN{FS=OFS="\t"}{if(NR==1){print "Pathway","Pathway_definition","Number_of_protein","Protein_ko_list","KO2acc"}else{print $1,$2,$3,$4,$7}}' pathways/pathway_table.xls >tmp.pathway_table.xls
awk 'BEGIN{FS=OFS="\t"}{if(NR>1){print $1,$2,$4}}' ${kegg_p} >kegg_layer.txt
${bin}/kegg_brite.pl tmp.pathway_table.xls kegg_layer.txt
${bin}/kegg_pathway_top20bars.pl -t tmp.pathway_table.xls -o pathway
"""
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         bin_v=self.bin_v,
                         uniprot2kegg=self.uniprot2kegg,
                         orgC=self.orgC,
                         KEGG_dir=self.KEGG_dir,
                         kegg_p=kegg_p,
                         kegg_identity=self.kegg_identity
                           )
        params = dict(
            cmd=cmd,
            node=3,
            memory=10
        )
        self.kegg_result.set_params(params)
        self.kegg_result.run()

    def run_cog_result(self):
        vscog = os.path.join(self.cog_result.work_dir, '', 'cog.list.xls')
        try:
            os.link(os.path.join(self.annot_diomand.work_dir, 'diomand_results', 'tmp_out', 'cog.list.xls'), vscog)
            os.link(os.path.join(self.raw_data, 'exp.list'), os.path.join(self.cog_result.work_dir, 'exp.list'))
        except:
            if not os.path.exists(vscog):
                super(V10Labelfree, self).end(normal=1, out='COG报错，diomand没有正常完成')
        if self.cog == 'COG':
            cmd = r"""awk 'BEGIN{FS=OFS="\t"}{if($1&&$2){print $1,$2}}' cog.list.xls>string.COG.list
${bin}/get_annot_list.pl exp.list ${uniprot2cog} pir.COG.list
${bin}/merge_2tab_file.py string.COG.list,pir.COG.list COG.list
${bin}/cog_annot.pl COG.list ${COGclassification} COG.annot.xls 
${bin}/cog_summary.pl -i COG.list
${bin}/cog_bar.pl -i pic_COG
${bin}/KOGinfo.py -type COG -a COG.annot.xls -c COG.class.catalog.xls -o COG.classification.xls
"""
        else:
            cmd = r"""awk 'BEGIN{FS=OFS="\t"}{if($1&&$3){print $1,$3}}' cog.list.xls>string.KOG.list
${bin}/get_annot_list.pl exp.list ${uniprot2kog} pir.KOG.list
{bin}/merge_2tab_file.py string.KOG.list,pir.KOG.list KOG.list
${bin}/kog_annot.pl KOG.list ${KOGclassification} KOG.annot.xls
${bin}/kog_summary.pl -i KOG.list
${bin}/kog_bar.pl -i pic_KOG
${bin}/KOGinfo.py -type KOG -a KOG.annot.xls -c KOG.class.catalog.xls -o KOG.classification.xls
"""
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         uniprot2cog=self.uniprot2cog,
                         uniprot2kog=self.uniprot2kog,
                         COGclassification=self.COGclassification,
                         KOGclassification=self.KOGclassification,
                         )
        params = dict(
            cmd=cmd,
            node=3,
            memory=10
        )
        self.cog_result.set_params(params)
        self.cog_result.run()

    def run_diff(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/diff/diff_labelfree.py'
        cmd += ' -exp ' + os.path.join(self.raw_data, 'exp.txt')
        cmd += ' -group_file ' + os.path.join(self.raw_data, 'group.txt')
        cmd += ' -control_file ' + os.path.join(self.raw_data, 'control.txt')
        cmd += ' -method ' + self.diff_method
        cmd += ' -alternative ' + self.alternative
        cmd += ' -p_correct ' + self.diff_pcorrect_method
        cmd += ' -up ' + self.up
        cmd += ' -down ' + self.down
        cmd += ' -pvalue ' + self.pvalue
        cmd += ' -cutoffs ' + self.cutoffs
        cmd += ' -average ' + self.average
        cmd += ' -out ' + os.path.join(self.diff.work_dir, 'diff_result')
        params = dict(
            cmd = cmd,
            node = 1,
            memory = 1
        )
        self.diff.set_params(params)
        self.diff.run()

    def run_cluster(self):
        tmp_list = os.listdir(os.path.join(self.diff.work_dir, 'diff_result'))
        if not tmp_list:
            super(V10Labelfree, self).end(normal=1, out='diff运行出错，请检查')
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/cluster/cluster_protein_labelfree.py'
        cmd += ' -exp ' + os.path.join(self.raw_data, 'exp.txt')
        cmd += ' -group_file ' + os.path.join(self.raw_data, 'group.txt')
        cmd += ' -diff_path ' + os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse')
        cmd += ' -clt ' + self.clt
        cmd += ' -mr ' + self.mb
        cmd += ' -ml ' + self.ml
        cmd += ' -out ' + os.path.join(self.cluster.work_dir, 'cluster_result')
        params = dict(
            cmd=cmd,
            node=1,
            memory=1
        )
        self.cluster.set_params(params)
        self.cluster.run()

    def run_annot_enrich(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/diff_annot_enrich/diff_annot_enrich.py'
        cmd += ' -go_list ' + os.path.join(self.go_result.work_dir, 'GO.list')
        cmd += ' -go_obo ' + self.go_obo
        cmd += ' -go_level ' + os.path.join(self.go_result.work_dir, 'GO.level.xls')
        cmd += ' -path_txt ' + os.path.join(self.kegg_result.work_dir, 'pathway.txt')
        cmd += ' -kegg_dir ' + self.KEGG_dir
        cmd += ' -pathway_table ' + os.path.join(self.kegg_result.work_dir, 'pathways', 'pathway_table.xls')
        cmd += ' -orgc ' + self.orgC
        cmd += ' -diff_path ' + os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse')
        cmd += ' -out ' + os.path.join(self.annot_enrich.work_dir, 'diff_annot_enrich_result')
        params = dict(
            cmd=cmd,
            node=1,
            memory=1
        )
        self.annot_enrich.set_params(params)
        self.annot_enrich.run()

    def run_venn(self):
        # with open(os.path.join(self.raw_data, 'exp.txt'), 'r') as e_r:
        #     header = e_r.readline().strip().split('\t')
        #     for line in e_r:
        #         if line.strip():
        #             line = line.strip().split('\t')
        #             for s in header[1:]:
        #                 n = header.index(s)
        #                 try:
        #                     e = float(line[n])
        #                     if e > 0:
        #                         with open(os.path.join(self.venn.work_dir, '%s.venn.list'%s), 'a') as fw:
        #                             fw.write(line[0] + '\n')
        #                 except:
        #                     pass
        venn_lists = glob.glob(os.path.join(self.diff.work_dir,  'diff_result_for_other_analyse', 'venn_pre', '*.venn.list'))
        if not venn_lists:
            super(V10Labelfree, self).end(normal=1, out='venn报错，diff没有正常生成venn前体文件')
        for venn in venn_lists:
            try:
                os.link(venn, os.path.join(self.venn.work_dir, os.path.basename(venn)))
            except:
                pass
        sconfig = os.path.join(self.venn.work_dir, 'sample.config')
        try:
            os.link(os.path.join(self.cluster.work_dir, 'cluster_result', 'group'), sconfig)
        except:
            pass
        if not os.path.exists(sconfig):
            super(V10Labelfree, self).end(normal=1, out='venn报错，cluster没有正常完成')
        cmd = '''${bin}/labelfree_venn.py ${sconfig}
'''
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         sconfig=sconfig
                         )
        params = dict(
            cmd=cmd,
            node=3,
            memory=6
        )
        self.venn.set_params(params)
        self.venn.run()

    def run_qc(self):
        cmd = '''cp ${cwd}/psm.xlsx ./
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
${bin}/Protein_Coverage_distribution.pl -i protein.xls -s "Accession,Coverage [%]"
${bin}/Protein_information.pl -i Protein_information.xls
'''
        cmd += r'''# cp ${cwd}/proteins.csv ./
# cp ${cwd}/protein-peptides.csv ./
# cp ${cwd}/Protein_information.xls ./
# 
# sed 's/,/\t/g' proteins.csv > proteins.xls.tmp
# sed 's/,/\t/g' protein-peptides.csv > protein-peptides.xls
# awk -F '\t' -vOFS='\t' '{if (NR==1){print $0"\tKDA"} else {print $0"\t"$(NF-1)/1000}}' proteins.xls.tmp >proteins.xls
# sed -i 's/\tCoverage (%)\t/\tCoverage\t/g' proteins.xls
# rm proteins.xls.tmp
# ${bin}/dMass_plot.pl -i DB_search_psm.xls -s "m/z,ppm"
# /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/pre_pipeline/bin/Protein_molecular_weight_distribution.pl -i proteins.xls -s "Accession,KDA" 
# ${bin}/Peptide_length_distribution.pl -i protein-peptides.xls -s "Peptide"
# ${bin}/Peptide_number_distribution.pl -i proteins.xls -s "Accession,#Peptides"
# ${bin}/Protein_information.pl -i Protein_information.xls
# ${bin}/Protein_Coverage_distribution.pl -i proteins.xls -s 'Accession,Coverage'
'''
        cmd = Template(cmd)
        cmd = cmd.render(bin=self.bin,
                         cwd=self.raw_data,
                         )
        params = dict(
            cmd=cmd,
            node=2,
            memory=6
        )
        self.qc.set_params(params)
        self.qc.run()

    def run_ppi(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/ppi/ppi_by_blast.py '
        cmd += ' -db_file ' + self.ppi_db
        cmd += ' -blast_file ' + os.path.join(self.annot_diomand.work_dir, 'diomand_results', 'exp.fasta_vs_string.xml')
        cmd += ' -diff_path ' + os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse')
        cmd += ' -identity ' + self.ppi_identity
        cmd += ' -out ' + os.path.join(self.ppi.work_dir, 'ppi_results')
        params = dict(
            cmd=cmd,
            node=1,
            memory=1
        )
        self.ppi.set_params(params)
        self.ppi.run()

    def run_ipath(self):
        pathways = glob.glob(os.path.join(self.annot_enrich.work_dir, 'diff_annot_enrich_result', '*pathway.txt'))
        if not pathways:
            super(V10Labelfree, self).end(normal=1, out='ipath报错，diff_annot_enrich没有正常完成')
        for file in pathways:
            if not u'.DE.' in file:
                try:
                    os.link(file, os.path.join(self.ipath.work_dir, os.path.basename(file)))
                except:
                    pass
        cmd = r'''
for i in `ls *.up.pathway.txt |awk -F'.' '{print $1}'`
do
    less $i.down.pathway.txt |awk 'BEGIN{FS=OFS="\t"}{print $2,"DOWN"}'> $i.tmp_d.list
    less $i.up.pathway.txt |awk 'BEGIN{FS=OFS="\t"}{print $2,"UP"}' > $i.tmp_u.list 
    /mnt/ilustre/users/bingxu.liu/workspace/tabletools_add.pl -i $i.tmp_d.list -t $i.tmp_u.list -n 1 -headi F -headt F |grep -v DOWN |awk 'BEGIN{FS=OFS="\t"}{print $1,"#FF0000","W10"}' > $i.up.Ipath 
    /mnt/ilustre/users/bingxu.liu/workspace/tabletools_add.pl -i $i.tmp_u.list -t $i.tmp_d.list -n 1 -headi F -headt F |grep -v UP |awk 'BEGIN{FS=OFS="\t"}{print $1,"#00FF00","W10"}' > $i.down.Ipath
    /mnt/ilustre/users/bingxu.liu/workspace/tabletools_add.pl -i $i.tmp_u.list -t $i.tmp_d.list -n 1 -headi F -headt F |grep UP |awk 'BEGIN{FS=OFS="\t"}{print $1,"#0000FF","W10"}' > $i.up_down.Ipath
    cat $i.up.Ipath $i.down.Ipath $i.up_down.Ipath > $i.DE.Ipath
done
'''
        params = dict(
            cmd=cmd,
            node=1,
            memory=1
        )
        self.ipath.set_params(params)
        self.ipath.run()

    def run_ipath_picture(self):
        pathway_table = os.path.join(self.kegg_result.work_dir, 'pathways', 'pathway_table.xls')
        diff_result = os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse')
        out = os.path.join(self.ipath_picture.work_dir, 'ipath_pictures')
        cmd = 'python /mnt/ilustre/centos7users/yitong.feng/script/ipath/ipath_protein_pipeline.py -diff_path %s -pathway_table %s -out %s' %(diff_result, pathway_table, out)
        params = dict(
            cmd=cmd,
            node=1,
            memory=1
        )
        self.ipath_picture.set_params(params)
        self.ipath_picture.run()

    def run_results_collect(self):
        des_file = os.path.join(self.raw_data, 'description.txt')
        if not os.path.exists(des_file):
            acc_lists = list()
            with open(os.path.join(self.annot_diomand.work_dir, 'diomand_results', 'exp.fasta_vs_nr.fast.xls'), 'r') as \
                fr, open(des_file, 'w') as fw:
                fw.write('Accession\tDescription\n')
                for line in fr:
                    if not line.strip():
                        continue
                    line = line.strip().split('\t')
                    acc, des = line[0], line[-1].split('>')[0]
                    if acc not in acc_lists:
                        fw.write(acc+'\t'+des+'\n')
                        acc_lists.append(acc)

        if self.single != 'T':
            MJ = 'MJ_PM_LABELFREE_autoreporter.py'
        else:
            MJ = 'MJ_PM_LABELFREE_nosingle_autoreporter.py'
        cmd = Template(filename='/mnt/ilustre/users/yitong.feng/scripts/labelfree_v10/result_collecter_labelfree.bash')
        cmd = cmd.render(outdir=self.output,
                         cog=self.cog,
                         rawdata=self.raw_data,
                         qc_dir=self.qc.work_dir,
                         go_result=self.go_result.work_dir,
                         kegg_result=self.kegg_result.work_dir,
                         cog_result=self.cog_result.work_dir,
                         diff_result=os.path.join(self.diff.work_dir, 'diff_result'),
                         diff_result_=os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse'),
                         ipath_result=self.ipath.work_dir,
                         ipath_picture=os.path.join(self.ipath_picture.work_dir, 'ipath_pictures'),
                         venn_result=self.venn.work_dir,
                         cluster_result=os.path.join(self.cluster.work_dir, 'cluster_result'),
                         annot_enrich_result=os.path.join(self.annot_enrich.work_dir, 'diff_annot_enrich_result'),
                         bin=self.bin,
                         up=self.up,
                         dup=self.dup,
                         MJ=MJ
                         )
        if self.run_ppi_:
            ppi_cmd = "mkdir -p 3.DiffExpAnalysis/3.7.ppi\n cp %s/* 3.DiffExpAnalysis/3.7.ppi" % os.path.join(self.ppi.work_dir, 'ppi_results')
            cmd = cmd.split('\n\n\n\n')
            cmd[1] += ppi_cmd
            cmd = '\n'.join(cmd)
        params = dict(
            cmd=cmd,
            node=3,
            memory=12
        )

        self.results_collect.set_params(params)
        self.results_collect.run()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_annot_diomand()
        self.run_diff()
        self.run_qc()
        super(V10Labelfree, self).fire()


if __name__ == '__main__':
    try:
        config = sys.argv[1]
    except:
        config = 'labelfree.ini'
    config = os.path.abspath(config)
    itraq = V10Labelfree(config)
    itraq.run()