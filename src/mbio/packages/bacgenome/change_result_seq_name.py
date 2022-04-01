__author__ = 'guanqing.zou'
import os
import pandas as pd


class ChangeResultSeqName():
    def __init__(self,map_dic):   #{"abc":'seq1'}
        self.map_dic = map_dic

    def change_fasta(self,infile):
        with open(infile) as f, open(infile+'_new','w') as fw:
            for line in f:
                if line.startswith('>'):
                    line=line.strip()
                    spline = line.split()
                    s_name = spline[0][1:]
                    if s_name in self.map_dic.keys():
                        spline[0] = '>'+self.map_dic[s_name]
                    fw.write(' '.join(spline)+'\n')
                else:
                    fw.write(line)

        os.remove(infile)
        os.rename(infile+'_new',infile)

    def change_fasta_spe(self,infile):
        with open(infile) as f, open(infile+'_new','w') as fw:
            for line in f:
                if line.startswith('>'):
                    line=line.strip()
                    spline = line.split()     #  >gene1 1882 3024 seq1_ORF_gene1 1882 3024
                    if len(spline) > 3:
                        c3 = spline[3].split("_")
                        if c3[0] in self.map_dic.keys():
                            c3[0] = self.map_dic[c3[0]]
                            spline[3] = '_'.join(c3)
                            fw.write(' '.join(spline)+'\n')
                        else:
                            fw.write(line+'\n')
                    else:
                        fw.write(line+'\n')

                else:
                    fw.write(line)

        os.remove(infile)
        os.rename(infile+'_new',infile)

    def change_csv(self,infile,head_name):
        data = pd.read_table(infile,sep="\t",header=0)
        ret = []
        for d in data[head_name]:
            if d in self.map_dic.keys():
                ret.append(self.map_dic[d])
            else:
                ret.append(d)
        data[head_name] = ret
        fw = open(infile,'w')
        data.to_csv(fw,sep='\t',index=False)

    def change_csv_split(self,infile,head_name):
        data = pd.read_table(infile,sep="\t",header=0)
        ret = []
        for d in data[head_name]:
            sd = d.split("_")
            sd1= sd[0]
            if sd1 in self.map_dic.keys():
                sd[0]=self.map_dic[sd1]
                ret.append('_'.join(sd))
            else:
                ret.append(d)
        data[head_name] = ret
        fw = open(infile,'w')
        data.to_csv(fw,sep='\t',index=False)


    def change_gff(self,infile):
        with open(infile) as f, open(infile+'_new','w') as fw:
            for line in f:
                if line.startswith('#'):
                    fw.write(line)
                else:
                    spline = line.split('\t',1)
                    if len(spline) < 2:
                        fw.write(line)
                    else:
                        if spline[0] in self.map_dic.keys():
                            spline[0] = self.map_dic[spline[0]]
                        fw.write('\t'.join(spline))
        os.remove(infile)
        os.rename(infile+'_new',infile)


    def change_ptt(self,infile):
        with open(infile) as f, open(infile+'_new','w') as fw:
            for line in f:
                line = line.strip()
                spline = line.split()
                if len(spline) ==1:
                    if spline[0] in self.map_dic.keys():
                        fw.write(self.map_dic[spline[0]]+'\n')
                        continue
                fw.write(line+'\n')

        os.remove(infile)
        os.rename(infile+'_new',infile)



    def change_pip(self,files_types):
        for e in files_types:
            if isinstance(e,list) :
                if len(e) > 1:
                    file = e[0]
                    type = e[1]
                else:
                    raise Exception("not >1")
                if type == 'csv' or type == "sp_column":
                    if  len(e) > 2:
                        column = e[2]
                    else:
                        raise ("csv must has column value")
            else:
                raise Exception('not list')

            if not os.path.exists(file):
                #print file
                continue

            #print 'deal: ' + file

            if type == 'fasta':
                self.change_fasta(file)
            elif type == 'gff':
                self.change_gff(file)
            elif type == 'ptt':
                self.change_ptt(file)
            elif type == 'csv':
                self.change_csv(file,column)
            elif type == 'sp_column':
                self.change_csv_split(file,column)
            elif type == "spe_fa":
                self.change_fasta_spe(file)


if __name__ == "__main__":
    # test = 0
    # if test == 1 :
    #     with open('test.txt','w') as fw:
    #         fw.write("a\tb\na1\tb1\na2\tb2\n")
    #     map_dic = {'b1':"c1"}
    #     C= ChangeResultSeqName(map_dic)
    #     C.change_csv('test.txt','b')
    def this_project_need_change_files(sample,pre_path):
        files = [
            ("%s.cds.ffn","spe_fa"), #spe_fa
            ("%s_summary.xls", "csv","Location"),
            #gene
            ("gene/%s.gene.gff","sp_column","Sequence id"), #sp_column
            ("gene/%s.trna.gff","sp_column","Sequence id"), #sp_column
            ("gene/%s.rrna.gff","sp_column","Sequence id"), #sp_column
            ("gene/%s.cds.ffn","spe_fa"),  #spe_fa
            ("gene/%s.cds.faa","spe_fa"),
            #annotation
            ("annotation/CAZy/%s_anno_cazy.xls","csv","Location"),
            ("annotation/COG/%s_cog_anno.xls","csv","Location"),
            ("annotation/GO/%s_go_anno.xls","csv","Location"),
            ("annotation/KEGG/%s_kegg_anno.xls","csv","Location"),
            ("annotation/NR/%s_anno_nr.xls","csv","Location"),
            ("annotation/NR/%s_whole_genome_anno_nr.xls","csv","Location"),
            ("annotation/Pfam/%s_anno_pfam.xls","csv","Location"),
            ("annotation/Pfam/%s_whole_genome_anno_pfam.xls","csv","Location"),
            ("annotation/Summary/%s_anno_summary.xls","csv","Location"),
            ("annotation/Swissprot/%s_anno_swissprot.xls","csv","Location"),
            ("annotation/Swissprot/%s_whole_genome_anno_swissprot.xls","csv","Location"),
            ("annotation/Two_component/%s.senser_regulator.xls","csv","Location"),
            #cgview
            #circos
            #gbk/analysis_gbk/
            #gbk/gbk/
            #gbk/
            ##("gbk/%s_seq1.gff","gff"),
            ##("gbk/%s_seq1.ptt","ptt"),
            # metabolic_system
            ("metabolic_system/antiSMASH/%s_gene_antismash.xls","csv","Location"),
            ("metabolic_system/antiSMASH/%s_antismash_anno.xls","sp_column","Cluster ID"),  #sp_column
            #mobile_elements/
            ("mobile_elements/CRISPR_Cas/%s_CRISPR_Cas_summary.xls","csv","Sample_chr"),
            ("mobile_elements/CRISPR_Cas/%s_CRISPR_Cas_detail.xls","csv","Sample_chr"),
            ("mobile_elements/Genomic_Islands/%s_GI_summary.xls","csv","Location"),
            ("mobile_elements/Genomic_Islands/%s_GI_detail.xls","csv","Location"),
            ("mobile_elements/prephage/%s_prephage_detail.xls","csv","Location"),
            ("mobile_elements/prephage/%s_prephage_summary.xls","csv","Location"),
            #pathogenic_system
            ("pathogenic_system/CARD/%s_card_anno.xls","csv","Location"),
            ("pathogenic_system/PHI/%s_phi_anno.xls","csv","Location"),
            ("pathogenic_system/SIGNALP/%s_Gram-_SignalP.txt","csv","Location"),
            ("pathogenic_system/SIGNALP/%s_Gram+_SignalP.txt","csv","Location"),
            ("pathogenic_system/TCDB/%s_whole_genome_tcdb_anno.xls","csv","Location"),
            ("pathogenic_system/TCDB/%s_tcdb_anno.xls","csv","Location"),
            ("pathogenic_system/TMHMM/%s_whole_genome_tmhmm_anno.xls","csv","Location"),
            ("pathogenic_system/TMHMM/%s_tmhmm_anno.xls","csv","Location"),
            ("pathogenic_system/VFDB/%s_vfdb_anno.xls","csv","Location"),
            ("structral_genome/promoter_predict/%s_promoter_result.xls","csv","Location")
            #tree
            ]

        new_files=[]
        for f in files:
            lf = list(f)
            lf[0]=pre_path + '/'+lf[0]%sample
            new_files.append(lf)
        return new_files

    #zouguanqing
    def change_back_result_files_name(pre,map,sample):
        pre_path =  pre+ '/' +sample
        deal_files = this_project_need_change_files(sample, pre_path)
        CR = ChangeResultSeqName(map)
        CR.change_pip(deal_files)

    map = {'AP012303.1':"seq1"}
    new_map = {v:k for k, v in map.items()}
    pre = '/mnt/ilustre/users/sanger-dev/workspace/20190706/Bacgenome_tsg_34741/output'
    change_back_result_files_name(pre,new_map,'GCA-000283755')


