# coding=utf-8
import sys
import argparse
import subprocess
import shlex
import pandas as pd
import json
import re
import os
__author__ = "gdq"
"""
Proteins are generally comprised of one or more functional regions, commonly termed domains. 
The presence of different domains in varying combinations in different proteins gives rise to 
the diverse repertoire of proteins found in nature. 
Identifying the domains present in a protein can provide insights into the function of that protein.
The Pfam database is a large collection of protein domain families. 
Each family is represented by multiple sequence alignments and a hidden Markov model (HMMs). 
Family:
    A collection of related protein regions
Domain:
    A structural unit
"""


def run_hmmscan(args):
    other_args = "".join(args.other_args)
    # format cmd and run cmd
    cmd = "{} -o {} --tblout {} --domtblout {} --seed 0 --notextw ".format(
        args.hmmscan, args.o, args.tblout, args.domtblout
    )
    cmd += "--cpu {} ".format(args.cpu)
    if "-T " not in args.other_args:
        cmd += "-E {} ".format(args.E)
        cmd += "--domE {} ".format(args.domE)
    if args.other_args:
        cmd += other_args
    cmd += " {hmm_db} {fasta}".format(hmm_db=args.hmmdb, fasta=args.seqfile)
    print(cmd)
    subprocess.check_call(shlex.split(cmd))
    # reformat domtblout
    output = "domain_predict.txt"
    with open(args.domtblout) as fr, open(output, 'w') as fw:
        for line in fr:
            if line.startswith("# target name"):
                col_names = [
                    "family",
                    "pfam_id",
                    "domain_len",
                    "query_id",
                    "query_accession",
                    "query_len",
                    "e_value",
                    "score",
                    "bias",
                    "domain_order",
                    "domain_num",
                    "c_evalue",
                    "i_evalue",
                    "i_score",
                    "i_bias",
                    "hmm_start",
                    "hmm_end",
                    "ali_start",
                    "ali_end",
                    "env_start",
                    "env_end",
                    "acc",
                    "description",
                ]
                fw.write("\t".join(col_names) + '\n')
            elif line.startswith("#") or line.startswith("Query file"):
                continue
            else:
                fw.write(re.sub("\s+", '\t', line, 22))
    # 处理转录本名AT4G22890.5被转换为AT4G22890X5的情况。
    # 后发现该情况实际是由于split_fasta导致的，split过程中把.转换为大写的X了。因此注释掉下面的代码
    # tmp_pd = pd.read_table(output, header=0)
    # query_id_list = tmp_pd['query_id']
    # with open(args.seqfile) as f:
    #     origin_query_id = list()
    #     for line in f:
    #         if line.startswith('>'):
    #             origin_query_id.append(line[1:].split()[0])
    # changed_query_id = [x.replace('.', 'X') for x in origin_query_id]
    # convert_dict = dict(zip(changed_query_id, origin_query_id))
    # back_query_id = [convert_dict[x] for x in query_id_list]
    # tmp_pd['query_id'] = back_query_id
    # tmp_pd.to_csv(output, sep='\t', header=True, index=False)
    return output


def get_domain_stat(domain_predict):
    domain_predict_pd = pd.read_table(domain_predict, sep='\t', header=0)
    domain2query = domain_predict_pd.loc[:, ['pfam_id', 'query_id']]
    domain2query.columns = ['pfam_id', 'query_id']
    domain_stat_dict = dict()
    group_domain = domain2query.groupby("query_id").groups
    for key in group_domain:
        tmp_pd = domain2query.iloc[group_domain[key], :]
        tmp_count = tmp_pd.groupby("pfam_id").count().to_dict()["query_id"]
        new_tmp_count = dict()
        for k, v in tmp_count.items():
            new_tmp_count[k.split(".")[0]] = v
        domain_stat_dict[key] = new_tmp_count
    return domain_stat_dict


def judge_plant_tf(domain_stat_dict, rule_dict_list):
    """
    :param domain_stat_dict:
    :param rule_dict_list:
    :return: {query_id:
                    {
                        family_name1: [d1,d2],
                        family_name2: [d2,d3],
                    },
            ......
            }
    对于没有参与判定转录因子家族的domain将不予记录
    对于没有被判定为转录因子的query将不予记录
    """
    query2family = dict()
    for query_id, domain_dict in domain_stat_dict.items():
        predict_domains = set(domain_dict.keys())
        query2family[query_id] = dict()
        for each_tf in rule_dict_list:
            judge_rule = each_tf["judge_rule"]
            binding_domains = set(judge_rule['binding'].keys())
            if not binding_domains:
                continue
            auxiliary_domains = set(judge_rule["auxiliary"])
            forbidden_domains = set(judge_rule["forbidden"])
            if len(predict_domains & forbidden_domains) != 0:
                continue
            if len(binding_domains - predict_domains) != 0:
                continue
            range_match = [1 for x in binding_domains if domain_dict[x] in judge_rule['binding'][x]]
            if sum(range_match) == len(binding_domains):
                if len(auxiliary_domains & predict_domains) >= 1:
                    if each_tf['Family'] in query2family[query_id]:
                        query2family[query_id].pop(each_tf['Family'])
                    based_domains = binding_domains | (auxiliary_domains & predict_domains)
                    query2family[query_id][each_tf['SubFamily']] = list(based_domains)
                else:
                    if not auxiliary_domains:
                        query2family[query_id][each_tf['SubFamily']] = list(binding_domains)
                    else:
                        continue
                        # 此处不满足条件
                        # if each_tf['SubFamily'] not in query2family[query_id]:
                        #     query2family[query_id][each_tf['Family']] = list(binding_domains)
            else:
                continue
        else:
            if not query2family[query_id]:
                _ = "{} was not judged as any kind of TF".format(query_id)
                query2family.pop(query_id)

    return query2family


def judge_animal_tf(domain_stat_dict, rule_dict):
    query2family = dict()
    for query_id, domain_dict in domain_stat_dict.items():
        query2family[query_id] = dict()
        for domain in domain_dict:
            family = rule_dict.get(domain)
            if family:
                query2family[query_id][family] = [domain, ]
        else:
            if not query2family[query_id]:
                _ = "{} was not judged as any kind of TF".format(query_id)
                query2family.pop(query_id)
    else:
        return query2family


def get_predicted_tf_fasta(in_file, target_ids, out_file, exclude_mode=False):
    with open(in_file) as f1, open(out_file, 'w') as f2:
        for line in f1:
            if line.startswith('#'):
                continue
            if line.startswith('>'):
                id_ = line.lstrip('>').split()[0]
                if exclude_mode:
                    if (id_ not in target_ids):
                        f2.write(line)
                else:
                    if (id_ in target_ids):
                        f2.write(line)
            else:
                if exclude_mode:
                    if (id_ not in target_ids):
                        f2.write(line)
                else:
                    if (id_ in target_ids):
                        f2.write(line)


def parse_plant_tf_judge_rules(raw_rules):
    """
    最后依据judge_rule_dict的转录因子判定规则:
        1. key包含的所有pfam_id都要被蛋白包含，且每个pfam_id的个数在指定的列表范围，满足此条件后才可进入2，3判读
        2. 对于auxiliary列表包含的pfam_id，如某蛋白含有其中任何一个pfam_id且满足1条件，判定具体到subfamily，否则为母类。
        3. 对于forbidden列表包含的pfam_id， 如某蛋白含有其中任何一个pfam_id，则不能认定是转录因子，即使1，2都符合。
    :param raw_rules:
    :return:
    """
    rule_txt = raw_rules
    rule_list = rule_txt.strip().split('\n')
    header = rule_list[0].split('\t')
    domain_num_limit = 10
    rule_dict_list = list()
    # Family	SubFamily	DNA-binding domain	Auxiliary domain	Forbidden domain
    for each in rule_list[1:]:
        each = each.replace("self-build", "PF0000x")  # 方便后续正则匹配，得出正确的判定规则
        tmp_dict = dict(zip(header, each.split('\t')))
        # sub_family = tmp_dict['SubFamily']
        binding_domain = tmp_dict["DNA-binding domain"]
        auxiliary_domain = tmp_dict['Auxiliary domain']
        forbidden_domain = tmp_dict['Forbidden domain']
        judge_rule_dict = {
            "binding": dict(),
            "auxiliary": list(),
            "forbidden": list(),
        }
        # 这里针对binding_domain 添加判定规则, 只考虑'>'和'<', 只考虑and关系，不考虑or关系，因为目前不存在这种状况
        match_result = re.findall(r"[^()]+?\s+\(([<>=]*)(\d+)\)\s+\((P[^()]+?)\)", binding_domain)
        if match_result:
            for sign, num, pf_id in match_result:
                if not sign:
                    judge_rule_dict['binding'][pf_id] = [int(num)]
                elif ">" in sign:
                    judge_rule_dict['binding'][pf_id] = range(int(num), domain_num_limit)
                elif "<" in sign:
                    judge_rule_dict['binding'][pf_id] = range(1, int(num)+1)
                else:
                    pass
        match_result = re.findall(r"[^()]+?\s+\((P[^()]+?)\)", binding_domain)
        if match_result:
            for pf_id in match_result:
                judge_rule_dict['binding'][pf_id] = range(1, domain_num_limit)
        #  这里针对auxiliary_domain, 只考虑or关系， 不考虑domain数量的要求情况，因为目前其他情况不存在
        match_result = re.findall(r"[^()]+?\s+\((P[^()]+?)\)", auxiliary_domain)
        if match_result:
            for pf_id in match_result:
                judge_rule_dict['auxiliary'].append(pf_id)
        #  这里针对forbidden_domain, 只考虑or关系， 不考虑domain数量的要求情况，因为目前其他情况不存在
        match_result = re.findall(r"[^()]+?\s+\((P[^()]+?)\)", forbidden_domain)
        if match_result:
            for pf_id in match_result:
                judge_rule_dict['forbidden'].append(pf_id)
        # 如果只有binding domain，那么就认为family和subfamily要保持一致
        if not judge_rule_dict['auxiliary'] and not judge_rule_dict['forbidden']:
            if not tmp_dict['SubFamily']: # 增加亚家族的判断，修正家族判断错误的问题
                tmp_dict['SubFamily'] = tmp_dict['Family']
        # save rule dict
        tmp_dict['judge_rule'] = judge_rule_dict
        rule_dict_list.append(tmp_dict)
    else:
        return rule_dict_list


def parse_animal_tf_judge_rules(raw_rules):
    rule_list = raw_rules.strip().split("\n")
    header = rule_list[0]
    family2pf_id = {x.split('\t')[2]: x.split('\t')[0] for x in rule_list}
    return family2pf_id


def run_blast(target, query, top=5, out="diamond.out.txt", out_format=6,
              sensitive="sensitive", p=12, evalue=1e-3, blast_type='blastp',
              diamond="/mnt/ilustre/users/sanger-dev/app/bioinfo/align/diamond-0.8.35/diamond"):
    """
    比对的目的是为了能够分配一个已知的TF_id，后续靶基因预测可以根据这个已知的TF_id获得转录因子靶向Motif.
    :param target: file for build database file
    :param query: input query file
    :param top: report alignments within this percentage range of top alignment score (overrides --max-target-seqs)
    :param out: output file
    :param format: output format
    :param sensitive:  sensitive, (default: fast), more-sensitive
    :param p: number of CPU threads
    :param evalue: maximum e-value to report alignments
    :param blast_type: blastp or blastx
    :param diamond: where is diamond
    :return:
    """
    # make db cmd
    if not target.endswith(".dmnd"):
        build_db_cmd = "{diamond} ".format(diamond=diamond)
        build_db_cmd += "makedb --in {} ".format(target)
        build_db_cmd += "-d {} ".format('seqdb')
        print(build_db_cmd)
        subprocess.check_call(shlex.split(build_db_cmd))
        target = 'seqdb'
    # blast cmd
    blast_cmd = "{diamond} ".format(diamond=diamond)
    blast_cmd += "{blast_type} ".format(blast_type=blast_type)
    blast_cmd += "-d {db} ".format(db=target)
    blast_cmd += "-q {query} ".format(query=query)
    blast_cmd += "-o {out} ".format(out=out)
    blast_cmd += "-k {hit_num} ".format(hit_num=top)
    blast_cmd += "-p {threads} ".format(threads=p)
    blast_cmd += "-f {format} ".format(format=out_format)
    blast_cmd += "-e {evalue} ".format(evalue=evalue)
    blast_cmd += "--{} ".format(sensitive)
    # run cmd
    print(blast_cmd)
    subprocess.check_call(shlex.split(blast_cmd))


def parse_args():
    parser = argparse.ArgumentParser(description="""
        A simple Wrapper for hmmscan, and its result will be used to predict TFs.
        hmmscan is used to search protein sequences against collections of pro-
        tein  profiles. For each sequence in <seqfile>, use that query sequence
        to search the target database of profiles in <hmmdb>, and output ranked
        lists  of  the  profiles  with  the  most  significant  matches  to the
        sequence.""", )
    parser.add_argument('-s', metavar="species", default="plant", help="plant or animal")
    # args for selecting peps to build database for blast
    parser.add_argument('-organism', default="unknown", help="organism name")
    parser.add_argument('-blast_all', default='yes', help="if 'all', blast without organism specified")
    parser.add_argument('-hmmscan', help="Path of hmmscan", default="/mnt/ilustre/users/sanger-dev/sg-users/litangjian/hmm/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan")
    parser.add_argument("-hmmdb", default="/mnt/ilustre/users/sanger-dev/app/database/pfam_31/Pfam-A.hmm", help="""
        The <hmmdb> needs to be  press'ed  using  hmmpress  before  it  can  be
        searched  with  hmmscan.   This  creates  four  binary  files, suffixed
        .h3{fimp}.
        """)
    parser.add_argument("-seqfile", required=True, help="""
         The <seqfile> may contain more than one query sequence. It  can  be  in
           FASTA  format,  or several other common sequence file formats (genbank,
           embl, and uniprot, among others), or in alignment file formats  (stock-
           holm,  aligned  fasta, and others). See the --qformat option for a com-
           plete list.
        """)
    parser.add_argument("-o", default="hmmscan_raw_result.txt", help="""
        Direct the main human-readable output to a file <f>  instead  of
        the default stdout.
        """)
    parser.add_argument("-tblout", default="tblout.txt", help="""
        Save  a  simple  tabular  (space-delimited) file summarizing the
        per-target output, with one  data  line  per  homologous  target
        model found.
        """)
    parser.add_argument("-domtblout", default="domtblout.tmp.txt", help="""
        Save  a  simple  tabular  (space-delimited) file summarizing the
        per-domain output, with one  data  line  per  homologous  domain
        detected in a query sequence for each homologous model.
        """)
    parser.add_argument("-E", default=0.001, type=float, help="""
        In the per-target output, report target profiles with an E-value
        of  <= <x>.  The default is 10.0, meaning that on average, about
        10 false positives will be reported per query, so  you  can  see
        the  top  of  the  noise  and decide for yourself if it's really
        noise.
        """)
    parser.add_argument("-domE", default=0.0001, type=float, help="""
        In  the per-domain output, for target profiles that have already
        satisfied the per-profile reporting threshold, report individual
        domains  with  a  conditional E-value of <= <x>.  The default is
        10.0.  A conditional E-value means the expected number of  addi-
        tional  false  positive  domains  in the smaller search space of
        those comparisons that already satisfied the per-profile report-
        ing threshold (and thus must have at least one homologous domain
        already).
        """)
    parser.add_argument("-cpu", default=12, type=int, help="""
        Set  the  number of parallel worker threads to <n>.  By default,
        HMMER sets this to the number of CPU cores it  detects  in  your
        machine  -  that is, it tries to maximize the use of your avail-
        able processor cores. Setting <n>  higher  than  the  number  of
        available  cores  is of little if any value, but you may want to
        set it to something less. You can also control  this  number  by
        setting an environment variable, HMMER_NCPU.
          This  option  is only available if HMMER was compiled with POSIX
        threads support. This is the  default,  but  it  may  have  been
        turned off for your site or machine for some reason.
        """)
    parser.add_argument("-other_args", nargs="*", default='', help="""
    All other optional arguments of hmmscan are supported. But, you have to quote
    """)
    # arg for diamond blast
    parser.add_argument('-tfdb', default="/mnt/ilustre/users/sanger-dev/app/database/TFDB/")
    parser.add_argument('-diamond', default="/mnt/ilustre/users/sanger-dev/app/bioinfo/align/diamond-0.8.35/diamond")
    parser.add_argument('-evalue', default=0.0001, help="diamond evalue")
    return parser.parse_args()


# tf_factor_rule
plant_tf_family_assignment_rules = """
Family	SubFamily	DNA-binding domain	Auxiliary domain	Forbidden domain
AP2/ERF	AP2	AP2 (>=2) (PF00847)	-	-
AP2/ERF	ERF	AP2 (1) (PF00847)	-	-
AP2/ERF	RAV	AP2 (PF00847) and B3 (PF02362)	-	-
B3 superfamily	ARF	B3 (PF02362)	Auxin_resp (PF06507)	-
B3 superfamily	B3	B3 (PF02362)	-	-
BBR-BPC	BBR-BPC	GAGA_bind (PF06217)	-	-
BES1	BES1	DUF822 (PF05687)	-	-
bHLH	bHLH	HLH (PF00010)	-	-
bZIP	bZIP	bZIP_1 (PF00170)	-	-
C2C2	CO-like	zf-B_box (PF00643)	CCT (PF06203)	-
C2C2	Dof	Zf-Dof (PF02701)	-	-
C2C2	GATA	GATA-zf (PF00320)	-	-
C2C2	LSD	Zf-LSD1 (PF06943)	-	Peptidase_C14 (PF00656)
C2C2	YABBY	YABBY (PF04690)	-	-
C2H2	C2H2	zf-C2H2 (PF00096)	-	RNase_T (PF00929)
C3H	C3H	Zf-CCCH (PF00642)	-	RRM_1 (PF00076) or Helicase_C (PF00271)
CAMTA	CAMTA	CG1 (PF03859)	-	-
CPP	CPP	TCR (PF03638)	-	-
DBB	DBB	zf-B_box (>=2) (PF00643)	-	-
E2F/DP	E2F/DP	E2F_TDP (PF02319)	-	-
EIL	EIL	EIN3 (PF04873)	-	-
FAR1	FAR1	FAR1 (PF03101)	-	-
GARP	ARR-B	G2-like (self-build)	Response_reg (PF00072)	-
GARP	G2-like	G2-like (self-build)	-	-
GeBP	GeBP	DUF573 (PF04504)	-	-
GRAS	GRAS	GRAS (PF03514)	-	-
GRF	GRF	WRC (PF08879)	QLQ (PF08880)	-
HB	HD-ZIP	Homeobox (PF00046)	HD-ZIP_I/II (self-build) or SMART (PF01852)	-
HB	TALE	Homeobox (PF00046)	BELL (self-build)or ELK (PF03789)	-
HB	WOX	homeobox (PF00046)	Wus type homeobox (self-build)	-
HB	HB-PHD	homeobox (PF00046)	PHD (PF00628)	-
HB	HB-other	homeobox (PF00046)	-	-
HRT-like	HRT-like	HRT-like (self-build)	-	-
HSF	HSF	HSF_dna_bind (PF00447)	-	-
LBD (AS2/LOB)	LBD (AS2/LOB)	DUF260 (PF03195)	-	-
LFY	LFY	FLO_LFY (PF01698)	-	-
MADS	M_type	SRF-TF (PF00319)	-	-
MADS	MIKC	SRF-TF (PF00319)	K-box (PF01486)	-
MYB superfamily	MYB	Myb_dna_bind (>=2) (PF00249)	-	SWIRM (PF04433)
MYB superfamily	MYB_related	Myb_dna_bind (1) (PF00249)	-	SWIRM (PF04433)
NAC	NAC	NAM (PF02365)	-	-
NF-X1	NF-X1	Zf-NF-X1 (PF01422)	-	-
NF-Y	NF-YA	CBFB_NFYA (PF02045)	-	-
NF-Y	NF-YB	NF-YB (self-build)	-	-
NF-Y	NF-YC	NF-YC (self-build)	-	-
Nin-like	Nin-like	RWP-RK (PF02042)	-	-
NZZ/SPL	NZZ/SPL	NOZZLE (PF08744)	-	-
S1Fa-like	S1Fa-like	S1FA (PF04689)	-	-
SAP	SAP	SAP (self-build)	-	-
SBP	SBP	SBP (PF03110)	-	-
SRS	SRS	DUF702 (PF05142)	-	-
STAT	STAT	STAT (self-build)	-	-
TCP	TCP	TCP (PF03634)	-	-
Trihelix	Trihelix	Trihelix (self-build)	-	-
VOZ	VOZ	VOZ (self-build)	-	-
Whirly	Whirly	Whirly (PF08536)	-	-
WRKY	WRKY	WRKY (PF03106)	-	-
ZF-HD	ZF-HD	ZF-HD_dimer (PF04770)	-	-
"""

animal_tf_family_assignment_rules = """
Family	DNA-binding domain	Pfam ID
AF-4	AF-4	PF05110
ARID	ARID	PF01388
bHLH	HLH	PF00010
CBF	CBF_alpha	PF02312
CEP-1	CEP1-DNA_bind	PF09287
CSL	BTD	PF09270
NF-YA	CBFB_NFYA	PF02045
CG-1	CG-1	PF03859
CP2	CP2	PF04516
CSD	CSD	PF00313
E2F	E2F_TDP	PF02319
ETS	Ets	PF00178
Fork head	Fork_head	PF00250
GCM	GCM	PF03615
GTF2I	GTF2I	PF02946
HMG	HMG_box	PF00505
HSF	HSF_DNA-bind	PF00447
HTH	HTH_psq	PF05225
IRF	IRF	PF00605
MYB	Myb_DNA-bd	PF00249
MBD	MBD	PF01429
NCU-G1	NCU-G1	PF15065
NDT80/PhoG	NDT80_PhoG	PF05224
Nrf1	Nrf1_DNA-bind	PF10491
PC4	PC4	PF02229
P53	P53	PF00870
PAX	PAX	PF00292
HPD	HPD	PF05044
RFX	RFX	PF02257
RHD	RHD	PF00554
Runt	Runt	PF00853
SAND	SAND	PF01342
SRF	SRF	PF00319
STAT	STAT_bind	PF02864
T-box	T-box	PF00907
TEA	TEA	PF01285
TSC22	TSC22	PF01166
Tub	Tub	PF01167
CTF/NFI	MH1	PF00859
MH1	MH1	PF03165
Homeobox	Homeobox	PF00046
Pou	Homeobox, Pou	PF00157
CUT	Homeobox, CUT	PF02376
TF_Otx	Homeobox, TF_Otx	PF03529
zf-C2HC	zf-C2HC	PF01530
zf-GAGA	zf-GAGA	PF09237
zf-BED	zf-BED	PF02892
ZBTB	zf-C2H2	PF00651
zf-C2H2	zf-C2H2	PF00096
DM	DM	PF00751
zf-GATA	zf-GATA	PF00320
zf-LITAF-like	zf-LITAF-like	PF10601
zf-MIZ	zf-MIZ	PF02891
zf-NF-X1	zf-NF-X1	PF01422
THAP	THAP	PF05485
"""


if __name__ == '__main__':
    # run hmmscan
    args = parse_args()
    domtblout = run_hmmscan(args)
    # tidy result
    domain_stat = get_domain_stat(domtblout)
    # print(domain_stat)
    if args.s == "plant":
        plant_tf_rule = parse_plant_tf_judge_rules(plant_tf_family_assignment_rules)
        query2family = judge_plant_tf(domain_stat, plant_tf_rule)
    else:
        animal_tf_rule = parse_animal_tf_judge_rules(animal_tf_family_assignment_rules)
        query2family = judge_animal_tf(domain_stat, animal_tf_rule)
    print("-----------Finish Hmmscan------------")

    # prepare input of blast
    if not query2family:
        print("NO TF was predicted and exit Normally")
        sys.exit(0)
    with open("predicted_tf.json", 'w') as f:
        json.dump(query2family, f, indent=4)
    get_predicted_tf_fasta(args.seqfile, query2family.keys(), "predicted_TFs.fa")
    if args.s == "plant":
        species_short_name = args.tfdb + '/plant_tf_pep/species_short_name.list'
        short_name_dict = dict()
        with open(species_short_name) as f:
            for line in f:
                if not line.strip():
                    continue
                short_name, sp = line.strip().split('\t')
                short_name_dict[sp.lower()] = short_name
    if args.s == "plant":
        pep_list = args.tfdb + '/plant_tf_pep/pep.list'
    else:
        pep_list = args.tfdb + '/animal_tf_pep/pep.list'
    pep_pd = pd.read_table(pep_list, header=0)

    # select sequence database for blast
    if args.blast_all == 'yes':
        top = 25 if args.organism != "unknown" else 1
        if args.s == "plant":
            target_seq = args.tfdb + '/plant_tf_pep/all_pep.dmnd'
        else:
            target_seq = args.tfdb + '/animal_tf_pep/all_pep.dmnd'
    else:
        top = 1
        if args.s == "plant":
            if args.organism == "unknown":
                target_seq = args.tfdb + '/plant_tf_pep/all_pep.dmnd'
            else:
                if args.organism.startswith("Ext-"):
                    if args.organism not in short_name_dict.values():
                        raise Exception("please provide correct short name of the plant species!")
                    target_seq = args.tfdb + '/plant_tf_pep/{}.pep.fa.dmnd'.format(args.organism)
                else:
                    if args.organism.capitalize() not in short_name_dict.values():
                        raise Exception("please provide correct short name of the plant species!")
                    target_seq = args.tfdb + '/plant_tf_pep/{}.pep.fa.dmnd'.format(args.organism.capitalize())
        else:
            if args.organism == "unknown":
                target_seq = args.tfdb + '/animal_tf_pep/all_pep.dmnd'
            else:
                if args.organism.lower().replace("_", ' ').capitalize() not in list(pep_pd['species']):
                    raise Exception("please provide correct name of the animal species!")
                target_seq = args.tfdb + '/animal_tf_pep/{}_transcription_factors.fasta.dmnd'.format(args.organism.capitalize())
    run_blast(
        target=target_seq,
        query="predicted_TFs.fa",
        top=top,
        out="diamond.out.txt",
        out_format=6,
        sensitive="sensitive",
        p=args.cpu,
        evalue=args.evalue,
        blast_type='blastp',
        diamond=args.diamond,
    )
    print("---------------finish diamond blast-----------------------")
    # select predicted target and format final result
    hmmscan_result = "domain_predict.txt"
    hmmscan_pd = pd.read_table(hmmscan_result, header=0)
    domain_score = hmmscan_pd.loc[:, ['query_id', 'pfam_id', 'e_value', 'score', "description"]].drop_duplicates()
    domain_score.index = [(x.split(".")[0]+'_'+y) for x, y in zip(domain_score['pfam_id'], domain_score['query_id'])]
    pep2species = dict(zip(pep_pd['pep_id'], pep_pd['species']))
    pep2family = dict(zip(pep_pd['pep_id'], pep_pd['family']))
    if args.s != "plant":
        pep2gene = dict(zip(pep_pd['pep_id'], pep_pd['gene_id']))
    blast_result = "diamond.out.txt"
    # blast_header = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    if os.path.getsize(blast_result) >= 2:  # 判断文件内容是否为空
        blast_pd = pd.read_table(blast_result, header=None, dtype={1: str})  # 第二列指定为字符串,以防有数字作为id且以0开头
        if top == 1:
            query2blast_hit = dict(zip(blast_pd[0], blast_pd[1]))
            query2blast_pident = dict(zip(blast_pd[0], blast_pd[2]))
            query2blast_evalue = dict(zip(blast_pd[0], blast_pd[10]))
        else:
            blast_stat = dict()
            if args.s != "plant":
                species_col = [pep2species[x].lower().replace(' ', '_') for x in blast_pd[1]]
            else:
                species_col = [pep2species[x].lower() for x in blast_pd[1]]
            blast_pd['species'] = species_col
            for _, row in blast_pd.iterrows():
                blast_stat.setdefault(row[0], list())
                blast_stat[row[0]].append(list(row)[1:])
            query2blast_hit = dict()
            query2blast_pident = dict()
            query2blast_evalue = dict()
            for query_id, hit_info in blast_stat.items():
                hit_species = [x[-1] for x in hit_info]
                if args.s == "plant":
                    hit_species = [short_name_dict[x].lower() for x in hit_species]
                if args.organism.lower() in hit_species:
                    tar_ind = hit_species.index(args.organism.lower())
                else:
                    tar_ind = 0
                query2blast_hit[query_id] = hit_info[tar_ind][0]
                query2blast_pident[query_id] = hit_info[tar_ind][1]
                query2blast_evalue[query_id] = hit_info[tar_ind][-3]
    else:
        print('本次diamond分析结果为空')
        query2blast_hit = dict()
        query2blast_pident = dict()
        query2blast_evalue = dict()

    # format final result
    final_result = list()
    for query_id, family_info in query2family.items():
        for family, domains in family_info.items():
            for each_domain in domains:
                if len(domain_score.loc[each_domain + '_' + query_id].shape) == 2:
                    # 由于通过pd.lc[x,y]可能返回多个结果，因此加入这样的判断
                    tmp_dict = dict(
                        query_id=query_id,
                        family=family,
                        domain=each_domain,
                        domain_link="http://pfam.xfam.org/family/{}".format(each_domain),
                        description=domain_score.loc[each_domain + '_' + query_id, "description"][0],
                        e_value=domain_score.loc[each_domain + '_' + query_id, "e_value"][0],
                        score=domain_score.loc[each_domain + '_' + query_id, "score"][0])
                else:
                    tmp_dict = dict(
                        query_id=query_id,
                        family=family,
                        domain=each_domain,
                        domain_link="http://pfam.xfam.org/family/{}".format(each_domain),
                        description=domain_score.loc[each_domain + '_' + query_id, "description"],
                        e_value=domain_score.loc[each_domain + '_' + query_id, "e_value"],
                        score=domain_score.loc[each_domain + '_' + query_id, "score"]
                    )
                if query_id in query2blast_hit:
                    hit_id = query2blast_hit[query_id]
                    if args.s == "plant":
                        short_name = short_name_dict[pep2species[hit_id].lower()]
                        if short_name.startswith("Ext-"):
                            ref_link = "http://planttfdb.gao-lab.org/tf_ext.php?sp={}&did={}".format(short_name.split("-")[1], hit_id)
                        else:
                            ref_link = "http://planttfdb.gao-lab.org/tf.php?sp={}&did={}".format(short_name, hit_id)
                    else:
                        gene_id = pep2gene[hit_id]
                        ref_link = "http://bioinfo.life.hust.edu.cn/AnimalTFDB#!/tf_gene_info?tf={}".format(gene_id)
                    tmp_dict.update(
                        blast_hit=query2blast_hit[query_id],
                        hit_link=ref_link,
                        hit_family=pep2family[hit_id],
                        hit_pident=query2blast_pident[query_id],
                        hit_evalue=query2blast_evalue[query_id],
                    )
                else:
                    tmp_dict.update(
                        blast_hit='',
                        hit_family='',
                        hit_pident='',
                        hit_link='',
                        hit_evalue='',
                    )
                final_result.append(tmp_dict)
    print("---------------------finish----------------------")
    # print(final_result)
    final_result_pd = pd.DataFrame(final_result)
    column_order = ['query_id', 'family', 'domain', 'domain_link', 'description', 'e_value', 'score']
    column_order += ['blast_hit', 'hit_link', 'hit_family', 'hit_pident', 'hit_evalue']
    final_result_pd = final_result_pd.loc[:, column_order]
    final_result_pd.to_csv("final_tf_predict.xls", header=True, index=False, sep='\t')

"""
< diomond output >
qseqid Query Seq - id 
qlen Query sequence length 
sseqid Subject Seq - id ----------------------------------->blast_hit
sallseqid All subject Seq - id(s), separated by a ’;’ 
slen Subject sequence length 
qstart Start of alignment in query 
qend End of alignment in query 
sstart Start of alignment in subject 
send End of alignment in subject 
qseq Aligned part of query sequence 
sseq Aligned part of subject sequence 
evalue Expect value ------------------------> hit_evalue
bitscore Bit score 
score Raw score 
length Alignment length 
pident Percentage of identical matches ------> hit_pident
nident Number of identical matches 
mismatch Number of mismatches 
positive Number of positive - scoring matches 
gapopen Number of gap openings 
gaps Total number of gaps 
ppos Percentage of positive - scoring matches 
qframe Query frame 
btop Blast traceback operations(BTOP) 
stitle Subject Title salltitles All Subject Title(s), separated by a ’<>’ 
qcovhsp Query Coverage Per HSP 
qtitle Query title 
Bydefault,thereare 12 preconﬁgured ﬁelds: 
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore. 
"""
