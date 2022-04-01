import sys
import re

def get_info(infile,index_1,index_2,head=False):
    info = []
    with open(infile) as f:
        if head:
            f.readline()
        for line in f:
            line = line.strip()
            if line == '':
                continue
            spline = line.split('\t')
            info_1 = spline[index_1]
            info_2 = spline[index_2]
            info.append([info_1,info_2])
    return info

def get_pfam_info(infile,head=False):
    info = []

    with open(infile) as f:
        if head:
            f.readline()
        for line in f:
            line = line.rstrip('\n')
            if line == '':
                continue
            spline = line.split('\t',2)
            gene = spline[0]
            pfam = spline[1].split('.')[0]

            info.append([gene,pfam, pfam + '\t' + spline[2]])
    return info

def search_senser_regulator(gene_pfam, senser_regulator_data, outfile,has_score):
    data = {}
    ret_dic = {}
    for i in senser_regulator_data:
        data[i[0]]= i[1]
    k_list = data.keys()
    gene_list = []
    ret_res = {}
    for j in gene_pfam:
        pfam = j[1]
        gene = j[0]
        des = j[2]
        pfam=pfam.split('.')[0]
        if pfam in k_list:
            if gene not in ret_dic.keys():
                ret_dic[gene] = data[pfam]
                gene_list.append(gene)
                ret_res[gene] = {"type":data[pfam],"desc":[des]}
            else:
                if ret_dic[gene] == 'hybrid':
                    continue
                elif ret_dic[gene] != data[pfam]:
                    ret_dic[gene] = 'hybrid'
                    ret_res[gene]['type'] = 'hybrid'
                    ret_res[gene]['desc'].append(des)


    with open(outfile, 'w') as fw:
        if has_score:
            fw.write("Gene\tType\tPfam_id\tDomain\tDomainDescription\tStart\tEnd\tPfamStart\tPfamEnd\tDomainE-Value\tscore\tLocation\tGene Description\n")
        else:
            fw.write("Gene\tType\tPfam_id\tDomain\tDomainDescription\tStart\tEnd\tPfamStart\tPfamEnd\tDomainE-Value\tLocation\tGene Description\n")
        for k in gene_list:
            if has_score:
                des_list = [[],[],[],[],[],[],[],[],[],[],[]]
            else:
                des_list = [[],[],[],[],[],[],[],[],[],[]]

            for des in ret_res[k]["desc"]:
                for id, e in enumerate(des.split('\t')):
                    des_list[id].append(e)
            des_list2 = []
            for d in des_list[:-2]:
                des_list2.append(';'.join(d))

            des_list2.append(des_list[-2][0])
            des_list2.append(des_list[-1][0])
            des = '\t'.join(des_list2)
            fw.write(k + '\t' + ret_res[k]['type']+'\t'+ des +'\n')
    return ret_dic

def do_stat(ret_dic,outfile):
    value_list = ret_dic.values()
    hybrid_n = value_list.count('hybrid')
    senser_n = value_list.count('sensor')
    regulator_n = value_list.count('regulator')
    with open(outfile,'w') as fw:
        fw.write('senser\tregulator\thybrid\n')
        fw.write('\t'.join([str(i) for i in [senser_n,regulator_n,hybrid_n]])+'\n')


if __name__ == '__main__':
    senser_regulator_file = sys.argv[1]
    database = get_info(senser_regulator_file,0,2)
    pfam_anno = sys.argv[2]
    with open(pfam_anno) as f:
        line = f.readline()
        spline = line.split('\t')
        if 'score' in spline:
            has_score = True
        else:
            has_score = False

    pfam_data  = get_pfam_info(pfam_anno)
    ret = search_senser_regulator(pfam_data, database, 'senser_regulator.xls',has_score)
    do_stat(ret,'senser_regulator.stat')
        
            

