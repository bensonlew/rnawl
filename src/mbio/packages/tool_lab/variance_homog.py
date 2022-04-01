from scipy.stats import levene, bartlett, fligner

fun_dic = {
    'levene' :levene,
    'bartlett' : bartlett,
    'ligner-killeen' : fligner
}

def deal_group(group_file):
    rdic = {}
    with open(group_file) as f:
        f.readline()
        for line in f:
            spl = line.strip().split('\t')
            if spl[1] in rdic :
                rdic[spl[1]].append(spl[0])
            else:
                rdic[spl[1]] = [spl[0]]
    return rdic

def get_sample_index(header):
    sample_index = dict()
    for id,h in enumerate(header[1:],1):
        sample_index[h] = id
    return sample_index

def get_share(sample_index, gdic):
    new_g_dic = dict()
    for g in sorted(gdic.keys()):
        tmp_ids = []
        for s in gdic[g]:
            if s in sample_index:
                sid = sample_index[s]
                tmp_ids.append(sid)
        if len(tmp_ids) > 0:
            new_g_dic[g] = tmp_ids

    return new_g_dic


def get_data_by_group(line,g_sample_id):
    ret_data = []

    for g in sorted(g_sample_id.keys()):
        tmp = []
        for sid in g_sample_id[g]:
            tmp.append(float(line[sid]))
        ret_data.append(tmp)
    return ret_data


def test_fun(infile, group, method, outfile):
    group_dic = deal_group(group)
    cur_fun = fun_dic[method]
    with open(infile) as fr,open(outfile,'w') as fw:
        fw.write('name\tstatistic\tpvalue\n')
        sph =fr.readline().strip().split('\t')
        sample_index = get_sample_index(sph)
        g_sample_id = get_share(sample_index, group_dic)
        for line in fr:
            spline = line.strip().split('\t')
            gdata = get_data_by_group(spline, g_sample_id)
            statistic ,pvalue = cur_fun(*gdata)
            fw.write(spline[0]+'\t'+str(statistic)+'\t'+str(pvalue)+'\n')


if __name__ == '__main__':
    import sys
    infile, group, method, outfile = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
    test_fun(infile, group, method, outfile)