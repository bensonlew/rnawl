#zouguanqing
import sys

def get_code(code_file):
    with open(code_file) as f:
        line = f.readline()
        code_data = eval(line.strip())
    return code_data


def trans_fun(seq,indic):
    faa = ''
    for i in range(0,len(seq),3):
        nn = seq[i:i+3]
        if nn not in indic.keys():
            print nn
            aa = '*'
        else:
            aa = indic[nn]['single']
        faa += aa
    if faa != '':
        faa += '\n'
    return faa


def main(code_file,trans_code,ffn, out_faa):
    code_data = get_code(code_file)
    indic = code_data[trans_code]
    seq = ''
    title = ''
    seq_list = []
    tmp = []
    with open(ffn) as f:
        for line in f:
            if line.startswith('>'):
                if len(seq)%3 != 0:
                    print title
                    raise Exception('seq length is not 3*')
                tmp.append(title)
                tmp.append(seq)
                seq_list.append(tmp)
                tmp = []
                title = line
                seq = ''
            else:
                seq += line.strip()

        seq_list.append([title,seq])

    fw = open(out_faa,'w')
    for t,seq in seq_list:
        seq_up = seq.upper()
        aa_seq = trans_fun(seq_up, indic)
        fw.write(t)
        fw.write(aa_seq)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('python trans_table.py code_file trans_table ffn out_faa')
        exit()
    code_file = sys.argv[1]
    trans_code = sys.argv[2]
    ffn = sys.argv[3]
    out_faa = sys.argv[4]
    main(code_file,trans_code,ffn, out_faa)



