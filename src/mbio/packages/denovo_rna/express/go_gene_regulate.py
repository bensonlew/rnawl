# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# __version__ = 'v1.0'
# __last_modified__ = '20170629'


import re
import os

def up_down_express_list(fp):
    """
    """
    group_name = os.path.basename(fp).split("_edgr_")[0]
    with open(fp) as f:
        header = f.readline().strip()
        if not header.lower().endswith('pvalue\tpadjust\tsignificant\tregulate\tncbi'):
            raise Exception('错误的表头:{}，表头应与基因差异表达分析结果格式一致'.format(header))
        down = []
        up = []
        for line in f:
            line_sp = line.strip().split('\t')
            if line_sp[-3] == 'yes':
                if line_sp[-2] == 'up':
                    up.append(line_sp[0])
                elif line_sp[-2] == 'down':
                    down.append(line_sp[0])
                elif line_sp[-2] == 'no change':
                    continue
                elif line_sp[-2] == 'undone':
                    raise Exception('文件中检查到‘undone’，表明文件没有上下调控信息')
                else:
                    print line_sp
                    raise Exception('未知的上下调说明类型:{}，必须为up或者为down'.format(line_sp[-1]))
        return up, down, group_name


def get_level_2_info(fp):
    with open(fp) as f:
        header = f.readline().strip()
        if header != 'term_type\tterm\tGO\tnumber\tpercent\tsequence':
            raise Exception('2层级GO注释统计表格式错误:{}'.format(header))
        all_gene_list = set()
        level2info = []
        for line in f:
            line_sp = line.strip().split('\t')
            gene_list = [re.sub(r'\(GO:.+\)', '', i) for i in line_sp[-1].split(';')]
            all_gene_list.update(gene_list)
            # print gene_list
            level2info.append([line_sp[0], line_sp[1], line_sp[2], gene_list])
    return level2info


def get_level_2_up_down(level2fp, up, down, outfile='./GO_up_down.xls', name=""):
    with open(outfile, 'wb') as w:
        all_genes = set()  # 保存注释和上下调同时存在基因
        all_records = []  # 保存上下调基因
        head = "term_type\tterm\tGO\t{}_up num\t{}_up percent\t{}_up genes\t{}_down num\t{}_down percent\t{}_down genes\n".format(name, name, name, name, name, name)
        w.write(head)
        go_info = get_level_2_info(level2fp)  # level2的go统计表的信息 二维列表
        # print go_info
        for record in go_info:
            # print record
            # break
            one_up = []
            one_down = []
            for one in record[3]:
                if one in up:
                    one_up.append(one)
                elif one in down:
                    one_down.append(one)
            all_records.append([one_up, one_down])
            all_genes.update(one_up)
            all_genes.update(one_down)
        genes_num = float(len(all_genes))

        for i in range(len(all_records)):
            up_num = len(all_records[i][0])
            down_num = len(all_records[i][1])
            up_percent = up_num / genes_num
            down_percent = down_num / genes_num
            up_genes = ';'.join(all_records[i][0])
            down_genes = ';'.join(all_records[i][1])
            new_line = '\t'.join(go_info[i][:3] + [str(up_num), str(up_percent), up_genes, str(down_num), str(down_percent), down_genes])
            w.write(new_line + '\n')


def GO_level_2_regulate(diff_express_fp, go_level_2_stat_fp, outfile):
    up, down, group_name = up_down_express_list(diff_express_fp)
    file_name = os.path.splitext(os.path.basename(diff_express_fp))[0].split("_edgr_stat")[0]
    # print up
    # print down
    get_level_2_up_down(go_level_2_stat_fp, up, down, outfile, group_name)


if __name__ == '__main__':
    # up_down_express_list('C:\Users\sheng.he.MAJORBIO\Desktop\P7_2_vs_E20_2_edgr_stat.xls')
    # get_level_2_up_down('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\go2level(2).xls',
    # ["c373_g1_i1", "c426_g1_i1", "c527_g1_i1", "c1571_g1_i1", "c1227_g1_i1", "c5101_g1_i1", "c1998_g1_i1", "c296_g1_i1", "c3383_g1_i1"],
    # ["c1546_g1_i1","c376_g1_i1","c1909_g1_i1","c2964_g1_i1"],
    # 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\temp.xls')
    GO_level_2_regulate("/mnt/ilustre/users/sanger-test/workspace/20170630/Refrna_demo1/Express/output/diff/genes_diff/diff_stat_dir/B1_vs_Z1_edgr_stat.xls","/mnt/ilustre/users/sanger-test/workspace/20170630/Refrna_demo1/MergeAnnot1/go2level.xls", "output")