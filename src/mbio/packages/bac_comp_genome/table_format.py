# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import pandas as pd
import os


def _get_header(functype):
    """不同注释文件中的表头
    :return : 返回对应注释的表头
    """
    headers = {
        'KEGG': [
            'Pathway', 'sample', 'KO', 'Definition', 'Enzyme',
            'Module', 'Hyperlink', 'Level1', 'Level2', 'Level3'
        ],
        'COG': [
            'Function', 'sample', 'NOG', 'NOG_description',
            'Function_description', 'Category'
        ]
    }
    if functype in headers:
        return headers[functype]


def _flat(line, functype, retval):
    """单行内的多个注释是展开成多行
    :params line: 一行注释信息
    :params functype: 可能存在个注释的注释功能类型，可按情况预定义
    :papams retval: 展开成多行的注释存入此列表中
    """
    lm = lambda x: line[x].split(';')

    if functype == 'KEGG':
        kos = lm('Pathway')
        l1 = lm('Level1')
        l2 = lm('Level2')
        l3 = lm('Level3')
        if len(kos) == len(l1):
            for i in range(len(kos)):
                newline = [
                    line['Gene ID'], kos[i], line['sample'],
                    line['KO'], line['Definition'], line['Enzyme'],
                    line['Module'], line['Hyperlink'],
                    l1[i], l2[i], l3[i]
                ]
                retval.append(newline)
    elif functype == 'COG':
        funcs = lm('Function')
        fun_des = lm('Fun_description')
        cat = lm('Category')
        for i in range(len(funcs)):
            newline = [
                line['Gene ID'], funcs[i], line['sample'],
                line['NOG'], line['NOG_description'],
                fun_des[i], cat[i]
            ]
            retval.append(newline)


def annot_format(annot, functype, level=None):
    """将以基因为索引的注释表转化为以选定功能层级为索引
    :params annot: 以基因为第一列的注释表
    :params functype: 以上文件注释的功能类型
    :params level: 列表，只返回包含指定注释层级的结果
    :return: 返回格式化后的pandas dataframe
    """
    functype = functype.upper()
    annot = pd.read_csv(annot, sep='\t')
    #annot = annot.rename(columns={'#Query': 'Gene ID'})
    if 'sample' not in annot:
        annot['sample'] = 'sample'
    if _get_header(functype):
        header = ['Gene ID', ] + _get_header(functype)
        tmp = []
        annot.apply(_flat, axis=1, args=(functype, tmp))
        annot = pd.DataFrame(data=tmp, columns=header)
    # exit()
    if level:
        s_level = []  # 用于前端传参的大小写兼容
        for l in level:
            for c in annot.columns:
                if c.lower() == l.lower():
                    s_level.append(c)
                    continue
        annot = annot[annot[s_level[0]] != '-']
        new_header = ['Gene ID', 'sample'] + s_level
        annot = annot[new_header]

    index = list(set(annot.columns) - set(['Gene ID', ]))
    annot = annot.groupby(index).agg(lambda x: ','.join(set(x)))

    return annot.reset_index().rename(columns={'Gene ID': 'genelist'})


def reshape(df, col, v, abundance='F', result_type='none'):
    """选定两列对pandas dataframe进行重塑
    :params col: 以col列作为列名进行重塑
    :params v: 以v列作为值进行重塑
    :return: 返回重塑后的dataframe
    """
    index = list(set(df.columns) - set([v, ]))
    fna = 0
    if result_type not in ['mean', 'median']:
        df = df.groupby(index).agg(lambda x: ','.join(set(x)))

        # 使用 set() 去除 groupby 后可能的重复
        if abundance == 'T':
            df[v] = map(lambda x: len(set(x.split(','))), df[v])
        else:
            fna = '-'
            df[v] = map(lambda x: ','.join(set(x.split(','))), df[v])
    else:
        df[v] = map(lambda x: len(set(x.split(','))), df[v])
        if result_type == 'mean':
            df = df.groupby(index).agg(lambda x: x.mean())
        else:
            df = df.groupby(index).agg(lambda x: x.median())

    df = df.unstack(col).fillna(fna).reset_index()
    columns = [b or a for a, b in df.columns]
    df.columns = columns

    return df


def prefix_name(fp, prefix):
    fp = os.path.split(os.path.abspath(fp))
    prefix_fp = prefix + '_' + fp[1]
    return os.path.join(fp[0], prefix_fp)


if __name__ == '__main__':
    import sys
    l = None
    if len(sys.argv) == 4:
        l = sys.argv[3].split(',')
    t = annot_format(sys.argv[1], sys.argv[2], l)
    o = prefix_name(sys.argv[1], 'format')
    print o
    t.to_csv(o, sep='\t', index=None)
