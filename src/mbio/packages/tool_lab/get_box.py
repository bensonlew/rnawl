# -*- coding:utf-8  -*-

import pandas as pd


def get_box_value_row(exp_path,group_path, out_file='box_out.xls',method=None):

    exp_data = pd.read_table(exp_path, sep='\t',index_col=0)
    if method == "T":
        exp_data = exp_data.apply(lambda x: (x*100) / x.sum())
    final_data = []
    group = pd.read_table(group_path, sep='\t')
    group.rename(columns={group.columns[0]:'sample',group.columns[1]:'group'}, inplace=True)
    g_map_s = dict(list(group.groupby(by='group')))
    for g in  g_map_s:
        tmp_data = exp_data[g_map_s[g]['sample']]
        box = tmp_data.apply(lambda x: x.describe(),axis=1)
        delta_q = box['75%'] - box['25%']
        tmp_data['up'] = box['75%'] + delta_q*1.5
        tmp_data['down'] = box['25%'] - delta_q*1.5

        sub_min = []
        sub_max = []
        abnormal = []
        ##找异常点
        for index in tmp_data.index:
            up = tmp_data['up'][index]
            down = tmp_data['down'][index]
            row = tmp_data.loc[index]
            row.pop('up')
            row.pop('down')
            ab = row[row.map(lambda x: True if x>up or x<down else False)]
            nor = row[row.map(lambda x: False if x>up or x<down else True)]
            sub_min.append(nor.min())
            sub_max.append(nor.max())
            abnormal.append(';'.join([str(s[0])+':'+str(s[1]) for s in zip(ab.index.tolist(),ab.tolist())]))

        box.drop(['count'],axis=1,inplace=True)
        box[g+'--sub-min'] = sub_min
        box[g+'--sub-max'] = sub_max
        box[g+'--abnormal'] = abnormal
        box.rename(columns={'mean':g+'--mean','std':g+'--std','min':g+'--min','max':g+'--max','25%':g+'--25%','75%':g+'--75%', '50%':g+'--50%'},inplace=True)
        final_data.append(box)

    final = pd.concat(final_data,axis=1)
    final.to_csv(out_file,sep='\t')


def get_box_value_column(result_path,group_path, method, out_file='box_out.xls'):

    data_list = []
    data = pd.read_table(result_path, sep='\t',index_col=0)
    data.fillna(0,inplace=True)

    group = pd.read_table(group_path, sep='\t')

    group.rename(columns={group.columns[0]:'sample',group.columns[1]:'group'}, inplace=True)
    g_map_s = dict(list(group.groupby(by='group')))
    s_map_g = dict(list(group.groupby(by='sample')))

    if method in ['mean', 'sum', 'median']:
        for g in g_map_s.keys():
            if method=='mean':
                data[g+'_merge'] = data[g_map_s[g]['sample'].tolist()].apply(lambda x: x.mean(),axis=1)
            elif method =='sum':
                data[g+'_merge'] = data[g_map_s[g]['sample'].tolist()].apply(lambda x: x.sum(),axis=1)
            else:
                data[g+'_merge'] = data[g_map_s[g]['sample'].tolist()].apply(lambda x: x.median(),axis=1)

        for g in g_map_s.keys():
            res = data[g+'_merge'].describe()
            delq= res['75%'] - res['25%']

            ##异常点
            up_out_line = res['75%'] + delq*1.5
            down_out_line = res['25%'] - delq*1.5
            error_dot = data[g+'_merge'][data[g+'_merge'].map(lambda x: True if x > up_out_line or x < down_out_line else False)]
            normal_dot = data[g+'_merge'][data[g+'_merge'].map(lambda x: False if x > up_out_line or x < down_out_line else True)]
            nor_min = normal_dot.min()
            nor_max = normal_dot.max()

            ##异常点
            scatter = []
            for index in error_dot.index:
                scatter.append(str(index)+'||'+ str(error_dot[g+'_merge'][index]))
            scatter = ';'.join(scatter)

            insert_data = {
                "name" : g,
                "q1" : res['25%'],
                "q3" : res['75%'],
                "median" : res['50%'],
                "all_min" : res['min'],
                "all_max" : res['max'],
                'min' : nor_min,
                'max' : nor_max,
                'abnormal' : scatter
            }

            data_list.append(insert_data)
    else:
        for s in s_map_g.keys():
            res = data[s].describe()
            delq= res['75%'] - res['25%']

            ##异常点
            up_out_line = res['75%'] + delq*1.5
            down_out_line = res['25%'] - delq*1.5
            error_dot = data[s][data[s].map(lambda x: True if x > up_out_line or x < down_out_line else False)]
            normal_dot = data[s][data[s].map(lambda x: False if x > up_out_line or x < down_out_line else True)]
            nor_min = normal_dot.min()
            nor_max = normal_dot.max()

            ##异常点
            scatter = []
            for index in error_dot.index:
                scatter.append(str(index)+'||'+ str(error_dot[s][index]))
            scatter = ';'.join(scatter)

            insert_data = {
                "name" : s,
                'category': s_map_g[s],
                "q1" : res['25%'],
                "q3" : res['75%'],
                "median" : res['50%'],
                "all_min" : res['min'],
                "all_max" : res['max'],
                'min' : nor_min,
                'max' : nor_max,
                'abnormal' : scatter
            }

            data_list.append(insert_data)

    out = pd.DataFrame(data_list)
    out.to_csv(out_file,sep='\t',index=False)

if __name__ == '__main__':
    import sys
    abund_file = sys.argv[1]
    group_file = sys.argv[2]
    method = sys.argv[3]
