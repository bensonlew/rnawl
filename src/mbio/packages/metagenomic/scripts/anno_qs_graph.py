#! /usr/bin/python
import sys
import pandas as pd
from pandas.core.frame import DataFrame

data=pd.read_table(sys.argv[1],sep='\t',header=0,index_col="#Class")
data =data.T
data['class']=data.index
data2=pd.read_table(sys.argv[2],sep='\t')
pp=['class','group']
data2.columns=pp 
data3 = data[data['class'].isin(list(data2['class']))]
names=list(data3.columns)
del names[-1]
da=pd.merge(data3,data2)
del da['class']
aa=dict(list(da.groupby(['group'])))
with open ("anno_qs_graph.xls",'w') as f:
    for i in aa.keys():
        for name in names:
            temp=dict(DataFrame(aa[i]).describe()[name])
            list1= list(aa[i][name])
            IQR =float(temp['75%'])-float(temp['25%'])
            list1=[float(x) for x in list1] 
            outlier_step = 1.5 * IQR
            min=float(temp['25%'])- outlier_step
            max=float(temp['75%']) + outlier_step
            new_list=[]
            for x in list1:
               if float(x) >=min and float(x) <=max:
                  new_list.append(x)
            outlier_list_col=[]
            for x in list1:
               if float(x) <min or float(x) > max:
                  outlier_list_col.append(x)
            new_list.sort()
            duu = []
            if float(new_list[0]) <float(temp['25%']):
                duu.append(float(new_list[0]))
            else:
                duu.append(float(temp['25%']))  
            if float(new_list[-1]) >float(temp['25%']):
                duu.append(float(new_list[-1]))
            else:
                duu.append(float(temp['75%']))  
            for j in temp.keys():
                if j == '25%' or j == '50%' or j == '75%':
                    duu.append(float(temp[j]))
            duu.sort()
            duu_new = [str(x) for x in duu]
            de=','.join(duu_new)
            duu2_new = [str(x) for x in outlier_list_col]
            if len(duu2_new) >=1:
                de2=','.join(duu2_new)
            else:
                de2="-"
            des= i + "\t" + name + "\t" + de + "\t" + de2
            f.write(des + "\n")
