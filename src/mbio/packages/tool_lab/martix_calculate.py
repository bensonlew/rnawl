# /usr/bin/python

import os,sys,re
import pandas as pd
from pandas.core.frame import DataFrame

data = pd.read_table(sys.argv[1],sep='\t', header=0)
data = data.set_index(['sample'],drop=True)
samples =data.index
data = data.replace('-', 0) 
ids = len(data.columns)
df1 = []
for i in samples:
   for j in samples:
      if i ==j:
         df1.append([i,j,0.0])
      else:
         num=0
         for d in zip(data.ix[i],data.ix[j]):
            if d[0] != d[1] and d[0] != 0 and d[1] != 0:
               num +=1
            elif d[0] != d[1] and d[0] ==0 or d[1] ==0:
               num +=0.5
         de = num/float(ids)
         df1.append([i,j,de])   
data3 = DataFrame(df1)
data4 = data3.groupby([0, 1])[2].apply(lambda x:float(x)).unstack()
des = len(samples)
data4.to_csv(sys.argv[2], sep='\t', header=False, index=True)
os.system("sed -i '1i\ {}' {}".format(" "+str(des),sys.argv[2]))
