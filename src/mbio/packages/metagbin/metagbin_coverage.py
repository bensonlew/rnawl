##! /usr/bin/python

import re,os,sys

files =os.listdir(sys.argv[1])

dict ={}
dict2 = {}
list =[]
list2 = []
n =0
for file in files:
    if re.search('.txt',file):
        with open (sys.argv[1] + '/' + file,'r') as f:
            n +=1
            lines =f.readlines()
            names =lines[0].strip('\n').split("\t")
            for name in names:
                if n ==1:
                    list.append(name)
                else:
                    if names.index(name) >=3:
                        list.append(name)
            for line in lines[1:]:
                lin =line.strip('\n').split("\t")
                if n ==1:
                    s = lin[0] + "\t" + lin[1]
                    list2.append(s)
                des ="\t".join(lin[3:])
                ss = lin[0] + "\t" + lin[1]
                if ss in dict:
                    dict[ss] +="\t" +  str(des)
                else:
                    dict[ss] = str(des)
                if ss in dict2:
                    dict2[ss] +=float(lin[2])
                else:
                    dict2[ss] =float(lin[2])

with open (sys.argv[2],'w') as g:
   g.write("\t".join(list)+"\n")
   for key in list2:
      g.write(key+ "\t"+ str(dict2[key]) + "\t" + dict[key]+"\n")