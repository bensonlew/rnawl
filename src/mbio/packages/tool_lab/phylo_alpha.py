import os
import sys
method = sys.argv[1]
tree = sys.argv[2]
otu = sys.argv[3]
rscript = sys.argv[4] #/mnt/ilustre/users/sanger-dev/app/program/miniconda3/envs/R/bin/Rscript
r = '''
library(picante)
tree = read.tree("{}")
otu = read.table("{}", header=TRUE, row.names=1,sep=\"\t\")
otu = t(otu)
dis = prune.sample(otu,tree)
#otu = otu[,colnames(dis)]
otu = otu[,dis$tip.label]
'''.format(tree,otu)

if method in ['mpd', 'nri']:
    r += '''
res = ses.mpd(otu, cophenetic(dis),abundance.weighted=TRUE,null.model="taxa.labels",runs=99)
'''
    if method == 'mpd':
        value_label = 'mpd.obs'
        r += '''
out_data = data.frame(sample=rownames(res), MPD=res$mpd.obs)
'''
    else:
        r += '''
out_data = data.frame(sample=rownames(res), NRI=res$mpd.obs.z*-1)
'''


elif method in ['mntd', 'nti']:


    r += '''
res = ses.mntd(otu, cophenetic(dis),abundance.weighted=TRUE,null.model="taxa.labels",runs=99)

'''
    if method == 'mntd':
         r += '''
out_data = data.frame(sample=rownames(res), MNTD=res$mntd.obs)
'''
    else:
         r += '''
out_data = data.frame(sample=rownames(res), NTI=res$mntd.obs.z*-1)

'''

elif method in ['pd']:
    r += '''
res = ses.pd(otu, dis, include.root=FALSE, null.model="taxa.labels", runs=99)
out_data = data.frame(sample=rownames(res), PD=res$pd.obs)

'''

r += '''
write.table(out_data, '{}.xls', sep='\t',quote=FALSE, row.names =FALSE)
'''.format(method)

with open('phylo_alpha.r', 'w') as fw:
    fw.write(r)


os.system('%s phylo_alpha.r'%(rscript))


