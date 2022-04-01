import os
import sys


new_gene_xls = sys.argv[1]
new_gene_xls_log = sys.argv[2]
png_out = sys.argv[3]
r_script = sys.argv[4]

exp_f = os.popen('grep "newgenes(g)" %s|cut -d "=" -f 2|sed "s/g/x/g"'%(new_gene_xls_log))
exp_str = exp_f.read()
exp_f.close()
print(exp_str)

if exp_str :
    exp_r = 'curve(%s,1,n,add=T,col="red")'%(exp_str)
else:
    exp_r = ''


s =  '''
data = read.table('%s', sep='\t', header=TRUE)

c2= c()
c1 = c()
c_ave = c()
c_ave_lab = c()
n=0
for (i in names(data)){
    n = n+1;
    for(j in data[i]){
        c2 = c(c2, j);
        }

    c1 = c(c1, rep(n, lengths(data[i])))
    c_ave = c(c_ave, mean(unlist(data[i])))
    c_ave_lab = c(c_ave_lab, n)
}

new_data = data.frame(genome=c1, genes=c2)
#new_data

ave_data = data.frame(genome=c_ave_lab, genes=c_ave)
#ave_data


png(file='%s')
barplot(ave_data$genes,names.arg=ave_data$genome, xlab='Genome Number', ylab='New Gene Cluster Number', main='New Genes',col='orange')
%s
dev.off()
'''%(new_gene_xls , png_out ,exp_r)


with open('new_gene_plot.R' ,'w') as fw:
    fw.write(s)

os.system('{} {}'.format(r_script, 'new_gene_plot.R'))
