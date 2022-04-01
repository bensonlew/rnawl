# import cPickle as cpickle
import sqlite3
# nr2des_dict = {}
import sys

acc_file = sys.argv[1]
f=open(acc_file , 'r')
conn = sqlite3.connect('gene2accession.db')
cursor = conn.cursor()
cursor.execute('create table gene2accession( \
tax_id int, \
gene_id int, \
status varchar(20), \
nucl_acc_ver varchar(10), \
nucl_gi varchar(20), \
prot_acc_ver varchar(10), \
prot_gi varchar(20), \
ge_nucl_ver varchar(20), \
ge_nucl_gi varchar(20), \
ge_start varchar(20), \
ge_end varchar(20), \
orient varchar(2), \
assembly varchar(20), \
mat_pep_ver varchar(20), \
mat_pep_gi varchar(20), \
symbol varchar(20))')

f.readline()

for i in f:
    cols = i.strip().split("\t")
    try:
        insert_sql = 'insert into gene2accession (tax_id, gene_id, status, nucl_acc_ver, nucl_gi, prot_acc_ver, prot_gi, ge_nucl_ver, ge_nucl_gi, ge_start, ge_end, orient, assembly, mat_pep_ver, mat_pep_gi, symbol) values ( {}, {}, \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \'{}\', \"{}\")'.format(
            str(cols[0]),
            str(cols[1]),
            str(cols[2]),
            str(cols[3]),
            str(cols[4]),
            str(cols[5]),
            str(cols[6]),
            str(cols[7]),
            str(cols[8]),
            str(cols[9]),
            str(cols[10]),
            str(cols[11]),
            str(cols[12]),
            str(cols[13]),
            str(cols[14]),
            str(cols[15]),
        )
        cursor.execute(insert_sql)
    except:
        print i
cursor.close()
conn.commit()
conn.close()
f.close()
