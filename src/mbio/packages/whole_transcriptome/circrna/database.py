import pandas as pd

def database_Homo(detail,database,organism_name,detail_new):
    if organism_name == 'Homo_sapiens':
        circbase_True_False = list()
        mine_circrna = pd.read_csv(detail,sep = '\t')
        circbase = pd.read_csv(database,sep = '\t')
        for circrna_id in mine_circrna['circrna_id']:
            if circrna_id in circbase['chr_start_end'].values:
                circbase_True_False.append("yes")
            else:
                circbase_True_False.append("no")
        mine_circrna['circbase']=circbase_True_False
        mine_circrna.to_csv(detail_new,index=False,sep = '\t')
    else:
        circbase_True_False = list()
        mine_circrna = pd.read_csv(detail,sep = '\t')
        circbase = pd.read_csv(database,sep = '\t')
        for start_end in mine_circrna['start_end']:
            if start_end in circbase['start_end'].values:
                circbase_True_False.append("yes")
            else:
                circbase_True_False.append("no")
        mine_circrna['circBase']=circbase_True_False
        mine_circrna.to_csv(detail_new,index=False,sep = '\t')

def main(args):
    database_Homo(args.detail,args.database,args.organism_name,args.detail_new)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='circBase')
    parser.add_argument('-i', action='store', required=True,
                        help='detail', dest='detail')
    parser.add_argument('-d', action='store', required=True,
                        help='database', dest='database')
    parser.add_argument('-r', action='store', required=True,
                        help='organism_name ', dest='organism_name')
    parser.add_argument('-o', action='store', required=True,
                        help='output ', dest='detail_new')
    args = parser.parse_args()

    main(args)
