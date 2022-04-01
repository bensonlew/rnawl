# -*- coding: utf-8 -*-
# __author__ = 'qiuping'


class cog_stat(object):
    def __init__(self):
        self.func_type = {
            'INFORMATION STORAGE AND PROCESSING': sorted(['J', 'A', 'K', 'L', 'B']),
            'CELLULAR PROCESSES AND SIGNALING': sorted(['D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O']),
            'METABOLISM': sorted(['C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q']),
            'POORLY CHARACTERIZED': sorted(['R', 'S']),
        }
        self.func_decs = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification', 'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure', 'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'C': 'Energy production and conversion', 'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only', 'S': 'Function unknown'
        }
        self.totalseq = 0
        self.gene_list = []

    def get_gene_list(self, gene_file):
        """传入基因序列文件，返回装有基因序列的列表"""
        with open(gene_file, 'rb') as f:
            for line in f:
                self.gene_list.append(line.strip('\n'))

    def get_gene_cog_list(self, cog_list, gene_list, outpath, trinity_mode=True):
        """
        将string2cog注释的结果文件筛选只包含基因的结果信息
        cog_list:string2cog tool运行得到的cog_list.xls结果文件；
        gene_list: 只包含基因序列名字的列表
        trinity_mode用于在新生成的xml的queryID是去除结尾的_i(数字) 的
        return: gene_cog_list.xls
        """
        with open(cog_list, 'rb') as c, open(outpath, 'wb') as w:
            head = c.readline()
            w.write(head)
            gene_name = []
            for line in c:
                line = line.strip('\n').split('\t')
                name = line[0]
                if name in gene_list:
                    if trinity_mode:
                        name = name.split('_i')[0]
                    w.write('{}\t{}\t{}\n'.format(name, line[1], line[2]))
                    if name not in gene_name:
                        self.totalseq += 1
                    gene_name.append(name)

    def get_gene_summary(self, cog_table, gene_list, out_dir, trinity_mode=True):
        """
        将string2cog注释的结果文件筛选只包含基因的结果信息:cog_table.xls, cog_summary.xls
        cog_table:string2cog tooly运行得到的cog_table.xls结果文件；
        gene_file: 只包含基因序列名字的列表
        trinity_mode用于在新生成的xml的queryID是去除结尾的_i(数字) 的
        return: gene_cog_table.xls, gene_cog_summary.xls
        """
        fun_stat = {'COG': {}, 'NOG': {}}  # {COG: {A:2, B:3,...,Z:3}, NOG: {A:2, B:3,...,Z:3}}
        cog_list = []
        gene_cog_table = out_dir + '/gene_cog_table.xls'
        gene_cog_summary = out_dir + '/gene_cog_summary.xls'
        with open(cog_table, 'rb') as c, open(gene_cog_table, 'wb') as t, open(gene_cog_summary, 'wb') as s:
            head = c.readline()
            t.write(head)
            first_gene = None
            for l in c:
                line = l.strip('\n').split('\t')
                gene = line[0]
                if not first_gene:
                    first_gene = gene
                if gene in gene_list:
                    w_line = line
                    if trinity_mode:
                        tmp = gene.split('_i')[0]
                        w_line[0] = tmp
                    w_line = '\t'.join(w_line)
                    t.write(w_line + '\n')
                    # deal with cog_summary.xls: get fun_stat
                    cog_id = line[-8]
                    if len(line[-6]) > 1:
                        func_cat = ','.split(line[-6])  # fun type: [H,Z] or [H]
                    else:
                        func_cat = [line[-6]]
                    if gene == first_gene:
                        if cog_id not in cog_list:
                            if cog_id.startswith('COG'):
                                for i in func_cat:
                                    if i in fun_stat['COG']:
                                        fun_stat['COG'][i] += 1
                                    else:
                                        fun_stat['COG'][i] = 1
                            else:
                                for i in func_cat:
                                    if i in fun_stat['NOG']:
                                        fun_stat['NOG'][i] += 1
                                    else:
                                        fun_stat['NOG'][i] = 1
                            cog_list.append(cog_id)
                    else:
                        if cog_id not in cog_list:
                            if cog_id.startswith('COG'):
                                for i in func_cat:
                                    if i in fun_stat['COG']:
                                        fun_stat['COG'][i] += 1
                                    else:
                                        fun_stat['COG'][i] = 1
                            else:
                                for i in func_cat:
                                    if i in fun_stat['NOG']:
                                        fun_stat['NOG'][i] += 1
                                    else:
                                        fun_stat['NOG'][i] = 1
                            cog_list.append(cog_id)
                        else:
                            if cog_id.startswith('COG'):
                                for i in func_cat:
                                    fun_stat['COG'][i] += 1
                            else:
                                for i in func_cat:
                                    fun_stat['NOG'][i] += 1
                        first_gene = gene
                    # get fun_stat end
            # write whole gene_cog_summary.xls
            s.write('#Total seqs with COG/NOG:' + str(self.totalseq) + '\n')
            s.write('#Tpye\tfunctional_categories\tCOG\tNOG\n')
            for thekey in ['INFORMATION STORAGE AND PROCESSING', 'CELLULAR PROCESSES AND SIGNALING', 'METABOLISM', 'POORLY CHARACTERIZED']:
                for g in self.func_type[thekey]:
                    detail = self.func_decs[g]
                    category = '[' + g + ']' + ' ' + detail
                    try:
                        cogcount = fun_stat['COG'][g]
                    except KeyError:
                        cogcount = 0
                    try:
                        nogcount = fun_stat['NOG'][g]
                    except KeyError:
                        nogcount = 0
                    s.write(thekey + '\t' + category + '\t' + str(cogcount) + '\t' + str(nogcount) + '\n')

    def stats(self, cog_list, gene_list, gene_cog_list, cog_table, out_dir):
        self.get_gene_cog_list(cog_list=cog_list, gene_list=gene_list, outpath=gene_cog_list)
        self.get_gene_summary(cog_table=cog_table, gene_list=gene_list, out_dir=out_dir)
