# -*- coding: utf-8 -*-
# __author__ = 'shenghe'


def format_human_result(f_path, outfile, rich_type="Pathway"):
    with open(f_path) as f, open(outfile, 'w') as w:
        head = f.readline().strip().split('\t')[2:]
        head = [i.replace('-hit-keg-mpm-cop-nul-nve-nve-xpe', '') for i in head]
        # richness = None
        for i in f:
            if i[: 8] == "Richness":
                # richness = i.strip().split('\t')[2:]
                break
        else:
            raise Exception("错误：humann注释结果没有Richness行")
        w.write(rich_type + '\t' + rich_type + '_Desc' + '\t' + '\t'.join(head) + '\n')
        # w.write('Richness\tRichness\t' + '\t'.join(richness) + '\n')
        for i in f:
            w.write(i)


if __name__ == '__main__':
    myfile = 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\human_profile\\04a-hit-keg-mpm-cop-nul-nve-nve-xpe.txt'
    outf = 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\human_profile\\Shenghe.modue.profile.txt'
    format_human_result(myfile, outf, rich_type="Module")
