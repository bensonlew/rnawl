# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import json
import re

parser = OptionParser(description='Filter miRNA.site.json by given mature.fa')
parser.add_option('-m', '--mature', dest='mature' ,help='given mature.fa used for filtering site in json file')
parser.add_option('-s', '--site', dest='site', help='json file of site relationship between hairpin and mature')
parser.add_option('-o', '--output', dest='output', help='path of result json file')
(opts, args) = parser.parse_args()

def main(mature_fa, mirna_site_json, output_json):
    mature_name_set = extract_mature_name(mature_fa)
    with open(mirna_site_json) as f:
        hairpin_mature_site = json.load(f)
    output_dict = dict()
    for hairpin, mature_site in hairpin_mature_site.iteritems():
        for mature in mature_site:
            if mature.lower() in mature_name_set:
                output_dict[hairpin] = mature_site
    with open(output_json, 'w') as f:
        json.dump(output_dict, f)

def extract_mature_name(mature_fa):
    mature_name_list = list()
    with open(mature_fa) as f:
        for line in f:
            if line[0] == '>':
                mature_name_line = line.strip()
                m = re.match(r'>(\S+)', mature_name_line)
                if m:
                    mature_name = m.group(1)
                    mature_name_list.append(mature_name.lower())
    return set(mature_name_list)

if __name__ == '__main__':
    if opts.mature and opts.site and opts.output:
        main(opts.mature, opts.site, opts.output)
    else:
        parser.print_help()
