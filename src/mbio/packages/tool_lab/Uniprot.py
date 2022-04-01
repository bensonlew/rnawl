import requests
import time
import random
import urllib.request
import argparse

class uniprot_subcellular(object):
    def __init__(self, uniprot_id, output, method, uniprot_trembl):
        self.uniprot_id = uniprot_id
        self.output = output
        self.method = method
        self.uniprot_trembl = uniprot_trembl

    def request_uniprot(self):
        time.sleep(3 + random.random())
        uniprot_id_list = self.uniprot_id.strip().split(';')
        new_uniprot_id_list = list()
        for i in uniprot_id_list:
            new_id = 'id:' + i
            new_uniprot_id_list.append(new_id)
        new_uniprot_id_str = '%20OR%20'.join(new_uniprot_id_list)
        uniprot_api_url = "https://www.uniprot.org/uniprot/?query={}&format=xlsx&columns=id,protein%20names,genes," \
                          "organism,feature(INTRAMEMBRANE),comment(SUBCELLULAR%20LOCATION),feature(TOPOLOGICAL%20DOMAIN)," \
                          "feature(TRANSMEMBRANE)&sort=score".format(new_uniprot_id_str)
        opener = urllib.request.build_opener()
        opener.addheaders = [('User-agent',
                              'Mozilla/5.0 (X11; Linux x86_64)'
                              ' AppleWebKit/537.36 (KHTML, like Gecko)'
                              ' Chrome/68.0.3440.106 Safari/537.36')]
        urllib.request.install_opener(opener)
        urllib.request.urlretrieve(uniprot_api_url, self.output)

    def run(self):
        if self.method == 'uniprot_api':
            self.request_uniprot()
        elif self.method == 'local':
            pass

if __name__ == '__main__':
    # def __init__(self, specie, list_path, vsstring, identity, useblast):
    parser = argparse.ArgumentParser(description="use the uniprot net api to get the SUBCELLULAR LOCATION")
    parser.add_argument("-uniprot_id", type=str, help="uniprot list")
    parser.add_argument("-output", type=str, help="the output file")
    parser.add_argument("-method", type=str, help="the method")
    parser.add_argument("-uniprot_trembl", type=str, help='the filter para to filter the string blast result')


    args = parser.parse_args()
    STRING = uniprot_subcellular(args.uniprot_id, args.output, args.method, args.uniprot_trembl)
    STRING.run()



