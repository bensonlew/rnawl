import requests, sys

gene_id = sys.argv[1]
server = "https://rest.ensembl.org"
ext = "/xrefs/id/{}?all_levels=1;external_db=EntrezGene".format(gene_id)

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
if 'primary_id' in decoded[0]:
    return decoded[0]['primary_id']
else:

