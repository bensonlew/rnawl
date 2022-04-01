#/bin/sh
species=${1}
database=${2}
web=${3}
wget -O $species\_gene.txt 'http://'$web'.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "'$species'_gene_ensembl" interface = "default" >
                <Attribute name = "ensembl_gene_id" />
                <Attribute name = "ensembl_gene_id_version" />
                <Attribute name = "ensembl_transcript_id" />
                <Attribute name = "ensembl_transcript_id_version" />
                <Attribute name = "3utr" />
        </Dataset>
</Query>'
