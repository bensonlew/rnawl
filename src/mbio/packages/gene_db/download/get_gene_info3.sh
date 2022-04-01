#/bin/sh
species=${1}
database=${2}
web=${3}
wget -O $species\_gene.txt 'http://'$web'.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "'$database'" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1"><Dataset name = "'$species'" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "ensembl_peptide_id" /><Attribute name = "description" /><Attribute name = "chromosome_name" /><Attribute name = "start_position" /><Attribute name = "end_position" /><Attribute name = "strand" /><Attribute name = "transcript_start" /><Attribute name = "transcript_end" /><Attribute name = "transcription_start_site" /><Attribute name = "transcript_length" /><Attribute name = "gene_biotype" /><Attribute name = "transcript_biotype" /></Dataset></Query>'
