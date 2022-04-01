mkdir ${outdir}/Result_Reporter
cd ${outdir}/Result_Reporter
mkdir -p 0.Data 1.QualityControl
mkdir -p 2.Annotation/{2.1.GO/,2.2.KEGG/,2.3.${cog}/}
mkdir -p 3.DiffExpAnalysis/{3.1.Statistics/{Venn/,Volcano/,Scatter/},3.2.Cluster/{All_samples/,Grouped/},3.3.GO/{Annotation/,Enrichment/},3.4.KEGG/{Annotation/,Enrichment/},3.5.Network,3.6.Ipath}

cp ${rawdata}/*.xlsx 0.Data
cd ${outdir}/Result_Reporter/0.Data
python ${bin}/add_info_to_proteinxlsx.py protein.xlsx ${go_result}/GO.level.xls ${kegg_result}/pathways/pathway_table.xls ${diff_result} ${rawdata}/description.txt protein_merge_info.xls
cd ${outdir}/Result_Reporter

cp ${qc_dir}/dMass.pdf 1.QualityControl
cp ${qc_dir}/Peptide_length_distribution.* 1.QualityControl
cp ${qc_dir}/Peptide_number_distribution.* 1.QualityControl
cp ${qc_dir}/Protein_coverage_distribution.* 1.QualityControl
cp ${qc_dir}/Protein_information.* 1.QualityControl
cp ${qc_dir}/Protein_molecular_weight_distribution.* 1.QualityControl

cp ${go_result}/GO.list 2.Annotation/2.1.GO
cp ${go_result}/GO.level.xls 2.Annotation/2.1.GO
cp ${go_result}/GO.level2.bar.pdf 2.Annotation/2.1.GO
cp ${go_result}/GO.level2.pie.pdf 2.Annotation/2.1.GO

cp ${kegg_result}/pathway.txt 2.Annotation/2.2.KEGG
cp ${kegg_result}/pathway.top20.pdf 2.Annotation/2.2.KEGG
cp ${kegg_result}/KEGG_brite.* 2.Annotation/2.2.KEGG
cp -r ${kegg_result}/pathways 2.Annotation/2.2.KEGG

cp ${cog_result}/${cog}.list 2.Annotation/2.3.${cog}
cp ${cog_result}/${cog}.classification.xls 2.Annotation/2.3.${cog}
cp ${cog_result}/${cog}.class.catalog.xls 2.Annotation/2.3.${cog}
cp ${cog_result}/${cog}.class.catalog.pdf 2.Annotation/2.3.${cog}

cp ${diff_result}/all_diff_up_down.xls 3.DiffExpAnalysis/3.1.Statistics
cp ${diff_result}/{*.diff.exp.xls,*volcano.pdf,*.up.list,*.down.list,*.DE.list} 3.DiffExpAnalysis/3.1.Statistics/Volcano
cp ${diff_result}/*scatter.pdf 3.DiffExpAnalysis/3.1.Statistics/Scatter

cp -r ${venn_result}/*combinations 3.DiffExpAnalysis/3.1.Statistics/Venn/
rm 3.DiffExpAnalysis/3.1.Statistics/Venn/*combinations/*log
cp -r ${cluster_result}/All_samples/{*Heatmap.pdf,*trendlines*.pdf,subclusters_*} 3.DiffExpAnalysis/3.2.Cluster/All_samples
cp -r ${cluster_result}/Grouped/{*Heatmap.pdf,*trendlines*.pdf,subclusters_*} 3.DiffExpAnalysis/3.2.Cluster/Grouped

cp ${annot_enrich_result}/{*.up-down.pdf,*.GO.level.xls} 3.DiffExpAnalysis/3.3.GO/Annotation
cp -r ${annot_enrich_result}/{*.go_enrichment*pdf,*.go_enrichment.xls} 3.DiffExpAnalysis/3.3.GO/Enrichment
cp -r ${annot_enrich_result}/*.GOplot.* 3.DiffExpAnalysis/3.3.GO/Enrichment
rm 3.DiffExpAnalysis/3.3.GO/Enrichment/plot*
rm 3.DiffExpAnalysis/3.3.GO/Enrichment/*.r

cp -r ${annot_enrich_result}/*.paths 3.DiffExpAnalysis/3.4.KEGG/Annotation
cp -r ${annot_enrich_result}/{*.kegg_enrichment*pdf,*.kegg_enrichment.xls} 3.DiffExpAnalysis/3.4.KEGG/Enrichment
cp -r ${annot_enrich_result}/*.KEGGplot.* 3.DiffExpAnalysis/3.4.KEGG/Enrichment
rm 3.DiffExpAnalysis/3.4.KEGG/Enrichment/plot*

cp ${diff_result}/*.DE.list 3.DiffExpAnalysis/3.5.Network

cp ${ipath_result}/*.Ipath 3.DiffExpAnalysis/3.6.Ipath
cp ${ipath_picture}/* 3.DiffExpAnalysis/3.6.Ipath

cd 3.DiffExpAnalysis/3.1.Statistics/Volcano/
for i in `ls *.diff.exp.xls |awk -F '.diff.exp.xls' '{print $1}'`
do
    ${bin}/get_diff_table_detail.py -i $i.diff.exp.xls -d '${go_result}/GO.list;${kegg_result}/pathway.txt;${cog_result}/${cog}.list;${rawdata}/description.txt' -l 'GO;KEGG;${cog};Description' -o $i.diff.exp.detail.xls

done
cd ${outdir}/Result_Reporter







cd ${outdir}
mkdir Reporter
cd Reporter
cp -r ${outdir}/Result_Reporter ./

cd ${outdir}/Reporter/Result_Reporter/3.DiffExpAnalysis/3.3.GO/Enrichment
for i in *go_enrichment.xls
do
    ${bin}/get_extra_FDR.py $i p_value $i.tmp
    mv $i.tmp $i
done

cd ${outdir}/Reporter/Result_Reporter/3.DiffExpAnalysis/3.4.KEGG/Enrichment
for i in *kegg_enrichment.xls
do
    ${bin}/get_extra_FDR.py $i p_value $i.tmp
    mv $i.tmp $i
done
cd ${outdir}/Reporter/
cp ${rawdata}/M1.docx ./
${bin}/${MJ} -m M1.docx -f Result_Reporter -org ${cog} -up ${up} -dump ${dup}


cd ${outdir}
mkdir Customer
cd Customer
cp -r ${outdir}/Reporter/* ./
cd ${outdir}/Customer/Result_Reporter/2.Annotation/2.2.KEGG
${bin}/merge_sheets_2_xlsx.py -f 'pathway.txt,KEGG_brite.A.xls,KEGG_brite.B.xls,KEGG_brite.C.xls' -l 'pathway,KEGG_brite.A,KEGG_brite.B,KEGG_brite.C' -o KEGG
#rm pathway.txt *xls
cd ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.1.Statistics/Volcano
for i in `ls *.diff.exp.xls |awk -F '.diff.exp.xls' '{print $1}'`
do
    ${bin}/merge_sheets_2_xlsx.py -f "$i.diff.exp.detail.xls,$i.diff.exp.xls,$i.DE.list,$i.up.list,$i.down.list" -l "T-TEST-DETAIL,T-TEST,DE-LIST,UP-LIST,DOWN-LIST" -o $i
    mv $i.xlsx ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.1.Statistics/$i.xlsx
done
rm *.list *exp.xls

cd ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/Grouped/
for i in subclusters_*
do
    cd $i
    ${bin}/merge_uncertainsheets_2_xlsx.py -f ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/Grouped/$i -p xls -o $i
    mv $i.xlsx ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/Grouped/
    cd ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/Grouped/
    rm -rf $i

done

cd ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/All_samples/
for i in subclusters_*
do
    cd $i
    ${bin}/merge_uncertainsheets_2_xlsx.py -f ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/All_samples/$i -p xls -o $i
    mv $i.xlsx ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/All_samples/
    cd ${outdir}/Customer/Result_Reporter/3.DiffExpAnalysis/3.2.Cluster/All_samples/
    rm -rf $i

done

cd ${outdir}