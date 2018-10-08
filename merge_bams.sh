#!/bin/bash

STAR_OUT="/data/ndcn-rnaseq/ndcn0569/RNAseq/star_out"

mkdir /data/ndcn-rnaseq/ndcn0569/RNAseq/star_out/merged_bams

MERGED_OUT="/data/ndcn-rnaseq/ndcn0569/RNAseq/star_out/merged_bams"



for FILE in $STAR_OUT/*_201106Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_136.txt
done

for FILE in $STAR_OUT/*_201191Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_109.txt
done

for FILE in $STAR_OUT/*_202118Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0136.txt
done

for FILE in $STAR_OUT/*_202179Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0109.txt
done

for FILE in $STAR_OUT/*_203130Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_151.txt
done

for FILE in $STAR_OUT/*_203167Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_181.txt
done

for FILE in $STAR_OUT/*_204142Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0151.txt
done

for FILE in $STAR_OUT/*_204155Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0181.txt
done

for FILE in $STAR_OUT/*_205143Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_169.txt
done

for FILE in $STAR_OUT/*_205154Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_134.txt
done

for FILE in $STAR_OUT/*_206131Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0169.txt
done

for FILE in $STAR_OUT/*_206166Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0134.txt
done

for FILE in $STAR_OUT/*_207119Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_168.txt
done

for FILE in $STAR_OUT/*_207178Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_147.txt
done

for FILE in $STAR_OUT/*_208107Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0168.txt
done

for FILE in $STAR_OUT/*_208190Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0147.txt
done

for FILE in $STAR_OUT/*_209190Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_192.txt
done

for FILE in $STAR_OUT/*_210178Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0192.txt
done

for FILE in $STAR_OUT/*_257101Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_191.txt
done

for FILE in $STAR_OUT/*_257196Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_131.txt
done

for FILE in $STAR_OUT/*_258113Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0191.txt
done

for FILE in $STAR_OUT/*_258184Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0131.txt
done

for FILE in $STAR_OUT/*_259125Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_43.txt
done

for FILE in $STAR_OUT/*_259172Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_189.txt
done

for FILE in $STAR_OUT/*_260137Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_4300.txt
done

for FILE in $STAR_OUT/*_260160Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0189.txt
done

for FILE in $STAR_OUT/*_261148Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_162.txt
done

for FILE in $STAR_OUT/*_261149Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_36.txt
done

for FILE in $STAR_OUT/*_262136Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0162.txt
done

for FILE in $STAR_OUT/*_262161Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_3600.txt
done

for FILE in $STAR_OUT/*_263124Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_175.txt
done

for FILE in $STAR_OUT/*_263173Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_167.txt
done

for FILE in $STAR_OUT/*_264112Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0175.txt
done

for FILE in $STAR_OUT/*_264185Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0167.txt
done

for FILE in $STAR_OUT/*_265102Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_178.txt
done

for FILE in $STAR_OUT/*_265195Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_160.txt
done

for FILE in $STAR_OUT/*_266114Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0178.txt
done

for FILE in $STAR_OUT/*_266183Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0160.txt
done

for FILE in $STAR_OUT/*_267126Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_188.txt
done

for FILE in $STAR_OUT/*_267171Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_164.txt
done

for FILE in $STAR_OUT/*_268138Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0188.txt
done

for FILE in $STAR_OUT/*_268159Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0164.txt
done

for FILE in $STAR_OUT/*_269147Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_143.txt
done

for FILE in $STAR_OUT/*_269150Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_59.txt
done

for FILE in $STAR_OUT/*_270135Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0143.txt
done

for FILE in $STAR_OUT/*_270162Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_5900.txt
done

for FILE in $STAR_OUT/*_271123Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_106.txt
done

for FILE in $STAR_OUT/*_271174Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_149.txt
done

for FILE in $STAR_OUT/*_272111Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0106.txt
done

for FILE in $STAR_OUT/*_272186Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0149.txt
done

for FILE in $STAR_OUT/*_273103Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_133.txt
done

for FILE in $STAR_OUT/*_273194 Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_119.txt
done

for FILE in $STAR_OUT/*_274115Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0133.txt
done

for FILE in $STAR_OUT/*_274182Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0119.txt
done

for FILE in $STAR_OUT/*_275127Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_159.txt
done

for FILE in $STAR_OUT/*_275170Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_190.txt
done

for FILE in $STAR_OUT/*_276139Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0159.txt
done

for FILE in $STAR_OUT/*_276158Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0190.txt
done

for FILE in $STAR_OUT/*_277146Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_166.txt
done

for FILE in $STAR_OUT/*_277151Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_128.txt
done

for FILE in $STAR_OUT/*_278134Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0166.txt
done

for FILE in $STAR_OUT/*_278163Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0128.txt
done

for FILE in $STAR_OUT/*_279122Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_163.txt
done

for FILE in $STAR_OUT/*_279175Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_148.txt
done

for FILE in $STAR_OUT/*_280110Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0163.txt
done

for FILE in $STAR_OUT/*_280187Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0148.txt
done

for FILE in $STAR_OUT/*_281104Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_155.txt
done

for FILE in $STAR_OUT/*_281193Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_142.txt
done

for FILE in $STAR_OUT/*_282116Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0155.txt
done

for FILE in $STAR_OUT/*_282181Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0142.txt
done

for FILE in $STAR_OUT/*_283128Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_79.txt
done

for FILE in $STAR_OUT/*_283169Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_170.txt
done

for FILE in $STAR_OUT/*_284140Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_7900.txt
done

for FILE in $STAR_OUT/*_284157Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0170.txt
done

for FILE in $STAR_OUT/*_285145Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_121.txt
done

for FILE in $STAR_OUT/*_285152Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_176.txt
done

for FILE in $STAR_OUT/*_286133Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0121.txt
done

for FILE in $STAR_OUT/*_286164Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0176.txt
done

for FILE in $STAR_OUT/*_287121Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_141.txt
done

for FILE in $STAR_OUT/*_287176Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_172.txt
done

for FILE in $STAR_OUT/*_288109Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0141.txt
done

for FILE in $STAR_OUT/*_288188Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0172.txt
done

for FILE in $STAR_OUT/*_289105Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_137.txt
done

for FILE in $STAR_OUT/*_289192Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_177.txt
done

for FILE in $STAR_OUT/*_290117Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0137.txt
done

for FILE in $STAR_OUT/*_290180Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0177.txt
done

for FILE in $STAR_OUT/*_291129Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_111.txt
done

for FILE in $STAR_OUT/*_291168Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_130.txt
done

for FILE in $STAR_OUT/*_292141Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0111.txt
done

for FILE in $STAR_OUT/*_292156Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0130.txt
done

for FILE in $STAR_OUT/*_293144Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_182.txt
done

for FILE in $STAR_OUT/*_293153Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_146.txt
done

for FILE in $STAR_OUT/*_294132Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0182.txt
done

for FILE in $STAR_OUT/*_294165Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0146.txt
done

for FILE in $STAR_OUT/*_295120Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_138.txt
done

for FILE in $STAR_OUT/*_295177Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_186.txt
done

for FILE in $STAR_OUT/*_296108Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0138.txt
done

for FILE in $STAR_OUT/*_296189Aligned.sortedByCoord.out.bam
do
echo ${FILE} >> $MERGED_OUT/CTS_0186.txt
done

for FILE in $MERGED_OUT/*.txt
do
echo "Merging replicates for" ${FILE}
cat $FILE
BAM_NAME=$(echo $FILE | rev | cut -c 5- | rev)
echo "Creating" $BAM_NAME.bam
/home/ndcn-rnaseq/ndcn0569/bin/samtools merge -b ${FILE} -u $BAM_NAME.bam
echo "Created" $BAM_NAME.bam
/home/ndcn-rnaseq/ndcn0569/bin/samtools sort -n $BAM_NAME.bam -o $BAM_NAME.sorted.bam -O bam -@ 15
echo "Removing unsorted" $BAM_NAME.bam
rm $BAM_NAME.bam
#/home/ndcn-rnaseq/ndcn0569/bin/samtools index $BAM_NAME.sorted.bam
/home/ndcn-rnaseq/ndcn0569/bin/samtools flagstat $BAM_NAME.sorted.bam > $BAM_NAME.flagstats.txt
done


 
