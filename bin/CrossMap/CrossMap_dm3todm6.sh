


##download the liftover files:

wget https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz

#convert to dm6:
for sample in input roX1 roX2
do
CrossMap.py bigwig dm3ToDm6.over.chain.gz dm3_MEL-${sample}.bw CrossMap_dm6_${sample}.bw
done

