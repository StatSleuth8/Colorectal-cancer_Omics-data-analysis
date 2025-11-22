###############################################part01_data_preprocessing###############################################
### Amplicon sequencing upstream analysis (QIIME2-2021.10 + vsearch)
### Author: shuai yang
### Purpose: 16S amplicon raw data → ASV table
### Required: QIIME2 ≥ 2021.10, SILVA 138.1 classifier
###for questions, please contact ShuaiYang (shuai_yang8@163.com)


mydir=/home/yang
threads=32

CRCmap=${mydir}/CRCmapping.txt             # sample list
raw_CRC=${mydir}/rawdata/CRC               # raw fastq
clean_CRC=${mydir}/cleandata/CRC           # clean fastq
temp_CRC=${mydir}/temp/CRC                 # tmp folder
result_CRC=${mydir}/result/CRC             # final results

mkdir -p ${clean_CRC} ${temp_CRC} ${result_CRC}

###############################################
### step 1 fastp quality trimming (paired-end)
###############################################
echo "### Step 1: fastp trimming"

time for i in `tail -n +2 ${CRCmap}`; do
fastp \
  -i  ${raw_CRC}/${i}_1.fastq \
  -I  ${raw_CRC}/${i}_2.fastq \
  -o  ${clean_CRC}/${i}_1.clean.fastq \
  -O  ${clean_CRC}/${i}_2.clean.fastq \
  -h  ${clean_CRC}/${i}.html \
  -j  ${clean_CRC}/${i}.json \
  -f 26 -F 26 \
  -t 0  -T 40 \
  -l 20 \
  -5 20 -3 20 \
  -M 20 -W 5 \
  --thread ${threads}
done


###############################################
### step 2 vsearch merge paired-end reads
###############################################
echo "### Step 2: vsearch merging"

mkdir -p ${temp_CRC}/merged

for i in `tail -n +2 ${CRCmap}`; do
vsearch \
  --fastq_mergepairs  ${clean_CRC}/${i}_1.clean.fastq \
  --reverse           ${clean_CRC}/${i}_2.clean.fastq \
  --fastq_allowmergestagger \
  --fastqout          ${temp_CRC}/merged/${i}.merged.fastq \
  --relabel           ${i}.
done


###############################################
### step 3 concatenate all merged reads
###############################################
echo "### Step 3: concatenate merged reads"

cat ${temp_CRC}/merged/*.merged.fastq > ${temp_CRC}/all.merged.fastq


###############################################
### step 4 filtering
###############################################
echo "### Step 4: quality filtering"

vsearch \
  --fastx_filter ${temp_CRC}/all.merged.fastq \
  --fastq_maxee_rate 0.01 \
  --fastaout ${temp_CRC}/filtered.fa


###############################################
### step 5 dereplication
###############################################
echo "### Step 5: dereplication"

vsearch \
  --derep_fulllength ${temp_CRC}/filtered.fa \
  --minuniquesize 3 \
  --sizeout \
  --output ${temp_CRC}/uniques.fa


###############################################
### step 6 UNOISE denoising to generate ASVs
###############################################
echo "### Step 6: ASV inference (UNOISE3)"

vsearch \
  --cluster_unoise ${temp_CRC}/uniques.fa \
  --centroids ${temp_CRC}/ASVs.fa \
  --sizeout \
  --relabel ASV_


###############################################
### step 7 chimera removal
###############################################
echo "### Step 7: chimera removal"

vsearch \
  --uchime3_denovo ${temp_CRC}/ASVs.fa \
  --nonchimeras ${result_CRC}/ASVs.nonchimeric.fa


###############################################
### step 8 OTU table (ASV table)
###############################################
echo "### Step 8: ASV/OTU table construction"

vsearch \
  --usearch_global ${temp_CRC}/filtered.fa \
  --db ${result_CRC}/ASVs.nonchimeric.fa \
  --id 0.97 \
  --threads ${threads} \
  --otutabout ${result_CRC}/otutab.txt


###############################################
### step 9 kraken2 taxonomy annotation
###############################################
echo "### Step 9: taxonomy annotation"

kraken2db=${mydir}/db/kraken2_silvadb
mkdir -p ${result_CRC}/abundance

kraken2 \
  --db ${kraken2db} \
  --threads ${threads} \
  --output ${result_CRC}/ASV.kraken2 \
  --report ${result_CRC}/ASV.kreport2 \
  ${result_CRC}/ASVs.nonchimeric.fa


# extract ASV–species mapping
awk -F "\t" '$1=="C" && $3!=0 {print $2"\t"$3}' \
  ${result_CRC}/ASV.kraken2 \
  | sed '1i\asvid\tspid' \
  > ${result_CRC}/asvclass.txt


###############################################
### step 10 compute relative abundance (python)
###############################################
echo "### Step 10: abundance calculation"

python3 ${mydir}/script/abundance.py \
  -asv   ${result_CRC}/otutab.txt \
  -class ${result_CRC}/asvclass.txt \
  -db    ${mydir}/db/dbpy.txt \
  -out   ${result_CRC}/abundance


