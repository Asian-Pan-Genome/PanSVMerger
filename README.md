# PanSVMerger-A-Pipeline-for-Merging-Multi-Allele-Structural-Variants
# Table of Contents

[About](#About)
   - [Pipeline Description](#Pipeline-Description)
   - [why PanSVMerge ?](#why-PanSVMerge)
   - [Requirements and Install](#Requirements-and-Install)

[Usage](#usage)
   - [Prepare non-overlapping multiallelic VCF file](#Prepare-non-overlapping-multiallelic-VCF-file)
   - [Merge SVs](#Merge-SVs)
   - [Construct benchmarking](#construct-benchmarking)
   - [Compare with benchmarking](#Compare-with-benchmarking)

[Support](#support)
   - [Contacts](#Contacts)
  
## About

PanSVMerger is a bioinformatics pipeline designed to merge multi-allele structural variants (SVs) located at the same genomic position based on their sequence similarity and length. This tool aims to provide a reliable and efficient method for consolidating SV data, which is crucial for accurate genomic analysis and interpretation.

## Pipeline-Description
- vcfwave:  Decomposes complex SVs from the original MC VCF file, resolving nested bubble structures by realigning reference and alternate alleles to parse out primitive alleles
- vcfcreatemulti: Reconstructs non-overlapping multi-allelic sites by removing nested alleles
- vsearch: Clusters alternative alleles based on sequence similarity, further divides clusters when sequence length differences exceed 50bp
- genotype updating: Reassigns genotypes based on cluster membership

## why-PanSVMerger
When working with Asian-Pan-Graph, we observed a high proportion of multi-allelic sites at single genomic loci. While existing tools like:
- [panpop](https://github.com/starskyzheng/panpop)
- [Truvari](https://github.com/ACEnglish/truvari)

could partially merge these alleles, they still left a significant number of multi-allelic sites unresolved. 

To address this, we developed **PanSVMerger** - a specialized method that:
1. Leverages the structure of APG-pan graphs
2. Implements sequence similarity-based clustering
3. Incorporates length difference metrics
4. Provides optimized allele merging for complex variation patterns

This approach significantly reduces fragmented multi-allelic representations while maintaining biological accuracy of the variation.

## Requirements-and-Install
- [vcflib](https://github.com/vcflib/vcflib)
- [pysam](https://github.com/pysam-developers/pysam) (stable version)
- vsearch

## Usage
### Prepare-non-overlapping-multiallelic-VCF-file
Due to the complexity of variants called from graph, particulally those complex regions with nested bubbles, VCF file genereated by MC still contains tricky and problematic sites although after "poping" using `vcfbubble`. [As Garrison et al. suggested](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md), we should first realign reference and alternate alleles to parse out the original "primitive" alleles into multiple records and next put them together again, which would be provided for `PanSVMerger` to consolidate SVs based on similarity and length (see below).
```
vcfwave -t16 MC.vcf.gz (-L 100000) > MC.wave.vcf ## you can ignore SVs longer than 100Kb
bcftools norm --threads 16 -m- MC.wave.vcf -Ov -o MC.wave.bi.vcf ## there are still few multiallelic sites
vcfcreatemulti MC.wave.bi.vcf > MC.wave.bi.createmulti.vcf
```

### Merge-SVs
After removing and reconstructing overlapping and nested sites, we would see multiple alleles in one complex SV record. However, they might be slightly different and should be considered as the same one, which would have great impacts on downstream analysis. Therefore, we'd like to consolidate and merge them to confident consensus. [See our paper for detailed methods.](https://github.com/tingting100/PanSVMerger#citation)
```
python SV.merge_variants.py MC.wave.bi.createmulti.vcf MC.wave.bi.createmulti.merge -t 16 ## if input compressed vcf file, please index first
bgzip -@16 MC.wave.bi.createmulti.merge.vcf ## if input compressed file, `SV.merge_variants.py` will output compressed merged vcf file
bcftools index -t --threads 16 MC.wave.bi.createmulti.merge.vcf.gz
```

### Construct-Benchmarking

This workflow integrates both assembly-based and long-read SV callers, followed by comprehensive merging to generate a high-confidence SV set.

#### 1. Assembly-based SV Calling

```bash
mkdir ./long_reads.SV/ && cd ./long_reads.SV/
# SVIM-asm
~/bin/minimap2/2.17/minimap2 -a -x asm5 --cs -r2k -t 4 genome.fa query.fa > alignments.sam
~/bin/samtools sort -m4G -@4 -o alignments.sorted.bam alignments.sam
~/bin/samtools index alignments.sorted.bam
~/bin/svim-asm diploid ./ alignments.sorted.bam genome.fa

# SyRI
~/bin/mummer-4.0.0/bin/nucmer --prefix=sample -l 100 genome.fa query.fa
~/bin/mummer-4.0.0/delta-filter -l 100 -i 90 -1 sample.delta > sample.delta.filter
~/bin/mummer-4.0.0/bin/show-coords -THrd sample.delta.filter > sample.filtered.coords

python3 ~/bin/syri-1.5/syri/bin/syri -c sample.filtered.coords -r genome.fa -q query.fa \
  -d sample.delta.filter -k --prefix sample -s ~/bin/show-snps --lf syri.log
```

#### 2. Long-read SV Calling
```
mkdir ./long_reads.SV/ && cd ./long_reads.SV/
# Sniffles
~/bin/samtools calmd --threads 6 -b HiFi_reads.sort.head.bam genome.fa > HiFi_reads.sort.baq.bam
~/bin/sniffles-core-1.0.12/sniffles -m HiFi_reads.sort.baq.bam -v HiFi_reads.sniffles.vcf \
  --report_BND --skip_parameter_estimation --min_support 3

# SVIM
~/bin/svim alignment ./ HiFi_reads.sort.head.bam genome.fa

# CuteSV
~/bin/cuteSV HiFi_reads.sort.head.bam genome.fa HiFi_reads.cuteSV.vcf ./ \
  --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
  --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 8 --genotype

# PBSV
~/bin/pbsv discover HiFi_reads.sort.head.bam HiFi_reads.svsig.gz
tabix -c '#' -s 3 -b 4 -e 4 HiFi_reads.svsig.gz
~/bin/pbsv call --ccs genome.fa HiFi_reads.svsig.gz HiFi_reads.pbsv.vcf
```

#### 3. Merging Results
```bash
# Combine VCFs from all callers
ls ./long_reads.SV/*.vcf > combine.txt
ls ./assembly_based.SV/*.vcf >> combine.txt

# Merge with SURVIVOR (max distance 1000bp, require 2 callers support)
~/bin/SURVIVOR merge combine.txt 1000 2 1 1 0 50 sampleN.truth.sv.vcf.gz
```

#### Output Files
```bash
sampleN.truth.sv.vcf.gz: Final merged structural variant calls for each individual

./long_reads.SV/: Individual long-reads-based SV VCFs

./assembly_based.SV/: Individual Assembly-based SV VCFs
```
### Compare-with-benchmarking
```
# 01_subsample_svs.sh
vcftools --gzvcf MC.APGpangraph.vcf.gz --indv sampleN --recode --stdout |bgzip -c > sampleN.sv.vcf.gz
tabix -p vcf sampleN.sv.vcf.gz

# 02_filter_large_svs.sh
vcfbub -l 0 -a 100000 -i sampleN.vcf.gz -o sampleN.filtered.vcf.gz

# 03_decompose_alleles.sh
vcfwave -I 1000 -i sampleN.filtered.vcf.gz -o sampleN.wave.vcf.gz

# 04_split_biallelic.sh
bcftools norm -m -any -i sampleN.wave.vcf.gz -o sampleN.biallelic.vcf.gz

# 05_run_truvari.sh
truvari bench -b sampleN.truth.sv.vcf.gz -c sampleN.biallelic.vcf.gz -o truvari_results
```


### Contacts
(yangting@genomics.cn)
(yangchentao@genomics.cn)
(quanyu_chen@outlook.com)
