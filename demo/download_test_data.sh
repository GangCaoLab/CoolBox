%%bash

mkdir data

# example .cool file
wget -O data/K562_MbolI_5kb.cool ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-K562-MboI-allreps-filtered.5kb.cool
# .looplist file
wget -O data/K562_MbolI_looplist.txt.gz ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_HiCCUPS_looplist.txt.gz

# bigwig files (ChIP-Seq and RNA-Seq data)
wget -O data/K562_H3K27ac.bigWig  ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733656/suppl/GSM733656_hg19_wgEncodeBroadHistoneK562H3k27acStdSig.bigWig
wget -O data/K562_H3K27me3.bigWig ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733658/suppl/GSM733658_hg19_wgEncodeBroadHistoneK562H3k27me3StdSig.bigWig
wget -O data/K562_H3K4me3.bigWig  ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM945nnn/GSM945297/suppl/GSM945297_hg19_wgEncodeUwHistoneK562H3k04me3StdZnf2c10c5RawRep1.bigWig
wget -O data/K562_RNASeq.bigWig ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM765nnn/GSM765393/suppl/GSM765393_wgEncodeCshlLongRnaSeqK562NucleolusTotalMinusRawSigRep4.bigWig
# ChIA-PET data
wget -O data/K562_chiapet_interaction.txt.gz ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM832nnn/GSM832465/suppl/GSM832465_CHK019M_L2_lane24.ChromatinInteractions.bed.gz
wget -O data/K562_chiapet.bigWig ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM832nnn/GSM832465/suppl/GSM832465_CHK019M_L2_lane24.bw

# RefSeq
wget -O data/refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz