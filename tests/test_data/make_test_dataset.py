import os
from os import system


# config
CHR = "chr9"
START = 4000000
END = 6000000
REGION = f"{CHR}:{START}-{END}"

REGION_MARK = REGION.replace("-", "_").replace(":", "_")
r_ = REGION[3:]


def extract_bigwig(source, target, chr, start, end):
    # https://bioinfocore.com/blogs/extract-subset-bigwig-file-for-a-given-genomic-region/
    import pyBigWig
    bw_in = pyBigWig.open(source)
    bw_out = pyBigWig.open(target, 'w')
    bw_out.addHeader([(chr, bw_in.chroms()[chr])])

    for x in bw_in.intervals(chr, start, end):
        bw_out.addEntries([chr], [x[0]], ends=[x[1]], values=[x[2]])

    bw_in.close()
    bw_out.close()


data_path = "/mnt/d/Data/bigwig"
files = os.listdir(data_path)
for f in files:
    name, type = f.split(".")
    target = f"{type}_{REGION_MARK}_{name}.{type}"
    print(target)
    extract_bigwig(f"{data_path}/{f}", target, CHR, START, END)


if False:
    import pysam
    # source files
    sources = {
        "cool": "/home/nanguage/DATA1/LD/differentiation/HiC/Ra-MseI-4sample-combin.multi.cool::resolutions/5000",
        "bam": "/home/nanguage/DATA2/public/ENCODE/ENCFF046HDL.sorted.bam",
        "chr_size": "/home/nanguage/S/Bioinfo/bedtools2/genomes/human.hg19.genome",
        "bed": "/home/nanguage/Test/BioInfo/CoolBox/demo/data/preprocessed/refGene.sorted.bed",
        "arcs": "/home/nanguage/Test/BioInfo/CoolBox/demo/data/preprocessed/K562_MbolI_looplist.arcs",
        "gtf": "~/DATA1/Genomes/Homo_sapiens/hg19/Homo_sapiens.GRCh37.75.chr.gtf.gz",
    }

    # chr size
    system(f"cp {sources['chr_size']} .")

    # cooler
    system(f"coolclip {sources['cool']} cool_{REGION_MARK}.cool {r_}")
    system(f"cooler zoomify --balance cool_{REGION_MARK}.cool -o cool_{REGION_MARK}.mcool")
    system(f"rm cool_{REGION_MARK}.cool")

    # bigwig and bam
    out_bam = f"bam_{REGION_MARK}.bam"
    with pysam.AlignmentFile(sources['bam'], 'rb') as bam_in:
        with pysam.AlignmentFile(out_bam, 'wb', template=bam_in) as bam_out:
            for read in bam_in.fetch(CHR, START, END):
                bam_out.write(read)
    bw = f"bigwig_{REGION_MARK}.bw"
    system(f"samtools index {out_bam}")
    system(f"bamCoverage --bam {out_bam} -o {bw} --binSize 1000")

    # bed
    system(f"cp {sources['bed']} .")
    bed_source = os.path.basename(sources['bed'])
    system(f"bgzip {bed_source}")
    system(f"tabix -p bed {bed_source}.gz")
    system(f"tabix {bed_source}.gz {REGION} > bed_{REGION_MARK}.bed")
    system(f"rm {bed_source}.gz {bed_source}.gz.tbi")

    # arcs
    system(f"cp {sources['arcs']} .")
    arcs_source = os.path.basename(sources['arcs'])
    system(f"sortBed -i {arcs_source} > {arcs_source}.sorted")
    system(f"bgzip {arcs_source}.sorted")
    system(f"tabix -p bed {arcs_source}.sorted.gz")
    system(f"tabix {arcs_source}.sorted.gz {r_} > arcs_{REGION_MARK}.arcs")
    system(f"rm {arcs_source} {arcs_source}.sorted.gz {arcs_source}.sorted.gz.tbi")

    # gtf
    system(f"tabix {sources['gtf']} {REGION} > gtf_{REGION_MARK}.gtf")
