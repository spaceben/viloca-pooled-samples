import subprocess
from pathlib import Path
import glob
from os import path
import pysam


def main(fin, coverage, out_dir, out_bams, merged_bam):
    print(out_bams)
    print(merged_bam)

    for i in glob.glob(path.join(fin, "*.fasta")):
        with open(i, "r") as f:
            freq = float(f.readline().split(":")[1])
        subprocess.run(f"art_illumina -sam --paired -m 320 -s 10 --rndSeed 20 -i {i} \
                    -f {str(round(coverage * freq))} -l 200 -na -o {i.split('.fasta')[0]}", shell=True)

    all_generated_sams = glob.glob(path.join(fin, "*.sam"))
    grps = [int(i.split("sample")[-1].split("-")[0]) for i in all_generated_sams]
    for i in set(grps):
        sams = [j for j, k in zip(all_generated_sams, grps) if k == i]
        fout = path.join(out_dir, 'sample'+str(i)+'.bam')

        pysam.merge("-f", "-o", fout, *sams)
        pysam.index(fout)

        # TODO [W::sam_parse1] unrecognized reference name "haplo11"; treated as unmapped

    pysam.merge("-f", "-o", merged_bam, *out_bams)
    pysam.index(merged_bam)

if __name__ == "__main__":
    main(
        Path(snakemake.input.fin),
        snakemake.params.coverage,
        Path(snakemake.output.out),
        snakemake.output.bam,
        snakemake.output.merged_bam,
    )