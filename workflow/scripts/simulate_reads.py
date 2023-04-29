import subprocess
from pathlib import Path
import glob
from os import path
import pysam
import fileinput
import re


def main(fin, coverage, out_dir, out_bams, merged_bam):

    for i in glob.glob(path.join(fin, "*.fasta")):
        with open(i, "r") as f:
            freq = float(f.readline().split(":")[1])
        subprocess.run(f"art_illumina -sam --paired -m 320 -s 10 --rndSeed 20 \
                        -i {i} -f {str(round(coverage * freq))} -l 200 -na \
                        -o {i.split('.fasta')[0]} -q", shell=True, check=True,
                        capture_output=True)

    all_generated_sams = glob.glob(path.join(fin, "*.sam"))
    for i in all_generated_sams:
        with pysam.AlignmentFile(i) as af:
            haplotype_name = af.get_reference_name(0).split(" ")[0]
        with fileinput.FileInput(i, inplace=True, backup=".bak") as fd:
            for line in fd:
                if "@SQ" in line:
                    # replace reference name in header
                    print(
                        re.sub(
                            rf"SN:.*?\t", f"SN:MasterSequence\t", line
                        ),
                        end="",
                    )
                else:
                    print(
                        re.sub(rf"\t{haplotype_name}.*?\t", f"\tMasterSequence\t", line),
                        end="",
                    )
        pysam.sort("-o", i, i) # sort in place

    grps = [int(i.split("sample")[-1].split("-")[0]) for i in all_generated_sams]
    for i in set(grps):
        sams = [j for j, k in zip(all_generated_sams, grps) if k == i]
        fout = path.join(out_dir, 'sample'+str(i)+'.bam')

        pysam.merge("-f", "-o", fout, *sams)
        pysam.index(fout)

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