import subprocess

def main(bam, merged_bam, ref, mode):
    if mode == "classic":
        bam = merged_bam
    else:
        bam = " ".join(bam)

    subprocess.run(f"shorah shotgun -b {bam} -f {ref} \
                   --sampler use_quality_scores --alpha 0.0001 \
                   --n_max_haplotypes 100", #--min_windows_coverage", # TODO min_windows_coverage
                   cwd=f"results/viloca_{mode}", shell=True, check=True,
                   capture_output=True)


if __name__ == "__main__":
    main(
        [i.replace("results", "..") for i in snakemake.input.bam],
        snakemake.input.merged_bam.replace("results", ".."),
        snakemake.input.ref.replace("results", ".."), # TODO better solution than replace?
        snakemake.wildcards.mode
    )