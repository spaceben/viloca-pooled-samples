import subprocess

def main(bam, ref, run):
    l = run.split("_")
    if len(l) == 2:
        mode, id = l
    else:
        mode = l[0]
    if mode == "pooled":
        bam = " ".join(bam)
    else:
        bam = bam[0].replace("0", str(id)) # TODO dirty solution

    print(bam)
    subprocess.run(f"shorah shotgun -b {bam} -f {ref} \
                    --sampler shorah --alpha 0.0001",
                    #--n_max_haplotypes 100", #--min_windows_coverage", # TODO min_windows_coverage
                    cwd=f"results/viloca_{run}", shell=True,
                    check=True, capture_output=True)


if __name__ == "__main__":
    main(
        [i.replace("results", "..") for i in snakemake.input.bam],
        snakemake.input.ref.replace("results", ".."), # TODO better solution than replace?
        snakemake.wildcards.run
    )