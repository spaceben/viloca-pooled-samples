# https://www.bioconductor.org/packages/devel/bioc/vignettes/QSutils/inst/doc/QSutils-Simulation.html
# https://www.bioconductor.org/packages/release/bioc/manuals/QSutils/man/QSutils.pdf
library(QSutils)
library(seqinr)
library(extraDistr)

len <- snakemake@params[["len"]]
set.seed(23)
dir <- "results/haplos/"

m1 <- GetRandomSeq(len)
write.fasta(as.list(m1), names="MasterSequence", as.string=FALSE,
            file.out=paste0("results/MasterSequence.fasta"))
n <- 5
v1 <- GenerateVars(m1, n-1, round(0.2*len),
                   ddnorm(1:round(0.2*len), round(0.1*len), round(0.03*len)))
w1 <- fn.ab(n,fn="pcf",r=4)

to_file <- function(haplos, freq, id) {
    stopifnot(length(haplos) == length(freq))
    names <- c()
    for (i in 1:length(haplos)) {
        names[i] = paste0("haplo", i-1, " freq:", freq[i] / sum(freq))
        write.fasta(as.list(haplos[i]), names=as.list(names[i]), as.string=FALSE,
                file.out=paste0(dir,id,"-",i,".fasta"), open="w")
    }
    # write.fasta(as.list(haplos), names=names, as.string=FALSE,
    #             file.out=paste0(dir,id,".fasta")) # TODO
}

to_file(c(m1, v1), w1, "sample0")

create_another_generation <- function(haplos, freq) {
    n2 = 5
    idx <- sample(1:length(haplos), 1, prob=freq/sum(freq))
    p2 <- Diverge(round(0.15*len):(round(0.15*len)+n2-1), haplos[idx]) # TODO
    w2 <- fn.ab(n2,fn="pcf",r=4)
    return(list(c(haplos, p2), c(freq, w2)))
}

s2 <- create_another_generation(c(m1, v1), w1)
to_file(s2[[1]], s2[[2]], "sample1")

s3 <- create_another_generation(s2[[1]], s2[[2]])
to_file(s3[[1]], s3[[2]], "sample2")

#data.frame(Hpl=DottedAlignment(s3[[1]]),Freq=s3[[2]])

create_ground_truth <- function(gen, name) {
    df <- setNames(data.frame(matrix(nrow = 0, ncol = 6)),
                c("", "type", "position", "variant", "haplotype", "frequency"))
    counter <- 1
    i_idx <- 0
    for (i in DottedAlignment(gen[[1]])) {
        if (i_idx == 0) {
            i_idx <- i_idx + 1
            next
        }
        j_idx <- 0
        for (j in strsplit(i, "")[[1]]) {
            if (j != ".") {
                df[counter,] <- c(counter-1, "mutation", j_idx, j,
                                paste0("haplo", i_idx),
                                gen[[2]][i_idx+1] / sum(gen[[2]]))
                counter <- counter + 1
            }
            j_idx <- j_idx + 1
        }
        i_idx <- i_idx + 1
    }

    write.csv(df, paste0("results/ground_truth", name, ".csv"), quote = FALSE, row.names = FALSE)
}

create_ground_truth(list(c(m1, v1), w1), "_sample0")
create_ground_truth(s2, "_sample1")
create_ground_truth(s3, "_sample2")