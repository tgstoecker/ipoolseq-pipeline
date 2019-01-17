library(data.table)
library(rtracklayer)
library(gwpcR)

# Simulate UMI count tables
simulate.trumicount.output <- function(knockouts, neutral,
                                       abundance.in=10**(runif(length(knockouts), min=0, max=6)),
                                       foldchange=2**ifelse(neutral,
                                                            rep(0, length(knockouts)),
                                                            runif(length(knockouts), min=-4, max=-1)),
                                       scale=0.4, dispersion=0.05,
                                       threshold.in=1, threshold.out=2,
                                       eff.in.5p=0.75, eff.in.3p=0.25, depth.in.5p=1.1, depth.in.3p=0.8,
                                       eff.out.5p=0.15, eff.out.3p=0.23, depth.out.5p=4.6, depth.out.3p=5.5,
                                       eff.sd=0.05, depth.sd=0.1)
{
  # Beta-distributed random variables with the given mean and variance
  rbeta.moments <- function(n, mean, var) {
    var <- min(var, mean*(1-mean) / (1 + .Machine$double.eps))
    alpha <- mean * (mean * (1 - mean) / var - 1)
    beta <- (1 - mean) * (mean * (1 - mean) / var - 1)
    rbeta(n, shape1=alpha, shape2=beta)
  }

  # Beta-distributed random variables with the given mean and variance
  rgamma.moments <- function(n, mean, var) {
    theta <- var / mean
    k <- mean / theta
    rgamma(n, shape=k, scale=theta)
  }
  
  # Simulation parameters
  stopifnot(length(knockouts) == length(neutral))
  stopifnot(length(knockouts) == length(abundance.in))
  stopifnot(length(knockouts) == length(foldchange))
  
  # Simulate PCR & sequencing parameters
  message("Simulating PCR & sequencing parameters (efficiency and depth)")
  n <- length(knockouts)
  random.params <- function(values, sd, type) {
    if (length(values) == length(knockouts))
      return(values)
    else if (length(values) == 1)
      switch(type,
             `[0,1]`=rbeta.moments(n, mean=values, var=sd**2),
             `[0,infty)`=rgamma.moments(n, mean=values, var=sd**2))
    else
      stop("need either a single value or as many as there are knockouts")
  }
  eff.in.5p <- random.params(eff.in.5p, eff.sd, type="[0,1]")
  eff.in.3p <- random.params(eff.in.3p, eff.sd, type="[0,1]")
  depth.in.5p <- random.params(depth.in.5p, depth.sd, type="[0,infty)")
  depth.in.3p <- random.params(depth.in.3p, depth.sd, type="[0,infty)")
  eff.out.5p <- random.params(eff.out.5p, eff.sd, type="[0,1]")
  eff.out.3p <- random.params(eff.out.3p, eff.sd, type="[0,1]")
  depth.out.5p <- random.params(depth.out.5p, depth.sd, type="[0,infty)")
  depth.out.3p <- random.params(depth.out.3p, depth.sd, type="[0,infty)")

  # Simulate the true input abundance (uniformly distributed),
  # and the true output abundance (gamma distributed)
  message("Simulating input & output abundances")
  a.in <- 10**(runif(n, min=0, max=6))
  a.out <- a.in * scale * foldchange * rgamma.moments(length(a.in), mean=1, var=dispersion)

  # Create TRUmiCount-compatible output
  generate.counts <- function(data) {
    data[, {
      loss <- pgwpcrpois(threshold-1, efficiency=efficiency, lambda0=depth)
      n.umis <- mapply(abundance, loss, FUN=function(a,l) { rpois(1, a*(1-l)) }, SIMPLIFY=TRUE)
      list(gene=paste0(knockout, ':', flank),
           n.umis=n.umis,
           n.tot=n.umis / (1 - loss),
           efficiency=efficiency,
           depth=depth,
           loss=loss,
           n.obs=n.umis)
    }]
  }
  
  # Simulate UMI counts as being Poissin distributed around the true abundances times (1-loss)
  message("Simulating count tables")
  table.in <- generate.counts(data.table(knockout=rep(knockouts, 2),
                                         flank=rep(c('5p', '3p'), each=n),
                                         abundance=c(a.in, a.in),
                                         efficiency=c(eff.in.5p, eff.in.3p),
                                         depth=c(depth.in.5p, depth.in.3p),
                                         threshold=threshold.in))
  table.out <- generate.counts(data.table(knockout=rep(knockouts, 2),
                                          flank=rep(c('5p', '3p'), each=n),
                                          abundance=c(a.out, a.out),
                                          efficiency=c(eff.out.5p, eff.out.3p),
                                          depth=c(depth.out.5p, depth.out.3p),
                                          threshold=threshold.out))
  
  # Table of true input abundances and foldchanges
  table.truth <- data.table(knockout=knockouts,
                            abundance.in=abundance.in,
                            foldchange=foldchange)
  
  return(list(`in`=table.in, `out`=table.out, `truth`=table.truth))
}

# Read list of knockouts
message("Reading knockouts (", snakemake@input$gff, ")")
knockouts <- readGFF(snakemake@input$gff, version=2)
# Simulate count table, assuming every second knockout is neutral
r <- simulate.trumicount.output(knockouts=knockouts$Name,
                                neutral=((knockouts$Neutral==1) |
                                          c(TRUE, FALSE)[1 + 1:nrow(knockouts) %% 2]))

message("Writing output files (", snakemake@output$pool_in, ", ", snakemake@output$pool_out, ", ", snakemake@output$truth, ")")
write.table(r$`in`, snakemake@output$pool_in, sep="\t", col.names=TRUE)
write.table(r$`out`, snakemake@output$pool_out, sep="\t", col.names=TRUE)
write.table(r$`truth`, snakemake@output$truth, sep="\t", col.names=TRUE)
