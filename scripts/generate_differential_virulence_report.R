library(rmarkdown)

if (exists("SNAKEMAKE.STANDIN") && get("SNAKEMAKE.STANDIN")) {
  setClass("SnakemakeStandin", representation(input = "list", output = "list"))
  
  snakemake <- new("SnakemakeStandin",
    input=list(`in`="data/Simulation/sim-in.count.tab",
               `out`="data/Simulation/sim-out.count.tab",
               `knockouts`="cfg/Uhse_et_al.2018/knockouts.gff"),
    output=list(`table`="data/Simulation/sim.da.tab",
                `html`="data/Simulation/sim.report.html")
  )
}

rmarkdown::render(input="scripts/differential_virulence.Rmd",
                  output_file=basename(snakemake@output$`html`),
                  output_dir=normalizePath(dirname(snakemake@output$`html`)),
                  output_format="html_document",
                  runtime="static",
                  intermediates_dir=tempdir(),
                  knit_root_dir=getwd(),
                  clean=TRUE)
