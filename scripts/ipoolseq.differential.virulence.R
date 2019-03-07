# ipoolseq.differential.virulence.R, Copyright 2018, 2019 Florian G. Pflug
# 
# This file is part of the iPool-Seq Analysis Pipeline
#
# The iPool-Seq Analysis Pipeline is free software: you can redistribute it
# and/or modify it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# The iPool-Seq Analysis Pipeline is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with the iPool-Seq Analysis Pipeline.  If not, see
# <http://www.gnu.org/licenses/

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
