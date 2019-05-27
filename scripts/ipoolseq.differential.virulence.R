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

rmarkdown::render(input="scripts/ipoolseq.differential.virulence.Rmd",
                  output_file=basename(snakemake@output$`html`),
                  output_dir=normalizePath(dirname(snakemake@output$`html`)),
                  runtime="static",
                  intermediates_dir=tempdir(),
                  knit_root_dir=getwd(),
                  clean=TRUE)
