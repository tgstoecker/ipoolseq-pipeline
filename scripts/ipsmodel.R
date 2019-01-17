# model.R, Copyright 2017, 2019 Florian G. Pflug
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

dipsmodel <- function(n.out, scale, disp, n.in, l.out, l.in, log=FALSE) {
  dnbinom(n.out, mu=scale*n.in*(1-l.out)/(1-l.in), size=n.in/(1+disp*n.in), log=log)
}

pipsmodel <- function(n.out, scale, disp, n.in, l.out, l.in, log.p=FALSE) {
  pnbinom(n.out, mu=scale*n.in*(1-l.out)/(1-l.in), size=n.in/(1+disp*n.in), log.p=log.p)
}

qipsmodel <- function(p, scale, disp, n.in, l.out, l.in, log.p=FALSE) {
  qnbinom(p, mu=scale*n.in*(1-l.out)/(1-l.in), size=n.in/(1+disp*n.in), log.p=log.p)
}

ipsmodel.fit <- function(n.out, n.in, l.out, l.in) {
  if ((length(n.out) != length(n.in)) || (length(n.in) != length(l.out)) || (length(l.out) != length(l.in)))
    stop("n.out, n.in, l.out, l.in must have the same length")
  
  
  p.init <- c(scale=mean(((n.out * (1-l.out)) / (n.in * (1-l.in)))[n.in > 0], na.rm=TRUE),
              disp=1)
  ll.fun <- function(p) {
    if (all(p > 0))
      -sum(dipsmodel(n.out, p['scale'], p['disp'], n.in, l.out, l.in, log=TRUE))
    else NA
  }
  r <- optim(par=p.init, fn=ll.fun)
  if (r$convergence == 0)
    as.list(r$par)
  else
    stop("numerical optimazation failed to converge")
}

if (FALSE) {
  # Test the model
  n.in.true <- outer(c(1,3), 10^c(0,0,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
  l.in <- c(0.66, 0.33, 0.5, 0.16)[1 + 0:(length(n.in.true)-1) %% 4]
  n.in <- mapply(function (ni, li) rpois(1, ni*(1-li)), ni=n.in.true, li=l.in)
  
  d <- c(0.6, 1.2, 0.8, 1.4, 1)[1 + 0:(length(n.in.true)-1) %% 4]

  scale <- 10
  n.out.true <- scale * n.in.true * d
  l.out <- c(0.2, 0.6, 0.4, 0.8)[1 + 0:(length(n.out.true)-1) %% 4]
  n.out <- mapply(function (no, lo) rpois(1, no*(1-lo)), no=n.out.true, lo=l.out)
  
  print(ipsmodel.fit(n.out[n.in > 0], n.in[n.in > 0], l.out[n.in > 0], l.in[n.in > 0]))
}
