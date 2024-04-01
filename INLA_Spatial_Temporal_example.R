library(INLA)


data('burkitt', package = 'splancs')

t(sapply(burkitt[, 1:3], summary))

k <- 6
tknots <- seq(min(burkitt$t), max(burkitt$t), length = k)
mesh.t <- inla.mesh.1d(tknots)
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(burbdy)),
                                          '0')))

mesh.s <- inla.mesh.2d(burpts,
                       boundary = inla.sp2segment(domainSP), 
                       max.edge = c(10, 25), cutoff = 5) # a crude mesh

spde <- inla.spde2.pcmatern(mesh = mesh.s,
                            prior.range = c(5, 0.01), # P(practic.range < 5) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
m <- spde$n.spde

Ast <- inla.spde.make.A(mesh = mesh.s, loc = burpts,
                        n.group = length(mesh.t$n), group = burkitt$t,
                        group.mesh = mesh.t)

idx <- inla.spde.make.index('s', spde$n.spde, n.group = mesh.t$n)

dmesh <- book.mesh.dual(mesh.s)

library(rgeos)
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i,], domainSP))
    return(gArea(gIntersection(dmesh[i,], domainSP)))
  else return(0)
})

st.vol <- rep(w, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), m)
n <- nrow(burkitt)
y <- rep(0:1, c(k * m, n))
expected <- c(st.vol, rep(0, n))
stk <- inla.stack(
  data = list(y = y, expect = expected), 
  A = list(rbind(Diagonal(n = k * m), Ast), 1), 
  effects = list(idx, list(a0 = rep(1, k * m + n))))


pcrho <- list(prior = 'pc.cor1', param = c(0.7, 0.7))
form <- y ~ 0 + a0 + f(s, model = spde, group = s.group, 
                       control.group = list(model = 'ar1',
                                            hyper = list(rho = pcrho)))

burk.res <- inla(form, family = 'poisson', 
                 data = inla.stack.data(stk), E = expect,
                 control.predictor = list(A = inla.stack.A(stk)),
                 control.inla = list(strategy = 'adaptive'))


eta.at.integration.points <- burk.res$summary.fix[1,1] +
  burk.res$summary.ran$s$mean
c(n = n, 'E(n)' = sum(st.vol * exp(eta.at.integration.points)))



r0 <- diff(range(burbdy[, 1])) / diff(range(burbdy[, 2]))
prj <- inla.mesh.projector(mesh.s, xlim = range(burbdy[, 1]),
                           ylim = range(burbdy[, 2]), dims = c(100, 100 / r0)) 
ov <- over(SpatialPoints(prj$lattice$loc), domainSP)
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         burk.res$summary.ran$s$mean[1:m + (j - 1) * m])
  r[is.na(ov)] <- NA
  return(r) 
})