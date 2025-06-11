library(INLA)
library(fields)
library(ggplot2)
library(viridisLite)

rm(list=ls())
set.seed(201803)
inla.seed = sample.int(n=1E6, size=1)
options(width=70, digits=3)

local.plot.field = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
  stopifnot(length(field) == mesh$n)
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  n.col = 20
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1, ...)
}

sigma.u = 1.5
# - the marginal standard deviation of the spatial field
range = 2
# - the correlation range of the spatial field
#kappa = sqrt(8)/range
# - kappa was the parameter we used in previous years
# - not used in this code, but you may encounter it


fake.locations = matrix(c(0,0,10,10, 0, 10, 10, 0), nrow = 4, byrow = T)
mesh.sim = inla.mesh.2d(loc = fake.locations, max.edge=c(0.4, 1))

plot(mesh.sim)
axis(1); axis(2)

spde = inla.spde2.pcmatern(mesh.sim, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

Qu = inla.spde.precision(spde, theta=c(log(range), log(sigma.u)))
u = inla.qsample(n=1, Q=Qu, seed = inla.seed)
u = u[ ,1]


local.plot.field(u, mesh.sim)
len = range
# - the true range
arrows(5-0.5*len, 6, 5+0.5*len, 6, length=0.05, angle=90, code=3, lwd=3)


n = 5*1E3
# - number of measurement locations
# - Don't do more than 5000 if you want the code to be quick to run (1 minute)
# - can do 100*1E3

loc.data = matrix(runif(2*n), n)*10 # coordinates

A = inla.spde.make.A(mesh=mesh.sim, loc=loc.data)
u = drop(A %*% u)

quilt.plot(x=loc.data[, 1],y=loc.data[, 2],z=u,nx=80,ny=80, 
           col = plasma(101), main="Field projected to data locations", 
           zlim = range(u))


sigma.iid = 0.3

x = runif(n)-0.5
# - mean 0 to not affect intercept
beta = c(1, 2)

lin.pred = beta[1] + beta[2]*x + u

y = lin.pred + sigma.iid*rnorm(n)

quilt.plot(x=loc.data[, 1],y=loc.data[, 2],z=y,nx=80,ny=80, 
           col = plasma(101), main="Observed data", 
           zlim = range(y))

df = data.frame(y=y, locx=loc.data[ ,1], locy=loc.data[ ,2], x = x)
summary(df)

mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(0.4, 1))
mesh$n

A = inla.spde.make.A(mesh=mesh, loc=data.matrix(df[ , c('locx', 'locy')]))
dim(A)

par(mar=c(1,1,1,1))
plot(mesh, asp=1)
points(df[ , c('locx', 'locy')], col='red', lwd=.1)


prior.median.sd = 1; prior.median.range = 5
spde = inla.spde2.pcmatern(mesh, prior.range = c(prior.median.range, .5), prior.sigma = c(prior.median.sd, .5))

stack = inla.stack(tag='est',
                   # - Name (nametag) of the stack
                   # - Here: est for estimating
                   data=list(y=df$y),
                   effects=list(
                     # - The Model Components
                     s=1:spde$n.spde, 
                     # - The first is 's' (for spatial)
                     data.frame(intercept=1, x=df$x)),
                   # - The second is all fixed effects
                   A=list(A, 1)
                   # - First projector matrix is for 's'
                   # - second is for 'fixed effects'
)


family = "gaussian"
prior.median.sd.g = 0.5 # prior median for sigma.epsilon
control.family = list(hyper = list(prec = list(
  prior = "pc.prec", param =
    c(prior.median.sd.g,0.5))))


formula = y ~ -1 + intercept + x + f(s, model=spde)
# - Remove standard intercept (always when using inla.stack)
# - Fixed effects + random effects
# - s belongs to the mesh
# - A-matrix will tell inla() how to go from mesh to data

initial.theta = c(2.35, 0.79, 0.46)
# - the first time you run this, set it to NULL
# - after running, set it to res$internal.summary.hyperpar$mean
# - and run the code again
# - Reason: Both faster and better inference

res = inla(formula, data=inla.stack.data(stack),
           family = family,
           control.family = control.family,
           control.predictor=list(A = inla.stack.A(stack)),
           quantiles=c(0.5, 0.025, 0.975, 0.1, 0.9, 0.25, 0.75),
           #control.compute = list(config=T, dic=T, cpo=T, waic=T), 
           # - Model comparisons
           #control.inla = list(int.strategy='grid'),
           # - More accurate integration over hyper-parameters
           control.mode = list(restart = T, theta = initial.theta))

summary(res)
