cor2cov <- function(R, S) {
 sweep(sweep(R, 1, S, "*"), 2, S, "*")
 }

pgnorm <-
function(quantile, mean=0, sd=1, tail="upper")
{
xmin <--4
xmax <- 4
x <- seq(xmin,xmax,length=5000)*sd + mean
hx <- dnorm(x,mean,sd)

if(length(quantile) == 2){
  quantile2 <- quantile[2]
  quantile <- quantile[1]
}

if(tail=="upper")
{
ub=quantile
prob <- round(pnorm(quantile,mean,sd,lower.tail=F),5)
cvupper=round(ub,5)
plot(x, hx, type="n", 
     xlab=expression("Y|X"), 
     ylab="Densité",
     main = bquote("Distribution Conditionnelle:  " ~ mu[y.x] == .(mean) ~ "et" ~ sigma[y.x] == .(sd)),
     cex.main = 1.4,
     cex.lab = 1.2,
     cex.axis = 1.2,
     axes=TRUE)
j <- x >=ub
lines(x, hx)
abline(h=0)
abline(v=mean, lwd = 2, lty = 2, col = "darkgreen")
abline(v=ub,col="red")
polygon(c(ub,x[j],xmax), c(0,hx[j],0), col="steelblue1")

result <- paste("P(X \u2265", cvupper, ") =", round(prob, 4), "\n Z = ", round(quantile / sd, 3))
# mtext(result, 3, cex = 1.4, col = "red")
legend("topright", legend =  result, cex = 1, text.col = "red",
       bty = "n")
 }

else if(tail=="lower")
{
lb=quantile
cvlower=round(lb,5)
prob <- round(pnorm(quantile,mean,sd,lower.tail=T),5)
plot(x, hx, type="n", 
     xlab=expression("Y|X"),  
     ylab="Densité",
     main = bquote("Distribution Conditionnelle:  " ~ mu[y.x] == .(mean) ~ "et" ~ sigma[y.x] == .(sd)),
     cex.main = 1.5,
     cex.lab = 1.2,
     cex.axis = 1.2,
     axes=T)
i <- x <= lb
lines(x, hx)
abline(h=0)
polygon(c(xmin,x[i],lb), c(0,hx[i],0), col="steelblue1")
abline(v=mean, lwd = 2, lty = 2, col = "darkgreen")
abline(v=lb,col="red")

result <- paste("P(X \u2264",cvlower,") =", round(prob, 4), "\n Z = ", round(quantile / sd, 3))
# mtext(result, 3, cex = 1.4, col = "red")
legend("topleft", legend =  result, cex = 1, text.col = "red",
       bty = "n")
}
else if(tail=="two")
{
  lb = quantile
  ub = quantile2
  cvlower=round(lb, 5)
  cvupper=round(ub, 5)
  prob <- round(pnorm(quantile2, mean, sd, lower.tail=T) - pnorm(quantile, mean, sd, lower.tail=T), 5)
  plot(x, hx, type="n", 
       xlab=expression("Y|X"), 
       ylab="Densité",
       main = paste("Distribution Conditionnelle, avec \u03BC =",mean,", \u03C3 =",sd),
       cex.main = 1.8,
       cex.lab = 1.5,
       cex.axis = 1.2,
       axes=T)
  i <- x <= lb
  j <- x >= ub
  lines(x, hx)
  abline(h=0)
  abline(v=mean, lwd = 2, lty = 2, col = "darkgreen")
  abline(v=c(lb, ub), col="red")
  #polygon(c(xmin,x[i],lb), c(0,hx[i],0), col="steelblue1")
  #polygon(c(ub,x[j],xmax), c(0,hx[j],0), col="steelblue1")
  polygon(c(lb, x[x>lb & x<ub], ub), c(0, hx[x>lb & x<ub], 0), col="steelblue1")
  
  result <- paste("P(", cvlower, "\u2264 X \u2264", cvupper,") =", round(prob, 4), "\n Z = [", round(lb / sd, 3), ",  ", round(ub / sd, 3), "]")
  # mtext(result, 3, cex = 1.4, col = "red")
  legend("topleft", legend =  result, cex = 1, text.col = "red",
         bty = "n")
}
}


# From HH library
bivariateNormal <-
  function(rho=0, layout=c(3,3), lwd=.2,
           angle=c(22.5, 67.5, 112.5, 337.5, 157.5, 292.5, 247.5, 202.5),
           col.regions=trellis.par.get("regions")$col, ...)
{
  x <- seq(-2, 2, length=33)
  y <- x
  fxy <- 1/(2*pi*sqrt(1-rho^2)) *
    exp(-.5/(1-rho^2) * outer(x^2, y^2, "+") - 2*rho*outer(x,y,"*"))
  n.angle <- length(angle)
  Angle <- rep(angle, rep(length(fxy), n.angle))
  Viewing.Angle <- ordered(Angle, angle)
  wireframe(rep(fxy, n.angle) ~ rep(x[row(fxy)], n.angle) * rep(y[col(fxy)], n.angle) |
            Viewing.Angle, r=rho, layout=layout, lwd=lwd, ...,
            panel = function(x, y, subscripts, z, angle, ...)
            {
              w <- unique(angle[subscripts])
              panel.wireframe(x=x, y=y, subscripts=subscripts, z=z,
                              screen = list(z = w, x = -60, y = 0), ...)
            },
            angle = Angle, ## this is how to pass down external element
            strip = strip.custom(strip.names = TRUE, style = "1", sep=": "),
            skip = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
            drape = TRUE, distance = 0.3,
            main = as.expression(substitute("Bivariate Normal, " * rho == r,
                                             c(alist(rho=rho), list(r=rho)))),
            xlab = list("x", cex = 0.6),
            ylab = list("y", cex = 0.6),
            zlab = list("f(x,y)", cex = 0.6))
}
