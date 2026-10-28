
intercept <- 0
diameter <- 4.6
pitch <- 5.4
trans <- 1.5
x0 <- diameter/8
y0 = 3

abline(h=y, v=x)

pep <- function(n = 1, intercept = 0, diameter = 4.6, pitch = 5.4,
                trans = 1.5, x0 = diameter/8, y0 = 0) {
  ymax <- max(n)*trans+y0
  y <- seq(from = ymax, by = -trans, length.out = length(n))
  x1 <- (y[1]-intercept)*(diameter/pitch)
  x <- seq(from = x1, by = -diameter/perturn, length.out = length(n))
  #x <- (y-intercept)*(diameter/pitch)
  x <- ifelse(x > diameter, x-(x%/%diameter)*diameter, x)
  xymat <- lapply(seq_along(n), function(i) {
    if (diameter - x[i] < x0 && diameter - x[i] > 0) {
      # Close to right border
      matrix(c(x[i], x[i]-diameter, rep(y[i], 2), rep(i, 2)),
             ncol = 3)
    } else if(x[i]-intercept < x0) {
      # Close to left border
      matrix(c(x[i]+diameter, x[i], rep(y[i], 2), rep(i, 2)),
             ncol = 3)
    } else {
      matrix(c(x[i], y[i], i), ncol = 3)
    }
  })
  
  #x2 <- ifelse(x > diameter, x-(x%/%diameter)*diameter, x)
  #matrix(c(x2, y), ncol=2)
  do.call(rbind, xymat)
}
pep(1:5)

plot(0, xlim=c(-5, 10), ylim=c(-5, 40), type = "n")
abline(v=c(0, diameter))
for (i in -1:7) abline(pitch*i, pitch/diameter)
pps <- pep(1:24, y0 = 0)
text(pps, labels = pps[,3])
points(pps, pch = 19)

#abline(h=pitch+trans)
