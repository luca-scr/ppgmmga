# Miscellaneous functions

projection <- function(B, W = diag(nrow(B)))
{
  # Returns the projection operator onto the subspace spanned by the cols of B
  # with respect to W inner product (by default the usual inner product).
  B <- as.matrix(B)
  B %*% MASS::ginv(t(B) %*% W %*% B) %*% t(B) %*% W
}

normalize <- function(x) 
{
  # Normalize the vector x to have unit length
  x <- as.vector(x)
  x <- x/drop(sqrt(crossprod(x)))
  return(x)
}

radian2degree <- function(x)
{
  # Convert a radian measure x to a degree measure.
  x*180/pi
}

Angle <- function (A, B, degree = FALSE)
{
  # Returns the sin of the maximal angle between two subspaces spanned by the 
  # (n x p) matrices A and B:
  #   sin(theta) = ||P_S(A) - P_S(B)||
  # where ||.|| is the spectral Euclidean norm, i.e. the maximum singular value,
  # and P_S(A) and P_S(B) are the projection matrices onto subspaces S(A) and
  # S(B). If degree = TRUE returns the angle in degree.
  A <- as.matrix(A)
  B <- as.matrix(B)
  if(any((dim(A) != dim(B))))
  { stop("matrices A and B do not have the same dimensions!") }
  A <- apply(A, 2, normalize)
  B <- apply(B, 2, normalize)
  th <- max(svd(projection(A) - projection(B))$d[-1])
  if(degree) th <- radian2degree(asin(th))
  return(th)
}

tr <- function(x)
{
  return(sum(diag(x)))
}

ppgmmga.table.angle <- function(...)
{
  models <- list(...)

  p <- models[[1]]$GMM$d
  d <- models[[1]]$d

  N <- length(models)
  angle <- matrix(data = NA,nrow = N,ncol = N)
  names <- rep(NA, N)

  for(j in 1:N)
  {
    for(i in j:N)
    {
      angle[j,i] <- Angle(A = models[[j]]$basis,
                          B = models[[i]]$basis,
                          degree = TRUE)
    }
    names[j] <- models[[j]]$approx
  }

  rownames(angle) <- names
  colnames(angle) <- names
  return(angle)
}


multiple_ggplot <- function(...)
{
  stopifnot(require(gridExtra))
  grid.arrange(...)
}

multiple_ggplot_sharedLegend <- function(..., ncol = length(list(...)), nrow  = 1, position = c("bottom", "right"))
{
  stopifnot(require(gridExtra))
  grid_arrange_shared_legend(..., ncol = ncol, nrow = nrow, position = position)
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) 
{
  stopifnot(require(gridExtra))
  stopifnot(require(grid))
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- arrangeGrob(do.call(arrangeGrob, gl),
                          legend,
                          ncol = ifelse(position == "bottom", 1, 2), 
                          heights = unit.c(unit(1, "npc") - lheight, lheight))
  grid.newpage()
  cobined <- grid.draw(combined)
  return(combined)

}

