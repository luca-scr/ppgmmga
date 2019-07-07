#################################################################################
#                    PLOT PPGMMGA                                               #
#################################################################################

plot.ppgmmga <- function(x, class = NULL,
                         dim = seq(x$d),
                         drawAxis = TRUE,
                         bins = nclass.Sturges,
                         ...)
{
  stopifnot(inherits(x, "ppgmmga"))
  d <- length(dim)
  d <- ifelse(d > 2, 3, length(dim))
  Zpp <- data.frame(x$Z[,dim,drop=FALSE])
  znames <- colnames(Zpp)
  if(is.null(colnames(x$data)))
    colnames(x$data) <- paste0("X", seq(x$GMM$d))
  if(!is.null(class))
  { 
    class.name <- deparse(substitute(class))
    if(!is.factor(class)) class <- factor(class)
  } else
  { 
    class <- rep(1, nrow(Zpp))
    class.name <- NULL
  }

  switch(d,
  "1" = # 1D case ----
  { 
    gg <- ggplot(Zpp, aes_string(x = znames[1])) +
      geom_rug(alpha = 1/2)
    if(!is.null(class.name))
    {
      gg <- gg + 
        geom_histogram(aes_string(y = "..density..",
                                  fill = "class", colour = "class"),
                       position = "identity",
                       alpha = 0.6,
                       bins = if(is.function(bins)) bins(Zpp[,1])
                              else as.numeric(bins)) +
        scale_fill_tableau("Classic 10") +
        scale_color_tableau("Classic 10") +
        labs(fill = class.name, col = class.name)
    } else
    {
      gg <- gg +
        geom_histogram(aes_string(y = "..density.."),
                       position = "identity",
                       alpha = 0.6,
                       col = "white",
                       bins = if(is.function(bins)) bins(Zpp[,1])
                              else as.numeric(bins))
    }
    gg <- gg + theme_bw()
     
    return(gg)
  },
  "2" = # 2D case ----
  { 
    gg <- ggplot(Zpp, aes_string(x = znames[1], y = znames[2]))
    if(!is.null(class.name))
    {
      gg <- gg +
        geom_point(cex = 1, aes_string(shape = "class", colour = "class")) +
        scale_colour_tableau("Classic 10") +
        labs(shape = class.name, colour = class.name)
    } else 
    {
      gg <- gg + geom_point(cex = 1, shape = 20, colour = "black")
    }

    if(drawAxis)
    { 
      df2 <- data.frame(varnames = abbreviate(colnames(x$data), 
			                                        minlength = 8),
                        x = x$basis[,dim[1]],
                        y = x$basis[,dim[2]],
                        stringsAsFactors = FALSE)
      mult <- min((max(Zpp[,1]) - min(Zpp[,1])/(max(df2$x)-min(df2$x))),
                  (max(Zpp[,2]) - min(Zpp[,2])/(max(df2$y)-min(df2$y))) )
      df2 <- transform(df2,
                       x = .7 * mult * df2$x,
                       y = .7 * mult * df2$y)
      gg <-  gg + 
        geom_segment(data = df2, 
                     aes_string(x = 0, y = 0, xend = "x", yend = "y"),
                     arrow = arrow(length=unit(1/2,"picas")),
                     alpha = 0.5, color = "gray30") +
        geom_text(data = df2, aes_string(x = "x", y = "y", 
                                         label = "varnames"),
                  nudge_x = 0.1*diff(range(df2$x))*sign(df2$x)[1],
                  nudge_y = 0.1*diff(range(df2$y))*sign(df2$y)[1],
                  alpha = 0.5, color = "gray30")
    }

    if(nlevels(class) > 6)
    { gg <- gg + scale_shape_manual(values=1:nlevels(class)) }
           
    gg <- gg + theme_bw()
    return(gg)
  },
  { # d > 2 case ----
    if(is.null(class.name))
    { 
      class <- rep(1,nrow(Zpp))
      col <- "black"
    } else
    { 
      col <- tableau_color_pal("Classic 10")(nlevels(class)) 
    }
    clPairs(data = Zpp, classification = class,
            colors = col, ...)
  }
  )
}

