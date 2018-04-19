
#################################################################################
#                    PLOT PPGMMGA                                               #
#################################################################################



plot.ppgmmga <- function(object, class = NULL,
                         dim = 1:object$d,
                         drawAxis = TRUE, # if 2D
                         bins = nclass.Sturges, # if 1D
                         ...)
{

  d <- length(dim)
  d <- ifelse(d > 2,3,length(dim))
  Zpp <- data.frame(object$Z[,dim,drop=FALSE])
  znames <- colnames(Zpp)
  if(is.null(colnames(object$data)))
    colnames(object$data) <- paste0("X", seq(object$GMM$d))

  if(!is.null(class))
  { 
    class.name <- deparse(substitute(class))
    if(!is.factor(class)) class <- factor(class)
    poi <- geom_point(cex = 1, aes(shape = class, colour = class))
  }
  else
  { poi <- geom_point(cex = 1, shape = 20, colour = "black") }

  switch(d,
         "1" = # 1D case
         { 
           gg <- ggplot(Zpp, aes_string(x = znames[1], fill = class)) +
             geom_histogram(aes(y=..density..),
                            position = "identity",
                            alpha = 0.6,
                            bins = if(is.function(bins)) bins(Zpp[,1])
                                   else as.numeric(bins)) +
             scale_fill_tableau("tableau10")
           if(!is.null(class))
             gg <- gg + labs(fill = class.name)
           gg <- gg + geom_rug(alpha = 1/2)
           return(gg)
         },
         "2" = # 2D case
         { 
           gg <- ggplot(Zpp, aes_string(x = znames[1], y = znames[2])) + 
             poi + scale_colour_tableau("tableau10")
           if(!is.null(class))
             gg <- gg + labs(shape = class.name, colour = class.name)

           if(drawAxis == TRUE)
           { df2 <- data.frame(varnames = abbreviate(colnames(object$data), 5),
                               x = object$basis[,1],
                               y = object$basis[,2],
                               stringsAsFactors = FALSE)
             mult <- min((max(Zpp[,1]) - min(Zpp[,1])/(max(df2$x)-min(df2$x))),
                         (max(Zpp[,2]) - min(Zpp[,2])/(max(df2$y)-min(df2$y))) )
             df2 <- transform(df2,
                              x = .7 * mult * df2$x,
                              y = .7 * mult * df2$y)
             gg <-  gg + geom_segment(data = df2, aes(x=0, y=0, xend=x, yend=y),
                                      arrow = arrow(length=unit(1/2,"picas")),
                                      alpha = 0.5, color = "gray30") +
               geom_text(data = df2, aes(x=x, y=y, label=varnames),
                         nudge_x = 0.1*sign(df2$x),
                         nudge_y = 0.1 * sign(df2$y),
                         alpha = 0.5, color="gray30")
           }

           if(nlevels(class) > 6)
           { gg <- gg + scale_shape_manual(values=1:nlevels(class)) }
           
           return(gg)
         },
         { # d > 2 case
           if(is.null(class))
           { class <- rep(1,nrow(Zpp))
             col <- "black"
           } else
           { col <- tableau_color_pal()(nlevels(class)) }
           clPairs(data = Zpp, classification = class,
                   CEX = 1, colors = col, ...)
         }
  )
}


