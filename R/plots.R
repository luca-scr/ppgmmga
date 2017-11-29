
#################################################################################
#                    PLOT PPGMMGA                                               #
#################################################################################



plot.ppgmmga <- function(object, class = NULL,
                         dim = 1:object$d,
                         drawAxis = TRUE, ...)
{

  if(!is.null(class) & !is.factor(class))
  { class <- factor(class) }
  if(is.null(class))
  { poi <- geom_point(cex = 1,
                      shape = 20,
                      colour = "black") }
  else
  { poi <- geom_point(cex = 1,
                      aes(shape = class, colour = class)) }

  d <- length(dim)
  d <- ifelse(d > 2,3,length(dim))
  Zpp <- data.frame(object$Z[,dim,drop=FALSE])
  znames <- colnames(Zpp)

  if(is.null(colnames(object$data))){
    colnames(object$data) <- paste0("V", seq(object$GMM$d))
  }

  switch(d,
         "1" =  { if(is.null(class))
         { gg <- ggplot(Zpp, aes_string(x = znames[1], fill = NULL)) +
           geom_histogram(aes(y = ..density..),
                          position = "identity",
                          alpha = 0.6,
                          bins = nclass.FD(Zpp[,1]),
                          ...)
         }
           else
           { gg <- ggplot(Zpp, aes_string(x = znames[1], fill = class)) +
             geom_histogram(aes(y=..density..),
                            position = "identity",
                            alpha = 0.6,
                            bins = nclass.FD(Zpp[,1]),
                            ...) +
             scale_fill_tableau('tableau10') +
             theme(legend.title=element_blank())
           }
           gg <- gg + geom_rug(alpha = 1/2)
           return(gg)
         },
         "2" = { gg <- ggplot(Zpp, aes_string(x = znames[1], y = znames[2])) +
           poi +
           scale_colour_tableau('tableau10') +
           theme(legend.title=element_blank())

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
         { col <- NULL
         if(is.null(class))
         { class <- rep(1,nrow(Zpp))
         col <- "black"
         clPairs(data = Zpp, classification = class,
                 CEX = 1, colors = "black", ...)}
         else
         { col <- tableau_color_pal()(nlevels(class)) }
         clPairs(data = Zpp, classification = class,
                 CEX = 1, colors = col, ...)
         }
  )
}





#################################################################################
#                    PLOT PPGMMGAAD                                              #
#################################################################################



plot.ppgmmgaAD <- function(object,
                           dim = 1:object$d,
                           drawRegion = FALSE,
                           drawAxis = FALSE,
                           class = NULL)

{

  if(is.null(class)){class <- object$Class}
  gg <- plot.ppgmmga(object, drawAxis = drawAxis, class = class,dim = dim)

  if(length(dim) != 2) {drawRegion = FALSE}


  if(drawRegion == TRUE & object$d == 2)
  {

    data <- gg$data
    switch(attributes(object$GMM$hypvol)$method,
           "ConvHull" = {
             dat <- data.convHull(coord = attributes(object$GMM$hypvol)$coord, data = data)
             gg <- gg + geom_line(data = dat,aes(x = x,y = y, group = class))},

           "box" = {
             dat <- data.box(data = data)
             gg <- gg + geom_line(data = dat,aes(x = x,y = y, group = class))})

  }else if (drawRegion == TRUE & length(dim) == 2){

    data <- gg$data
    switch(attributes(object$GMM$hypvol)$method,
           "ConvHull" = {
             vol <- volume(data = data,method = "ConvHull")
             dat <- data.convHull(coord = vol$coord, data = data)
             gg <- gg + geom_line(data = dat,aes(x = x,y = y, group = class))},

           "box" = {
             dat <- data.box(data = data)
             gg <- gg + geom_line(data = dat,aes(x = x,y = y, group = class))})


  }
  gg
}




data.convHull <- function(coord, data)
{

  convHull <- coord
  row <- nrow(convHull)
  convHull <- as.vector(t(convHull))
  dat <- data.frame(matrix(NA,ncol = 3,
                           nrow = length(convHull),
                           dimnames = list(NULL,c("x","y","class"))))

  dat[,1:2] <- data[convHull,]
  dat[,3] <- rep(1:row,each = 2)

  return(dat)
}


data.box <- function(data)
{

  xcoord <- c(maxx = max(data[,1]),
              minx = min(data[,1]))
  ycoord <- c(maxy = max(data[,2]),
              miny = min(data[,2]))
  points <- as.data.frame(expand.grid(x = xcoord,y = ycoord))
  colnames(points) <- colnames(data)
  vol <- volume(data = points,method = "ConvHull")
  #V <- vol$volume
  convHull <- vol$coord
  row <- nrow(convHull)
  convHull <- as.vector(t(convHull))
  dat <- data.frame(matrix(NA,ncol = 3,
                           nrow = length(convHull),
                           dimnames = list(NULL,c("x","y","class"))))

  #dat[,1:2] <- attributes(vol)$data[convHull,]
  dat[,1:2] <- points[convHull,]
  dat[,3] <- rep(1:row,each = 2)

  return(dat)
}

