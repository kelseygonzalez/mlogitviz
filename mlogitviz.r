# CLustering predicted probabilities
# Kelsey Gonzalez
# 2019-08-27

rm(list = ls()) # clean up the environment



##############################################################################################################
# main function
##############################################################################################################


#load data variables for ease now
ml <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")


runmodel = function(my, variables, mdata){
  # loading packages (and install first)
  packages = c("foreign", "nnet", "ggplot2", "reshape2")
  for( i in packages ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages(i, dependencies = TRUE )
      require(i, character.only = TRUE )
    }
  }
  # run multinomial logit model
  model = multinom(as.formula(paste(my, "~", paste(variables, collapse ="+"))), data = mdata)
  # create empty matrix to fill with results. set width by outcome variable length. 
  predictions = matrix(nrow=0, ncol=length(levels(mdata[,get("my")])))
  colnames(predictions) = levels(mdata[,get("my")])
  # gather predicted probabilities for each variable level for factor variables
  for (var in variables){
    if (class(mdata[,get("var")]) == "factor") {
      for (one_level in levels(mdata[,get("var")])){
        predictions = rbind(colMeans(fitted(model)[mdata[,get("var")] == one_level,]), predictions)
        rownames(predictions)[1] = paste0(var,"; ",one_level)
      }
    }
  # for numerical variables, run predicted probabilities with min, median and max. 
    if (class(mdata[,get("var")]) == "numeric"){
      # Min
      predictions = rbind(colMeans(fitted(model)[mdata[,get("var")] == (fivenum(mdata[,get("var")])[1]),]), predictions)
      rownames(predictions)[1] = paste0(var,"; Min")
      
      # Median
      predictions = rbind(colMeans(fitted(model)[mdata[,get("var")] == (fivenum(mdata[,get("var")])[3]),]), predictions)
      rownames(predictions)[1] = paste0(var,"; Median")
      
      # Max
      predictions = rbind(colMeans(fitted(model)[mdata[,get("var")] == (fivenum(mdata[,get("var")])[5]),]), predictions)
      rownames(predictions)[1] = paste0(var,"; Max")
    }
  }
  predictions
  
  
  #define euclidean distance function
  euc.distance <- function(x,y) {
    dist <- sum((x - y)^2)
    dist <- sqrt(dist)
    return(dist)
  }
  
  #create distance matrix (AE) with scaled matrix by column
  AE <- matrix(0, nrow = nrow(scale(predictions)), ncol = nrow(scale(predictions)))
  rownames(AE) <- colnames(AE) <- rownames(scale(predictions))
  for (i in 1:(nrow(scale(predictions))-1)) {
    for (j in (i + 1):nrow(scale(predictions))) {
      AE[j,i] <- AE[i,j] <- euc.distance(scale(predictions)[i,], scale(predictions)[j,])
    }}
  
  #perform hierarchical clustering on distance matrix (AE)
  hc <- hclust(d = as.dist(AE), method = "complete")
  
  #dendro.plot <- ggdendrogram(data = (as.dendrogram(hclust(d = as.dist(AE)))), rotate = TRUE) + 
  #  theme(axis.text.y = element_text(size = 9)) +
  #  geom_hline(yintercept=quantile(AE, .75), color="black", size=1)
  
  #Use that line (75% quartile of distance measures) to create cutoff point for number of groups using the average distance
  mycl <- cutree(hc, h=quantile(AE, .75))
  
  #Let's sort our scaled matrix by hc ordering
  scaledmatrix.sort <- (scale(predictions))[hc$order,]
  
  # Let's develop a sorted list of cluster memberships. The cluster memberships are from  mycl <- cutree(hc, h=quantile(AE, .75)) of foo:
  clusmemb <-  mycl[hc$order]
  
  # Here are the unique cluster numbers IN THE SAME ORDER as the permuted matrix:
  # unique.clusmemb <- unique(clusmemb)
  
  # PRESERVING THIS ORDER, let's find the number of members of each cluster:
  cluster.sizes <- rep(0, length(unique(clusmemb)))
  for (k in 1:length(cluster.sizes)) {
    cluster.sizes[k] <- sum(clusmemb == (unique(clusmemb))[k])
  }
  
  # Now draw the rectangles that go with these cluster sizes: The tricky part is to remember that cell (1,1) is in the
  # bottom left. We want to draw cumulative sum of cluster lengths from 1. Also, we add .5 because the "lines" function wants to 
  # draw the line in the *middle* of each square. The following web pages is helpful on this: url =
  # https://stackoverflow.com/questions/54717536/how-to-keep-abline-within-boundaries-of-the-corrplot
  cluster.heights <- cumsum(cluster.sizes)+ 0.5
  
  # Compare the cluster sizes, the cumulative cluster sizes, and the cluster heights:
  # cumsums <- cumsum(cluster.sizes)
  # cbind(cluster.sizes, cumsums, cluster.heights)
  
  # To map this onto our future heatmap, we need to map the different line's starting x and y coordinates
  # and ending x and y coordinates. We have Left vertical, Right verticle, minimal line, and then the cut points. 
  
  line = data.frame (
    Y = c(0.5, 0.5, 0.5, cluster.heights),
    Yend = c((nrow(scaledmatrix.sort)+0.5),(nrow(scaledmatrix.sort)+0.5), 0.5, cluster.heights),
    X = c(0.5, (ncol(scaledmatrix.sort)+0.5), rep(0.5,length(cluster.heights)+1)),
    Xend = c(0.5, (ncol(scaledmatrix.sort)+0.5), rep((ncol(scaledmatrix.sort)+0.5), length(cluster.heights)+1))
  )
  
  #convert wide form matrix to long form
  scaledmatrix.sort.long <- melt(scaledmatrix.sort, id = rownames(scaledmatrix.sort))
  #renaming columns for ease
  colnames(scaledmatrix.sort.long)[colnames(scaledmatrix.sort.long)=="Var1"] <- "Predictors"
  colnames(scaledmatrix.sort.long)[colnames(scaledmatrix.sort.long)=="Var2"] <- "Outcomes"
  # scaledmatrix.sort.long
  
  #create heatmap with rows as outcomns, columns as categorical variables, and colors as high to low from scaling
  heatmap.plot <- ggplot(data = scaledmatrix.sort.long, aes(x = Outcomes, y = Predictors)) +
    scale_fill_gradient2(low = "red4", mid = "white",
                         high = "forestgreen", midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
    geom_tile(aes(fill = value)) +  theme(legend.position = "top") +
    geom_segment(data=line, aes(x=X, y=Y, xend=Xend, yend=Yend)) 
  
  print(heatmap.plot)
  #google how to remove extra grey area in ggplot (padding?)
}

# works
runmodel("prog", c("ses","write","female","schtyp","awards"), ml)

# broken, why? need to debug
runmodel("prog", c("ses","write","math", "science", "female","ses","schtyp","awards"), ml)
runmodel("prog", c("ses","write","math", "science", "female","schtyp","awards"), ml)


