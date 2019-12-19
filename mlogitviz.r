# CLustering predicted probabilities
# Kelsey Gonzalez
# 2019-05-27
# last update: 2019-12-05


##############################################################################################################
# main function
##############################################################################################################

runmodel = function(my, variables, mdata){
  # loading packages (and install first)
  packages = c("foreign", "nnet", "ggplot2", "reshape2", "gridExtra", "svglite","sjPlot", "sjmisc", "sjlabelled")
  for( i in packages ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages(i, dependencies = TRUE )
      require(i, character.only = TRUE)
    }
  }
  # run multinomial logit model
  model = multinom(as.formula(paste(my, "~", paste(variables, collapse ="+"))), data = mdata)
  tab_model(model, p.style = "asterisk", string.resp = TRUE, show.ci = FALSE, file = "output.html")

  # create empty matrix to fill with results. set width by outcome variable length. 
  pred = matrix(nrow=0, ncol=length(levels(mdata[,get("my")])))
  colnames(pred) = levels(mdata[,get("my")])
  
  # gather predicted probabilities for each variable level for factor variables
  for (var in variables){
    if (class(mdata[,get("var")]) == "factor") {
      for (one_level in levels(mdata[,get("var")])){
        pred = rbind(colMeans(fitted(model)[mdata[,get("var")] == one_level,]), pred)
        rownames(pred)[1] = paste0(var,"; ",one_level)
      }
    }
  
  # for numerical variables, run predicted probabilities with min, median and max. 
    if (class(mdata[,get("var")]) == "numeric"){
      
      #get descriptive statistics for the variable
      five_num = fivenum(mdata[,get("var")])
      
      # Min
      fitted_1 = fitted(model)[mdata[,get("var")] == (five_num[1]),]
      if (class(fitted_1) == "numeric") {fitted_1 = t(as.matrix(fitted_1))}
      col_Means = colMeans(fitted_1)      
      pred = rbind(col_Means, pred)
      rownames(pred)[1] = paste0(var,"; Min")
      
      # Median
      fitted_1 = fitted(model)[mdata[,get("var")] == (five_num[3]),]
      if (class(fitted_1) == "numeric") {fitted_1 = t(as.matrix(fitted_1))}
      col_Means = colMeans(fitted_1)      
      pred = rbind(col_Means, pred)
      rownames(pred)[1] = paste0(var,"; Median")
      
      # Max
      fitted_1 = fitted(model)[mdata[,get("var")] == (five_num[5]),]
      if (class(fitted_1) == "numeric") {fitted_1 = t(as.matrix(fitted_1))}
      col_Means = colMeans(fitted_1)      
      pred = rbind(col_Means, pred)
      rownames(pred)[1] = paste0(var,"; Max")
      
      if (length(unique(mdata[,get("var")])) > 5) {
        # 1st Quart
        fitted_1 = fitted(model)[mdata[,get("var")] == (five_num[2]),]
        if (class(fitted_1) == "numeric") {fitted_1 = t(as.matrix(fitted_1))}
        col_Means = colMeans(fitted_1)      
        pred = rbind(col_Means, pred)
        rownames(pred)[1] = paste0(var,"; 1stquart")
        
        # 3rd Quart
        fitted_1 = fitted(model)[mdata[,get("var")] == (five_num[4]),]
        if (class(fitted_1) == "numeric") {fitted_1 = t(as.matrix(fitted_1))}
        col_Means = colMeans(fitted_1)      
        pred = rbind(col_Means, pred)
        rownames(pred)[1] = paste0(var,"; 3rdquart")
      }
    }
  }
  # until I fixed the fact that fivenum may give values no representative in the dataset, 
  # I will remove variables that are NA
  pred <- na.omit(pred) 
  write.csv(pred, file = "1_unsorted_pred.csv")
  
  
  # perform hierarchical clustering on a euclidean distance matrix
  hc <- hclust(d = dist(pred, method = "euclidean"), method = "ward.D2")
  plot(hc, cex = 0.6, hang = -1, main="Hierarchal Cluster of Predicted Probabilities")
  abline(h=quantile(dist(pred, method = "euclidean"), .75))
  rect.hclust(hc, h=quantile(dist(pred, method = "euclidean"), .75), border = "#c67d08")

    # Use that line (75% quartile of distance measures) to create cutoff point for number of groups using the average distance
  mycl <- cutree(hc, h=quantile(dist(pred, method = "euclidean"), .75))
  mycl
  print(mycl)

  # Let's sort our scaled matrix by hc ordering
  matrix.sort <- (pred)[hc$order,]
  write.csv(matrix.sort, file = "2_sorted_pred.csv")
  
  scaledmatrix.sort <- (scale(matrix.sort))
  write.csv(scaledmatrix.sort, file = "3_scaled_sorted_pred.csv")
  
  pred.long <- melt(pred[hc$order,], id = rownames(pred))
  colnames(pred.long)[colnames(pred.long)=="value"] <- "originalvalue"

  # Let's develop a sorted list of cluster memberships. The cluster memberships are from  mycl <- cutree(hc, h=quantile(distance_matrix, .75)) of foo:
  clusmemb <-  mycl[hc$order]
  write.csv(cbind(scaledmatrix.sort,clusmemb), file = "scaledmatrix_sort_long.csv")
  
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
  
  # To map this onto our future heatmap, we need to map the different line's starting x and y coordinates
  # and ending x and y coordinates. We have Left vertical, Right verticle, minimal line, and then the cut points. 
  
  line = data.frame (
    Y = c(0.5, 0.5, 0.5, cluster.heights),
    Yend = c((nrow(scaledmatrix.sort)+0.5),(nrow(scaledmatrix.sort)+0.5), 0.5, cluster.heights),
    X = c(0.5, (ncol(scaledmatrix.sort)+0.5), rep(0.5,length(cluster.heights)+1)),
    Xend = c(0.5, (ncol(scaledmatrix.sort)+0.5), rep((ncol(scaledmatrix.sort)+0.5), length(cluster.heights)+1))
  )
  
  #convert wide form matrix to long form, required for heatmap plot
  scaledmatrix.sort.long <- melt(scaledmatrix.sort, id = rownames(scaledmatrix.sort))
  
  #renaming columns for ggplot ease of call
  colnames(scaledmatrix.sort.long)[colnames(scaledmatrix.sort.long)=="Var1"] <- "Predictors"
  colnames(scaledmatrix.sort.long)[colnames(scaledmatrix.sort.long)=="Var2"] <- "Outcomes"
  colnames(scaledmatrix.sort.long)[colnames(scaledmatrix.sort.long)=="value"] <- "Column\nZ-Score"
  
  # allows original values to label the cells instead of scaled values
  scaledmatrix.sort.long <- cbind(scaledmatrix.sort.long, pred.long)
  
  #create heatmap with rows as outcomns, columns as categorical variables, and colors as high to low from scaling
  heatmap.plot <- ggplot(data = scaledmatrix.sort.long, 
                         aes(x = Outcomes, 
                             y = Predictors)) +
    scale_fill_gradient2(low = "#fac877", 
                         mid = "white",
                         high = "#799250", 
                         midpoint = 0, 
                         space = "Lab",
                         na.value = "grey50", 
                         guide = "colourbar", 
                         aesthetics = "fill") +
    geom_tile(aes(fill = `Column\nZ-Score`)) +  
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggtitle("Heatmap of Predicted Probabilities") +
    theme(legend.position = "right", axis.text.x = element_text(angle = 90)) +
    geom_segment(data=line, aes(x=X, y=Y, xend=Xend, yend=Yend)) +
    geom_text(size = 4, aes(label = format(paste0(round(originalvalue, 2)*100,"%"))))
  print(heatmap.plot)
  ggsave(file="finished.svg", plot=heatmap.plot, width=6.56, height=7.49)

}
