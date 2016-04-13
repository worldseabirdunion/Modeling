##########################################################################
# PROGRAM: rf.modelSel (FOR CLASSIFICATION OR REGRESSION)
# USE: RANDOM FOREST MODEL SELECTION USING SCALED IMPORTANCE VALUES
# REQUIRES: DATAFRAME WITH Y-VAR AND X-VAR IN SEQUENTIAL COLS
#           R 2.13.0
#           randomForest 4.6-2
#
# ARGUMENTS: 
#       ydata      Y Data for model 
#       xdata      X Data for model
#       imp.scale  Type of scaling for importance values (mir or se), default is mir
#       r          Vector of importance percentiles to test i.e., c(0.1, 0.3, 0.5, 0.7, 0.9)
#       final      Run final model with selected variables (TRUE/FALSE)
#       plot.imp   Plot variable importance (TRUE/FALSE)
#       ...        Arguments to pass to randomForest (e.g., ntree=1000, 
#                  replace=TRUE, proximity=TRUE)
#
# VALUE:
#     A LIST WITH THE FOLLOWING OBJECTS
#         rf.final - FINAL RF MODEL (LIST OBJECT) USING SELECTED VARIABLES (IF final=TRUE)
#         SELVARS - LIST OF FINAL SELECTED VARIABLES
#         TEST - VALIDATION USED IN MODEL SELECTION
#         IMPORTANCE - IMPORTANCE VALUES FROM SELECTED MODEL
#
# NOTES: 
#        IF YOU WANT TO RUN CLASSIFICATION MAKE SURE Y IS A FACTOR
#        OTHERWISE RF RUNS IN REGRESSION MODE
#
#        The mir scale option perfroms a row standardization and the se option performs
#        normalizes using The “standard errors” of the permutation-based importance measure.
#        Both options result in a 0-1 range but se summs to 1.  
#                  mir = i/max(i)
#                  se = (i / se) / ( sum(i) / se) 
#
#        IMPORTANCE CANNONT BE FALSE AND IS SET IN THE FUNCTION, SO DON'T USE IMPORTANCE FLAG
#        
#        For regression the model selection criteria is; largest %variation 
#        explained, smallest MSE, and fewest parameters. 
#         
#        For classification; Smallest OOB error, smallest maximum within 
#        class error, and fewest parameters. 
#
# REFERENCES:
#    Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species 
#      Using Random Forest. Landscape Ecology 5:673-683.
#
#    Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas 
#      connectivity in Yellowstone National Park with landscape genetics. 
#      Ecology 91:252-261
#
#    Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species 
#      distribution and change using Random Forests CH.8 in Predictive Modeling in 
#      Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
# 
# EXAMPLES: 
# # CLASSIFICATION
#     data(iris)
#     iris$Species <- as.factor(iris$Species) 
#     ( rf.class <- rf.modelSel(iris[,1:4], iris[,"Species"], imp.scale="mir") )                            
# # REGRESSION
#     data(airquality)
#     airquality <- na.omit(airquality)
#     ( rf.regress <- rf.modelSel(airquality[,2:6], airquality[,1], imp.scale="se") )
#                               
# CONTACT: 
#     Jeffrey S. Evans 
#     Senior Landscape Ecologist 
#     The Nature Conservancy - Central Science
#     Adjunct Faculty
#     University of Wyoming
#     Laramie, WY
#     (970)672-6766
#     jeffrey_evans@tnc.org
##########################################################################
rf.modelSel <- function(xdata, ydata, imp.scale="mir", r=c(0.25, 0.50, 0.75),  
                        final=TRUE, plot.imp=TRUE, ...) 
  {
 if (!require (randomForest)) stop("randomForest PACKAGE MISSING")
 # IMPORTANCE SCALE FUNCTION
 rf.ImpScale <- function (x, scale="mir") { 
  if (!inherits(x, "randomForest")) 
       stop(deparse(substitute(x)), " Must be a randomForest object")
  if (x$type == "regression") {
   if (is.null(x$importanceSD) == TRUE | "%IncMSE" %in% 
       names(as.data.frame(x$importance)) == FALSE)
        stop("OBJECT DOES NOT CONTAIN PERMUTATED IMPORTANCE, PLEASE RUN 
       randomForest WITH importance=TRUE")   
    if (scale == "mir") {
      i <- ( as.data.frame(x$importance[,"%IncMSE"]) / 
		      max(x$importance[,"%IncMSE"]) )
		}	  
    if (scale == "se") {
	  i <- ( as.data.frame(x$importance[,"%IncMSE"]) / 
             as.data.frame(x$importanceSD) ) / 
             sum(as.data.frame(x$importance[,"%IncMSE"]) / 
             as.data.frame(x$importanceSD) )
        }
	 }
  if (x$type == "classification" | x$type == "unsupervised") {
   if (is.null(x$importanceSD) == TRUE | "MeanDecreaseAccuracy" %in% 
       names(as.data.frame(x$importance)) == FALSE)
        stop("OBJECT DOES NOT CONTAIN PERMUTATED IMPORTANCE, PLEASE RUN 
       randomForest WITH importance=TRUE")   
    if (scale == "mir") {
      i <- ( as.data.frame(x$importance[,"MeanDecreaseAccuracy"]) / 
		      max(x$importance[,"MeanDecreaseAccuracy"]) )
		}	  
    if (scale == "se") {
	  i <- ( as.data.frame(x$importance[,"MeanDecreaseAccuracy"]) / 
             as.data.frame(x$importanceSD[,"MeanDecreaseAccuracy"]) ) / 
             sum( as.data.frame(x$importance[,"MeanDecreaseAccuracy"]) / 
             as.data.frame(x$importanceSD[,"MeanDecreaseAccuracy"]) )
        }
	 }
    names(i) <- "importance"
   return( i )            
 }
RFtype <- is.factor(ydata) #TEST FOR FACTOR IN Y 
##CLASSIFICATION##
if (RFtype == "TRUE") {
  rf.all <- randomForest(x=xdata, y=ydata, importance=TRUE, ...)                    
       class.errors <- as.data.frame(rf.all$err.rate)
        class.errors <- na.omit(class.errors)  
         class.errors[class.errors == NaN] <- 0
          class.errors[class.errors == Inf] <- 1         
        i <- as.vector(array(0, dim=c((0),(1))))
       for ( l in 2:nlevels(as.factor(names(class.errors))) )
          {
          x.bar <- mean(class.errors[,l])              
            i <- as.vector(append(i, x.bar, after=length(i) ))
            }        
         max.error = max(i[2:length(i)] ) 
	imp <- rf.ImpScale(rf.all, scale=imp.scale) 
    results <- as.data.frame(array(0, dim=c( 0, 4 )))
      x <- c(0, (median(rf.all$err.rate[,"OOB"]) * 100), max.error * 100, dim(xdata)[2] )
    results <- rbind(results, x) 	 	 
     for (p in 1:length(r) ) {
		 t = quantile(imp[,1], probs=r[p])
         sel.imp <- subset(imp, importance >= t)
           sel.vars <- rownames(sel.imp)
     if (length( sel.vars ) > 1) {                             
         xdata.sub <- xdata[,sel.vars]       
      rf.model <- randomForest(x=xdata.sub, y=ydata, importance=TRUE, ...)          
           class.errors <- as.data.frame(rf.model$err.rate)
            class.errors <- na.omit(class.errors)  
             class.errors[class.errors == NaN] <- 0
              class.errors[class.errors == Inf] <- 1      
        i <- as.vector(array(0, dim=c((0),(1))))
       for ( l in 2:nlevels(as.factor(names(class.errors))) )
          {
          x.bar <- mean(class.errors[,l])              
            i <- as.vector(append(i, x.bar, after=length(i) ))
            }        
         max.error = max(i[2:length(i)] )     
         x <- c(t, median(rf.model$err.rate[,1]) * 100, max.error * 100, length(sel.vars) )
         results <- rbind(results, x)
         }
        }
  names(results) <- c("THRESHOLD", "OOBERROR", "CLASS.ERROR", "NPARAMETERS")
  results <- results[order(results$CLASS.ERROR, results$OOBERROR, results$NPARAMETERS),]  
    t <- as.vector(results[,"THRESHOLD"])[1]     
        sel.imp <- subset(imp, importance >= t)    
        sel.vars <- rownames(sel.imp)	                              
  } # END OF CLASSIFICATION    
##REGRESSION## 
if (RFtype == "FALSE") {      
 rf.all <- randomForest(x=xdata, y=ydata, importance=TRUE, ...) 
	imp <- rf.ImpScale(rf.all, scale=imp.scale) 
     results <- as.data.frame(array(0, dim=c( 0, 4 )))
     x <- c(0, (median(rf.all$rsq)), mean(rf.all$mse), dim(xdata)[2] )
     results <- rbind(results, x)     
   for (p in 1:length(r) ) {
        t = quantile(imp[,1], probs=r[p])		 
        sel.vars <- rownames(subset(imp, importance >= t))  
     if (length( sel.vars ) > 1) {                             
      xdata.sub <- as.data.frame(xdata[,sel.vars]) 
      rf.model <- randomForest(x=xdata.sub, y=ydata, importance=TRUE, ...)          
      x <- c(t, (median(rf.model$rsq)), mean(rf.model$mse), length(sel.vars) )
      results <- rbind(results, x)
      }
    }
   names(results) <- c("THRESHOLD", "VAREXP", "MSE", "NPARAMETERS")
   results <- results[order(-results$VAREXP, results$MSE, results$NPARAMETERS),]  
   t <- as.vector(results[,"THRESHOLD"])[1] 
      sel.imp <- subset(imp, importance >= t)    
      sel.vars <- rownames(sel.imp)	  
   } # END OF REGRESSION 
    if (plot.imp == TRUE) { 
	  p <- as.matrix(subset(imp, importance >= t))    
       ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]])  
       dotchart(p[ord,1], main="Model Improvment Ratio", pch=19)
    }
    if (final == TRUE) {
       sub.xdata <- xdata[,sel.vars]  #SUBSET VARIABLES
        rf.final <- randomForest(x=sub.xdata, y=ydata, importance=TRUE, ...)           
      ( list(MODEL=rf.final, SELVARS=sel.vars, TEST=results, IMPORTANCE=sel.imp) )      
         } else {
      ( list(SELVARS=sel.vars, TEST=results, IMPORTANCE=sel.imp) ) 
    }     
 }
