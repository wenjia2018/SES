printind = function(object){ 
  
  varnames <- object$var.names 
  
  cat("\n")
  cat("DIRECTION DEPENDENCE ANALYSIS: Independence Properties", "\n", "\n") 
  # ------------------------------------------------------------------------------------------- Print Target Model:
  cat(paste("Target Model:", varnames[2], "->", varnames[1], sep = " "), "\n", "\n")
  cat("Homoscedasticity tests:", "\n") 
  
  sigtests1.yx <- rbind( c(object$breusch_pagan[[1]]$statistic, object$breusch_pagan[[1]]$parameter, object$breusch_pagan[[1]]$p.value),  
                         c(object$breusch_pagan[[2]]$statistic, object$breusch_pagan[[2]]$parameter, object$breusch_pagan[[2]]$p.value)
  )
  sigtests1.yx <- round(sigtests1.yx, 4)
  rownames(sigtests1.yx) <- c("BP-test", "Robust BP-test")
  colnames(sigtests1.yx) <- c("X-squared", "df", "p-value")
  print.default(format( sigtests1.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
  cat("\n")
  
  if( !is.null(object$nlfun) ){
    
    sigtests2.yx <- rbind(object$nlcor.yx$t1, object$nlcor.yx$t2, object$nlcor.yx$t3)
    sigtests2.yx <- round(sigtests2.yx, 4)		 
    
    if(is.na(suppressWarnings(as.numeric(object$nlcor.yx$func)))){
      cat(paste("Non-linear correlation tests:", object$nlcor.yx$func, "transformation"))
      
      rownames(sigtests2.yx) <- c(paste("Cor[", object$nlcor.yx$func, "(", "r_", varnames[1], "), ", varnames[2],"]", sep=""),
                                  paste("Cor[", "r_", varnames[1], ", ", object$nlcor.yx$func, "(", varnames[2], ")]", sep=""),
                                  paste("Cor[", object$nlcor.yx$func, "(", "r_", varnames[1], "), ", object$nlcor.yx$func, "(", varnames[2],")]", sep="")
      )
    } else{
      cat(paste("Non-linear correlation tests: Power transformation using", object$nlcor.yx$func))
      
      rownames(sigtests2.yx) <- c(paste("Cor[", "r_", varnames[1], "^", object$nlcor.yx$func, ", ", varnames[2],"]", sep=""),
                                  paste("Cor[", "r_", varnames[1], ", ", varnames[2], "^", object$nlcor.yx$func, "]", sep=""),
                                  paste("Cor[", "r_", varnames[1], "^", object$nlcor.yx$func, ", ", varnames[2], "^", object$nlcor.yx$func, "]", sep="")
      )
    }
    cat("\n")
    
    colnames(sigtests2.yx) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    print.default(format( sigtests2.yx, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    
    cat("\n")
  }	 
  
  if( !is.null(object$hsic.method) ){
    
    if(object$hsic.method[1] == "gamma") cat("Hilbert-Schmidt Independence Criterion: Gamma-Approximation", "\n")
    if(object$hsic.method[1] == "boot") cat(paste("Hilbert-Schmidt Independence Criterion: Bootstrap-Approximation", " (", object$hsic.method[2], " resamples)", sep = ""), "\n")
    
    cat(paste("HSIC = ", round(object$hsic.yx$statistic, 4), " p-value = ", round(object$hsic.yx$p.value, 4), sep = ""))
    cat("\n", "\n")
  } 
  
  # ------------------------------------------------------------------------------------------- Print Alternative Model:
  
  cat(paste("Alternative Model:", varnames[1], "->", varnames[2], sep = " "), "\n", "\n")
  
  cat("Homoscedasticity tests:", "\n") 
  
  sigtests1.xy <- rbind( c(object$breusch_pagan[[3]]$statistic, object$breusch_pagan[[3]]$parameter, object$breusch_pagan[[3]]$p.value),  
                         c(object$breusch_pagan[[4]]$statistic, object$breusch_pagan[[4]]$parameter, object$breusch_pagan[[4]]$p.value)
  )
  sigtests1.xy <- round(sigtests1.xy, 4)
  rownames(sigtests1.xy) <- c("BP-test", "Robust BP-test")
  colnames(sigtests1.xy) <- c("X-squared", "df", "p-value")
  print.default(format( sigtests1.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
  cat("\n")
  
  if(!is.null(object$nlfun)){
    
    sigtests2.xy <- rbind(object$nlcor.xy$t1, object$nlcor.xy$t2, object$nlcor.xy$t3)
    sigtests2.xy <- round(sigtests2.xy, 4)
    
    if(is.na(suppressWarnings(as.numeric(object$nlcor.xy$func)))){
      cat(paste("Non-linear correlation tests:", object$nlcor.xy$func, "transformation"))
      
      rownames(sigtests2.xy) <- c(paste("Cor[", object$nlcor.xy$func, "(", "r_", varnames[2], "), ", varnames[1],"]", sep=""),
                                  paste("Cor[", "r_", varnames[2], ", ", object$nlcor.xy$func, "(", varnames[1], ")]", sep=""),
                                  paste("Cor[", object$nlcor.xy$func, "(", "r_", varnames[2], "), ", object$nlcor.xy$func, "(", varnames[1],")]", sep="")
      )
    } else{
      cat(paste("Non-linear correlation tests: Power transformation using", object$nlcor.xy$func))
      
      rownames(sigtests2.xy) <- c(paste("Cor[", "r_", varnames[2], "^", object$nlcor.xy$func, ", ", varnames[1],"]", sep=""),
                                  paste("Cor[", "r_", varnames[2], ", ", varnames[1], "^", object$nlcor.xy$func, "]", sep=""),
                                  paste("Cor[", "r_", varnames[2], "^", object$nlcor.xy$func, ", ", varnames[1], "^", object$nlcor.xy$func, "]", sep="")
      )
    }
    cat("\n")
    
    colnames(sigtests2.xy) <- c("estimate", "t-value", "df", "Pr(>|t|)")
    print.default(format( sigtests2.xy, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
    
    cat("\n")
  }
  
  if(!is.null(object$hsic.method)){
    
    if(object$hsic.method[1] == "gamma") cat("Hilbert-Schmidt Independence Criterion: Gamma-Approximation", "\n")
    if(object$hsic.method[1] == "boot") cat(paste("Hilbert-Schmidt Independence Criterion: Bootstrap-Approximation", " (", object$hsic.method[2], " resamples)", sep = ""), "\n")
    
    cat(paste("HSIC = ", round(object$hsic.xy$statistic, 4), " p-value = ", round(object$hsic.xy$p.value, 4), sep = ""))
    cat("\n", "\n")
  }
  
}

printresid = function(object) {
  {
    
    varnames <- object$var.names 
    
    cat("\n")
    cat("DIRECTION DEPENDENCE ANALYSIS: Residual Distributions", "\n", "\n") 
    cat("Skewness and kurtosis tests:", "\n") 
    
    sigtests <- rbind( c(object[[1]]$target$statistic, object[[1]]$target$p.value, object[[1]]$alternative$statistic, object[[1]]$alternative$p.value),  
                       c(object[[2]]$target$statistic, object[[2]]$target$p.value, object[[2]]$alternative$statistic, object[[2]]$alternative$p.value)  
    )
    sigtests <- round(sigtests, 4)
    rownames(sigtests) <- c("Skewness", "Kurtosis")
    colnames(sigtests) <- c("target", "z-value", "Pr(>|z|)", "alternative", "z-value", "Pr(>|z|)")
    print.default(format( sigtests, digits = max(3L, getOption("digits") - 3L), scientific = NA, scipen = 999), print.gap = 2L, quote = FALSE)
    
    if(is.null(object$boot.args)){
      cat("\n")
      cat("Skewness and kurtosis difference tests:", "\n") 
      
      citests <- rbind(object$skewdiff, object$kurtdiff)
      citests <- round(citests, 4)
      rownames(citests) <- c("Skewness", "Kurtosis")
      colnames(citests) <- c("diff", "z-value", "Pr(>|z|)")
      print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
      cat("\n")
    }
    
    if(!is.null(object$boot.args)){
      ci.level <- as.numeric(object$boot.args[2]) * 100
      cat("\n")
      cat("Skewness and kurtosis difference tests:", "\n") 
      
      citests <- rbind(object$skewdiff, object$kurtdiff)
      citests <- round(citests, 4)
      rownames(citests) <- c("Skewness", "Kurtosis")
      colnames(citests) <- c("diff", "z-value", "Pr(>|z|)", "lower", "upper")
      print.default(format( citests, digits = max(3L, getOption("digits") - 3L)), print.gap = 2L, quote = FALSE)
      cat("\n")
      cat(paste("Number of resamples:", object$boot.args[3]))
      cat("\n")
      if(object$boot.args[1] == "bca") cat(ci.level, "% ", "BCa bootstrap CIs reported", "\n", "\n", sep = "") 
      if(object$boot.args[1] == "perc") cat(ci.level, "% ", "Percentile bootstrap CIs reported", "\n", "\n", sep = "") 
    }
    
    cat("---")
    cat("\n")
    cat(paste("Note: Target is", varnames[2], "->", varnames[1], sep = " "))
    cat("\n")
    cat(paste("      Alternative is", varnames[1], "->", varnames[2], sep = " "))
    cat("\n")
  }
}
