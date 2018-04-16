#' Draw rank-abundance curves
#'
#' @description This function draws a rank-abundance curve for community data.
#' If you wish to draw multiple rank-abundance curves for selected samples use \link{racurves}.
#'
#' @param matrix Community data, a matrix-like object with samples in rows.
#' @param main The main title.
#' @param ylog If set on \code{TRUE} the y-axis is displayed on a log-scale.
#' @section Details:
#' Rank abundance curves or Whittaker plots (see \cite{Whittaker 1965}) are used to display relative species abundance as biodiversity component.
#' They are a means to visualise species richness and species evenness.
#' @examples
#' ## Draw simple rank-abundance curve
#' racurve(meadows)
#'
#' ## Draw simple rank-abundance curve with log-scaled axis
#' racurve(meadows, ylog=T)
#' @seealso \code{\link{racurves}}
#' @references Whittaker, R. H. (1965). Dominance and Diversity in Land Plant Communities: Numerical relations of species express the importance of competition in community function and evolution. \emph{Science} \strong{147 :} 250-260.
#' @author Friedemann Goral \email{fgoral@gwdg.de}
#' @export

racurve <-  function(matrix, main = "Rank-abundance diagram", nlab = 0, ylog = FALSE, frequency = FALSE, ylim = NULL, xlim = NULL) {
  if(!is.data.frame(matrix)) {
    matrix <- data.frame(matrix)
  }
  
  freq <- sort(apply(matrix>0, 2, sum), decreasing = T)
  
  abund <- sort(apply(matrix, 2, sum), decreasing = T)
  sum(abund)
  rel.abund <- sort((abund / sum(abund)), decreasing = T)
  
  labels <- names(head(rel.abund, n = nlab))
  
  if(frequency == FALSE) {
    
    if(ylog == TRUE) {
      plot(rel.abund, xlab="Abundance Rank", ylab="Relative Abundance",
           main=main, log="y", ylim = ylim, xlim = xlim)
      lines(rel.abund)
      if(nlab != 0) {
        text(head(rel.abund, n = nlab), labels = labels, pos = 4, cex = 0.7)
      }
      
    } else {
      plot(rel.abund, xlab="Abundance Rank", ylab="Relative Abundance",
           main=main, , ylim = ylim, xlim = xlim)
      lines(rel.abund)
      if(nlab != 0) {
        text(head(rel.abund, n = nlab), labels = labels, pos = 4, cex = 0.7)
      }
    }
  } else if(frequency == TRUE) {
    
    if(ylog == TRUE) {
      plot(freq, xlab="Frequency Rank", ylab="Frequency",
           main=main, log="y", ylim = ylim, xlim = xlim)
      lines(freq)
      if(nlab != 0) {
        text(head(freq, n = nlab), labels = labels, pos = 4, cex = 0.7)
      }
      
    } else {
      plot(freq, xlab="Frequency Rank", ylab="Frequency",
           main=main, ylim = ylim, xlim = xlim)
      lines(freq)
      if(nlab != 0) {
        text(head(freq, n = nlab), labels = labels, pos = 4, cex = 0.7)
      }
    }
  }
  
  out <- list(abund = abund, rel.abund = rel.abund, freq = freq)
  invisible(out)
}

