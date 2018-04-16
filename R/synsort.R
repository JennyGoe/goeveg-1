#' Sorting functions for synoptic tables
#' 
#' @description Synoptic tables are a tool for interpretation of cluster species composition.
#' This functions provide sorting options for synoptic tables, sorting criteria can be either 
#' values in synoptic tables, such as frequencies, as well as combined criteria with considering
#' differential character, too.
#' Sorting algorithm aims to sort species in given cluster column order to blocked structure. 
#' Thereby, species with high frequencies and/or differential character are displayed blocked 
#' for each cluster or several neighbouring clusters.
#' 
#' @param syn1 Input synoptic table 1 (as dataframe) with priority entries for sorting. 
#' Usually dataframe from \code{link[goeveg] {syntable}} function output, 
#' but function should work with every synoptioc table input, as long as formats are 
#' appropriate. The values of this table will be displayed in the final output table. 
#' @param syn2 Optional second input table with additional sorting criteria. Note that 
#' values of this second input table will be considered in sorting, but not be displayed
#' in final synoptic table with method="allspec".
#' @param cluster Integer vector with classification cluster identity. Ensure matching order 
#' of cluster identity and samples in dataframe for correct allocation of cluster numbers to samples. 
#' @param method Sorting algorithm (method=c("allspec", "p_diff", "n_diff", pn_diff", "accspec", 
#' "all_diff")). See Details.
#' @param min1 Treshold minimum value for considering species of syn1 in ordering algorithm. 
#' Species below that minimum will neither be considered in algorithm nor displayed in final 
#' synoptic table, but will be listed in the "others" vector.  
#' @param min2 Treshold minimum value for considering species of syn2 in ordering algorithm. 
#' Species below that minimum will neither be considered in algorithm nor displayed in final 
#' synoptic table, but will be listed in the "others" vector.  
#' 
#' @details 
#' Sorting methods can consider numeric values of input synoptic tables only (total or percentage 
#' frequencies, method="allspec"), as well as combined sorting of frequencies and 
#' differential character (method=c("p_diff", "accspec", all_diff").
#' Sorting criteria can be either given by one input table, as well as two input tables with 
#' showing values of only the first, but considering the second in sorting algorithm.
#' 
#' Method="allspec" sorts all species exceeding minimum tresholds (min1, min2). The algorithm first 
#' detects maximum values for each species frequencies among clusters. After that, species are sorted 
#' by descending frequencies. With combining two sorting criteria, both 
#' input table values are considered, but sorting is prior to the values of syn1 (see examples).
#'     
#' With including differential species character as sorting criterion (method=("p_diff", 
#' "n_diff", pn_diff", "accspec", "all_diff")), input table1 must be numeric, the second one with 
#' information on differential character (outoput from \code{link[goeveg] {syntable}} function with 
#' type="diffspec"). Again, algorithm detects highest cluster values of species in syn1 as base for sorting, 
#' but will sort them considering differentiating character criterion (from second 
#' input table syn2). Species with high values in syn1 AND differential character will then be listed 
#' on the top of a species block. Within differentiating species, prevalence of diagnostic character 
#' is considered by favouring positive and/or cluster-specific differential character. 
#'  
#' 
#' @section Output
#' Output is a list of result components, including sorting method ($output), species sorting criteria 
#' ($species), sample sizes in clusters ($samplesize), sorted final synoptic tables with values from 
#' input table 1 ($syntable) and other species, that failed to be included in the final table due to 
#' treshold values given by min1 and min2.
#' In case of combined sorting with considering differential species character, a table with differential 
#' character of species is additionally contained in result list ($differential).
#' 
#' 
#' @references 
#' Tsiripidis, I., Bergmeier, E., Fotiadis, G. & Dimopoulos, P. 2009: A new algorithm for the determination of differential taxa. - Journal of Vegetation Science 20: 233-240.
#' 
#' @author Jenny Schellenberg (\email{jschell@gwdg.de})
#' 
#' @examples 
#' # Synoptic table of Scheden vegetation data: 
#' library(cluster)
#' pam1 <- pam(schedenveg, 4)
#' 
#' # 1) creating unordered synoptic tables
#' unordered <- syntable(schedenveg, pam1$clustering, abund = "perc", type = "percfreq")   # for unordered synoptiv percentage frequency table
#' differential <- syntable(schedenveg, pam1$clustering, abund = "perc", type = "diffspec") # differential species analysis
#' phitable <- syntable(schedenveg, pam1$clustering, abund = "perc", type = "phi") # for fidelity phi


#' # 2) sort synoptic tables
#' # Common synoptic table: sort by percentage frequency, show only species with minimum frequency of 25% in table
#' sorted <- synsort(syn1 = unordered$syntable, cluster = pam1$clustering, method = "allspec", min1 = 25)
#' 
#' # Sorted synoptic table, with only positive differentiating species with minimum 25% frequency in table  
#' positive <- synsort(syn1 = unordered$syntable, syn2 = differential$syntable, cluster = pam1$clustering, 
#' method = "p_diff", min1 = 25)
#' 
#' # Sorted complete synoptic table with percentage frequency (only species >25%) and differential character. 
#' complete <- synsort(syn1 = unordered$syntable, syn2 = differential$syntable, cluster = pam1$clustering, 
#' method = "all_diff", min1 = 25)
#' 
#' # Sorted synoptic table of species with minimum phi-value of 0.3, show percentage frequency
#' phi_complete <- synsort(syn1 = unordered$syntable, syn2 = phitable$syntable, cluster = pam1$clustering, method = "allspec", min1 = 25, min2 = 0.3)
#' # sort table by phi values, show only species with phi>0.3
#' phi_table <- synsort(syn1 = phitable$syntable, cluster = pam1$clustering, method = "allspec", min1 = 0.3)
#' 
#' # Synoptic table showing phi-values of species with differential character
#' phi_diff <- synsort(syn1 = phitable$syntable, syn2 = differential$syntable, cluster = pam1$clustering, 
#' method = "all_diff", min1 = 0.3)










synsort <- function(syn1, syn2 = syn1 , cluster, method, min1, min2) {

library(Hmisc)
  
if (method == "allspec") {
      if (all(syn2 == syn1)) {
      frames <- list()
      for (i in 1:length(unique(cluster))) {
        all <- syn1[apply(syn1,1,max) >= min1,]
        frames[[i]] <- assign(paste0("frame",i), all[apply(all,1,max) == all[,i],]) 
        frames[[i]] <- frames[[i]][order(frames[[i]][,i], decreasing=TRUE),] }
      for ( i in 2:max(cluster)) {
      duprows <- rownames(frames[[i]]) %in% rownames(frames[[1]])
      frames[[1]] <- rbind(frames[[1]], frames[[i]][!duprows,]) }
      allspec <- frames[[1]]
      results <- list("output" = "synoptic table sorted by values of one input table",
                      "species" = paste0("species with value in input table Minimum", sep=" ", min1, ", others listet below"),
                      "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                      "syntable" = allspec,
                      "others" = sort(rownames(syn1[apply(syn1,1,max)<=min1,])))
      print(results)
      return(results)
      } else {
        frames <- list()
        frames2 <- list() 
        frames3 <- list()
        for (i in 1:length(unique(cluster))) {
        all <- syn1[apply(syn1,1,max)>=min1,]
        frames[[i]] = assign(paste0("frame",i), all[apply(all,1,max) == all[,i],]) 
        
        all2 = syn2[apply(syn2,1,max)>=min2,] 
        frames2[[i]] = assign(paste0("frame2_",i), all2[apply(all2,1,max) == all2[,i],])   
        
        frames3[[i]] = assign(paste0("frame3_",i), merge(frames[[i]], frames2[[i]], by = "row.names"))
        rownames(frames3[[i]]) = frames3[[i]][,1]
        frames3[[i]] = frames3[[i]][,-1]
        frames3[[i]] = frames3[[i]][,1:length(unique(cluster))]
        names(frames3[[i]]) = c(sort(unique(cluster))) 
        frames3[[i]] = frames3[[i]][order (frames3[[i]][,i], decreasing=TRUE),]    }
        allspec <- unique(do.call("rbind", frames3))
        results <- list("output" = "synoptic table sorted by values of two input tables, values of first input table shown",
                      "species" = paste0("species with minimim value of", sep=" ", min1, " in table 1 and Minimum of", sep=" ", min2, 
                                         " in input table 2, others listet below"),
                      "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                      "syntable" = allspec,
                      "others" = sort(rownames(syn1)[rownames(syn1) %nin% rownames(allspec)]))
        print(results)
        return(results)}

} else if (method =="p_diff" | method == "accspec" | method == "all_diff") {
      frames <- list()
      for (i in 1:length(unique(cluster))) {
      all = syn1[apply(syn1,1,max) >= min1,]
      frames[[i]] = assign(paste0("frame",i), all)}
      
      synord <- syn2
      for ( i in 1:length(syn2[1,])) {
        for ( k in 1:length(syn2[,1])) {
          if (synord[k,i] == "p" & syn1[k,i] >= min1) {synord[k,i] <- 1}
          else {synord[k,i] <- 0}   }}
      synord <- data.frame(lapply(synord, as.numeric))
      rownames(synord) <- rownames(syn2)  
      names(synord) <- seq(1,length(unique(cluster)),1)
      pos_synord <- synord
      posidiff <- apply(pos_synord,1,sum)
      
      
      di <- c("n", "pn", "-")
      for ( l in 1:3) {
       synord <- syn2
       for ( i in 1:length(syn2[1,])) {
         for ( k in 1:length(syn2[,1])) {
           if (synord[k,i] == di[l]) {synord[k,i] <- 1}
           else {synord[k,i] <- 0}   }}
       synord <- data.frame(lapply(synord, as.numeric))
       rownames(synord) <- rownames(syn2)  
       names(synord) <- seq(1,length(unique(cluster)),1)
      if ( l == 1) {
        neg_synord <- synord
        negdiff <- apply(neg_synord,1,sum)
      } else if ( l == 2 ) {
        posneg_synord <- synord
        posnegdiff <- apply(posneg_synord,1,sum)
      } else {
        acc_synord <- synord
        accdiff <- apply(acc_synord,1,sum)
      }     }
      
      if (method == "accspec" | method == "all_diff") {
          frames5 <- list()
               pos_in5 <- synord[accdiff == max(accdiff),]
               in5 <- pos_in5
               if ( length(pos_in5[,1]) == 0 ) {
                 } else { diff_in5 <- pos_in5
                     for (i in 1:length(diff_in5[,1])) {
                     diff_in5[i,] <- in5[rownames(in5) == rownames(diff_in5[i,]),]} 
               }      # evtl anpassen an groessere datensaetze
                
    
          if (method == "all_diff") {
              synacc <- pos_in5
                  (if ( length(diff_in5[,1]) == 0 ) {diffspec <- diff_in5}
                  else { diffspec <- diff_in5
                  for ( i in 1:length(diffspec[1,])) {
                  for ( k in 1:length(diffspec[,1])) {
                       if (diffspec[k,i]=="1") {diffspec[k,i] <- "-"}
                       else {diffspec[k,i] <- "NA"}   }} }  )
              diffspecacc <- diffspec 
              accdiff <- diff_in5
          } else {
              syntab   <- pos_in5
              diffspec <- diff_in5
              for ( i in 1:length(diffspec[1,])) {
              for ( k in 1:length(diffspec[,1])) {
                  if (diffspec[k,i]=="1") {diffspec[k,i] <- "-"}
                  else {diffspec[k,i] <- "NA"}   }}
              results <- list("output" = "synoptic table sorted by values of numeric input table and non-differential, accompanying species",
                    "species" = paste0("species with minimim value of", sep=" ", min1, " in input numeric table and that are ACCOMPANYING"),
                    "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                    "syntable" = syntab,
                    "differential" = diffspec,
                    "others" = sort(rownames(syn1)[rownames(syn1) %nin% rownames(syntab)]))
              print(results)
              return(results)
          }
      } 
} else {stop("Sorting of synoptic table failed: wrong method entry.")}


vektor <- list(posidiff, negdiff, posnegdiff)
Rahmen <- list(pos_synord, neg_synord, posneg_synord)

for ( m in 1:3)  {
    posdiff <- vektor[[m]]
    rahmen <- Rahmen[[m]] 
    
    frames1 <- list()
    in1 <- rahmen[posdiff == 1,] 
    for (i in 1:length(unique(cluster))) {
      frames1[[i]] = assign(paste0("frame1_",i), in1[apply(in1,1,max) == in1[,i],])
      frames1[[i]] = assign(paste0("frame1_",i), merge(frames[[i]], frames1[[i]], by = "row.names"))
      rownames(frames1[[i]]) = frames1[[i]][,1]
      frames1[[i]] = frames1[[i]][,-1]
      frames1[[i]] = frames1[[i]][,1:length(unique(cluster))]
      names(frames1[[i]]) = c(sort(unique(cluster)))
      frames1[[i]] = frames1[[i]][order(frames1[[i]][,i], decreasing=TRUE),]  }
    pos_in1 <- do.call("rbind", frames1)
        (if (length(pos_in1[,1])==0) {diff_in1 <- pos_in1
         } else {
           diff_in1 <- pos_in1
           for (i in 1:length(diff_in1[,1])) {
               diff_in1[i,] <- in1[rownames(in1) == rownames(diff_in1[i,]),]}   })
    
    frames2 <- list()
    in2 <- rahmen[posdiff == 2,]
    for (i in 1:length(unique(cluster))) {
    frames2[[i]] = assign(paste0("frame2_",i), in2[apply(in2,1,max) == in2[,i],])
    frames2[[i]] = assign(paste0("frame2_",i), merge(frames[[i]], frames2[[i]], by = "row.names"))
    rownames(frames2[[i]]) = frames2[[i]][,1]
    frames2[[i]] = frames2[[i]][,-1]
    frames2[[i]] = frames2[[i]][,1:length(unique(cluster))]
    names(frames2[[i]]) = c(sort(unique(cluster)))
    frames2[[i]] = frames2[[i]][order(frames2[[i]][,i], decreasing=TRUE),]  }
    pos_in2 <- unique(do.call("rbind", frames2))
    (if (length(pos_in2[,1])==0) {diff_in2 <- pos_in2 }
     else {diff_in2 <- pos_in2
     for (i in 1:length(diff_in2[,1])) {
      diff_in2[i,] <- in2[rownames(in2) == rownames(diff_in2[i,]),]}  })
  
    frames3 <- list()
    in3 <- rahmen[posdiff == 3,]
    for (i in 1:length(unique(cluster))) {
    frames3[[i]] = assign(paste0("frame3_",i), in3[apply(in3,1,max) == in3[,i],])
    frames3[[i]] = assign(paste0("frame3_",i), merge(frames[[i]], frames3[[i]], by = "row.names"))
    rownames(frames3[[i]]) = frames3[[i]][,1]
    frames3[[i]] = frames3[[i]][,-1]
    frames3[[i]] = frames3[[i]][,1:length(unique(cluster))]
    names(frames3[[i]]) = c(sort(unique(cluster))) 
    frames3[[i]] = frames3[[i]][order(frames3[[i]][,i], decreasing=TRUE),]  }
    pos_in3 <- unique(do.call("rbind", frames3))
    (if (length(pos_in3[,1])==0) {diff_in3 <- pos_in3 }
     else { diff_in3 <- pos_in3
     for (i in 1:length(diff_in3[,1])) {
      diff_in3[i,] <- in3[rownames(in3) == rownames(diff_in3[i,]),]}   })
  
    frames4 <- list()
    in4 <- rahmen[posdiff == 4,]
    for (i in 1:length(unique(cluster))) {
    frames4[[i]] = assign(paste0("frame4_",i), in4[apply(in4,1,max) == in4[,i],])
    frames4[[i]] = assign(paste0("frame4_",i), merge(frames[[i]], frames4[[i]], by = "row.names"))
    rownames(frames4[[i]]) = frames4[[i]][,1]
    frames4[[i]] = frames4[[i]][,-1]
    frames4[[i]] = frames4[[i]][,1:length(unique(cluster))]
    names(frames4[[i]]) = c(sort(unique(cluster))) 
    frames4[[i]] = frames4[[i]][order(frames4[[i]][,i], decreasing=TRUE),]  }
    pos_in4 <- unique(do.call("rbind", frames4))
    (if (length(pos_in4[,1])==0) {diff_in4 <- pos_in4 }
     else {diff_in4 <- pos_in4
    for (i in 1:length(diff_in4[,1])) {
      diff_in4[i,] <- in4[rownames(in4) == rownames(diff_in4[i,]),]}   })
  
    frames5 <- list()
    in5 <- rahmen[posdiff == 5,]
    for (i in 1:length(unique(cluster))) {
    frames5[[i]] = assign(paste0("frame5_",i), in5[apply(in5,1,max) == in5[,i],])
    frames5[[i]] = assign(paste0("frame5_",i), merge(frames[[i]], frames5[[i]], by = "row.names"))
    rownames(frames5[[i]]) = frames5[[i]][,1]
    frames5[[i]] = frames5[[i]][,-1]
    frames5[[i]] = frames5[[i]][,1:length(unique(cluster))]
    names(frames5[[i]]) = c(sort(unique(cluster)))
    frames5[[i]] = frames5[[i]][order(frames5[[i]][,i], decreasing=TRUE),]  }
    pos_in5 <- unique(do.call("rbind", frames5))
    (if (length(pos_in5[,1])==0) {diff_in5 <- pos_in5 }
    else {diff_in5 <- pos_in5
    for (i in 1:length(diff_in5[,1])) {
      diff_in5[i,] <- in5[rownames(in5) == rownames(diff_in5[i,]),]} })
  
    frames6 <- list()
    in6 <- rahmen[posdiff == 6,]
    for (i in 1:length(unique(cluster))) {
    frames6[[i]] = assign(paste0("frame6_",i), in6[apply(in6,1,max) == in6[,i],])
    frames6[[i]] = assign(paste0("frame6_",i), merge(frames[[i]], frames6[[i]], by = "row.names"))
    rownames(frames6[[i]]) = frames6[[i]][,1]
    frames6[[i]] = frames6[[i]][,-1]
    frames6[[i]] = frames6[[i]][,1:length(unique(cluster))]
    names(frames6[[i]]) = c(sort(unique(cluster)))
    frames6[[i]] = frames6[[i]][order(frames6[[i]][,i], decreasing=TRUE),]  }
    pos_in6 <- unique(do.call("rbind", frames6))
    (if (length(pos_in6[,1])==0) {diff_in6 <- pos_in6 }
     else {diff_in6 <- pos_in6
     for (i in 1:length(diff_in6[,1])) {
      diff_in6[i,] <- in6[rownames(in6) == rownames(diff_in6[i,]),]}  })
  
    frames7 <- list()
    in7 <- rahmen[posdiff == 7,]
    for (i in 1:length(unique(cluster))) {
    frames7[[i]] = assign(paste0("frame7_",i), in7[apply(in7,1,max) == in7[,i],])
    frames7[[i]] = assign(paste0("frame7_",i), merge(frames[[i]], frames7[[i]], by = "row.names"))
    rownames(frames7[[i]]) = frames7[[i]][,1]
    frames7[[i]] = frames7[[i]][,-1]
    frames7[[i]] = frames7[[i]][,1:length(unique(cluster))]
    names(frames7[[i]]) = c(sort(unique(cluster)))
    frames7[[i]] = frames7[[i]][order(frames7[[i]][,i], decreasing=TRUE),]  }
    pos_in7 <- unique(do.call("rbind", frames7))
    (if (length(pos_in7[,1])==0) { diff_in7 <- pos_in7 }
     else {diff_in7 <- pos_in7
    for (i in 1:length(diff_in7[,1])) {
      diff_in7[i,] <- in7[rownames(in7) == rownames(diff_in7[i,]),]}  })
    
    frames8 <- list()
    in8 <- rahmen[posdiff >= 8,]
    for (i in 1:length(unique(cluster))) {
      frames8[[i]] = assign(paste0("frame8_",i), in8[apply(in8,1,max) == in8[,i],])
      frames8[[i]] = assign(paste0("frame8_",i), merge(frames[[i]], frames8[[i]], by = "row.names"))
      rownames(frames8[[i]]) = frames8[[i]][,1]
      frames8[[i]] = frames8[[i]][,-1]
      frames8[[i]] = frames8[[i]][,1:length(unique(cluster))]
      names(frames8[[i]]) = c(sort(unique(cluster)))
      frames8[[i]] = frames8[[i]][order(frames8[[i]][,i], decreasing=TRUE),]  }
    pos_in8 <- unique(do.call("rbind", frames8))
    (if (length(pos_in8[,1])==0) { diff_in8 <- pos_in8 }
      else {diff_in8 <- pos_in8
      for (i in 1:length(diff_in8[,1])) {
        diff_in8[i,] <- in8[rownames(in8) == rownames(diff_in8[i,]),]}  })
 

    if (method == "p_diff") {
     syntab   <- rbind(pos_in1, pos_in2, pos_in3, pos_in4, pos_in5, pos_in6, pos_in7, pos_in8)
     diffspec <- rbind(diff_in1, diff_in2, diff_in3, diff_in4, diff_in5, diff_in6, diff_in7, diff_in8)
     for ( i in 1:length(diffspec[1,])) {
         for ( k in 1:length(diffspec[,1])) {
           if (diffspec[k,i]=="1") {diffspec[k,i] <- "p"}
           else {diffspec[k,i] <- "-"}   }}
       results <- list("output" = "synoptic table sorted by values of numeric input table and differential species and synoptic table of POSITIVE differentiating species",
                    "species" = paste0("species with minimim value of", sep=" ", min1," in input numeric table and that are POSITIVE differential species (algorithm of Tsiripidis et al.) in one ore max. 7 clusters"),
                    "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                    "syntable" = syntab,
                    "differential" = diffspec,
                    "others" = sort(rownames(syn1)[rownames(syn1) %nin% rownames(syntab)]))
       print(results)  
       return(results)
       break
         
    } else if (method == "all_diff") {
          if (identical(pos_synord, Rahmen[[m]])) {
          synpos <- rbind(pos_in1, pos_in2, pos_in3, pos_in4, pos_in5, pos_in6, pos_in7, pos_in8)
          diffspec <- rbind(diff_in1, diff_in2, diff_in3, diff_in4, diff_in5, diff_in6, diff_in7, diff_in8)
          for ( i in 1:length(diffspec[1,])) {
            for ( k in 1:length(diffspec[,1])) {
              if (diffspec[k,i]=="1") {diffspec[k,i] <- "p"}
              else {diffspec[k,i] <- "-"}   }} 
          diffspecpos <- diffspec
          
          } else if ( identical(neg_synord, Rahmen[[m]])) {
          synneg <- rbind(pos_in1, pos_in2, pos_in3, pos_in4, pos_in5, pos_in6, pos_in7, pos_in8)
          diffspec <- rbind(diff_in1, diff_in2, diff_in3, diff_in4, diff_in5, diff_in6, diff_in7, diff_in8)
               if (length(diffspec[,1]) == 0) { diffspecneg <- diffspec
               } else {
                  for ( i in 1:length(diffspec[1,])) {
                    for ( k in 1:length(diffspec[,1])) {
                        if (diffspec[k,i]=="1") {diffspec[k,i] <- "n"}
                        else {diffspec[k,i] <- "-"}   }} }
          diffspecneg <- diffspec
    
          } else if ( identical(posneg_synord, Rahmen[[m]]))  {
          synposneg <- rbind(pos_in1, pos_in2, pos_in3, pos_in4, pos_in5, pos_in6, pos_in7, pos_in8)
          diffspec <- rbind(diff_in1, diff_in2, diff_in3, diff_in4, diff_in5, diff_in6, diff_in7, diff_in8)
               if (length(diffspec[,1]) == 0) { diffspecposneg <- diffspec 
               } else {
                 for ( i in 1:length(diffspec[1,])) {
                  for ( k in 1:length(diffspec[,1])) {
                    if (diffspec[k,i]=="1") {diffspec[k,i] <- "pn"}
                    else {diffspec[k,i] <- "-"}   }}   }
            diffspecposneg <- diffspec
          } else {stop("Failure merging frames in sorting algorithm")}
    } else {}
}

if (method == "all_diff") {
 m1 <- merge(diffspecpos, diffspecneg, by= "row.names", all.x = TRUE, sort=F)
 rownames(m1) <- m1[,1]
 m1 <- m1[,-1]
 m1[is.na(m1)] <- "NA"
 for ( i in 1:length(unique(cluster))) {
    for ( l in 1:length(m1[,1])) {
        if (m1[l,i] == "-" & m1[l,(i+length(unique(cluster)))] == "n") {m1[l,i] = "n"}
        else if (m1[l,i] == "p") {m1[l,i] = "p"}
        else {m1[l,i] = "NA"}
      } } 
 m1 <- m1[,1:length(unique(cluster))]
    
 m2 <- merge(m1, diffspecposneg, by= "row.names", all.x = TRUE, sort=F)
 rownames(m2) <- m2[,1]
 m2 <- m2[,-1]
 m2[is.na(m2)] <- "NA"
 for ( i in 1:length(unique(cluster))) {
      for ( l in 1:length(m2[,1])) {
        if (m2[l,i] == "NA" & m2[l,(i+length(unique(cluster)))] == "pn") {m2[l,i] = "pn"}
        else if (m2[l,i] == "p") {m2[l,i] = "p"}
        else if (m2[l,i] == "n") {m2[l,i] = "n"}
        else {m2[l,i] = "-"}
      } } 
 m2 <- m2[,1:length(unique(cluster))]
 names(m2) <- seq(1:length(unique(cluster)))
 m3 <- rbind(m2, diffspecacc)
  
 completetable <- merge(m3, all, all.x=TRUE, by= "row.names", sort=F)
 rownames(completetable) <- completetable[,1]
 completetable <- completetable[,-1]
      name <- c("")
      for(i in 1:length(unique(cluster))) {
      name[i] = paste0("diff ", sort(unique(cluster))[i])
      name[i+length(unique(cluster))] = paste0("cluster ", sort(unique(cluster))[i])   }   
      names(completetable) <- name    
      
      frames <- list()
      for (i in 1:length(unique(cluster))) {
      frames[[i]] <- assign(paste0("frame",i), completetable[
        apply(completetable[,(min(sort(unique(cluster)))+length(unique(cluster))):
                              (max(sort(unique(cluster)))+length(unique(cluster)))],
              1,max) == completetable[,(i+length(unique(cluster)))],]) 
      frames[[i]] <- frames[[i]][order(frames[[i]][,c(i+length(unique(cluster)),i)], decreasing=TRUE),]
      frames[[i]] <- frames[[i]] [complete.cases(frames[[i]]),] }
    
      completetable <- unique(do.call("rbind", frames))
      results <- list("output" = "complete synoptic table, sorted by values of numeric input table and differential species character",
                    "species" = paste0("species with minimim value of", sep=" ", min1, " and their differentiating character"),
                    "samplesize" = tapply(rep(1,length(cluster)),cluster,sum),
                    "syntable" = completetable,
                    "others" = sort(rownames(syn1)[rownames(syn1) %nin% rownames(completetable)]))
      print(results)
} else {stop("Failed to create sorted synoptic table; check input and parameters of formula.")} 
  
}
