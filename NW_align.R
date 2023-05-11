# setwd('') #remember to define a working directory


#### Dependencies check ####

if(!require(plot.matrix)){ # this package is needed for output visualization
  install.packages("plot.matrix")
  library(plot.matrix)
}

library(parallel)    # this is already included in base R  installation



#### functions for matrices initialization ####

init_scorematrix<-function(seq1, seq2, gap) { # initialize score matrix function 
  
  mat <- matrix(0, nrow =nchar(seq2)+1 , ncol = nchar(seq1)+1) # matrix of 0s
  
  #naming rows/cols according to sequences
  names_seq1<- strsplit(seq1, '')[[1]]
  names_seq2<- strsplit(seq2, '')[[1]]
  
  colnames(mat)<- c('', names_seq1)
  rownames(mat) <-c('', names_seq2)
  
  # initialization values
  vecrow <- -seq(from = 0,
                 by = abs(gap),
                 length.out = length(mat[1,])) # horizontal gap in first row
  
  vecol <- -seq(from = 0,
                by = abs(gap),
                length.out = length(mat[,1])) # vertical gap in first column
  
  mat[1,] <- vecrow
  mat[,1] <- vecol
  
  return(mat)
}



init_trbmatrix <- function(seq1, seq2) { # initialize matrix of directions function
  
  mat <- matrix('', nrow =nchar(seq2)+1, ncol = nchar(seq1)+1) # matrix of empty strings 
  
  #naming rows/cols according to sequences
  names_seq1<- strsplit(seq1, '')[[1]]
  names_seq2<- strsplit(seq2, '')[[1]]
  
  colnames(mat)<- c('', names_seq1)
  rownames(mat) <-c('', names_seq2)
  
  # initialization values
  vecrow <- replicate(length(mat[1,]), 'H') # horizontal gap in first row
  mat[1,] <- vecrow
  
  vecol <- replicate(length(mat[,1]), 'V') # vertical gap in first column
  mat[,1] <- vecol
  
  mat[1,1] <- '0'
  
  return(mat)
}




#### Score and Movement matrices computing #####

Nw_scores <- function(seq1,
                      seq2, 
                      match, 
                      mismatch, 
                      gap) { # takes the 2 sequences to align and the scoring system
  
  
  # create 2 matrices
  score_matr <- init_scorematrix(seq1,seq2, gap)
  score_matr
  trb_matr<- init_trbmatrix(seq1,seq2)
  trb_matr
  
# loop preliminary operations  
  j <- 1
  i <- 1
  d_score <- 0
  h_score <- 0 
  v_score <- 0
  
# real loop 
  for (j in 1:nchar(seq2)) {
    for(i in 1:nchar(seq1)) { ## computes possible values of the cell
      
      if (rownames(score_matr)[j+1] == colnames(score_matr)[i+1]) { # diagonal (match/mismatch) score
        d_score <- score_matr[j,i] + match 
      } else {d_score<-score_matr[j,i]+ mismatch}  
      
      h_score <- score_matr[j+1,i] + gap   # horizontal gap score
      v_score <- score_matr[j,i+1] + gap   # vertical gap score
      
      # assign score
      score <- max(h_score, d_score, v_score)
      #score <- max(scorevec)  # value to put in the cell 
      
      score_matr[j+1,i+1] <- score # score matrix filling
      
      # Movement matrix computing
      if (score == d_score) {trb_matr[j+1, i+1] <- paste(trb_matr[j+1, i+1], 'D', sep="")}
      if (score == v_score) {trb_matr[j+1, i+1] <- paste(trb_matr[j+1, i+1], 'V', sep="")}
      if (score == h_score) {trb_matr[j+1, i+1] <- paste(trb_matr[j+1, i+1], 'H', sep="")}
      
    }
  }
  
  trb_matr <- gsub('NA','', trb_matr)
  
  return(list(score_matr, trb_matr)) # score and directions matrix
}



#### Traceback follow #####

threading <- function(traceback_matrix, 
                      i = length(traceback_matrix[,1]), 
                      j = length(traceback_matrix[1,]), 
                      aln_list = list('','','')) { # Function takes the directions matrix 
                                                   # and computes all possible paths  
                                                     # others are default values
                                                     # (do not modify them unless you know what you're doing) 
  
  
  
  if (traceback_matrix[i,j] == '0') { # If the recursion has arrived at first cell 
    
    # stop and store the path in a temporary file
    write.table(aln_list,file = 'tempFile_Aln.txt',
                append = T,
                row.names = F,
                col.names = F)   
    
    
    
  } else { # If still elsewhere inside the matrix
           ## store paths starting from the bottom left and going backwards
            ## (so paths will need to be reversed)
    
    k <- 1
    while (k <= nchar(traceback_matrix[i,j])) { #Iterate through each move in that cell
      move <- unlist(strsplit(traceback_matrix[i,j], ''))[[k]]
      
      
      if(k==1) { ## First character:
        
        if (move == "D") { #If a move is D...
          aln_list[[k]] <- paste(as.character(aln_list[[k]]), "D", sep="") #Append that to the path
          threading(traceback_matrix,i-1,j-1,aln_list[[k]]) #Send the function to the diagonal cell
        }
        
        if (move == "V") { #If a move is V...
          aln_list[[k]] <- paste(as.character(aln_list[[k]]), "V", sep="") #Append that to the path
          threading(traceback_matrix,i-1,j,aln_list[[k]]) #Send the function to the above cell
        }
        
        if (move == "H") { #If a move is H...
          aln_list[[k]] <- paste(as.character(aln_list[[k]]), "H", sep="") #Append that to the path
          threading(traceback_matrix,i,j-1,aln_list[[k]]) #Send the function to the left cell
        }
      }
      
      
      if(k==2) { ## Second character: (can only be H or V)
        
       # if (move == "D") { #If a move is D...
        #  aln_list[[k]] <- paste(substr(as.character(aln_list[[k-1]]),
         #                               start = 0,
          #                              stop = ifelse(nchar(as.character(aln_list[[k-1]]))-1>0,
           #                                           nchar(as.character(aln_list[[k-1]]))-1, 
            #                                          0)), 
             #                    "D",
              #                   sep="") #Append that to the path
          
        #  threading(traceback_matrix,i-1,j-1,aln_list[[k]]) #Send the function to the diagonal cell
        #}
        
        if (move == "V") { #If a move is V...
          aln_list[[k]] <- paste(substr(as.character(aln_list[[k-1]]),
                                        start = 0, 
                                        stop = ifelse(nchar(as.character(aln_list[[k-1]]))-1>0,
                                                      nchar(as.character(aln_list[[k-1]]))-1, 
                                                      0)), # previous path less last move
                                 "V", 
                                 sep="")  # Append that to the path
          
          threading(traceback_matrix,i-1,j,aln_list[[k]]) #Send the function to the above cell
        }
        
        if (move == "H") { #If a move is H...
          aln_list[[k]] <- paste(substr(as.character(aln_list[[k-1]]),
                                        start = 0, 
                                        stop = ifelse(nchar(as.character(aln_list[[k-1]]))-1>0,
                                                      nchar(as.character(aln_list[[k-1]]))-1, 
                                                      0)), # previous path less last move
                                 "H", 
                                 sep="") #Append that to the path
          
          threading(traceback_matrix,i,j-1,aln_list[[k]]) #Send the function to the left cell
        }
      }
      
      
      if(k==3) { ## Third character: (can only be H)
        
        #if (move == "D") { #If a move is D...
         # aln_list[[k]] <- paste(substr(as.character(aln_list[[k-1]]),
          #                              start = 0,
           #                             stop = ifelse(nchar(as.character(aln_list[[k-1]]))-1>0,
            #                                          nchar(as.character(aln_list[[k-1]]))-1, 
             #                                         0)),
              #                   "D",
               #                  sep="") #Append that to the path
          
          #threading(traceback_matrix,i-1,j-1,aln_list[[k]]) #Send the function to the diagonal cell
        #}
        
        #if (move == "V") { #If a move is V...
         # aln_list[[k]] <- paste(substr(as.character(aln_list[[k-1]]),
          #                              start = 0, 
           #                             stop = ifelse(nchar(as.character(aln_list[[k-1]]))-1>0,
            #                                          nchar(as.character(aln_list[[k-1]]))-1, 
             #                                         0)),
              #                   "V", 
               #                  sep="")  # Append that to the path
          
          #threading(traceback_matrix,i-1,j,aln_list[[k]]) #Send the function to the above cell
        #}
        
        if (move == "H") { #If a move is H...
          aln_list[[k]] <- paste(substr(as.character(aln_list[[k-1]]),
                                        start = 0, 
                                        stop = ifelse(nchar(as.character(aln_list[[k-1]]))-1>0,
                                                      nchar(as.character(aln_list[[k-1]]))-1, 
                                                      0)), # previous path less last move
                                 "H", 
                                 sep="") #Append that to the path
          threading(traceback_matrix,i,j-1,aln_list[[k]]) #Send the function to the left cell
        }
      }
      k<-k+1
    }
  }
}


## Traceback computing ####

nw_traceback <- function(matrix_path, 
                         traceback_matrix) { # takes the generated file of moves  and
                                             # translates them into aligned sequences
  
  temp <- nchar(matrix_path)
  
  i <- 1 # index of traceback matrix columns
  j <- 1 # index of traceback matrix rows
  c <- 0 # index for steps
  line1 <- '' # sequence 1
  line2 <- '' # sequence 2
  
  
  for (c in 1:temp) { # revese the path and loop each move
    step <- rev(unlist(strsplit(matrix_path, '')))[[c]]
    
    # translation
    if(step == 'D') { #
      line1 <- paste(line1, colnames(traceback_matrix)[i])
      line2 <- paste(line2, rownames(traceback_matrix)[j])
      i <- i+1
      j <- j+1 
    }
    
    if(step == 'V') { # vertical gap doesn't change column
      line1 <- paste(line1,'-')
      line2 <- paste(line2, rownames(traceback_matrix)[j])
      j <- j+1
    }
    
    if(step == 'H') { # horizontal gap doesn't change row
      line1 <- paste(line1, colnames(traceback_matrix)[i])
      line2 <- paste(line2,'-')
      i<- i+1
    }
  }
  
  
  
  suppressWarnings(write.table(matrix(c(line1, line2, ''),          # adds the pairs of aligned sequences 
               nrow = 3,                           # to a predefined file
               dimnames = list(c('seq1', 'seq2','______________________________________________'),  
                               'Possible alignment')),                
        file = 'Test_output_final.txt', 
        row.names = T,
        col.names = T,
        append = T))
  
  
  return(matrix(c(line1, line2),
                nrow = 2, 
                dimnames = list(c('Sequence 1:', 'Sequence 2:'),'Possible alignment'))) # console return
                                                        #ex:
                                                          #   Possible alignment
                                                          # seq1: 'AA-TC[...]'
                                                          # seq2: 'A-GCC[...]'
}


#### Definitive ####

G_align <- function(seq1, seq2, match, mismatch,gap) {
  
  res <- Nw_scores(seq1, seq2, match, mismatch, gap) # returns score and movement matrices 
  
  if (file.exists('tempFile_Aln.txt')) { # delete (permanetly) the file if exists
    file.remove('tempFile_Aln.txt')
  }
  moves <- threading(res[[2]])
  directions <-list(t(read.table('tempFile_Aln.txt')))
  
  global_score <- res[[1]][length(res[[1]][,1]), length(res[[1]][1,])] # last cell of score matrix
  pox_aligns <- length(directions[[1]])
  
  write.table(matrix(c(length(res[[1]][1,])-1,  # initialize output file 
                       length(res[[1]][,1])-1,  # aligned sequences will then be appended here
                       global_score,
                       pox_aligns, 
                       ''),
                     nrow = 5,                                        
                     dimnames = list(c('Sequence 1 length:',
                                       'Sequence 2 length:',
                                       'Global Alignment Score:',
                                       'Possible alignments:',
                                       '___________________________________________'),            
                                     'SUMMMARY')),          
              file = 'Test_output_final.txt',
              row.names = T,
              col.names = T,
              append = F)
  
  
  return(mclapply(directions[[1]],res[[2]][-1,-1],
           FUN = 'nw_traceback'))
  
}


G_align(str1, str2, 1,-1,-1)






