# Collection of functions as a library
# function for checking small deletion
var.del <- function(y,i) {
  del.ref <- as.character(y$Reference)
  del.alt <- as.character(y$Alternate)
  if (nchar(del.ref) > nchar(del.alt)) {
    alt.seq <- substr(del.ref,1+nchar(del.alt),nchar(del.ref))
    if (itv$Ref[i]==alt.seq) {
      return(y)
    } else {
      return(NULL)
    }           # outdt <- rbind(outdt,eachrc)
  }
}

# function for checking small insersion
var.ins <- function(x,i) {
  del.ref <- as.character(x$Reference)
  del.alt <- as.character(x$Alternate)
  if (nchar(del.ref) < nchar(del.alt)) {
    alt.seq <- substr(del.alt,1+nchar(del.ref),nchar(del.alt))
    if (itv$Alt[i]==alt.seq) {
      return(x)
    } else {
      return(NULL)
    }
  }       #outdt <- rbind(outdt,eachrc)
}

# function for SNV checking SNV
var.snv <- function (z,i) {
  del.ref <- as.character(z$Reference)
  del.alt <- as.character(z$Alternate)
  if ((nchar(del.ref) == 1) && (nchar(del.alt) == 1)) {
    if (itv$Alt[i]==z$Alternate && itv$Ref[i]==z$Reference)
      return(z)
    else
      return(NULL)
  }
}

# Main codes
var.match <- function(itv, gndt) {
  rc <- dim(itv)[1]
  outdt <- NULL
  for (i in 1:rc) {
    if (itv$Alt[i]=="-") { # small deletion variants
      gpos <- itv$Start[i]-1 # relocate the position for small deletion variants
      deldt <- gndt[gndt$Position==gpos,] # genomic location matching
      checkpoint <- dim(deldt)[1]
      if (checkpoint==1) { # check one or multiple records matched
        eachrc <- var.del(deldt,i)
        outdt <- rbind(outdt,eachrc) 
      } else if (checkpoint>1) {
        for (j in 1:checkpoint) {
          eachrc <- var.del(deldt[j,],i)
          outdt <- rbind(outdt,eachrc)
        }
      }
    } else if (itv$Ref[i]=="-") { # insersion
      deldt <- gndt[gndt$Position==itv$Start[i],]
      checkpoint <- dim(deldt)[1]
      if (checkpoint==1) {
        eachrc <- var.ins(deldt,i)
        outdt <- rbind(outdt,eachrc)
      } else if (checkpoint>1) {
        for (k in 1:checkpoint) {
          eachrc <- var.ins(deldt[k,],i)
          outdt <- rbind(outdt,eachrc)
        }
      }
    } else { #SNV
      deldt <- gndt[gndt$Position==itv$Start[i],]
      checkpoint <- dim(deldt)[1]
      if (checkpoint==1) {
        eachrc <- var.snv(deldt,i)
        outdt <- rbind(outdt,eachrc)
      } else if (checkpoint>1) {
        for (l in 1:checkpoint) {
          eachrc <- var.snv(deldt[l,],i)
          outdt <- rbind(outdt,eachrc)
        }
      }
    }
  }
  return(outdt)
}
