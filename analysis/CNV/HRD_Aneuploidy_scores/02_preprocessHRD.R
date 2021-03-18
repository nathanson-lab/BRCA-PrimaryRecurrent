#' Define additional properties of the sequencing data necessary to compute HRD
#' 
#' @param seq.dat the raw sequencing data
#' @param ref the reference genome to use (grch37, grch38)
#' @details define chromosome size, arms, telomere position, and centromere position from reference data; 
#' compute allelic imbalance and segment size
#' 
#' @return seq.dat
#' 
#'@examples
#'seq.dat <- sub01.segments[ sub01.segments$chromosome == "chr1",]
#'
#'seq.dat <- preprocessHRD( seq.dat )
#'seq.dat
#'@export

# ------------------------------- preprocessHRD ---------------------------------------- #
# define allelic imbalance (AI), telomere positions, segment length, cross arms
# these are used in the various HRD scores
preprocessHRD <- function( seq.dat,  ref )
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   seq.dat (data.frame), the sequencing data with more factors added
{
  
  if( ref == "grch37" )
  {
    ref.dat <- grch37.ref.dat
  } else
    if( ref == "grch38" )
    {
      ref.dat <- grch38.ref.dat
    } else
    {
      print(paste(ref, "is not a valid reference genome.", sep = " "))
      stop("select one of: grch37, grch38")
    }
  
  # check the sequencing data input to make sure it has the data we need
  if( any(!(seq.cols.needed %in% colnames(seq.dat))))
  {
    print(paste("column", 
                seq.cols.needed[which(!(seq.cols.needed %in% colnames(seq.dat)))],
                "missing in the seq data.", sep=" "))
    print(paste("Looking for columns:", seq.cols.needed, sep=" "))
    print(paste("Found:", seq.cols.needed[seq.cols.needed %in% colnames(seq.dat)]))
    stop("Column mismatch. Exiting.")
  } 
  
  # remove rows w/ missing entries
  if( any( is.na(seq.dat$CNt) ))
  {
    print(paste("Row: ", which(is.na(seq.dat$CNt)), " contains NA; removing.", sep = ""))
    seq.dat <- subset(seq.dat, !is.na(seq.dat$CNt))
  }
  
  levels(seq.dat$chromosome) <- levels(ref.dat$chromosome)
  
  seq.dat$frac.chr <- (seq.dat$end.pos - seq.dat$start.pos) / ref.dat$chr.size[match(seq.dat$chromosome, ref.dat$chromosome)]
  
  # segment length
  seq.dat$brk.len <- 0
  seq.dat <- combineSeg(seq.dat, 3e06)
  
  # matches reference data to the correct chromosome in the subject data
  key = match(seq.dat$chromosome, ref.dat$chromosome)
  
  seq.dat$seg.len <- seq.dat$end.pos - seq.dat$start.pos
  seq.dat$chr.size <- ref.dat$chr.size[match(seq.dat$chromosome, ref.dat$chromosome)]
  
  seq.dat$AI <- (seq.dat$A > seq.dat$B) & (seq.dat$A != 1) & (seq.dat$B != 0)
  
  seq.dat$start.arm <- seq.dat$end.arm <- rep("NA", dim(seq.dat)[1])
  
  ind = seq.dat$start.pos - ref.dat$centromere.start[key] > 1
  seq.dat$start.arm[ind]  = "p"
  seq.dat$start.arm[!ind] = "q"
  
  ind = seq.dat$end.pos - ref.dat$centromere.end[key] > 1
  seq.dat$end.arm[ind]  = "p"
  seq.dat$end.arm[!ind] = "q"
  
  seq.dat$cross.arm <- seq.dat$start.arm != seq.dat$end.arm
  
  
  # define end telomeres
  seq.dat$post.telomere <- (seq.dat$start.pos - ref.dat$p.telomere.end[key] <= 1000)
  seq.dat$pre.telomere  <- (ref.dat$q.telomere.start[key] - seq.dat$end.pos <= 1000)
  # ---
  
  return(seq.dat)
}
# ------------------------------------------------------------------------------------- #