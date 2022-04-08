# load packages
library(GenomicRanges)

# Generates a random GRanges object
#
# You can use this function to get a random interval list for testing.
# Don't change this function, it's just for you to use while
# you work on the other functions
#' @examples
# gr1 = randomGRanges()
randomGRanges = function(){
	starts = sample(1:1e7, 500)
	sizes = sample(1:1e4, 500)
	ends = starts + sizes
	ret = GenomicRanges::GRanges(seqnames="chrX", ranges=IRanges::IRanges(starts,ends))
	GenomicRanges::reduce(ret)
}

# Counts overlaps between query and database
#
# This function counts the number of intervals from a 
# database set of intervals that overlap with a single
# query interval. 
#
#' @param query  
#'   Interval to check for overlaps. Should be a list 
#'   object with 'start' and 'end' values.
#' @param database
#'   A data.frame, with rows corresponding to intervals,
#'   with 'starts' and 'ends' as columns.
#' @return
#'   Number of overlaps counted (as a numeric)
#' @examples
#'   query = list(start=50, end=60)
#'   database = data.frame(starts=c(35, 45, 55), ends=c(45, 65, 85))
#'   countOverlapsSimple(query, database)  # returns 2
countOverlapsSimple = function(query, database) {

  ## count overlap by filtering database
  # d = database[(database$starts<= query$end & database$ends >= query$start),]
  # n_overlaps = nrow(d)
  
  ## count overlap with sequential search, loop throgh the database
	# initiate number of overlaps
  n_overlaps = 0
  # check if query interval is well-formed
  if (query$start > query$end){
    stop("query interval is not well-formed: start > end")
  }
  # loop through database
  for (i in 1:nrow(database)){
    # check if intervals are well-formed
    if (database[i,]$starts > database[i,]$ends){
      stop(paste0("interval #", i, " is not well-formed: start > end"))
    }
    # check overlap
    if( database[i,]$starts <= query$end & database[i,]$ends >= query$start){
      n_overlaps = n_overlaps + 1
      }
  }
  
    return (n_overlaps)
}

# Measure the Jaccard similarity between two interval sets
#
# This function calculates the Jaccard Score between two
# given interval sets, provided as GenomicRanges::GRanges 
# objects. The Jaccard score is computed as the intersection
# between the two sets divided by the union of the two sets.
# @param gr1  First GRanges object
# @param gr2  Second GRanges object
# @return The Jaccard score (as numeric)
calculateJaccardScore = function(gr1, gr2){

	# get the intersection
  intersection = intersect(gr1, gr2)
  # get the union
  union = union(gr1, gr2)
  # jaccard_score = size of intersection / size of union
  jaccard_score = length(intersection)/length(union)
  # return 
  return (jaccard_score)
}

# Calculate pairwise Jaccard similarity among several interval sets
#
# This function makes use of \code{calculateJaccardScore}. It simply
# loops through each pairwise comparison and calculates.
#' Round the result to 3 significant figures using \code{signif}
#' @param lst A base R list of GRanges objects to compare
#' @return
#'   A matrix of size n-by-n, where n is the number of elements in lst.
#' @examples
#' lst = replicate(10, randomGRanges())
#' pairwiseJaccard(lst)
pairwiseJaccard = function(lst) {

  # get the number of interval sets
  n = length(lst)
  # initiate matrix
  matrix = matrix(data=NA, nrow=n, ncol=n)
  # fill in the matrix
  for(i in 1:n){
    for(j in 1:n){
      # calculate jaccard score
      jaccard_score = calculateJaccardScore(lst[[i]], lst[[j]])
      # Round the result to 3 significant figures
      matrix[i,j] = signif(jaccard_score, digits = 3)
    }
  }
  
  return (matrix)
}