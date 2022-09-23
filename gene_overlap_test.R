gene_overlap_test = function(list1, list2, background,verbose=TRUE){
  n_A = length(list1)
  n_B = length(list2)
  n_C = length(background)
  n_A_B = length(intersect(list1,list2))
  hyp = phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE)
  jac = n_A_B/(n_A+n_B-n_A_B)
  if(verbose){print(paste0(n_A_B," matching from lists of lengths ",n_A," and ",n_B,"; p=",round(hyp,4)))}
  res = list(hypergeom = round(hyp,4),
             jaccard = jac,
             intersect = n_A_B)
  return(res)
}