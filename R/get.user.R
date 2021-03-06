#############################################################
#########	MAKE USER-DEFINED EVALUATION GROUPS	#############
#############################################################

get.user <- function(occ.grp, bg.grp){
  if (class(occ.grp) == 'data.frame' | class(bg.grp) == 'data.frame'){ 
    warning('For user-defined partitioning, occ.grp and bg.grp must be supplied aas vectors.')
  }
  
	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)	
}
