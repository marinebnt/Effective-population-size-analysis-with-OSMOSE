
# Auxiliary functions for lists manipulation

# function to merge lists of dataframes/matrices into single dataframe/matrix
myf = function(x) do.call(cbind, x)

# function to merge lists of dataframes/matrices into single dataframe/matrix
myf2 = function(x) {
  for(i in seq_along(x)) {
    for(j in seq_along(x[[i]])) {
      check2 = any(sapply(x[[i]][[j]], is_tibble))
      if(check2) {
        for(k in seq_along(x[[i]][[j]])) {
          # print(x[[i]][[j]][[k]])
          x[[i]][[j]][[k]] = as.numeric(x[[i]][[j]][[k]][[2]])
          # print(x[[i]][[j]][[k]])
        }
      }
      check = is.matrix(x[[i]][[j]]) | is.data.frame(x[[i]][[j]])
      if(isTRUE(check)) {
        x[[i]][[j]] <- as.data.frame(x[[i]][[j]])
        x[[i]][[j]]$cohort = rownames(x[[i]][[j]])
        x[[i]][[j]]$variable = names(x[[i]])[[j]]
        x[[i]][[j]] = pivot_longer(x[[i]][[j]], !c(cohort, variable), names_to = "replicate", values_to = "value")
        x[[i]][[j]]$replicate = stringr::str_replace(x[[i]][[j]]$replicate, pattern = "rep_", replacement = "")
        # print(sprintf("Expanding grid %s (%d) %s (%d)", names(x)[[i]], i, names(x[[i]])[j], j)) 
      }
      check3 = any(sapply(x[[i]][[j]], is.matrix))
      if(isTRUE(check3)) {
        # print(sprintf("Skipping %s (%d) %s (%d)", names(x)[[i]], i, names(x[[i]])[j], j)) 
        next
      }
      tmp = try(myf(x[[i]][[j]]), silent = TRUE)
      if(!inherits(tmp, "try-error") & !check) {
        # print(sprintf("Merging %s (%d) %s (%d)", names(x)[[i]], i, names(x[[i]])[j], j))
        x[[i]][[j]] = tmp 
      }
    }
  }
  return(x)
}

# when the element of the list is not a vector but a matrix of cohort x time
myf3 = function(x) {
    Nb_H_l = expand.grid(cohort=(1:ncol(x))+(years-ncol(x)),replicate=as.numeric(gsub(names(x), pattern="rep_", replacement="")))
    Nb_H_l$value = as.vector(INIT2b[[SP]]$data_Nb$Nb)
    Nb_H_l$variable = "Nb"
}

# function to reshape list into a nested list by species and variable
reshape_list = function(x) {
  xxx = unlist(x, recursive = FALSE)
  nm = lapply(lapply(strsplit(names(xxx), split = ".", fixed = TRUE), FUN="[", i=c(2,1)), FUN=paste, collapse=".")
  names(xxx) = nm
  # xxx = xxx[sort(names(xxx))]
  
  sp = sapply(strsplit(names(xxx), split=".", fixed=TRUE), FUN = "[", i=1)
  var = sapply(strsplit(names(xxx), split=".", fixed=TRUE), FUN = "[", i=2)
  
  out = setNames(vector(length(unique(sp)), mode="list"), nm=unique(sp))
  for(i in seq_along(sp)) {
    out[[sp[i]]][[var[i]]] = xxx[[i]]
  }
  return(out)
}

# function to reshape list into a dataframe
reshape_data = function(x, varName="var_name", as.variable=FALSE) {
  if(!is.list(x)) return(x)
  nm = try(strsplit(names(x), split=".", fixed=TRUE), silent=TRUE)
  if(inherits(nm, "try-error")) return(x)
  ll = sapply(nm, length)
  exit = TRUE
  if(all(ll==1)) {
    dat = expand.grid(cohort=rownames(x), replicate=unlist(nm))
    exit=FALSE
  }
  if(all(ll==2)) {
    dat = data.frame(cohort=sapply(nm, FUN = "[", i=2), 
                     replicate=sapply(nm, FUN = "[", i=1))
    exit=FALSE
  }
  if(isTRUE(exit)) return(x)
  if(isTRUE(as.variable)) {
    dat$variable = varName
    varName = "value"
  }
  dat[[varName]] = unlist(x)
  dat$cohort = as.numeric(gsub(dat$cohort, pattern="cohort_", replacement=""))
  dat$replicate = as.numeric(gsub(dat$replicate, pattern="rep_", replacement=""))
  return(dat)
}

# function to reshape all elements of a nested list
reshape_all = function(x) {
  for(i in seq_along(x)) {
    for(j in seq_along(x[[i]])) {
      if(!is.list(x[[i]][[j]])) next
      x[[i]][[j]] = lapply(x[[i]][[j]], FUN = reshape_data)
      }
  }
  return(x)
}