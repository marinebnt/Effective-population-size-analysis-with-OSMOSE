#* beneat marine
#* 29/08/25
#* from osmose restart to a matrix of the number of mature individuals by age by time step


abund_mature_by_sp <- function(input_dir, outputDir_Ne=outputDir_Ne){

  conf = ns1
  
  # list the restart files
  path_nc <- paste0(input_dir, "/restart")
  nc_list <- list.files(path_nc, full.names = T)
  
  if(!length(nc_list)==0){
    # create output variable
    restart_name <- table(unlist(stringr::str_split(nc_list, pattern = "\\.nc."))[seq(1,2*length(nc_list), by=2)])
    nb_restarts_diff <- length(restart_name)
    nb_restarts_repl <- max(as.numeric(unlist(stringr::str_split(nc_list, pattern = "\\.nc."))[seq(2,2*length(nc_list), by=2)]))
    table = array(-99, dim = c(nb_restarts_diff, c(max_age+1), nb_sp, c(nb_restarts_repl+1)))
    dimnames(table)[[3]] <- as.character(sp_names)
    dimnames(table)[[4]] <- paste0("rep_", c(1:c(nb_restarts_repl+1)))
    
    i=0
    Time = c()
    k=0
    for (rest in names(restart_name)){
      dir <- nc_list[grep(pattern = (rest), nc_list)]
      k=k+1
      for (sp in species){
        for (rep in 1:nb_rep){
          i=i+1
          dir_rep <- dir[grep(pattern = (rep-1), dir)]
          
          colnames(table) <- c(0:max_age)
          Time[k] <- floor(as.numeric(stringr::str_split(stringr::str_split(dir_rep, "step")[[1]][2], paste0(".nc.", c(rep-1)))[[1]][1])/24)
          
          nc_tempo <- nc_open(dir_rep)
          mat      <- ncvar_get(nc_tempo, "maturity")
          spe      <- ncvar_get(nc_tempo, "species")
          abu      <- ncvar_get(nc_tempo, "abundance")
          age      <- floor(ncvar_get(nc_tempo, "age"))
          nc_close(nc_tempo)
          df_nc_restart <- data.frame("mat"=mat, "abu"=abu, "age"=age, "spe" = spe)
          df_nc_restart <- df_nc_restart[which(df_nc_restart$mat == 1),]
          namessp_i     <- length(species[1:i])-1
          df_nc_restart <- df_nc_restart[which(df_nc_restart$spe == namessp_i),]
          if(dim(df_nc_restart)[1]==0){next}
          abund_age     <- aggregate(df_nc_restart, abu~age, sum)
          abund_age     <- rbind(data.frame("age"=0:c(min(abund_age)-1), "abu"=rep(0, length(0:c(min(abund_age)-1)))), abund_age)
          missing       <- max_age-length(abund_age$age)+1
          abund_age_l30 <- c(abund_age$abu, rep(0, missing))
          
          table[k, , sp, paste0("rep_",rep)] <- abund_age_l30
        }
      }
    }
    dimnames(table)[[2]] <- 0:c(max_age)
    dimnames(table)[[1]] <- Time
    save(table, file=paste0(outputDir_Ne, "/abundance_mature.RData"))
  }
  return(table)
}

# 
# title = rep("", max_age)
# title[1] = c("Data coming from osmose restart : a matrix of the number of mature individuals by age (columns) by time step (rows)")
# table <- rbind(title,table)
# rownames(table) <- Time
# 
# write.csv()
