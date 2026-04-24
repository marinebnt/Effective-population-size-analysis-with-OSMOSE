library(ncdf4)
library(tidyr)

main = "New_gpu_conf.R"
source("osmose-ns-v4.4.1_clim/R/functions/source.R")

 # list_dirs <- list.dirs(path="osmose-ns-v4.4.1_clim/output/test", recursive=TRUE, full.names=TRUE)
list_restart <- list_dirs[grep(pattern="restart", x=list_dirs)]

# marine trophic index documentation 
#https://www.un.org/esa/sustdev/natlinfo/indicators/methodology_sheets/oceans_seas_coasts/marine_trophic_index.pdf
#https://d1wqtxts1xzle7.cloudfront.net/104639264/regional-libre.pdf?1690769155=&response-content-disposition=inline%3B+filename%3DThe_state_of_biodiversity_and_fisheries.pdf&Expires=1769972201&Signature=aGzH1Q1Dv~EiLAS~QK31xavSwlRZVAAqOyi3Ysk-b0j7VevkhyMJrzACxmOcfpltW2LZ~2WIrOP4z4t-HxiZj6CHug0RDrKYDuLhiUoYOd8ZEf5n1uXd6e8vZ-b1TPpfAt91spErhCofrHTeOnA28MtCzgihjfFaKz-JKCIyg9E84W~Xc1OEK2w-sTaZYbF2o8Ytsj8j2jiCzxjIdqnkYqahvrqCpoTgJgZfM~VonWi-S5A6av4100D5oCuEWuNwy8TXEwWnGMZ-w1Z~ChWFjkasDuh0LK3GVeE~T-4E3zR7XDNbro2dSWcptLTnxaDhB2c52a96opWpKJJyaQCO9g__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA#page=45
#https://pmc.ncbi.nlm.nih.gov/articles/PMC1569461/pdf/rstb20041597.pdf



for (re in list_restart[5]){
    # what restart file his in the corresponding restart folder ?
    files_re <- list.files(path=re, pattern="snapshot", full.names=TRUE) 
    name = stringr::str_split(files_re,"conf_")[[1]][2]
    name = stringr::str_split(name,"_sp")[[1]][1]

    for (rep in (0:19)[1]){
        # what file has the corresponding replicate number ?
        files_re_rep <- files_re[grep(pattern=paste0("nc.", rep), x=files_re)]

        #read .csv file corresponding to the abundance by TL
        df_abund_TL <- read.csv(file=paste0(re, "/../Indicators/ns_abundanceDistribByTL_Simu", rep, ".csv"), skip=c(1))
        df_abund_TL <- pivot_longer(df_abund_TL, cols= -c(Time, TL), names_to = "species_id", values_to = "abundance")
        df_abund_TL$species <- (as.numeric(df_abund_TL$species_id=="AtlanticMackerel")+
                                as.numeric(df_abund_TL$species_id=="Sandeel")*2+
                                as.numeric(df_abund_TL$species_id=="EuropeanSprat")*3+ 
                                as.numeric(df_abund_TL$species_id=="NorwayPout")*4+
                                as.numeric(df_abund_TL$species_id=="EuropeanPlaice")*5+
                                as.numeric(df_abund_TL$species_id=="NorwayPout")*6+
                                as.numeric(df_abund_TL$species_id=="CommonSole")*7+
                                as.numeric(df_abund_TL$species_id=="Saithe")*8+
                                as.numeric(df_abund_TL$species_id=="AtlanticCod")*9+
                                as.numeric(df_abund_TL$species_id=="Haddock")*10+
                                as.numeric(df_abund_TL$species_id=="HorseMackerel")*11+
                                as.numeric(df_abund_TL$species_id=="Whiting")*12+
                                as.numeric(df_abund_TL$species_id=="CommonDab")*13+
                                as.numeric(df_abund_TL$species_id=="GreyGurnard")*14+
                                as.numeric(df_abund_TL$species_id=="Hake")*15+
                                as.numeric(df_abund_TL$species_id=="Shrimp")*16)
        tab <- table(df_abund_TL$species, df_abund_TL$species_id)
        tab_sp <- data.frame("id"=rownames(which(tab>0, arr.ind=TRUE)), 
                             "species"=colnames(tab[,which(tab>0, arr.ind=TRUE)[,2]]))
        tab_sp$id <- as.numeric(tab_sp$id)
        tab_sp <- tab_sp[order(tab_sp[,1]),]
        for (fi in files_re_rep[1]){
            
            rep_i = ifelse(fi == files_re_rep[1], rep, rep+20)

            # what time step of osmose is this nc file ? 
            step <- stringr::str_split(stringr::str_split(fi, pattern="_step")[[1]][2], pattern="\\.nc")[[1]][1]
            step_div <- floor(as.numeric(step)/24)
            
            df_abund_TL_t <- df_abund_TL[df_abund_TL$Time == step_div,]
            
            # read the restart, compare the abundance by TL with the one in the .csv file and extract alo the biomass by TL
            nc_data <- nc_open(fi)
            vals_df <- data.frame("TL" = ncvar_get(nc_data, "trophiclevel"))
            vals_df$species <- ncvar_get(nc_data, "species")
            vals_df$weight <- ncvar_get(nc_data, "weight")
            vals_df$abundance <- ncvar_get(nc_data, "abundance")
            vals_df$mature <- ncvar_get(nc_data, "maturity")
            vals_df$age <- ncvar_get(nc_data, "age")
            vals_df$gweight <- ncvar_get(nc_data, "gonadWeight")
            nc_close(nc_data)
            vals_df$biomass <- vals_df$weight * vals_df$abundance
            vals_df <- tibble(vals_df)

            # check that abundance by TL are the same in both files :not working well
            # vals_df_summary <- vals_df
            # vals_df_summary$TL <- round(vals_df_summary$TL, 1)
            # vals_df_summary <- aggregate(abundance ~ TL+ species, data=vals_df_summary, sum)
            # merged_df <- merge(vals_df_summary, as.data.frame(df_abund_TL_t), by = c("TL", "species"))
            # merged_df$range_diff <-merged_df$abundance.x/merged_df$abundance.y
            # merged_df$range_diff[is.infinite(merged_df$range_diff)] <- NA
            # merged_df[which(merged_df$range_diff>1),] 

            eggw <- unlist(osmose::get_par(conf, "egg.weight"))
            names(eggw) <- stringr::str_remove(names(eggw), pattern="species.egg.weight.sp")
            eggw_df <- as.data.frame(eggw)
            eggw_df$species <- as.numeric(rownames(eggw_df))
            vals_df <- vals_df %>% left_join(eggw_df, by="species")
            vals_df$fec <- (vals_df$gweight)/(vals_df$eggw)


            vals_df_std <- vals_df %>% 
                group_by(species) %>%
                mutate(weight = weight/mean(weight, na.rm=TRUE),
                       fec = fec/mean(fec, na.rm=TRUE))

            ggplot( vals_df_std, aes(x=weight, y=fec)) + 
                geom_line() + xlab("weight") + 
                ylab("fecundity") + 
                ggtitle(paste0("File  ", re)) +
                geom_abline(intercept=0, slope=1, color="red")+
                facet_wrap(~species, labeller = labeller(species = setNames(tab_sp$species, tab_sp$id)))


            
    #         # estimate the mean TL by species
    #         sp_names <- tab_sp[which(tab_sp$id %in% unique(vals_df$species)), 2]
    #         n_sp <- length(sp_names)
    #         TL_all <-  matrix(-99, nrow=n_sp, ncol=40)
    #         dimnames(TL_all) <- list(sp_names, as.character(0:39))
    #         TL_adult <- matrix(-99, nrow=n_sp, ncol=40)
    #         dimnames(TL_adult) <- list(sp_names, as.character(0:39))
    #         TL_juv <- matrix(-99, nrow=n_sp, ncol=40)
    #         dimnames(TL_juv) <- list(sp_names, as.character(0:39))
    #         for (sp in vals_df$species){
    #             sp_nm <- tab_sp[which(tab_sp$id == sp), 2]
    #             df_sp <- vals_df[vals_df$species==sp,]
    #             df_sp <- df_sp[df_sp$age>1.3,] # remove larvae stage
    #             df_sp_juv <- df_sp[df_sp$mature==0,] # only immature individuals
    #             df_sp_adu <- df_sp[df_sp$mature==1,] # only mature individuals
    #             TL_all[which(rownames(TL_all) == sp_nm), which(colnames(TL_all) == rep_i)] <- sum(df_sp$TL * df_sp$biomass) / sum(df_sp$biomass)
    #             TL_adult[which(rownames(TL_adult) == sp_nm), which(colnames(TL_adult) == rep_i)] <- sum(df_sp_adu$TL * df_sp_adu$biomass) / sum(df_sp_adu$biomass)
    #             TL_juv[which(rownames(TL_juv) == sp_nm), which(colnames(TL_juv) == rep_i)] <- sum(df_sp_juv$TL * df_sp_juv$biomass) / sum(df_sp_juv$biomass)
            }
        }
    # }
    # if (!exists("TL_all_df")){
    #     TL_all_df <- data.frame(name=rowMeans(TL_all))
    #     TL_juv_df <-  data.frame(name=rowMeans(TL_juv))
    #     TL_adult_df <-  data.frame(name=rowMeans(TL_adult))
    #     colnames(TL_all_df) = name
    #     colnames(TL_juv_df) = name
    #     colnames(TL_adult_df) = name
    # } else {
    #     TL_all_df[, name] = rowMeans(TL_all)
    #     TL_juv_df[, name] = rowMeans(TL_juv)
    #     TL_adult_df[, name] = rowMeans(TL_adult)
    # }

}

write.csv(TL_all_df, file=paste0("osmose-ns-v4.4.1_clim/output/", "TL_all_species.csv"))
write.csv(TL_juv_df, file=paste0("osmose-ns-v4.4.1_clim/output/", "TL_juv_species.csv"))
write.csv(TL_adult_df, file=paste0("osmose-ns-v4.4.1_clim/output/", "TL_adult_species.csv"))