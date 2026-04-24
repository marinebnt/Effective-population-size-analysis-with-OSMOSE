library(dplyr)

# ---- SETTINGS ----
main_dir <- "osmose-ns-v4.4.1_clim/output/test"    # folder containing all the simulation folders
prefix   <- "ns_mortalityRate-"                   # csv filename prefix to detect
rows_to_add <- data.frame(                # the rows you want to append
  col1 = c("Predation (Mpred), Starvation (Mstarv), Additional mortality (Madd), Fishing (F) & Out-of-domain (Zout) mortality rates per time step of saving, except for Madd Eggs that is expressed in osmose time step. Z is the total mortality for migratory fish outside the simulation grid. To get annual mortality rates, sum the mortality rates within one year.", "Time", ""),																								
  col2 = c("", "Mpred", "Eggs"),
  col3 = c("", "Mpred", "Juvenil"),
  col4 = c("", "Mpred", "Adult"),
  col5 = c("", "Mstarv", "Eggs"),
  col6 = c("", "Mstarv", "Juvenil"),
  col7 = c("", "Mstarv", "Adult"),
  col8 = c("", "Madd", "Eggs"),
  col9 = c("", "Madd", "Juvenil"),
  col10 = c("", "Madd", "Adult"),
  col11 = c("", "F", "Eggs"),
  col12 = c("", "F", "Juvenil"),
  col13 = c("", "F", "Adult"),
  col14 = c("", "Zout", "Eggs"),
  col15 = c("", "Zout", "Juvenil"),
  col16 = c("", "Zout", "Adult"),
  col17 = c("", "Mfor", "Eggs"),
  col18 = c("", "Mfor", "Juvenil"),
  col19 = c("", "Mfor", "Adult"),
  col20 = c("", "Mdis", "Eggs"),
  col21 = c("", "Mdis", "Juvenil"),
  col22 = c("", "Mdis", "Adult"),
  col23 = c("", "Mage", "Eggs"),
  col24 = c("", "Mage", "Juvenil"),
  col25 = c("", "Mage", "Adult"),
  col25 = c("", "", "")
)
rows_to_add <- as.data.frame(rows_to_add)
rows_to_add$x <- 1:3
# ---- SCRIPT ----

# Get folders whose names contain ".R"
sim_folders <- list.dirs(main_dir, recursive = FALSE)
sim_folders <- sim_folders[grepl("\\.R", basename(sim_folders))]

for (folder in sim_folders) {
  
  # Find CSVs starting with the prefix
  csv_files <- list.files(folder, pattern = paste0("^", prefix), 
                          recursive = TRUE, full.names = TRUE)
  
  for (csv in csv_files) {
    
    # Read CSV
      # read.csv(csv, header = FALSE, stringsAsFactors = FALSE)
      # error = function(e) NULL

      df <- as.data.frame(read.csv(csv, header = FALSE, stringsAsFactors = FALSE))
      df$x <- 4:c(nrow(df)+3)
      # Extract the first element
      first_elem <- df[1, 1]
      
      # Check if numeric
      is_num <- suppressWarnings(!is.na(as.numeric(first_elem)))
      
      if (is_num) {
        # Append the two rows
        colnames(rows_to_add)<- names(df)
        updated <- rbind(rows_to_add, df)
        names(updated) <- NULL
        
        # Save the updated CSV (overwrite)
        write.csv(updated, csv, row.names = FALSE, )
        
        message("Updated: ", csv)
      }
    }
  }
