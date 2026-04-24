
function_INIT1 <- function(MULTIP, REPO, FORCE = FALSE){
  
  for (m in MULTIP){
    simulation = paste0(REPO, "/output_", multip_name, m)
    outputDir  = file.path(configDir, "output", simulation)
    if(!dir.exists(outputDir) || FORCE){
      cat(m, "\n")
      parameters = paste0("-Psimulation.nschool.multiplier=", m)
      options = paste0("-Xmx3g -Xms1g")
      run_osmose(input = configFile, output = outputDir, osmose = jarFile, version=version, options=options, parameters = parameters, force=T)
    }
    else{next}
  }
}
