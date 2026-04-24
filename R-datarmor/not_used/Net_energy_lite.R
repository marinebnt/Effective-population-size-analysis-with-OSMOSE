library(fields)
library(calibrar)
library(osmose)

source("osmose-ns-v4.4.1_clim/R/functions/functions_lite_Enet.R")
# get_par(conf, "species.name")
# # Atlantic Herring example ------------------------------------------------
# ##############
# # Atlantic Herring
# Imax = 14.01618
# e_crit = 0.57 # predation.efficiency.critical
# # maintainance 
# c_m = 1.49e+11
# e_m = 0.5956822
# # mobilisation
# e_M = 1.929569
# e_D = 1.92967 # greater than e_M
# Tp = 43.09747
# 
# beta = 0.75
# assimilation = 0.8
# 
# Tmin = -1
# Topt = 12.56059132
# Tmax = 26.21059132
# 
# nb=0
# q = list()
# w = 100
# T= seq(from=-8, to=28, by=0.005)
# T0= seq(from=-5, to=25, by=0.05)
# hs = seq(from=0, to=1, by=0.005)
# Hs = seq(from=0, to=1, by=0.05)
# qs = seq(from=0, to=1, by=0.005)
# for (k in 1){
#   nb = nb+1
#   Topt_i <- as.numeric(Topt[nb])
#   Tp_i <- as.numeric(Tp[nb])
#   c_m_i = as.numeric(c_m[nb])
#   e_m_i = as.numeric(e_m[nb])
#   e_M_i = as.numeric(e_M[nb])
#   e_D_i = as.numeric(e_D[nb])
#   Imax_i = as.numeric(Imax[nb])
#   q = phim(Topt_i, c_m_i, e_m_i)/Imax_i/assimilation # predation.efficiency.critical.sp0?
#   out = list()
#   for(ih in seq_along(hs)) {
#     out[[ih]] = Enet(h=hs[ih], T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
#   }
#   enet = do.call(rbind, sapply(out, FUN="[", i="Enet"))
#   xTopt = unlist(sapply(out, FUN="[", i="Topt"))
#   xTmin = unlist(sapply(out, FUN="[", i="Tmin"))
#   xTmax = unlist(sapply(out, FUN="[", i="Tmax"))
#   iN = which(unlist(sapply(out, FUN="[", i="inNiche")))
#   
#   h1.00 = Enet(h=1, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
#   h0.80 = Enet(h=0.80, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
#   h0.57 = Enet(h=0.57, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
#   h0.38 = Enet(h=0.38, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
#   
#   plot(h1.00, add=FALSE, full=TRUE, lwd=3)
#   text(x=-2, y=c(15, 14, 13)+0.5, 
#        labels=c("NET ENERGY", "GROSS ENERGY", "MAINTENANCE ENERGY"),
#        col=c("black", "blue", "red"), adj=0, font=2)
#   plot(h1.00, add=FALSE, lwd=2)
#   text(26, 5, expression(xi == 1.0), col="black", font=2)
#   plot(h0.80, col="blue", lwd=2)
#   text(26, 4, expression(xi == 0.8), col="blue", font=2)
#   plot(h0.57, col="green", lwd=2)
#   text(26, 3, expression(xi == 0.57), col="green", font=2)
#   plot(h0.38, col="orange", lwd=2)
#   text(26, 2, expression(xi == 0.38), col="orange", font=2)
#   
#   
#   H1.00 = Enetlite(h=1, T, q, Imax_i, min(xTmin[iN]), max(xTopt[iN]), max(xTmax[iN]), assimilation=0.8)
#   plot(H1.00, add=FALSE, lwd=2, ylim=c(-2, 7.5), las=1)
#   abline(h=max(H1.00$Enet), lty=3, col="red")
#   text(5, 0.35 + max(H1.00$Enet), expression(alpha*I[max](1-xi[crit])),font=2)
# }
# #############

# All north sea species
conf <- read_osmose("osmose-ns-v4.4.1_clim/New_param_Eco_Bioen_neutral_notFullEnet.R")
Imax <- get_par(conf, "predation.ingestion")
e_crit <- get_par(conf, "predation.efficiency.critical")
c_m <- get_par(conf, ".c_m.")
e_m <- get_par(conf, "e.maint")
e_m <- get_par(conf, "e.maint")
e_M <- get_par(conf, "e.Mob")
e_D <- get_par(conf, "species.bioen.mobilized.e.D.")
c_M <- get_par(conf, "c_M")
c_M <- get_par(conf, "c_M")
Tp <- get_par(conf, "mobilized.tp")
Tmin <- get_par(conf, "tmin")
Tmax <- get_par(conf, "tmax")
Topt <- get_par(conf, "topt")
sp <- get_par(conf, "species.name")
beta = 0.75
assimilation = 0.8

nb=0
q = list()
w = 100
T= seq(from=-8, to=28, by=0.005)
T0= seq(from=-5, to=25, by=0.05)
hs = seq(from=0, to=1, by=0.005)
Hs = seq(from=0, to=1, by=0.05)
qs = seq(from=0, to=1, by=0.005)
for (k in names(e_crit)){
  nb = nb+1
  Topt_i <- as.numeric(Topt[nb])
  Tp_i <- as.numeric(Tp[nb])
  c_m_i = as.numeric(c_m[nb])
  e_m_i = as.numeric(e_m[nb])
  e_M_i = as.numeric(e_M[nb])
  e_D_i = as.numeric(e_D[nb])
  Imax_i = as.numeric(Imax[nb])
  q[[k]] = phim(Topt_i, c_m_i, e_m_i)/Imax_i/assimilation # predation.efficiency.critical.sp0?
  out = list()
  for(ih in seq_along(hs)) {
    out[[ih]] = Enet(h=hs[ih], T, Imax = Imax_i, e_M = e_M_i, 
                      e_D = e_D_i, Tp = Tp_i, c_m = c_m_i, e_m = e_m_i)
  }
  enet = do.call(rbind, sapply(out, FUN="[", i="Enet"))
  xTopt = unlist(sapply(out, FUN="[", i="Topt"))
  xTmin = unlist(sapply(out, FUN="[", i="Tmin"))
  xTmax = unlist(sapply(out, FUN="[", i="Tmax"))
  iN = which(unlist(sapply(out, FUN="[", i="inNiche")))
  xh1.00 = Enet(h=1, T0, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
  
  out1 = list()
  par0 = c(0.3, min(xTmin[iN]), max(xTopt[iN]), max(xTmax[iN]))
  opt = optim2(par0, fn=.fit, lower=c(0, -10, 0, 0), upper=c(1,5,20,30), method="AHR-ES")
  heta_crit = opt$par[1]
  hTmin = opt$par[2]
  hTopt = opt$par[3]
  hTmax = opt$par[4]
  for(ih in seq_along(Hs)) {
    out1[[ih]] = Enet(h=Hs[ih], T0, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
  }
  
  enet1 = do.call(rbind, sapply(out1, FUN="[", i="Enet"))
  out2 = list()
  opt2 = optim2(par0, fn=.fit2, lower=c(0, -10, 0, 0), upper=c(1,5,20,30), method="AHR-ES")
  geta_crit = opt2$par[1]
  gTmin = opt2$par[2]
  gTopt = opt2$par[3]
  gTmax = opt2$par[4]
  for(ih in seq_along(hs)) {
    out2[[ih]] = Enetlite(h=hs[ih], T, eta_crit=geta_crit, 
                          Imax_i, Tmin=gTmin, Topt=gTopt, Tmax=gTmax)
  }
  
  assign(paste0("out_", k), value = list(out, out1, out2))
  
}
e_crit <- (as.matrix(q))
colnames(e_crit) <- 
write_osmose(e_crit, "osmose-ns-v4.4.1_clim/pred_eff_Enetlite_estimate.R", sep = "=")
# eta_crit_new = q

nb=0
for (k in names(e_crit)){ # loop for the plots
  nb=nb+1
  Topt_i <- as.numeric(Topt[nb])
  Tmin_i <- as.numeric(Tmin[nb])
  Tmax_i <- as.numeric(Tmax[nb])
  Tp_i <- as.numeric(Tp[nb])
  c_m_i = as.numeric(c_m[nb])
  e_m_i = as.numeric(e_m[nb])
  e_M_i = as.numeric(e_M[nb])
  e_D_i = as.numeric(e_D[nb])
  Imax_i = as.numeric(Imax[nb])
  q_i = q[nb]
  sp_i <- sp[nb]
  output = get(paste0("out_", k))
  
  enet = do.call(rbind, sapply(output[[1]], FUN="[", i="Enet"))
  xTopt = unlist(sapply(output[[1]], FUN="[", i="Topt"))
  xTmin = unlist(sapply(output[[1]], FUN="[", i="Tmin"))
  xTmax = unlist(sapply(output[[1]], FUN="[", i="Tmax"))
  iN = which(unlist(sapply(output[[1]], FUN="[", i="inNiche")))
  
  enet2 = do.call(rbind, sapply(output[[3]], FUN="[", i="Enet"))
  xTopt2 = unlist(sapply(output[[2]], FUN="[", i="Topt"))
  xTmin2 = unlist(sapply(output[[2]], FUN="[", i="Tmin"))
  xTmax2 = unlist(sapply(output[[2]], FUN="[", i="Tmax"))
  iN2 = which(unlist(sapply(output[[2]], FUN="[", i="inNiche")))
  
  H1.00 = Enetlite(h=1, T, q[[k]], Imax_i, min(xTmin[iN]), max(xTopt[iN]), max(xTmax[iN]), assimilation=0.8)
  h1.00 = Enet(h=1, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
  h0.80 = Enet(h=0.80, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
  h0.57 = Enet(h=0.57, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
  h0.38 = Enet(h=0.38, T, Imax_i, e_M_i, e_D_i, Tp_i, c_m_i, e_m_i)
  
  dir.create(paste0("osmose-ns-v4.4.1_clim/R/plots/", sp_i), recursive = T)
  
  # Figure 1 ----------------------------------------------------------------
  png(paste0("osmose-ns-v4.4.1_clim/R/plots/", sp_i, "/E_net.png"))
  plot(h1.00, add=FALSE, full=TRUE, lwd=3, main = k)
  text(x=-2, y=c(15, 14, 13)+0.5, 
       labels=c("NET ENERGY", "GROSS ENERGY", "MAINTENANCE ENERGY"),
       col=c("black", "blue", "red"), adj=0, font=2)
  dev.off()
  
  # Figure 2 ----------------------------------------------------------------
  png(paste0("osmose-ns-v4.4.1_clim/R/plots/", sp_i, "/E_net_compare_eff_crit.png"))
  plot(h1.00, add=FALSE, lwd=2, main = k)
  text(26, 5, expression(xi == 1.0), col="black", font=2)
  plot(h0.80, col="blue", lwd=2)
  text(26, 4, expression(xi == 0.8), col="blue", font=2)
  plot(h0.57, col="green", lwd=2)
  text(26, 3, expression(xi == 0.57), col="green", font=2)
  plot(h0.38, col="orange", lwd=2)
  text(26, 2, expression(xi == 0.38), col="orange", font=2)
  dev.off()
  
  # Figure 3 ----------------------------------------------------------------
  png(paste0("osmose-ns-v4.4.1_clim/R/plots/", sp_i, "/E_netlite.png"))
  plot(H1.00, add=FALSE, lwd=2, ylim=c(-2, 7.5), las=1)
  abline(h=max(H1.00$Enet), lty=3, col="red")
  text(5, 0.35 + max(H1.00$Enet), expression(alpha*I[max](1-xi[crit])),font=2)
  dev.off()
  
  # Figure 4 ----------------------------------------------------------------
  png(paste0("osmose-ns-v4.4.1_clim/R/plots/", sp_i, "/E_netlite_and_regular.png"))
  plot(H1.00, add=FALSE, lwd=2, ylim=c(-2, 6.5), las=1)
  plot(h1.00, lwd=2, col="blue")
  text(x=25, y=c(6,5.5), labels=c("FULL", "LITE"), col=c("black", "blue"), adj=0, font=2)
  dev.off()
  
  # # Figure 5 ---------------------------------------------------------------
  # zlim = range(c(enet, enet2))/Imax_i
  # zlim = range(pretty(zlim))*1.05
  # n = 512
  # col = colorful::divergencePalette(n=n, zlim=zlim, 
  #                                   col=c("firebrick3", "dodgerblue3"), symmetric = FALSE)
  # # Open PNG graphics device
  # png(paste0("osmose-ns-v4.4.1_clim/R/plots/", sp_i, "/combined_heatmaps.png"), 
  #     width = 1200, height = 800, res = 150)
  # # Multi-panel layout setup
  # par(mar = c(4, 4, 1, 1), oma = c(1, 1, 1, 1))
  # layout(mat = matrix(c(1, 3, 2, 3), nrow = 2), heights = c(5, 1))
  # # --- Panel 1 ---
  # image(hs, T, enet / Imax_i, col = col, las = 1, zlim = zlim, 
  #       xlab = "Predation Efficiency", ylab = "Temperature")
  # lines(hs[iN], xTopt[iN], lty = 3)
  # lines(c(rev(hs[iN]), hs[iN]), c(rev(xTmin[iN]), xTmax[iN]), col = "grey80")
  # # --- Panel 2 ---
  # image(hs, T, enet2 / Imax_i, col = col, las = 1, zlim = zlim, 
  #       xlab = "Predation Efficiency", ylab = "Temperature")
  # lines(hs[iN2], xTopt2[iN2], lty = 3)
  # lines(c(rev(hs[iN2]), hs[iN2]), c(rev(xTmin2[iN2]), xTmax2[iN2]), col = "grey80")
  # # --- Color Legend (Panel 3) ---
  # par(mar = c(2, 10, 3, 10))
  # zcol = matrix(seq(from = zlim[1], to = zlim[2], length.out = n), ncol = 1)
  # image(as.numeric(zcol), 1, zcol, zlim = zlim, col = col, axes = FALSE,
  #       ylab = "", xlab = "")
  # mtext("Relative Net Energy", 3, cex = 0.8)
  # axis(1)
  # box()
  # # Close graphics device
  # dev.off()
  
}

  
  # plot
  
    
  
  
  # TPC plot ----------------------------------------------------------------
  
  x = h1.00
  par(mar=c(2,3,1,1))
  plot(tail(x$T, -500), tail(x$Enet, -500), axes=FALSE, xlab="", ylab="", type="l",
       ylim=c(-0.5, 1.1*max(x$Enet)), xlim=c(-8, 30), col=4, lwd=3)
  arrows(x0=par("usr")[1]+2, y0=-0.5, y1=par("usr")[4], lty=1, lwd=3)
  segments(x0=x$Topt, y0=0, y1=max(x$Enet), lty=3)
  axis(2, at=0, tick = FALSE, las=1, line=-0.5, font=2)
  ptr = c(0.7*x$Tmin+0.3*x$Topt, 0.5*x$Tmax+0.5*x$Topt)
  text(x=mean(ptr), y=0.7, labels = "PTR", font=2, cex=1.5)
  polygon(x=c(ptr, rev(ptr)), y=c(0,0,0.3,0.3), angle=45, col=3)
  polygon(x=c(ptr, rev(ptr)), y=c(0,0,0.3,0.3), angle=45, col=1, density=10, lwd=2)
  arrows(x0=par("usr")[1], x1=par("usr")[2], y0=0, lty=1, lwd=3)
  abline(h=0, lty=1, lwd=3)
  ref = c(x$Tmin, x$Topt, 0.7*x$Tmax+0.3*x$Topt, x$Tmax)
  axis(1, at=ref, labels=c("Tmin", "Topt", "Tpk", "Tmax"), pos=0, font=2, lwd.ticks=3)
  axis(1, at=ref[4], labels="Tmax", pos=0, font=2)
  mtext("TRAIT PERFORMANCE", side=2, line=-1, font=2)
  mtext("TEMPERATURE", side=1, font=2)
  
  plot(qs, log10(output[[3]]))
  q = qs[which.min(out2)]
  
  H1.00 = Enetlite(h=1.00, T, eta_crit=xq, Imax_i, xTmin, xTopt, xTmax)
  H0.80 = Enetlite(h=0.80, T, eta_crit=xq, Imax, xTmin, xTopt, xTmax)
  H0.57 = Enetlite(h=0.57, T, eta_crit=xq, Imax, xTmin, xTopt, xTmax)
  H0.38 = Enetlite(h=0.38, T, eta_crit=xq, Imax, xTmin, xTopt, xTmax)
  
  m1.00 = Enetlite(h=1.00, T, q=xq, Imax, xTmin-3, max(Topt), xTmax+3)
  x1.00 = m1.00
  x1.00$Enet = ifelse(x1.00$T < x1.00$Tmin, -(abs(x1.00$Enet)+1)^0.75+1, x1.00$Enet)
  
  plot(x1.00, add=FALSE, lwd=2)
  
  plot(h1.00, col="red")
  plot(m1.00)
  
  plot(H1.00, add=FALSE, lwd=2)
  plot(H0.80, col="blue", lwd=2)
  plot(H0.57, col="green", lwd=2)
  plot(H0.38, col="orange", lwd=2)
  
  plot(h1.00)
  plot(h0.80, col="blue")
  plot(h0.57, col="green")
  plot(h0.38, col="orange")
  
}
















# Tmin, Topt and Tmax as function of predation efficiency (I/Imax). Classical OSMOSE assumes critical predation efficiency as 0.57, where ingestion only covers maintenance. 

tpc = phi(T, Tmin_i, Tmax_i, Topt_i)
plot(T, tpc, type="l")

