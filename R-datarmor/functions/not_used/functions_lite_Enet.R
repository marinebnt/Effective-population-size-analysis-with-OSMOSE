
phi = function(T, T1, T2, Topt) {
  numerator = (T - T1) * (T - T2)
  denominator = numerator - (T - Topt)^2
  phi = numerator / denominator
  # phi = ifelse(T>Tmax | T<Tmin, 0, phi)
  return(phi)
}

phim = function(T, c_m, e_m) {
  k = 8.62e-5
  Tk = T + 273.15
  return(c_m*exp(-(e_m/k)*(1/Tk)))
}

phiT = function(T, e_M, e_D, Tp) {
  # Recovers the temperature of the school cell
  k = 8.62e-5
  numerator = exp(-e_M / (k * (T + 273.15)))
  denominator = (1 + (e_M / (e_D - e_M)) * 
                   exp((e_D) / k * (1 / (Tp + 273.15) - 1 / (T + 273.15))))
  normalizer = exp(-e_M / (k * (Tp + 273.15))) /
    (1 + e_M / (e_D - e_M))
  
  output = numerator / denominator / normalizer
  
  return(output)
}

Enet = function(h, T, Imax, e_M, e_D, Tp, c_m, e_m, assimilation=0.8) {
  # T = seq(from=Tmin-2*dT, to=Tmax+dT/2, by=0.01)
  T = sort(T)
  dT = max(diff(T))
  Egross = Imax*assimilation*phiT(T, e_M, e_D, Tp)
  Emaint = phim(T, c_m, e_m)
  Enet = h*Egross - Emaint
  Topt = T[which.max(Enet)]
  Tmin = T[which.min(ifelse(T<=Topt, Enet, Inf)^2)]
  Tmax = T[which.min(ifelse(T<Topt, Inf, Enet)^2)]
  inNiche = any(Enet>= -dT)
  out = list(T=T, Enet=Enet, Egross=Egross, Emaint=Emaint, Tmin=Tmin, Topt=Topt, Tmax=Tmax, inNiche=inNiche)
  class(out) = "EnergyBudget"
  return(out)
}

Enetlite = function(h, T, eta_crit, Imax, Tmin, Topt, Tmax, assimilation=0.8) {
  # T = seq(from=Tmin-2*dT, to=Tmax+dT/2, by=0.01)
  T = sort(T)
  dT = max(diff(T))
  
  r = Tmax - Tmin
  hmin = (eta_crit/(1-eta_crit))*(Tmin-Topt)^2
  hmax = (eta_crit/(1-eta_crit))*(Tmax-Topt)^2
  H = (hmax-hmin)/r
  
  k1 = (-(r+H) + sqrt((r+H)^2 + 4*hmin))/2
  k2 = k1 + H
  T1 = Tmin - k1
  T2 = Tmax + k2
  
  Enet = assimilation*Imax*(h*phi(T, T1, T2, Topt) - eta_crit)
  
  Topt = T[which.max(Enet)]
  Tmin = T[which.min(ifelse(T<=Topt, Enet, Inf)^2)]
  Tmax = T[which.min(ifelse(T<Topt, Inf, Enet)^2)]
  inNiche = any(Enet>= -dT)
  out = list(T=T, Enet=Enet, Egross=NA, Emaint=NA, Tmin=Tmin, Topt=Topt, Tmax=Tmax, inNiche=inNiche)
  class(out) = "EnergyBudget"
  return(out)
}

plot.EnergyBudget = function(x, add=TRUE, full=FALSE, ylim=NULL, las=1, ...) {
  if(!isTRUE(add)) {
    if(isTRUE(full)) {
      ylim = range(x$Enet, x$Egross, x$Emaint, na.rm=TRUE)
    }
    plot(x$T, x$Enet, type="l", ylim=ylim, las=las, xlab="Temperature", ylab="Net Energy", ...)
    abline(h=0, lty=3, col="red")
  } else {
    lines(x$T, x$Enet, ...)
  }
  
  if(isTRUE(full)) {
    lines(x$T, x$Egross, col="blue", ...)
    lines(x$T, x$Emaint, col="red", ...)
    abline(v=c(x$Tmin, x$Topt, x$Tmax), lty=3)
  }
  
  return(invisible())
}


.fit = function(par) {
  eta_crit = par[1]
  Tmin = par[2]
  Topt = par[3]
  Tmax = par[4]
  Imax = Imax_i
  sim = Enetlite(h=1, T0, eta_crit=eta_crit, Imax=Imax, Tmin=Tmin, Topt=Topt, Tmax=Tmax)
  out = sum((sim$Enet-xh1.00$Enet)^2)
  return(out)
}


.fit2 = function(par) {
  eta_crit = par[1]
  Tmin = par[2]
  Topt = par[3]
  Tmax = par[4]
  Imax = Imax_i
  xout = list()
  for(ih in seq_along(Hs)) {
    xout[[ih]] = Enetlite(h=Hs[ih], T0, eta_crit=eta_crit, Imax, Tmin, Topt, Tmax)
  }
  sim = do.call(rbind, sapply(xout, FUN="[", i="Enet"))
  out = sum((sim - enet1)^2)
  return(out)
}
