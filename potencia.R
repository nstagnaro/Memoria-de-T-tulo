library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

fishercalc<-function(x, cdf){
  if(min(cdf)<0 || max(cdf)>1){ stop("invalid cumulative distribution function")}
  cdf<-sort(unique(cdf), decreasing = FALSE)
  if(max(cdf)<1){cdf<-c(cdf,1)}
  verify<-prod(x %in% cdf)
  if(verify==0){ stop("given p-values do not come from specified cumulative distribution function")}
  z<-rep(0,length(x))
  vals<-sort(unique(x), decreasing = FALSE) 
  m<-length(cdf)
  z[x==cdf[1]]<-2*(1-log(cdf[1]))
  for(i in 2:m){
    if(cdf[i]-cdf[i-1]>10^(-8)){z[x==cdf[i]]<-2*(1-(cdf[i]*log(cdf[i])- cdf[i-1]*log(cdf[i-1]))/(cdf[i]-cdf[i-1])) }
    else{z[x==cdf[i]]<--2*log((cdf[i]+cdf[i-1])/2)}
  }
  return(z) # Z grande de Fisher
}


edgingtoncalc<-function(x, cdf){
  if(min(cdf)<0 || max(cdf)>1){ stop("invalid cumulative distribution function")}
  cdf<-sort(unique(cdf), decreasing = FALSE)
  if(max(cdf)<1){cdf<-c(cdf,1)}
  verify<-prod(x %in% cdf)
  if(verify==0){ stop("given p-values do not come from specified cumulative distribution function")}
  z<-rep(0,length(x))
  vals<-sort(unique(x), decreasing = FALSE)
  m<-length(cdf)
  z[x==cdf[1]]<-cdf[1]/2
  for(i in 2:m){
    z[x==cdf[i]]<- (cdf[i]+cdf[i-1])/2
  }
  return(z)
}


pearsoncalc<-function(x, cdf){
  if(min(cdf)<0 || max(cdf)>1){ stop("invalid cumulative distribution function")}
  cdf<-sort(unique(cdf), decreasing = FALSE)
  if(max(cdf)<1){cdf<-c(cdf,1)}
  verify<-prod(x %in% cdf)
  if(verify==0){ stop("given p-values do not come from specified cumulative distribution function")}
  z<-rep(0,length(x))
  x[x==0]<-min(x[x>0]) 
  vals<-sort(unique(x), decreasing = FALSE)
  m<-length(cdf)
  z[x==cdf[1]]<- 2+2*log(1-cdf[1])*(1-cdf[1])/cdf[1]
  for(i in 2:(m-1)){
    if(cdf[i]-cdf[i-1]>10^(-8)){z[x==cdf[i]]<-2*(1+((1-cdf[i])*log(1-cdf[i])- (1-cdf[i-1])*log(1-cdf[i-1]))/(cdf[i]-cdf[i-1])) }
    else{z[x==cdf[i]]<--2*log(1-(cdf[i]+cdf[i-1])/2)}
  }
  z[x==cdf[m]]<-2*(1-( (1-cdf[m-1])*log(1-cdf[m-1]))/(1-cdf[m-1]))
  return(z)
}

#### Stouffer's method:  standard normal quantile of P
stouffercalc<-function(x, cdf){
  if(min(cdf)<0 || max(cdf)>1){ stop("invalid cumulative distribution function")}
  cdf<-sort(unique(cdf), decreasing = FALSE)
  if(max(cdf)<1){cdf<-c(cdf,1)}
  verify<-prod(x %in% cdf)
  if(verify==0){ stop("given p-values do not come from specified cumulative distribution function")}
  z<-rep(0,length(x))
  x[x==0]<-min(x[x>0])
  cdf<-cdf[cdf>0]
  vals<-sort(unique(x), decreasing = FALSE)
  m<-length(cdf)
  z[x==cdf[1]]<- (-exp(-qnorm(cdf[1])^2/2))/(cdf[1])
  for(i in 2:m){
    if(cdf[i]-cdf[i-1]>10^(-8)){z[x==cdf[i]]<-(-exp(-qnorm(cdf[i])^2/2)+exp(-qnorm(cdf[i-1])^2/2))/(cdf[i]-cdf[i-1]) }
    else{z[x==cdf[i]]<-pnorm((qnorm(cdf[i])+qnorm(cdf[i-1]))/2)}
  }
  return(z/sqrt(2*pi))
}

### Acá falta George que no se pone porque es (Pearson - Fisher)/2


##### varianzas
varfisher<-function(cdf){
  cdf<-sort(unique(cdf), decreasing = FALSE) 
  n<-length(cdf)-1
  c<-cdf*log(cdf)
  aux<-c(0,cdf[1:n]*log(cdf[1:n]))
  prob<-cdf-c(0,cdf[1:n])
  c<-c[prob>0]
  aux<-aux[prob>0]
  prob<-prob[prob>0]
  return(4*sum((c-aux)^2/prob))
}


varpearson<-function(cdf){
  cdf<-sort(unique(cdf), decreasing = FALSE) # Ensure the cdf is unique and sorted in increasing order
  cdf[cdf>1-10^(-10)]<-1-10^(-10)
  n<-length(cdf)-1
  c<-(1-cdf)*log(1-cdf)
  aux<-c(0,(1-cdf[1:n])*log(1-cdf[1:n]))
  prob<-cdf-c(0,cdf[1:n])
  c<-c[prob>0]
  aux<-aux[prob>0]
  prob<-prob[prob>0]
  return(4*sum((c-aux)^2/prob))
}


varedgington<-function(cdf){
  cdf<-sort(unique(cdf), decreasing = FALSE) # Ensure the cdf is unique and sorted in increasing order
  n<-length(cdf)-1
  aux<-c(0,cdf[1:n])
  prob<-cdf-aux
  return(sum(cdf*aux*prob/4))
}


varstouffer<-function(cdf){
  cdf<-sort(unique(cdf), decreasing = FALSE) # Ensure the cdf is unique and sorted in increasing order
  cdf<-cdf[cdf>10^(-80)]
  n<-length(cdf)-1
  aux<-c(0,cdf[1:n])
  prob<-cdf-aux
  return(sum(stouffercalc(cdf,cdf)^2*prob ))
}


vargeorge<-function(cdf){
  cdf<-sort(unique(cdf), decreasing = FALSE) # Ensure the cdf is unique and sorted in increasing order
  cdf[cdf>1-10^(-10)]<-1-10^(-10)
  cdf[cdf<10^(-10)]<-10^(-10)
  c<-cdf*log(cdf)
  c<-c(0,c)
  c[length(c)]<-0
  d<-(1-cdf)*log(1-cdf)
  d[length(d)]<-0
  d<-c(0,d)
  n<-length(cdf)-1
  aux<-c(0,cdf[1:n])
  prob<-cdf-aux
  aux2<-c[1:length(cdf)]+d[1:length(cdf)]
  aux3<-c[2:length(d)]+d[2:length(d)]
  aux2<-aux2[prob>0]
  aux3<-aux3[prob>0]
  prob<-prob[prob>0]
  return(sum((aux2-aux3)^2/prob)) 
}

# P-valores discretos por cola bajo H0 
pvals_from_pmf <- function(pmf, side=c("left","right","two.sided"),
                           two_sided=c("doubling","prob","midp")) {
  side <- match.arg(side)
  two_sided <- match.arg(two_sided)
  
  if (side=="left")  return(cumsum(pmf))
  if (side=="right") return(rev(cumsum(rev(pmf))))
  
  # Bilateral 
  # Necesitamos los p de una sola realización X sobre el SOPORTE: 
  # pmf está en el orden del soporte (x = x_min, ..., x_max)
  p_left  <- cumsum(pmf)
  p_right <- rev(cumsum(rev(pmf)))
  
  if (two_sided == "doubling") {
    # p_two = min(1, 2 * min{pL, pR})
    return(pmin(1, 2*pmin(p_left, p_right)))
  }
  
  if (two_sided == "midp") {
    # mid-p: pL + pR - pmf, y doblamos la parte estricta
    p_dbl <- pmin(1, 2*pmin(p_left - 0.5*pmf, p_right - 0.5*pmf))
    return(pmax(0, p_dbl))
  }
  
  if (two_sided == "prob") {
    o <- order(pmf)               # de más pequeño a más grande
    rk <- match(seq_along(pmf), o)
    cs <- cumsum(pmf[o])
    return(cs[rk])
  }
}

# Plots de las distribuciones de los p-vals que se usan en los metodos (pmf y cdf)

# Aca se calculan la pmf y cdf

build_pval_distributions <- function(support, pmf, digits = 12) {
  stopifnot(length(support) == length(pmf))
  stopifnot(abs(sum(pmf) - 1) < 1e-8)
  
  pv_left  <- pvals_from_pmf(pmf, "left")
  pv_right <- pvals_from_pmf(pmf, "right")
  pv_two   <- pvals_from_pmf(pmf, "two.sided")
  
  agg <- function(pv, pmf, tail_label) {
    tibble(p = round(pv, digits), prob = pmf) |>
      group_by(p) |>
      summarize(prob = sum(prob), .groups = "drop") |>
      arrange(p) |>
      mutate(tail = tail_label,
             cdf  = cumsum(prob))
  }
  
  out <- bind_rows(
    agg(pv_left,  pmf, "left"),
    agg(pv_right, pmf, "right"),
    agg(pv_two,   pmf, "two.sided")
  )
  
  out |> mutate(
    tail = factor(
      tail,
      levels = c("left","right","two.sided"),
      labels = c("Cola izquierda","Cola derecha","Bilateral")
    )
  )
}

.make_subtitle <- function(vars) {
  if (is.null(vars)) return(NULL)
  if (!is.list(vars)) vars <- as.list(vars)
  if (is.null(names(vars)) || any(names(vars) == "")) {
    stop("`subtitle_vars` debe ser una lista nombrada, p.ej. list(ncase=1000, ncontr=4000, nmuta=100)")
  }
  fmt_val <- function(x) if (is.numeric(x)) format(x, trim=TRUE, scientific=FALSE, digits=6) else as.character(x)
  vals <- vapply(vars, fmt_val, character(1))
  paste(sprintf("%s = %s", names(vars), vals), collapse = ", ")
}

# Plots en azul + subtítulo + orden Cola izq -> Cola der -> Bilateral
plot_pval_pmf <- function(dist_df, title = "Distribución de los p-valores (PMF)",
                          subtitle_vars = NULL, show_points = TRUE) {
  subtitle <- .make_subtitle(subtitle_vars)
  
  ggplot(dist_df, aes(x = p, y = prob, group = tail)) +
    geom_line(linewidth = 0.9, color = "blue") +
    { if (show_points) geom_point(size = 1.4, color = "blue") } +
    facet_wrap(~ tail, nrow = 1) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    labs(x = "p-valor", y = "Masa de probabilidad",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12)
}

plot_pval_cdf <- function(dist_df, title = "Distribución de los p-valores (CDF)",
                          subtitle_vars = NULL, show_points = FALSE) {
  subtitle <- .make_subtitle(subtitle_vars)
  
  ggplot(dist_df, aes(x = p, y = cdf, group = tail)) +
    geom_step(linewidth = 0.9, color = "blue") +
    { if (show_points) geom_point(size = 1.2, color = "blue") } +
    facet_wrap(~ tail, nrow = 1) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    labs(x = "p-valor", y = "F(p) = P(P ≤ p)",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12)
}

# Precomputación bajo H0 (w=1) para cortar a nivel alpha
prep_nc_experiment <- function(ncase, ncontr, nmuta, alpha, n_comb){
  N <- ncase + ncontr
  support <- max(0, ncase+nmuta-N):min(ncase, nmuta)
  pmf_H0  <- dhyper(support, nmuta, N-nmuta, ncase)
  stopifnot(abs(sum(pmf_H0)-1) < 1e-8)
  
  sides <- c("left","right","two.sided")
  prep_por_cola <- lapply(setNames(sides, sides), function(lado){
    pv    <- pvals_from_pmf(pmf_H0, lado)
    cdf_p <- sort(unique(pv))
    pos_of <- match(pv, cdf_p)
    
    zF <- fishercalc(cdf_p, cdf_p)
    zP <- pearsoncalc(cdf_p, cdf_p)
    zS <- stouffercalc(cdf_p, cdf_p)
    zE <- edgingtoncalc(cdf_p, cdf_p)
    zG <- (zP - zF)/2
    
    nuF <- varfisher(sort(pv))
    nuP <- varpearson(sort(pv))
    nuS <- varstouffer(sort(pv))
    nuE <- varedgington(sort(pv))
    nuG <- vargeorge(sort(pv))
    
    i <- 1:n_comb
    cutF <- qgamma(1-alpha, shape=4*i/nuF, scale=nuF/2)   # >
    cutP <- qgamma(alpha,   shape=4*i/nuP, scale=nuP/2)   # <
    cutS <- qnorm(alpha, mean=0,   sd=sqrt(i*nuS))        # <
    cutE <- qnorm(alpha, mean=i/2, sd=sqrt(i*nuE))        # <
    cutG <- qnorm(alpha, mean=0,   sd=sqrt(i*nuG))        # <
    
    list(lado=lado, support=support, pos_of=pos_of,
         zF=zF, zP=zP, zS=zS, zE=zE, zG=zG,
         cutF=cutF, cutP=cutP, cutS=cutS, cutE=cutE, cutG=cutG)
  })
  
  list(support=support, prep_por_cola=prep_por_cola)
}

# Preparación H0 (w=1) y ratios solo para "n=n_comb"
simulate_nc_power <- function(ncase, ncontr, nmuta,
                                     alpha,
                                     n_comb,
                                     iter,
                                     logw_left,         # vector numérico (p.ej. seq(-0.12, 0.08, by=0.02))
                                     logw_right,        # vector numérico
                                     logw_two,          # vector numérico
                                     print_every = 100) {
  
  # Pre-cálculos bajo H0 para cortes y mapas z 
  base <- prep_nc_experiment(ncase, ncontr, nmuta, alpha, n_comb)
  support <- base$support
  P_left  <- base$prep_por_cola[["left"]]
  P_right <- base$prep_por_cola[["right"]]
  P_two   <- base$prep_por_cola[["two.sided"]]
  
  # Helpers para una cola dada (simulación en una malla de w)
  .sim_one_side <- function(logw_grid, P, lado_label) {
    w_grid <- exp(logw_grid)
    out <- vector("list", length(w_grid))
    
    for (iw in seq_along(w_grid)) {
      w <- w_grid[iw]
      cat(sprintf("[%s %2d/%d] log(w)=%.5f  w=%.6f ... ",
                  lado_label, iw, length(w_grid), log(w), w)); flush.console()
      
      # Contadores SOLO en i = n_comb
      acc <- list(Fisher=0L, Pearson=0L, Stouffer=0L, Edgington=0L, George=0L)
      
      for (t in 1:iter) {
        # Generar n_comb obs no centrales (truco U^w)
        X_vec <- replicate(n_comb, {
          Z <- c(runif(ncontr), runif(ncase)^w)
          sum(order(Z)[1:nmuta] > ncontr)
        })
        idx <- X_vec - support[1] + 1L
        pos <- P$pos_of[idx]
        
        # Suma final (i = n_comb) de cada método y comparación con cortes H0
        sF <- sum(P$zF[pos]); sP <- sum(P$zP[pos]); sS <- sum(P$zS[pos])
        sE <- sum(P$zE[pos]); sG <- sum(P$zG[pos])
        cF <- P$cutF[n_comb];  cP <- P$cutP[n_comb];  cS <- P$cutS[n_comb]
        cE <- P$cutE[n_comb];  cG <- P$cutG[n_comb]
        
        acc$Fisher    <- acc$Fisher    + as.integer(sF > cF)
        acc$Pearson   <- acc$Pearson   + as.integer(sP < cP)
        acc$Stouffer  <- acc$Stouffer  + as.integer(sS < cS)
        acc$Edgington <- acc$Edgington + as.integer(sE < cE)
        acc$George    <- acc$George    + as.integer(sG < cG)
        
        if (t %% print_every == 0) { cat("."); flush.console() }
      }
      cat(" ok\n")
      
      out[[iw]] <- data.frame(
        logw = log(w), w = w, lado = lado_label, alpha = alpha, n_comb = n_comb,
        Fisher    = acc$Fisher    / iter,
        Pearson   = acc$Pearson   / iter,
        Stouffer  = acc$Stouffer  / iter,
        Edgington = acc$Edgington / iter,
        George    = acc$George    / iter
      )
    }
    dplyr::bind_rows(out)
  }
  
  # Caso especial: si left y right tienen EXACTAMENTE la misma malla, simulamos una sola vez por w y evaluamos ambos lados.
  .sim_left_right_same_grid <- function(logw_grid, P_left, P_right) {
    w_grid <- exp(logw_grid)
    out <- vector("list", 2 * length(w_grid))
    k <- 0L
    
    for (iw in seq_along(w_grid)) {
      w <- w_grid[iw]
      cat(sprintf("[L/R %2d/%d] log(w)=%.5f  w=%.6f ... ", iw, length(w_grid), log(w), w)); flush.console()
      
      accL <- list(Fisher=0L, Pearson=0L, Stouffer=0L, Edgington=0L, George=0L)
      accR <- list(Fisher=0L, Pearson=0L, Stouffer=0L, Edgington=0L, George=0L)
      
      for (t in 1:iter) {
        X_vec <- replicate(n_comb, {
          Z <- c(runif(ncontr), runif(ncase)^w)
          sum(order(Z)[1:nmuta] > ncontr)
        })
        idx  <- X_vec - support[1] + 1L
        
        # LEFT
        posL <- P_left$pos_of[idx]
        sFL <- sum(P_left$zF[posL]); sPL <- sum(P_left$zP[posL]); sSL <- sum(P_left$zS[posL])
        sEL <- sum(P_left$zE[posL]); sGL <- sum(P_left$zG[posL])
        cFL <- P_left$cutF[n_comb];  cPL <- P_left$cutP[n_comb];  cSL <- P_left$cutS[n_comb]
        cEL <- P_left$cutE[n_comb];  cGL <- P_left$cutG[n_comb]
        accL$Fisher    <- accL$Fisher    + as.integer(sFL > cFL)
        accL$Pearson   <- accL$Pearson   + as.integer(sPL < cPL)
        accL$Stouffer  <- accL$Stouffer  + as.integer(sSL < cSL)
        accL$Edgington <- accL$Edgington + as.integer(sEL < cEL)
        accL$George    <- accL$George    + as.integer(sGL < cGL)
        
        # RIGHT
        posR <- P_right$pos_of[idx]
        sFR <- sum(P_right$zF[posR]); sPR <- sum(P_right$zP[posR]); sSR <- sum(P_right$zS[posR])
        sER <- sum(P_right$zE[posR]); sGR <- sum(P_right$zG[posR])
        cFR <- P_right$cutF[n_comb];  cPR <- P_right$cutP[n_comb];  cSR <- P_right$cutS[n_comb]
        cER <- P_right$cutE[n_comb];  cGR <- P_right$cutG[n_comb]
        accR$Fisher    <- accR$Fisher    + as.integer(sFR > cFR)
        accR$Pearson   <- accR$Pearson   + as.integer(sPR < cPR)
        accR$Stouffer  <- accR$Stouffer  + as.integer(sSR < cSR)
        accR$Edgington <- accR$Edgington + as.integer(sER < cER)
        accR$George    <- accR$George    + as.integer(sGR < cGR)
        
        if (t %% print_every == 0) { cat("."); flush.console() }
      }
      cat(" ok\n")
      
      k <- k + 1L
      out[[k]] <- data.frame(
        logw = log(w), w = w, lado = "Cola izquierda", alpha = alpha, n_comb = n_comb,
        Fisher=accL$Fisher/iter, Pearson=accL$Pearson/iter, Stouffer=accL$Stouffer/iter,
        Edgington=accL$Edgington/iter, George=accL$George/iter
      )
      k <- k + 1L
      out[[k]] <- data.frame(
        logw = log(w), w = w, lado = "Cola derecha", alpha = alpha, n_comb = n_comb,
        Fisher=accR$Fisher/iter, Pearson=accR$Pearson/iter, Stouffer=accR$Stouffer/iter,
        Edgington=accR$Edgington/iter, George=accR$George/iter
      )
    }
    dplyr::bind_rows(out)
  }
  
  cat(sprintf("\n== Poder (n=%d, alpha=%g, iter=%d) con mallas arbitrarias ==\n", n_comb, alpha, iter))
  # Fase L/R
  if (length(logw_left) == length(logw_right) && all(abs(logw_left - logw_right) < 1e-12)) {
    tab_lr <- .sim_left_right_same_grid(logw_left, P_left, P_right)
  } else {
    tab_l  <- .sim_one_side(logw_left,  P_left,  "Cola izquierda")
    tab_r  <- .sim_one_side(logw_right, P_right, "Cola derecha")
    tab_lr <- dplyr::bind_rows(tab_l, tab_r)
  }
  # Fase bilateral
  tab_two <- .sim_one_side(logw_two, P_two, "Bilateral")
  
  dplyr::bind_rows(tab_lr, tab_two)
}

# Plot: curva de poder vs log(w) 
plot_power_all <- function(power_table, title = NULL, subtitle_vars = NULL,
                           facet_spacing_x = 1.5,
                           bilateral_x_range = NULL) {
  
  df <- power_table |>
    dplyr::mutate(
      lado = factor(lado,
                    levels = c("Cola izquierda","Cola derecha","Bilateral"))
    )
  
  # Pasamos a formato largo de una vez
  df_long <- df |>
    tidyr::pivot_longer(
      cols = c(Fisher, Pearson, Stouffer, Edgington, George),
      names_to = "Metodo", values_to = "Power"
    )
  
  # Rango especial para Bilateral (funciona aunque sea mayor o menor)
  if (!is.null(bilateral_x_range)) {
    stopifnot(length(bilateral_x_range) == 2)
    rng <- sort(bilateral_x_range)
    
    # Separamos bilateral del resto
    df_lr  <- dplyr::filter(df_long, lado != "Bilateral")
    df_bil <- dplyr::filter(df_long, lado == "Bilateral")
    
    # 1) si el rango es más angosto -> filtramos
    df_bil <- dplyr::filter(df_bil, logw >= rng[1], logw <= rng[2])
    
    # 2) si el rango es más ancho -> extendemos con puntos "fantasma" (Power = NA)
    if (nrow(df_bil) > 0) {
      # Usamos una sola categoría de 'Metodo' para forzar el eje (vale con NA)
      met0 <- levels(factor(df_long$Metodo))[1]
      # armamos filas mínimas con las columnas necesarias para el mapeo y facet
      pad <- data.frame(
        logw   = rng,
        Power  = NA_real_,
        Metodo = met0,
        lado   = factor("Bilateral",
                        levels = c("Cola izquierda","Cola derecha","Bilateral"))
      )
      df_bil <- dplyr::bind_rows(df_bil, pad)
    } else {
      # si no hay datos (raro), creamos solo el padding para mostrar el eje
      met0 <- levels(factor(df_long$Metodo))[1]
      df_bil <- data.frame(
        logw   = rng,
        Power  = NA_real_,
        Metodo = met0,
        lado   = factor("Bilateral",
                        levels = c("Cola izquierda","Cola derecha","Bilateral"))
      )
    }
    
    # recombinamos
    df_long <- dplyr::bind_rows(df_lr, df_bil)
  }
  
  if (is.null(title)) {
    title <- sprintf("Curvas de poder (n = %d, \u03B1 = %g)",
                     unique(df_long$n_comb), unique(df_long$alpha))
  }
  subtitle <- if (is.null(subtitle_vars)) NULL else {
    if (!is.list(subtitle_vars)) subtitle_vars <- as.list(subtitle_vars)
    paste(sprintf("%s = %s", names(subtitle_vars),
                  vapply(subtitle_vars,
                         function(x) if (is.numeric(x))
                           format(x, digits=6, trim=TRUE) else x, "")),
          collapse = ", ")
  }
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = logw, y = Power, color = Metodo, shape = Metodo)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) formatC(x, format = "f", digits = 2)
    ) +
    ggplot2::scale_x_continuous(name = "log(w)") +
    ggplot2::facet_wrap(~ lado, nrow = 1, scales = "free_x") +  # eje X libre por faceta
    ggplot2::labs(title = title, subtitle = subtitle, y = "Power") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title.position = "plot",
      panel.spacing.x = grid::unit(facet_spacing_x, "lines")
    )
}


### En este caso notemos que estamos interesados en el conteo de rechazos, pero asumiendo H1 V 
# Pues lo que se busca aca es medir la potencia del test basado en los 5 metodos discretizados
# Y sabemos que la potencia es la capacidad que tiene el test (con su metodo respectivo) de eliminar la basura
# Esto es de rechazar H0 cuando este es falso
# Así, mientras mas se rechace H0 para un metodo (dado que es falso), mejor

# Notar que en el code hacemos la construccion de los optimos de Wassertein para cada metodo
# Y la construccion de las dist continuas que aproximan a cada metodo discretizado igual que en el caso de error tipo I
# Pues esta construccion partia asumiendo H0 para tener la dist continua objetivo inicial de los 5 metodos

#### Ejemplo hipergeometrica no central

# Parámetros base
ncase  <- 20
ncontr <- 100 #Fijo
nmuta  <- 2

alpha  <- 0.01
n_comb <- 100
iter   <- 1000            # subir para curvas más suaves (p.ej. 10000)

# Mallas arbitrarias (no simétricas)
logw_left  <- seq(-0.25, 0.02, by = 27/700)
logw_right <- seq(-0.02, 0.25, by = 27/700)     # distinta a la de left
logw_two   <- seq(-0.3,  0.5, by = 0.16)


print_every <- 100        # mostrar "." cada 200 iter

set.seed(123)


# P0
P0 <- prep_nc_experiment(ncase, ncontr, nmuta, alpha, n_comb)

# pmf bajo H0 (w=1)
support <- P0$support
pmf_H0  <- dhyper(support, nmuta, ncase + ncontr - nmuta, ncase)

# Tabla PMF/CDF de p-valores (izq, der, bilateral)
dist_df <- build_pval_distributions(support, pmf_H0)

# Subtítulo con parámetros
subt <- list(ncase = ncase, ncontr = ncontr, nmuta = nmuta)

# Plots
p_pmf <- plot_pval_pmf(dist_df,
                       title = "Distribución de los p-valores (PMF)",
                       subtitle_vars = subt)

p_cdf <- plot_pval_cdf(dist_df,
                       title = "Distribución de los p-valores (CDF)",
                       subtitle_vars = subt)

print(p_pmf)
print(p_cdf)


tab_power <- simulate_nc_power(
  ncase = ncase, ncontr = ncontr, nmuta = nmuta,
  alpha = alpha,
  n_comb = n_comb,
  iter = iter,
  logw_left  = logw_left,
  logw_right = logw_right,
  logw_two   = logw_two,
  print_every = print_every
)

# Curvas de poder por cola:
subt <- list(ncase = ncase, ncontr = ncontr, nmuta = nmuta, iter = iter)
p_all <- plot_power_all(
  tab_power,
  subtitle_vars = subt,
  facet_spacing_x = 2
)
print(p_all)




# Parámetros base
ncase  <- 20
ncontr <- 100 #Fijo
nmuta  <- 2

alpha  <- 0.05
n_comb <- 100
iter   <- 1000            # subir para curvas más suaves (p.ej. 10000)

# Mallas arbitrarias (no simétricas)
logw_left  <- seq(-0.2, 0.02, by = 11/350)
logw_right <- seq(-0.02, 0.2, by = 11/350)     # distinta a la de left
logw_two   <- seq(-0.3,  0.5, by = 0.16)


print_every <- 100        # mostrar "." cada 200 iter

set.seed(123)


# P0
P0 <- prep_nc_experiment(ncase, ncontr, nmuta, alpha, n_comb)

# pmf bajo H0 (w=1)
support <- P0$support
pmf_H0  <- dhyper(support, nmuta, ncase + ncontr - nmuta, ncase)

# Tabla PMF/CDF de p-valores (izq, der, bilateral)
dist_df <- build_pval_distributions(support, pmf_H0)

# Subtítulo con parámetros
subt <- list(ncase = ncase, ncontr = ncontr, nmuta = nmuta)

# Plots
p_pmf <- plot_pval_pmf(dist_df,
                       title = "Distribución de los p-valores (PMF)",
                       subtitle_vars = subt)

p_cdf <- plot_pval_cdf(dist_df,
                       title = "Distribución de los p-valores (CDF)",
                       subtitle_vars = subt)

print(p_pmf)
print(p_cdf)


tab_power <- simulate_nc_power(
  ncase = ncase, ncontr = ncontr, nmuta = nmuta,
  alpha = alpha,
  n_comb = n_comb,
  iter = iter,
  logw_left  = logw_left,
  logw_right = logw_right,
  logw_two   = logw_two,
  print_every = print_every
)

# Curvas de poder por cola:
subt <- list(ncase = ncase, ncontr = ncontr, nmuta = nmuta, iter = iter)
p_all <- plot_power_all(
  tab_power,
  subtitle_vars = subt,
  facet_spacing_x = 2
)
print(p_all)