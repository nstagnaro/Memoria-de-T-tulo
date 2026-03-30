#####calculadoras de distribuciones
## los x son los p-valores discretos asociados a cada obs
## cdf real de los p-valores (por ejemplo si son binom pues es la cdf de una binom bajo H0)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)  # <— para combinar plots

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

# Hacemos una funcion que generaliza los cumsum dependiendo de si es izq, der o two

pvals_from_pmf <- function(pmf, side=c("left","right","two.sided")){
  side <- match.arg(side)
  if (side=="left")      return(cumsum(pmf))
  if (side=="right")     return(rev(cumsum(rev(pmf))))
  sapply(pmf, function(p) sum(pmf[pmf<=p]))  # two.sided
}

# Hacemos una funcion que calcula las tablas con ratios de rechazo para alpha y lado generalizados
tabla_rechazos_todos <- function(alpha,
                                 lado   = c("left","right","two.sided"),
                                 support, pmf,
                                 n_comb = 100, iter = 100000,
                                 sampler = NULL, verbose = TRUE){
  
  lado <- match.arg(lado)
  stopifnot(length(support)==length(pmf), abs(sum(pmf)-1) < 1e-8)
  
  if (is.null(sampler)) {
    sampler <- function(n) sample(support, size = n, replace = TRUE, prob = pmf)
  }
  
  # p-valores discretos según el lado
  pv <- pvals_from_pmf(pmf, lado)
  
  # Mapas z_i por método (vector del tamaño de 'support')
  zmap <- list(
    Fisher    = fishercalc(pv, pv),
    Pearson   = pearsoncalc(pv, pv),
    Stouffer  = stouffercalc(pv, pv),
    Edgington = edgingtoncalc(pv, pv),
    George    = (pearsoncalc(pv, pv) - fishercalc(pv, pv))/2
  )
  
  # Varianzas de las distribuciones aproximadas por metodo
  vars <- list(
    Fisher    = varfisher(sort(pv)),
    Pearson   = varpearson(sort(pv)),
    Stouffer  = varstouffer(sort(pv)),
    Edgington = varedgington(sort(pv)),
    George    = vargeorge(sort(pv))
  )
  
  methods <- c("Fisher","Pearson","Stouffer","Edgington","George")
  
  # Matrices empíricas para cada método
  emp <- lapply(methods, function(.m) matrix(0, nrow = n_comb, ncol = iter))
  names(emp) <- methods
  
  # Simulación
  for (w in seq_len(iter)){
    X   <- sampler(n_comb)
    idx <- match(X, support) # Generaliza el X+1-Support[1] para dist en general (debiese ser similar a esto para la hipergeom centralizada)
    if (anyNA(idx)) stop("El sampler generó valores fuera de 'support'.")
    for (.m in methods){
      emp[[.m]][, w] <- cumsum(zmap[[.m]][idx])
    }
    if (verbose && (w %% 100 == 0)) message("iter = ", w)
  }
  
  # Cuantiles críticos por método (los que se comparan con casa est de prueba)
  i <- 1:n_comb
  crit <- list(
    Fisher    = qgamma(1 - alpha, shape = 4*i/vars$Fisher,    scale = vars$Fisher/2),   # >
    Pearson   = qgamma(alpha,     shape = 4*i/vars$Pearson,   scale = vars$Pearson/2),  # <
    Stouffer  = qnorm(alpha, mean = 0,     sd = sqrt(i*vars$Stouffer)),                 # <
    Edgington = qnorm(alpha, mean = i/2,   sd = sqrt(i*vars$Edgington)),                # <
    George    = qnorm(alpha, mean = 0,     sd = sqrt(i*vars$George))                    # <
  )
  
  # Proporciones de rechazo por método (los que van en la tabla final)
  out <- data.frame(
    i         = i,
    Fisher    = rowMeans(emp$Fisher    > crit$Fisher),
    Pearson   = rowMeans(emp$Pearson   < crit$Pearson),
    Stouffer  = rowMeans(emp$Stouffer  < crit$Stouffer),
    Edgington = rowMeans(emp$Edgington < crit$Edgington),
    George    = rowMeans(emp$George    < crit$George)
  )
  out$alpha <- alpha
  out$lado  <- lado
  out
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

.shape_map <- c(Edgington=16, Fisher=17, George=24, Pearson=3, Stouffer=15)

.make_subtitle <- function(vars){
  if (is.null(vars)) return(NULL)
  if (!is.list(vars) || is.null(names(vars)) || any(names(vars)==""))
    stop("subtitle_vars debe ser lista NOMBRADA, p.ej. list(ncase=..., ncontr=..., nmuta=..., iter=...)")
  fmt <- function(x) if (is.numeric(x)) format(x, trim=TRUE, scientific=FALSE, digits=6) else as.character(x)
  paste(sprintf("%s = %s", names(vars), vapply(vars, fmt, character(1))), collapse = ", ")
}

# Esta función elimina breaks automáticos que estén "encimados" con alpha
.get_smart_breaks <- function(limits, alpha) {
  # Generamos breaks bonitos (pedimos 5 para tener buena referencia)
  br <- scales::pretty_breaks(n = 5)(limits)
  
  # Calculamos la distancia mínima aceptable (15% del rango del panel)
  # Esto evita que los números se toquen visualmente
  range_size <- diff(limits)
  tolerance <- range_size * 0.15
  
  # Filtramos: quitamos los breaks que están demasiado cerca de alpha
  br <- br[abs(br - alpha) > tolerance]
  
  # Retornamos los breaks limpios incluyendo siempre el alpha
  sort(unique(c(br, alpha)))
}

plot_ratios_grid <- function(df_left, df_right, df_two,
                             alpha,
                             subtitle_vars = NULL,
                             title = NULL,
                             x_max = 100,
                             point_every = 5,
                             y_decimals = NULL){
  
  df_all <- dplyr::bind_rows(
    df_left  |> mutate(lado_label = "Cola izquierda"),
    df_right |> mutate(lado_label = "Cola derecha"),
    df_two   |> mutate(lado_label = "Bilateral")
  ) |>
    mutate(lado_label = factor(lado_label, levels = c("Cola izquierda", "Cola derecha", "Bilateral"))) |>
    tidyr::pivot_longer(cols = any_of(c("Fisher", "Pearson", "Stouffer", "Edgington", "George")),
                        names_to = "Metodo", values_to = "ratio")
  
  pts_all <- df_all |> dplyr::filter(i %% point_every == 0)
  
  # Si no se pasan decimales, calculamos cuántos ceros tiene alpha para no perder precisión
  if(is.null(y_decimals)) {
    # Si alpha es muy pequeño (p.ej 0.0001), necesitamos al menos 4 o 5 decimales.
    # Usamos format.info para detectar la precisión necesaria.
    precision <- max(4, -floor(log10(alpha)) + 1) 
  } else {
    precision <- y_decimals
  }
  
  # Usamos una función que decide entre formato fijo o científico según el valor
  label_fn <- function(x) {
    # Si el valor es muy pequeño pero no es cero, usamos formato 'g' o 'f' con precisión
    scales::label_number(accuracy = 10^-precision, decimal.mark = ".")(x)
  }
  
  # Construcción del gráfico
  p <- ggplot(df_all, aes(x = i, y = ratio, color = Metodo, shape = Metodo)) +
    geom_hline(yintercept = alpha, color = "red", linewidth = 0.7) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_point(data = pts_all, size = 1.8, alpha = 0.9) +
    
    facet_wrap(~ lado_label, nrow = 1, scales = "free_y") +
    
    # Eje Y optimizado
    scale_y_continuous(
      breaks = function(limits) .get_smart_breaks(limits, alpha),
      labels = label_fn
    ) +
    
    scale_shape_manual(values = .shape_map, name = "Método") +
    scale_color_discrete(name = "Método") +
    scale_x_continuous(limits = c(0, x_max),
                       breaks = seq(0, x_max, length.out = 5),
                       labels = scales::label_number(accuracy = 1)) +
    
    labs(
      title    = if (is.null(title)) sprintf("P-valores hipergeométricos centralizados — Error Tipo I=%g", alpha) else title,
      subtitle = .make_subtitle(subtitle_vars),
      x = "P-valores combinados",
      y = "Ratio de rechazo"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title.position = "plot",
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(), 
          plot.subtitle = element_text(size = 11, margin = margin(b = 10)))
  
  return(p)
}

#### Ejemplo hipergeometrica

# Recordar que esta distribucion tiene los parametros NCHYG(ncase,ncontr,nmuta,w)
# w=1 (en un principio veremos la centralizada)

ncase <- 100;
ncontr <- 500;
nmuta <- 10
support <- max(0, ncase+nmuta-(ncase+ncontr)):min(ncase, nmuta)
pmf <- dhyper(support, nmuta, ncase+ncontr-nmuta, ncase) #f.d.p
# Esta dhyper asume el w=1, pues es la hipergeometrica centralizada

iter<-100000
n_comb<-100 # Cantidad de p-vals que se combinaran (se generaliza el 100)


# Creamos las obs aleatorias con dist hiper como en el ejemplo particular
sampler_hyg <- function(n) rhyper(n, nmuta, ncase+ncontr-nmuta, ncase)
set.seed(999)
pvcontccl0001 <- tabla_rechazos_todos(
  alpha   = 0.0001,
  lado    = "left",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontccl001 <- tabla_rechazos_todos(
  alpha   = 0.001,
  lado    = "left",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontccl01 <- tabla_rechazos_todos(
  alpha   = 0.01,
  lado    = "left",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontccr0001 <- tabla_rechazos_todos(
  alpha   = 0.0001,
  lado    = "right",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontccr001 <- tabla_rechazos_todos(
  alpha   = 0.001,
  lado    = "right",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontccr01 <- tabla_rechazos_todos(
  alpha   = 0.01,
  lado    = "right",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontcct0001 <- tabla_rechazos_todos(
  alpha   = 0.0001,
  lado    = "two.sided",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontcct001 <- tabla_rechazos_todos(
  alpha   = 0.001,
  lado    = "two.sided",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)
pvcontcct01 <- tabla_rechazos_todos(
  alpha   = 0.01,
  lado    = "two.sided",
  support = support,
  pmf     = pmf,
  n_comb  = n_comb,
  iter    = iter,
  sampler = sampler_hyg,
  verbose = TRUE
)

# Una pregunta natural es el pensar porque los errores tipo I estimados para cada una de las metodologias se calculan así:
# La respuesta es que sabemos que el error tipo I es la prob de rechazar H0 dado de que es V
# En este caso estamos contando los rechazos comparando con los cuantiles de las dist continuas approx
# Bajo los supuestos de H0 a la hora de construir los valores observados
# Es decir, rechazamos dado que H0 es V xd, literal la def de error tipo I
# La cuestion es que acá al igual que en el caso de potencia el contraste de hip es claramente:
# H0: \psi = 1  v/s  H1: \psi\neq 1 (< o > respectivamente)

# Lo que pasa es que para el mundo del control de error tipo I todo se hace combinando p-valores HG centralizados
# Porque se asume H0 es V

# Y es obvio que queremos que las estimaciones esten cerca del nivel nominal, así el test
# se aproxima mejor a lo que queremos a nivel de rechazos

###################################################################################

# Construir tablas de p-valores para izq, der y bilateral

dist_df <- build_pval_distributions(support, pmf, digits = 12)

subt <- list(ncase = ncase, ncontr = ncontr, nmuta = nmuta)

p_pmf <- plot_pval_pmf(dist_df,
                       title = "Distribución de los p-valores (PMF)",
                       subtitle_vars = subt)

p_cdf <- plot_pval_cdf(dist_df,
                       title = "Distribución de los p-valores",
                       subtitle_vars = subt)

print(p_pmf)
print(p_cdf)

# Construir los 3 graficos de ratios de rechazo

subt <- list(ncase = ncase, ncontr = ncontr, nmuta = nmuta, iter = iter)

# Todos los graficos alpha = 0.01
p_01 <- plot_ratios_grid(
  df_left  = pvcontccl01,   
  df_right = pvcontccr01,
  df_two   = pvcontcct01,
  alpha    = 0.01,
  subtitle_vars = subt
)
print(p_01)

# Todos los graficos alpha = 0.001
p_001 <- plot_ratios_grid(
  df_left  = pvcontccl001,   
  df_right = pvcontccr001,
  df_two   = pvcontcct001,
  alpha    = 0.001,
  subtitle_vars = subt
)
print(p_001)

# Todos los graficos alpha = 0.0001
p_0001 <- plot_ratios_grid(
  df_left  = pvcontccl0001,   
  df_right = pvcontccr0001,
  df_two   = pvcontcct0001,
  alpha    = 0.0001,
  subtitle_vars = subt
)
print(p_0001)


###################################################################################