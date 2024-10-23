#' Cálculo de Tamaño Muestral y Potencia para Estudios de Razón de Prevalencia
#' 
#' @param alfa Nivel de significación
#' @param potencia Potencia deseada (opcional si se proporciona n)
#' @param n Tamaño muestral (opcional si se proporciona potencia)
#' @param p0 Prevalencia en el grupo no expuesto (requerido si metodo = "rho")
#' @param RP Razón de Prevalencia (opcional si se proporciona p1)
#' @param p1 Prevalencia en el grupo expuesto (opcional si se proporciona RP)
#' @param IC_superior Límite superior del intervalo de confianza de la RP
#' @param IC_inferior Límite inferior del intervalo de confianza de la RP (opcional si se proporciona EE)
#' @param EE Error estándar del coeficiente (opcional si se proporcionan los IC)
#' @param n_previo Tamaño de muestra del estudio previo (requerido para método "ee")
#' @param r Razón de tamaños de grupo (expuestos/no expuestos), por defecto es 1
#' @param metodo Método de cálculo: "rho" o "ee" (error estándar)
#' @param valores_rho Vector de valores rho para calcular tamaños muestrales (si metodo = "rho")
#' @return Un data frame con resultados según el método elegido
#' @export
SamplePrevRatioMulti <- function(alfa, potencia = NULL, n = NULL,
                                p0 = NULL, RP = NULL, p1 = NULL,
                                IC_superior = NULL, IC_inferior = NULL,
                                EE = NULL, n_previo = NULL, r = 1, 
                                metodo = "rho",
                                valores_rho = seq(0, 0.9, by = 0.1)) {
  
  # Verificaciones según el método
  if (metodo == "rho") {
    if (is.null(p0)) {
      stop("Para método 'rho' se requiere p0")
    }
    if (is.null(RP) && is.null(p1)) {
      stop("Para método 'rho' se debe proporcionar RP o p1")
    }
    if (!is.null(RP) && !is.null(p1)) {
      stop("Solo se debe proporcionar RP o p1, no ambos")
    }
  } else if (metodo == "ee") {
    if (is.null(RP)) {
      stop("Para método 'ee' se requiere RP")
    }
    if (is.null(EE) && (is.null(IC_superior) || is.null(IC_inferior))) {
      stop("Para método 'ee' se requiere EE o ambos límites del IC")
    }
    if (is.null(n_previo)) {
      stop("Para método 'ee' se requiere n_previo")
    }
  } else {
    stop("El método debe ser 'rho' o 'ee'")
  }
  
  # Si se proporcionaron los IC, calcular el EE
  if (!is.null(IC_superior) && !is.null(IC_inferior)) {
    EE <- (log(IC_superior) - log(RP)) / qnorm(0.975)
    cat("Error Estándar calculado:", EE, "\n")
  }
  
  if (metodo == "rho") {
    # Calcular p1 si se proporcionó RP
    if (!is.null(RP)) {
      p1 <- RP * p0
    } else {
      RP <- p1/p0
    }
    
    q1 <- r / (1 + r)
    q0 <- 1 / (1 + r)
    z_alfa <- qnorm(1 - alfa/2)
    p_agrupada <- (q1 * p1 + q0 * p0)
    
    if (!is.null(potencia)) {
      z_beta <- qnorm(potencia)
      calcular_n <- function(rho) {
        A <- z_alfa * sqrt(p_agrupada * (1 - p_agrupada) * (1/q1 + 1/q0))
        B <- z_beta * sqrt(p1 * (1 - p1) / q1 + p0 * (1 - p0) / q0)
        C <- (p1 - p0)^2
        N <- ((A + B)^2 / C) / (1 - rho^2)
        return(ceiling(N))
      }
      resultados <- data.frame(
        rho = valores_rho,
        tamaño_total = sapply(valores_rho, calcular_n)
      )
      resultados$tamaño_expuestos <- ceiling(resultados$tamaño_total * q1)
      resultados$tamaño_no_expuestos <- ceiling(resultados$tamaño_total * q0)
    } else {
      calcular_potencia <- function(rho) {
        A <- z_alfa * sqrt(p_agrupada * (1 - p_agrupada) * (1/q1 + 1/q0))
        C <- (p1 - p0)^2
        B <- sqrt(n * C * (1 - rho^2)) - A
        potencia <- pnorm(B / sqrt(p1 * (1 - p1) / q1 + p0 * (1 - p0) / q0))
        return(potencia)
      }
      resultados <- data.frame(
        rho = valores_rho,
        potencia = sapply(valores_rho, calcular_potencia)
      )
    }
  } else {
    # Método basado en error estándar
    z_alfa <- qnorm(1 - alfa/2)
    if (!is.null(potencia)) {
      z_gamma <- qnorm(potencia)
      n <- ceiling((z_alfa + z_gamma)^2 * n_previo * EE^2 / (log(RP))^2)
      resultados <- data.frame(
        tamaño_total = n
      )
      resultados$tamaño_expuestos <- ceiling(n * r/(1 + r))
      resultados$tamaño_no_expuestos <- ceiling(n * 1/(1 + r))
    } else {
      z_gamma <- sqrt(n * (log(RP))^2 / (n_previo * EE^2)) - z_alfa
      potencia <- pnorm(z_gamma)
      resultados <- data.frame(
        potencia = potencia
      )
    }
  }
  return(resultados)
}

#' Logística del Estudio Transversal
#'
#' @param n_final Tamaño muestral final calculado
#' @param tasa_rechazo Tasa de rechazo esperada (proporción entre 0 y 1)
#' @param tasa_elegibilidad Tasa de elegibilidad esperada (proporción entre 0 y 1)
#' @param sujetos_por_dia Número de sujetos/registros que pueden ser procesados por día
#' @param dias_laborables_mes Número de días laborables por mes
#' @return Una lista con los cálculos logísticos del estudio
#' @export
logistica_estudio_transversal <- function(n_final, 
                                         tasa_rechazo, 
                                         tasa_elegibilidad, 
                                         sujetos_por_dia, 
                                         dias_laborables_mes) {
  
  # Validación de entradas
  if (tasa_rechazo < 0 || tasa_rechazo >= 1) 
    stop("La tasa de rechazo debe estar entre 0 y 1")
  if (tasa_elegibilidad <= 0 || tasa_elegibilidad > 1) 
    stop("La tasa de elegibilidad debe estar entre 0 y 1")
  if (sujetos_por_dia <= 0) 
    stop("El número de sujetos por día debe ser positivo")
  if (dias_laborables_mes <= 0) 
    stop("El número de días laborables por mes debe ser positivo")
  
  # Cálculos
  n_evaluar <- n_final / (1 - tasa_rechazo)
  n_invitar <- n_evaluar / tasa_elegibilidad
  total_dias <- n_invitar / sujetos_por_dia
  duracion_meses <- total_dias / dias_laborables_mes
  
  # Resultados
  resultados <- list(
    muestra_final = n_final,
    muestra_evaluar = ceiling(n_evaluar),
    muestra_invitar = ceiling(n_invitar),
    dias_reclutamiento = ceiling(total_dias),
    meses_reclutamiento = round(duracion_meses, 2)
  )
  
  # Imprimir resumen
  cat("
Resumen de logística del estudio:
")
  cat("----------------------------------------
")
  cat("Muestra final requerida:", n_final, "
")
  cat("Personas a evaluar:", ceiling(n_evaluar), 
      "(considerando", tasa_rechazo*100, "% de rechazo)
")
  cat("Personas a invitar:", ceiling(n_invitar), 
      "(considerando", tasa_elegibilidad*100, "% de elegibilidad)
")
  cat("Días necesarios:", ceiling(total_dias),
      "(evaluando", sujetos_por_dia, "personas por día)
")
  cat("Meses necesarios:", round(duracion_meses, 2),
      "(con", dias_laborables_mes, "días laborables por mes)
")
  cat("----------------------------------------
")
  
  return(invisible(resultados))
}
