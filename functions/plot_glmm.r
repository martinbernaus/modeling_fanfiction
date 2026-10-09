
plot_glmm <- function(xdata, model, xvar, yvar, xlabel, ylabel, family = "gaussian",
                      fit_fun = "auto", nb_link = NULL, offset_value = 1,
                      random_var = FALSE, random_slope = FALSE, random_formula = FALSE,
                      glmer_control = NULL, lmer_control = NULL, glmmTMB_control = NULL, REML = NULL,
                      ci_method = "wald", ci_level = 0.95,
                      nboots = 1000, boot_seed = NULL, boot_parallel = "no", boot_ncpus = 1,
                      aggregate_by = FALSE,
                      factor_var = FALSE, factor_label = NULL,
                      slice_var = FALSE, slice_label = NULL,
                      color_var = FALSE, size_var = FALSE, extra_var_label = NULL,
                      size_factor = 0.5, barplot_style = FALSE, point_style = "outline",
                      log_yaxis = "auto", log_xaxis = "auto", log_ybase = 2,
                      legend_location = "topright", legend_outside_width = "auto") {

  xvar_is_factor <- is.factor(xdata[[xvar]])
  crossed <- xvar_is_factor && !isFALSE(factor_var)

  if (xvar_is_factor && !isFALSE(slice_var))
    stop("xvar is already a factor; slice_var is not supported in this mode.")
  if (crossed && !is.factor(xdata[[factor_var]]))
    stop("factor_var must be a factor when crossed with a factor xvar.")
  if (!xvar_is_factor && !isFALSE(factor_var) && !isFALSE(slice_var))
    stop("Use either factor_var or slice_var, not both.")
  if (crossed && !isTRUE(barplot_style)) {
    warning("plot_glmm: crossed factor design (xvar + factor_var both factors) requires barplot_style -- forcing barplot_style = TRUE.")
    barplot_style <- TRUE
  }
  if (!identical(ci_method, "wald") && !identical(ci_method, "boot"))
    stop("ci_method must be 'wald' or 'boot'.")
  if (!identical(point_style, "outline") && !identical(point_style, "solid"))
    stop("point_style must be 'outline' or 'solid'.")

  # --- Legend positioning: "outside" places legends in a right margin sized to fit
  # the actual legend content (via strwidth), avoiding excess empty space. ---
  legend_outside <- identical(legend_location, "outside")

  compute_outside_margin <- function() {
    texts <- character(0)
    if (!isFALSE(factor_var)) {
      texts <- c(texts, factor_label, levels(xdata[[factor_var]]))
    } else if (!isFALSE(slice_var)) {
      texts <- c(texts, slice_label, "-1 SD", "Mean", "+1 SD")
    }
    if (!isFALSE(color_var)) {
      lbl <- if (!is.null(extra_var_label)) extra_var_label else color_var
      brks <- pretty(range(xdata[[color_var]], na.rm = TRUE), n = 4)
      texts <- c(texts, lbl, as.character(round(brks, 2)))
    }
    if (!isFALSE(size_var)) {
      lbl <- if (!is.null(extra_var_label)) extra_var_label else size_var
      brks <- pretty(range(xdata[[size_var]], na.rm = TRUE), n = 4)
      texts <- c(texts, lbl, as.character(round(brks, 2)))
    }
    if (length(texts) == 0) texts <- "legend"
    op <- par(no.readonly = TRUE); on.exit(par(op))
    par(family = "sans")
    max_w_in <- max(strwidth(texts, units = "inches", cex = 1)) * 1.05
    line_h_in <- par("csi")
    mar_lines <- max_w_in / line_h_in + 3
    min(max(mar_lines, 4), 16)
  }

  right_mar <- if (legend_outside) {
    if (identical(legend_outside_width, "auto")) compute_outside_margin() else legend_outside_width
  } else 0.5

  outside_next_y <- NULL
  draw_legend <- function(...) {
    if (legend_outside) {
      usr <- par("usr")
      pin <- par("pin")
      x_user_per_in <- (usr[2] - usr[1]) / pin[1]
      y_user_per_in <- (usr[4] - usr[3]) / pin[2]
      x_pos <- usr[2] + 0.15 * x_user_per_in
      if (is.null(outside_next_y)) outside_next_y <<- usr[4]
      lg <- legend(x = x_pos, y = outside_next_y, xpd = NA, xjust = 0, yjust = 1, ...)
      outside_next_y <<- lg$rect$top - lg$rect$h - 0.3 * y_user_per_in
      invisible(lg)
    } else {
      legend(x = legend_location, inset = 0.02, ...)
    }
  }

  # --- Shared small helpers (deduplicated across Branch A / Branch B) ---
  set_plot_mar <- function(mar, lwd = NULL) {
    args <- list(family = "sans", mar = mar, mgp = c(2.5, 0.5, 0), tcl = -0.3, las = 1,
                 cex.axis = 1.1, cex.lab = 1.2, font.lab = 2)
    if (!is.null(lwd)) args$lwd <- lwd
    do.call(par, args)
  }

  compute_log_ticks <- function(values, base, verbose_label = NULL) {
    if (!base %in% c(2, 10)) stop("log_ybase must be 2 or 10.")
    rng <- range(exp(values), na.rm = TRUE)
    y.at <- if (base == 2) 2^(floor(log2(rng[1])):ceiling(log2(rng[2])))
            else 10^(floor(log10(rng[1])):ceiling(log10(rng[2])))
    y.at <- round(y.at, 4)
    if (!is.null(verbose_label)) { cat("\n[", verbose_label, "] y-axis ticks:\n", sep = ""); print(y.at) }
    y.at
  }

  make_ci_df <- function(fit_link, se_link, crit, linkinv) {
    data.frame(fit = linkinv(fit_link),
               lwr = linkinv(fit_link - crit * se_link),
               upr = linkinv(fit_link + crit * se_link))
  }

  apply_wrapped_nuisance <- function(pred.data, xdata, wrapped_nuisance_vars) {
    for (v in wrapped_nuisance_vars) {
      if (is.factor(xdata[[v]])) {
        pred.data[[v]] <- factor(levels(xdata[[v]])[1], levels = levels(xdata[[v]]))
        cat("    -> wrapped nuisance factor:", v, "held at reference level:", levels(xdata[[v]])[1], "\n")
      } else {
        pred.data[[v]] <- mean(xdata[[v]], na.rm = TRUE)
        cat("    -> wrapped nuisance continuous:", v, "held at mean:", round(mean(xdata[[v]], na.rm = TRUE), 3), "\n")
      }
    }
    pred.data
  }

  compute_point_aesthetics <- function(xdata, color_var, size_var, size_factor) {
    out <- list(point_cols = NULL, cvals = NULL, rng_c = NULL, palette = NULL,
                point_cex = NULL, svals = NULL, rng_s = NULL, cex_range = c(0.5, 2.5))
    if (!isFALSE(color_var)) {
      out$cvals <- xdata[[color_var]]
      out$rng_c <- range(out$cvals, na.rm = TRUE)
      out$palette <- colorRampPalette(c("white", "black"))(100)
      out$point_cols <- out$palette[round(1 + 99 * (out$cvals - out$rng_c[1]) / diff(out$rng_c))]
      cat("\n[14] color_var:", color_var, "| range:", out$rng_c, "\n")
    } else cat("\n[14] No color_var requested.\n")
    if (!isFALSE(size_var)) {
      out$svals <- xdata[[size_var]]
      out$rng_s <- range(out$svals, na.rm = TRUE)
      norm <- (out$svals - out$rng_s[1]) / diff(out$rng_s)
      out$point_cex <- out$cex_range[1] + size_factor * diff(out$cex_range) * sqrt(norm)
      cat("[14] size_var:", size_var, "| range:", out$rng_s,
          "| size_factor:", size_factor, "| cex range:", range(out$point_cex, na.rm = TRUE), "\n")
    } else cat("[14] No size_var requested.\n")
    out
  }

  draw_aes_legends <- function(pa, color_var, size_var, extra_var_label, size_factor, draw_legend) {
    if (!isFALSE(color_var)) {
      lbl  <- if (!is.null(extra_var_label)) extra_var_label else color_var
      mid  <- pretty(pa$cvals, n = 4); mid <- mid[mid > pa$rng_c[1] & mid < pa$rng_c[2]]
      brks <- round(sort(unique(c(pa$rng_c[1], mid, pa$rng_c[2]))), 2)
      cols <- pa$palette[round(1 + 99 * (brks - pa$rng_c[1]) / diff(pa$rng_c))]
      cat("\n[17] color legend breaks (min/max forced in):\n"); print(brks)
      draw_legend(legend = brks, pt.bg = cols, pch = 21, col = "black",
                  pt.cex = 1.4, bty = "n", title = lbl, title.font = 2)
    } else cat("\n[17] No color legend needed.\n")
    if (!isFALSE(size_var)) {
      lbl  <- if (!is.null(extra_var_label)) extra_var_label else size_var
      mid  <- pretty(pa$svals, n = 4); mid <- mid[mid > pa$rng_s[1] & mid < pa$rng_s[2]]
      brks <- round(sort(unique(c(pa$rng_s[1], mid, pa$rng_s[2]))), 2)
      cexs <- pa$cex_range[1] + size_factor * diff(pa$cex_range) * sqrt((brks - pa$rng_s[1]) / diff(pa$rng_s))
      cat("[17] size legend breaks (min/max forced in):\n"); print(brks)
      draw_legend(legend = brks, pch = 21, pt.bg = "white", col = "black",
                  pt.cex = cexs, bty = "n", title = lbl, title.font = 2)
    } else cat("[17] No size legend needed.\n")
  }

  # Unified raw-point drawing. pch_shape: group-specific shape used ONLY when
  # point_style=="outline" and neither point_cols nor point_cex is supplied
  # (i.e. groups distinguished by shape rather than by fill color).
  draw_points <- function(x, y, base_col, point_style, point_cols = NULL, point_cex = NULL,
                           sel = NULL, pch_shape = NULL, default_bg = "white") {
    if (!is.null(sel)) { x <- x[sel]; y <- y[sel] }
    pc   <- if (!is.null(point_cols) && !is.null(sel)) point_cols[sel] else point_cols
    pcex <- if (!is.null(point_cex) && !is.null(sel)) point_cex[sel] else point_cex
    cex_use <- if (!is.null(pcex)) pcex else 1.1

    if (identical(point_style, "solid")) {
      bg  <- if (!is.null(pc)) pc else adjustcolor(base_col, alpha.f = 0.45)
      col <- if (!is.null(pc)) "black" else adjustcolor(base_col, alpha.f = 0.8)
      points(x, y, pch = 21, bg = bg, col = col, cex = cex_use)
    } else if (!is.null(pc) || !is.null(pcex) || is.null(pch_shape)) {
      points(x, y, pch = 21, bg = if (!is.null(pc)) pc else default_bg,
             col = if (is.null(pch_shape)) "black" else base_col, cex = cex_use)
    } else {
      points(x, y, pch = pch_shape, col = "black", cex = cex_use)
    }
  }

  is_merMod  <- inherits(model, "merMod")
  is_lmerMod <- inherits(model, "lmerMod")
  is_glmmTMB <- inherits(model, "glmmTMB")
  is_mixed_model <- is_merMod || is_glmmTMB

  fit_fun_lc <- if (is.character(fit_fun)) tolower(fit_fun) else fit_fun
  want_tmb  <- is_glmmTMB || identical(fit_fun_lc, "glmmtmb")
  want_glmm <- !want_tmb && (is_merMod || !isFALSE(random_var) || !isFALSE(random_formula) ||
               identical(fit_fun_lc, "glmer") || identical(fit_fun_lc, "glmer.nb") ||
               identical(fit_fun_lc, "lmer"))
  want_mixed <- want_tmb || want_glmm

  if (want_mixed && isFALSE(random_var) && isFALSE(random_formula) && !is_mixed_model)
    stop("GLMM/LMM requested but no random_var/random_formula given, and template `model` is not a merMod/glmmTMB to borrow random structure from.")
  if (want_mixed && !requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required (for nobars()/findbars() formula utilities) but is not installed.")
  if (want_tmb && !requireNamespace("glmmTMB", quietly = TRUE))
    stop("Package 'glmmTMB' is required for glmmTMB support but is not installed.")
  if (identical(ci_method, "boot") && !want_mixed) {
    warning("plot_glmm: ci_method = 'boot' is only supported for GLMMs/LMMs/glmmTMB models -- falling back to 'wald'.")
    ci_method <- "wald"
  }

  is_gaussian_family <- if (want_tmb) FALSE else if (identical(fit_fun_lc, "lmer")) TRUE else if (is.character(family)) {
    tolower(family) == "gaussian"
  } else {
    tryCatch({
      fam_obj <- if (is.function(family)) family() else family
      identical(fam_obj$family, "gaussian") && identical(fam_obj$link, "identity")
    }, error = function(e) FALSE)
  }
  if (is_merMod && is_lmerMod) is_gaussian_family <- TRUE

  cat("\n========== plot_glmm START ==========\n")
  cat("xvar:", xvar, "| yvar:", yvar, "| family:", if (is.character(family)) family else "<object>",
      "| fit_fun:", fit_fun, "| is_merMod:", is_merMod, "| is_lmerMod:", is_lmerMod,
      "| is_glmmTMB:", is_glmmTMB, "| want_tmb:", want_tmb,
      "| is_gaussian_family:", is_gaussian_family, "| want_glmm:", want_glmm, "| random_var:", random_var,
      "| ci_method:", ci_method, "| ci_level:", ci_level, "| aggregate_by:", aggregate_by,
      "| xvar_is_factor:", xvar_is_factor, "| crossed:", crossed,
      "| factor_var:", factor_var, "| slice_var:", slice_var,
      "| barplot_style:", barplot_style, "| point_style:", point_style,
      "| legend_location:", legend_location, "| right_mar:", right_mar,
      "| color_var:", color_var, "| size_var:", size_var, "\n")

  fe_formula_orig <- if (is_mixed_model) lme4::nobars(formula(model)) else formula(model)
  model_terms_obj <- terms(fe_formula_orig)
  term.labels <- attr(model_terms_obj, "term.labels")
  term_vars   <- lapply(term.labels, function(tl) all.vars(str2lang(tl)))
  n_term_vars <- sapply(term_vars, length)
  cat("\n[1] term.labels (fixed-effects terms from original model):\n"); print(term.labels)

  offset_pos <- attr(model_terms_obj, "offset")
  if (!is.null(offset_pos)) {
    vars_list    <- as.list(attr(model_terms_obj, "variables"))[-1]
    offset_exprs <- vars_list[offset_pos]
    offset_terms <- sapply(offset_exprs, function(e) paste(deparse(e), collapse = ""))
    offset_vars  <- unique(unlist(lapply(offset_exprs, all.vars)))
  } else {
    offset_exprs <- list()
    offset_terms <- character(0)
    offset_vars  <- character(0)
  }
  cat("\n[1b] offset_terms (preserved verbatim in refit formula):\n"); print(offset_terms)
  cat("[1b] offset_vars (held at offset_value =", offset_value, "in pred.data, interpreted as 'rate'):\n")
  print(offset_vars)

  yvar_vals <- xdata[[yvar]]
  yvar_is_matrix <- is.matrix(yvar_vals)
  if (yvar_is_matrix) {
    if (ncol(yvar_vals) != 2)
      stop("plot_glmm: yvar is a matrix with ", ncol(yvar_vals),
           " columns -- expected a 2-column cbind(success, failure) style binomial response.")
    y_plot <- yvar_vals[, 1] / rowSums(yvar_vals)
    cat("[1b] yvar is a 2-column matrix (cbind(success, failure) binomial response).\n")
    cat("[1b] y_plot computed as proportion = successes / (successes + failures), for raw-point/bar plotting.\n")
    if (length(offset_vars) > 0)
      warning("plot_glmm: yvar is matrix-valued; offset_value rescaling is skipped for y_plot (not meaningful for a proportion response).")
  } else if (length(offset_vars) > 0) {
    y_plot <- yvar_vals / xdata[[offset_vars[1]]] * offset_value
    cat("[1b] y_plot: yvar converted to rate (yvar /", offset_vars[1], "* offset_value) for plotting raw points.\n")
    cat("[1b] NOTE: y-axis now represents a RATE, not a raw count. Update ylabel accordingly if needed.\n")
  } else {
    y_plot <- yvar_vals
  }

  if (!isFALSE(aggregate_by)) {
    if (!(aggregate_by %in% names(xdata)))
      stop("aggregate_by: column '", aggregate_by, "' not found in xdata.")
    n_unique_xvar <- tapply(xdata[[xvar]], xdata[[aggregate_by]], function(v) length(unique(v)))
    if (any(n_unique_xvar > 1))
      warning("plot_glmm: aggregate_by = '", aggregate_by, "' has groups with more than one unique ",
              "value of xvar ('", xvar, "') -- aggregated points may not align to a single x-position.")
    y_plot <- ave(y_plot, xdata[[aggregate_by]], FUN = function(v) mean(v, na.rm = TRUE))
    cat("[1b] y_plot aggregated: mean(", yvar, ") computed within each level of '", aggregate_by,
        "', broadcast back to all rows (", length(unique(xdata[[aggregate_by]])), " groups).\n", sep = "")
  }

  xvar_scaled   <- if (xvar_is_factor) FALSE else grepl("^z\\.", xvar)
  xvar_prescale <- if (xvar_is_factor) xvar else sub("^z\\.", "", xvar)
  xvar_log      <- if (xvar_is_factor) FALSE else
                    if (identical(log_xaxis, "auto")) grepl("^log\\.", xvar_prescale) else isTRUE(log_xaxis)
  xvar_raw      <- if (xvar_log) sub("^log\\.", "", xvar_prescale) else xvar_prescale
  yvar_log <- if (identical(log_yaxis, "auto")) grepl("^log\\.", yvar) else isTRUE(log_yaxis)

  xvar_has_raw <- xvar_prescale %in% names(xdata)
  if (!xvar_is_factor && (xvar_scaled || xvar_log) && !xvar_has_raw) {
    warning("plot_glmm: raw column '", xvar_prescale, "' not found in xdata -- ",
            "x-axis will show z-scores/model-scale values as-is instead of back-transformed values.")
  }

  factor_label <- if (!isFALSE(factor_var) && is.null(factor_label)) factor_var else factor_label
  slice_prescale <- if (!isFALSE(slice_var)) sub("^z\\.", "", slice_var) else NULL
  slice_label <- if (!isFALSE(slice_var) && is.null(slice_label)) slice_prescale else slice_label
  cat("\n[2] xvar_scaled:", xvar_scaled, "| xvar_log:", xvar_log,
      "| xvar_raw:", xvar_raw, "| xvar_has_raw:", xvar_has_raw, "| yvar_log:", yvar_log,
      "| factor_label:", factor_label, "| slice_prescale:", slice_prescale,
      "| slice_label:", slice_label, "\n")

  slice_final <- FALSE
  if (!isFALSE(slice_var)) {
    if (slice_var %in% names(xdata)) {
      slice_final <- slice_var
    } else if (paste0("z.", slice_prescale) %in% names(xdata)) {
      slice_final <- paste0("z.", slice_prescale)
    } else if (slice_prescale %in% names(xdata)) {
      slice_final <- paste0("z.", slice_prescale)
      xdata[[slice_final]] <- as.vector(scale(xdata[[slice_prescale]]))
      cat("    -> created:", slice_final, "\n")
    } else {
      stop("slice_var: could not find or create a matching column for '", slice_var, "'.")
    }
  }
  cat("[2] slice_final (resolved column used for slicing):", slice_final, "\n")

  xvar_candidates   <- unique(c(xvar, xvar_prescale, xvar_raw))
  slice_candidates  <- if (!isFALSE(slice_final)) unique(c(slice_var, slice_prescale, slice_final)) else character(0)
  factor_candidates <- if (!isFALSE(factor_var)) factor_var else character(0)

  is_xvar_term   <- n_term_vars == 1 & sapply(term_vars, function(v) v[1] %in% xvar_candidates)
  is_factor_term <- n_term_vars == 1 & sapply(term_vars, function(v) v[1] %in% factor_candidates)
  is_slice_term  <- n_term_vars == 1 & sapply(term_vars, function(v) v[1] %in% slice_candidates)
  is_nuisance    <- n_term_vars == 1 & !is_xvar_term & !is_factor_term & !is_slice_term
  is_interaction <- n_term_vars > 1
  cat("\n[3] xvar term(s) matched:\n"); print(term.labels[is_xvar_term])
  cat("[3] factor term matched:\n"); print(term.labels[is_factor_term])
  cat("[3] slice term matched:\n"); print(term.labels[is_slice_term])
  cat("[3] nuisance terms:\n"); print(term.labels[is_nuisance])
  cat("[3] interaction terms detected (preserved automatically):\n"); print(term.labels[is_interaction])

  xvar_terms <- term.labels[is_xvar_term]
  if (length(xvar_terms) == 0) {
    xvar_terms_final <- xvar
  } else {
    xvar_syms <- sapply(term_vars[is_xvar_term], `[`, 1)
    xvar_terms_final <- mapply(function(term, sym) gsub(sym, xvar, term, fixed = TRUE),
                                xvar_terms, xvar_syms, USE.NAMES = FALSE)
  }
  cat("\n[4] xvar_terms_final (rebuilt, used in refit formula):\n"); print(xvar_terms_final)

  nuisance_terms <- term.labels[is_nuisance]
  nuisance_vars  <- if (length(nuisance_terms) > 0) sapply(term_vars[is_nuisance], `[`, 1) else character(0)
  is_wrapped     <- nuisance_terms != nuisance_vars
  simple_nuisance_vars   <- nuisance_vars[!is_wrapped]
  wrapped_nuisance_terms <- nuisance_terms[is_wrapped]
  wrapped_nuisance_vars  <- nuisance_vars[is_wrapped]
  cat("\n[5] simple nuisance vars (auto z-scale/center as before):\n"); print(simple_nuisance_vars)
  cat("[5] wrapped nuisance terms (reused verbatim, base var held at mean/reference):\n")
  print(wrapped_nuisance_terms)

  if (length(simple_nuisance_vars) > 0) {
    is_factor_simple  <- sapply(xdata[simple_nuisance_vars], is.factor)
    cont_predictors   <- simple_nuisance_vars[!is_factor_simple]
    factor_predictors <- simple_nuisance_vars[is_factor_simple]
  } else {
    cont_predictors <- character(0); factor_predictors <- character(0)
  }
  cat("\n[6] cont_predictors:\n"); print(cont_predictors)
  cat("[6] factor_predictors:\n"); print(factor_predictors)

  scaled     <- cont_predictors[grepl("^z\\.", cont_predictors)]
  unscaled   <- setdiff(cont_predictors, scaled)
  z_unscaled <- if (length(unscaled) > 0) paste0("z.", unscaled) else character(0)
  for (i in seq_along(unscaled)) {
    xdata[[z_unscaled[i]]] <- as.vector(scale(xdata[[unscaled[i]]]))
    cat("    -> created:", z_unscaled[i], "| mean:", round(mean(xdata[[unscaled[i]]], na.rm = TRUE), 3),
        "| sd:", round(sd(xdata[[unscaled[i]]], na.rm = TRUE), 3), "\n")
  }
  z_names <- c(scaled, z_unscaled)
  cat("[6] z_names:\n"); print(z_names)

  f_names <- character(0)
  if (length(factor_predictors) > 0) {
    f_names <- paste0(factor_predictors, ".code")
    for (i in seq_along(factor_predictors)) {
      f     <- factor_predictors[i]
      dummy <- as.numeric(xdata[[f]] == levels(xdata[[f]])[2])
      xdata[[f_names[i]]] <- dummy - mean(dummy)
      cat("    -> factor:", f, "| ref:", levels(xdata[[f]])[1],
          "| coded:", levels(xdata[[f]])[2], "| dummy mean:", round(mean(dummy), 3), "\n")
    }
  }
  cat("[6] f_names:\n"); print(f_names)

  raw_names   <- c(xvar_prescale, if (!isFALSE(factor_var)) factor_var,
                    if (!isFALSE(slice_final)) slice_prescale,
                    sub("^z\\.", "", scaled), unscaled, factor_predictors, wrapped_nuisance_vars)
  final_names <- c(xvar,         if (!isFALSE(factor_var)) factor_var,
                    if (!isFALSE(slice_final)) slice_final,
                    scaled,               z_unscaled, f_names,             wrapped_nuisance_vars)
  rename_map <- c(setNames(final_names, raw_names), setNames(final_names, final_names))
  cat("\n[7] rename_map:\n"); print(rename_map)

  interaction_terms <- term.labels[is_interaction]
  translate_term <- function(term) {
    vars_in_term <- all.vars(str2lang(term))
    for (v in vars_in_term) {
      if (v %in% names(rename_map)) term <- gsub(v, rename_map[[v]], term, fixed = TRUE)
    }
    term
  }
  interaction_terms_final <- if (length(interaction_terms) > 0)
    sapply(interaction_terms, translate_term, USE.NAMES = FALSE) else character(0)
  cat("\n[8] interaction_terms_final (translated, added to refit formula):\n")
  print(interaction_terms_final)
  if (!isFALSE(slice_final) && length(interaction_terms_final) > 0 &&
      !any(grepl(slice_final, interaction_terms_final, fixed = TRUE)))
    cat("[8] NOTE: slice_var does not appear in any detected interaction -- lines may end up parallel.\n")

  random_terms_final <- character(0)
  if (want_mixed) {
    if (!isFALSE(random_formula)) {
      random_terms_final <- random_formula
    } else if (!isFALSE(random_var)) {
      random_terms_final <- if (!isFALSE(random_slope))
        paste0("(1 + ", random_slope, " | ", random_var, ")")
      else paste0("(1 | ", random_var, ")")
    } else if (is_mixed_model) {
      bars <- lme4::findbars(formula(model))
      random_terms_final <- sapply(bars, function(b) paste0("(", paste(deparse(b), collapse = ""), ")"))
    }
    cat("\n[8b] random_terms_final (GLMM/LMM/glmmTMB random-effect structure):\n"); print(random_terms_final)
  }

  rhs <- c(xvar_terms_final, if (!isFALSE(factor_var)) factor_var, if (!isFALSE(slice_final)) slice_final,
           z_names, f_names, wrapped_nuisance_terms, interaction_terms_final, offset_terms, random_terms_final)
  new_formula <- as.formula(paste(yvar, "~", paste(unique(rhs), collapse = " + ")))
  cat("\n[9] Final formula:\n"); print(new_formula)

  family_str   <- if (is.character(family)) tolower(family) else NA
  is_nb_family <- !is.na(family_str) && family_str %in% c("nb", "negbin", "negative.binomial", "negative binomial")
  use_nb_fit   <- if (identical(fit_fun_lc, "auto")) is_nb_family
                  else if (want_glmm) identical(fit_fun_lc, "glmer.nb")
                  else identical(fit_fun_lc, "glm.nb")
  cat("[9] is_nb_family:", is_nb_family, "| use_nb_fit:", use_nb_fit,
      "| want_tmb:", want_tmb, "| want_glmm:", want_glmm, "\n")

  resolve_tmb_family <- function(fam) {
    if (!is.character(fam)) return(fam)
    switch(tolower(fam),
      "beta"     = glmmTMB::beta_family(),
      "nbinom2"  = ,
      "nb"       = glmmTMB::nbinom2(),
      "nbinom1"  = glmmTMB::nbinom1(),
      "tweedie"  = glmmTMB::tweedie(),
      "gaussian" = gaussian(),
      "binomial" = binomial(),
      "poisson"  = poisson(),
      "gamma"    = Gamma(link = "log"),
      stop("plot_glmm: unrecognized family string '", fam, "' for glmmTMB. ",
           "Pass a family object/function instead (e.g. glmmTMB::beta_family()).")
    )
  }

  if (want_tmb) {
    fam_resolved <- resolve_tmb_family(family)
    plot.res <- if (!is.null(glmmTMB_control))
      glmmTMB::glmmTMB(new_formula, data = xdata, family = fam_resolved, control = glmmTMB_control)
    else glmmTMB::glmmTMB(new_formula, data = xdata, family = fam_resolved)
    cat("[9] Model fit via glmmTMB::glmmTMB(), family:", if (is.character(family)) family else "<object>", "\n")

  } else if (want_glmm) {
    if (use_nb_fit) {
      if (!is.null(nb_link))
        warning("plot_glmm: nb_link is ignored for glmer.nb() (link is fixed internally).")
      plot.res <- if (!is.null(glmer_control))
        lme4::glmer.nb(new_formula, data = xdata, control = glmer_control)
      else lme4::glmer.nb(new_formula, data = xdata)
      cat("[9] Model fit via lme4::glmer.nb().\n")

    } else if (is_gaussian_family) {
      use_reml <- if (!is.null(REML)) REML
                  else if (is_lmerMod) lme4::isREML(model)
                  else FALSE
      plot.res <- if (!is.null(lmer_control))
        lme4::lmer(new_formula, data = xdata, REML = use_reml, control = lmer_control)
      else lme4::lmer(new_formula, data = xdata, REML = use_reml)
      cat("[9] Model fit via lme4::lmer(). REML:", use_reml,
          if (is.null(REML) && is_lmerMod) "(auto-detected from template model)"
          else if (is.null(REML)) "(default: FALSE)" else "(user-specified)", "\n")

    } else {
      plot.res <- if (!is.null(glmer_control))
        lme4::glmer(new_formula, data = xdata, family = family, control = glmer_control)
      else lme4::glmer(new_formula, data = xdata, family = family)
      cat("[9] Model fit via lme4::glmer(), family:", if (is.character(family)) family else "<object>", "\n")
    }
  } else {
    if (use_nb_fit) {
      if (!requireNamespace("MASS", quietly = TRUE))
        stop("Package 'MASS' is required for glm.nb() but is not installed.")
      plot.res <- if (is.null(nb_link)) MASS::glm.nb(new_formula, data = xdata)
                  else MASS::glm.nb(new_formula, data = xdata, link = nb_link)
      cat("[9] Model fit via MASS::glm.nb(). Estimated theta:", round(plot.res$theta, 4),
          "| SE(theta):", round(plot.res$SE.theta, 4), "\n")
    } else {
      plot.res <- glm(new_formula, data = xdata, family = family)
      cat("[9] Model fit via glm(), family:", if (is.character(family)) family else "<object>", "\n")
    }
  }

  cat("[9] Fixed-effect coefficients:\n")
  print(if (want_tmb) glmmTMB::fixef(plot.res)$cond
        else if (want_glmm) lme4::fixef(plot.res)
        else coef(plot.res))

  if (want_glmm && isTRUE(plot.res@optinfo$conv$opt != 0 ||
                          length(plot.res@optinfo$conv$lme4$messages) > 0)) {
    cat("[9] WARNING: convergence issue detected on refit (lme4). Consider passing glmer_control/lmer_control ",
        "matching the settings used to fit the original model.\n", sep = "")
  }
  if (want_tmb) {
    tmb_conv_ok <- isTRUE(plot.res$sdr$pdHess) &&
      (is.null(plot.res$fit$convergence) || identical(plot.res$fit$convergence, 0L))
    if (!tmb_conv_ok)
      cat("[9] WARNING: convergence issue detected on refit (glmmTMB) -- non-positive-definite Hessian ",
          "or optimizer non-convergence. Consider simplifying the random-effect structure (e.g. || instead ",
          "of |) or passing glmmTMB_control.\n", sep = "")
  }

  fam_obj  <- family(plot.res)
  fam_name <- fam_obj$family
  linkinv  <- if (!is.null(fam_obj$linkinv)) fam_obj$linkinv else make.link(fam_obj$link)$linkinv

  alpha <- 1 - ci_level
  if (want_mixed) {
    use_t <- FALSE
    crit  <- qnorm(1 - alpha / 2)
  } else {
    use_t <- !(fam_name %in% c("binomial", "poisson"))
    crit  <- if (use_t) qt(1 - alpha / 2, df = plot.res$df.residual) else qnorm(1 - alpha / 2)
  }
  cat("[9] fam_name:", fam_name, "| use_t (qt vs qnorm):", use_t, "| crit:", round(crit, 4),
      "| ci_level:", ci_level, "\n")

  grouping <- !isFALSE(factor_var) || !isFALSE(slice_final)

  compute_ci <- function(pred.data) {
    if (want_tmb) {
      if (identical(ci_method, "boot")) {
        cat("    [boot] Running manual parametric bootstrap (glmmTMB) with nboots =", nboots, "...\n")
        if (!is.null(boot_seed)) set.seed(boot_seed)
        t_start <- Sys.time()
        sims <- tryCatch(simulate(plot.res, nsim = nboots), error = function(e) NULL)
        if (is.null(sims)) stop("plot_glmm: simulate() failed for glmmTMB parametric bootstrap.")
        run_one <- function(ysim) {
          dat2 <- xdata
          dat2[[yvar]] <- ysim
          refit <- tryCatch(suppressWarnings(update(plot.res, data = dat2)), error = function(e) NULL)
          if (is.null(refit)) return(rep(NA_real_, nrow(pred.data)))
          tryCatch(predict(refit, newdata = pred.data, type = "response",
                            re.form = NA, allow.new.levels = TRUE),
                    error = function(e) rep(NA_real_, nrow(pred.data)))
        }
        boot_list <- if (identical(boot_parallel, "multicore") && boot_ncpus > 1 &&
                         requireNamespace("parallel", quietly = TRUE))
          parallel::mclapply(sims, run_one, mc.cores = boot_ncpus)
        else lapply(sims, run_one)
        boot_mat <- do.call(rbind, boot_list)
        cat("    [boot] Done in", round(difftime(Sys.time(), t_start, units = "secs"), 1), "sec.",
            "| failed replicates:", sum(apply(boot_mat, 1, function(r) any(is.na(r)))), "/", nboots, "\n")
        qs <- apply(boot_mat, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
        point_est <- predict(plot.res, newdata = pred.data, type = "response",
                              re.form = NA, allow.new.levels = TRUE)
        data.frame(fit = point_est, lwr = qs[1, ], upr = qs[2, ])
      } else {
        pl <- predict(plot.res, newdata = pred.data, type = "link", se.fit = TRUE,
                      re.form = NA, allow.new.levels = TRUE)
        make_ci_df(pl$fit, pl$se.fit, crit, linkinv)
      }
    } else if (want_glmm && identical(ci_method, "boot")) {
      cat("    [boot] Running bootMer with nboots =", nboots, "...\n")
      boot_FUN <- function(fit) {
        predict(fit, newdata = pred.data, type = "response", re.form = NA)
      }
      t_start <- Sys.time()
      boot_out <- lme4::bootMer(plot.res, FUN = boot_FUN, nsim = nboots,
                                 type = "parametric", use.u = FALSE,
                                 parallel = boot_parallel, ncpus = boot_ncpus,
                                 seed = boot_seed)
      cat("    [boot] Done in", round(difftime(Sys.time(), t_start, units = "secs"), 1), "sec.",
          "| failed replicates:", sum(apply(boot_out$t, 1, function(r) any(is.na(r)))), "/", nboots, "\n")
      qs <- apply(boot_out$t, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
      data.frame(fit = boot_out$t0, lwr = qs[1, ], upr = qs[2, ])
    } else if (want_glmm) {
      fe_terms <- delete.response(terms(lme4::nobars(formula(plot.res))))
      mm <- model.matrix(fe_terms, data = pred.data)
      beta <- lme4::fixef(plot.res)
      mm <- mm[, names(beta), drop = FALSE]
      fit_link <- as.vector(mm %*% beta)
      if (length(offset_exprs) > 0) {
        offset_contrib <- rowSums(sapply(offset_exprs, function(e) eval(e, envir = pred.data)))
        fit_link <- fit_link + offset_contrib
      }
      vc <- as.matrix(vcov(plot.res))
      se_link <- sqrt(rowSums((mm %*% vc) * mm))
      make_ci_df(fit_link, se_link, crit, linkinv)
    } else {
      pl <- predict(plot.res, newdata = pred.data, type = "link", se.fit = TRUE)
      make_ci_df(pl$fit, pl$se.fit, crit, linkinv)
    }
  }

  # ============================================================
  # BRANCH A: xvar is a FACTOR -- single or crossed categorical design
  # ============================================================
  if (xvar_is_factor) {

    group_cols <- c(xvar, if (crossed) factor_var)
    x_levels <- levels(xdata[[xvar]])
    if (crossed) {
      f_levels <- levels(xdata[[factor_var]])
      cells <- expand.grid(setNames(list(x_levels, f_levels), group_cols),
                            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    } else {
      cells <- setNames(data.frame(x_levels, stringsAsFactors = FALSE), xvar)
    }
    cells[[xvar]] <- factor(cells[[xvar]], levels = x_levels)
    if (crossed) cells[[factor_var]] <- factor(cells[[factor_var]], levels = f_levels)
    n_cells <- nrow(cells)
    cat("\n[10] cells (one row per combination):\n"); print(cells)

    pred.data <- cells
    for (z in z_names) pred.data[[z]] <- 0
    for (f in f_names) pred.data[[f]] <- 0
    pred.data <- apply_wrapped_nuisance(pred.data, xdata, wrapped_nuisance_vars)
    for (v in offset_vars) pred.data[[v]] <- offset_value
    cat("\n[11] pred.data:\n"); print(pred.data)

    ci.plot <- compute_ci(pred.data)
    cat("\n[12] ci.plot (back-transformed):\n"); print(ci.plot)

    if (isTRUE(barplot_style)) {
      if (!isFALSE(color_var) || !isFALSE(size_var))
        warning("plot_glmm: color_var/size_var are ignored in barplot_style (no individual points plotted).")

      agg <- aggregate(x = y_plot, by = xdata[group_cols], FUN = function(v) mean(v, na.rm = TRUE))
      cells_key <- do.call(paste, c(cells[group_cols], sep = "___"))
      agg_key   <- do.call(paste, c(agg[group_cols], sep = "___"))
      bar_heights <- agg$x[match(cells_key, agg_key)]
      cat("\n[13] bar_heights (aggregated raw rates, matched to cells order):\n"); print(bar_heights)

      # FIX: 0 is always included in the range so barplot()'s zero-baseline bars
      # (which it always draws from y=0 upward, regardless of ylim) are never clipped.
      default_ylim <- if (fam_name == "binomial") c(0, 1) else
        range(c(0, bar_heights, ci.plot$lwr, ci.plot$upr), na.rm = TRUE)

      set_plot_mar(c(if (crossed) 4.5 else 4, if (crossed) 6 else 4, 1.5, right_mar), lwd = 2)

      x.at <- barplot(bar_heights, col = "white", ylim = default_ylim,
                       ylab = ylabel, xlab = "",
                       yaxt = if (yvar_log) "n" else "s")
      x.at <- x.at[, 1]
      cat("[14] bar midpoint x-positions:\n"); print(x.at)

      mtext(text = as.character(cells[[xvar]]), side = 1, line = 0.2, at = x.at, cex = 0.9)
      if (crossed) {
        f_at <- sapply(f_levels, function(lv) mean(x.at[cells[[factor_var]] == lv]))
        mtext(text = f_levels, side = 1, line = 1.6, at = f_at, cex = 0.9, font = 2)
        mtext(text = xlabel, side = 1, line = 0.2, at = par("usr")[1], adj = 1.15,
              cex = 0.9, font = 2, xpd = NA)
        mtext(text = factor_label, side = 1, line = 3.0, cex = 1.2, font = 2)
      } else {
        mtext(text = xlabel, side = 1, line = 1.6, cex = 1.2, font = 2)
      }

      hll <- if (n_cells > 1) mean(diff(x.at)) / 4 else 0.3
      segments(x0 = x.at - hll, x1 = x.at + hll, y0 = ci.plot$fit, y1 = ci.plot$fit, lwd = 3)
      arrows(x0 = x.at, x1 = x.at, y0 = ci.plot$lwr, y1 = ci.plot$upr,
             code = 3, length = 0.1, angle = 90, lwd = 2)
      cat("[15] Fitted-value segments and CI arrows drawn.\n")

      if (yvar_log) axis(side = 2, at = log(compute_log_ticks(y_plot, log_ybase, "16")), las = 1,
                          labels = compute_log_ticks(y_plot, log_ybase))

    } else {
      xpos <- seq_len(n_cells)
      pa <- compute_point_aesthetics(xdata, color_var, size_var, size_factor)

      set_plot_mar(c(4, 4, 0.5, right_mar))

      plot(NA, xlim = c(0.5, n_cells + 0.5),
           ylim = range(c(y_plot, ci.plot$lwr, ci.plot$upr)),
           xlab = xlabel, ylab = ylabel, bty = "l", xaxt = "n",
           yaxt = if (yvar_log) "n" else "s")
      axis(side = 1, at = xpos, labels = x_levels, las = 1)

      arrows(x0 = xpos, y0 = ci.plot$lwr, x1 = xpos, y1 = ci.plot$upr,
             angle = 90, code = 3, length = 0.05, col = "steelblue4", lwd = 2)
      points(xpos, ci.plot$fit, pch = 19, col = "steelblue4", cex = 1.3)

      xjit <- jitter(as.numeric(xdata[[xvar]]), amount = 0.15)
      draw_points(xjit, y_plot, "steelblue4", point_style, pa$point_cols, pa$point_cex)

      draw_aes_legends(pa, color_var, size_var, extra_var_label, size_factor, draw_legend)

      if (yvar_log) axis(side = 2, at = log(compute_log_ticks(y_plot, log_ybase)), las = 1,
                          labels = compute_log_ticks(y_plot, log_ybase))
    }

    cat("\n========== plot_glmm END ==========\n\n")
    return(invisible(plot.res))
  }

  # ============================================================
  # BRANCH B: xvar is CONTINUOUS -- line/ribbon logic
  # ============================================================
  x.seq <- seq(min(xdata[[xvar]], na.rm = TRUE), max(xdata[[xvar]], na.rm = TRUE), length.out = 100)
  cat("\n[10] x.seq range:", range(x.seq), "\n")

  if (!isFALSE(factor_var)) {
    group_vals <- levels(xdata[[factor_var]]); n_levels <- length(group_vals)
    cat("\n[11] Grouping by factor_var. Levels:", group_vals, "\n")
  } else if (!isFALSE(slice_final)) {
    group_vals <- c(-1, 0, 1); n_levels <- 3
    cat("\n[11] Grouping by slice_var. Values (SD units):", group_vals, "\n")
  } else {
    cat("\n[11] No grouping (single line).\n")
  }

  if (!grouping) {
    pred.data <- data.frame(x = x.seq); names(pred.data) <- xvar
    for (z in setdiff(z_names, xvar)) pred.data[[z]] <- 0
    for (f in f_names) pred.data[[f]] <- 0
  } else {
    pred.data <- data.frame(x = rep(x.seq, times = n_levels)); names(pred.data) <- xvar
    for (z in setdiff(z_names, c(xvar, slice_final))) pred.data[[z]] <- 0
    for (f in f_names) pred.data[[f]] <- 0
    if (!isFALSE(factor_var)) {
      pred.data[[factor_var]] <- factor(rep(group_vals, each = length(x.seq)), levels = group_vals)
    } else {
      pred.data[[slice_final]] <- rep(group_vals, each = length(x.seq))
    }
  }
  pred.data <- apply_wrapped_nuisance(pred.data, xdata, wrapped_nuisance_vars)
  for (v in offset_vars) pred.data[[v]] <- offset_value
  cat("\n[12] pred.data (head):\n"); print(head(pred.data))
  cat("[12] pred.data xvar column range (should match x.seq, NOT be constant):\n")
  print(range(pred.data[[xvar]]))

  if (identical(ci_method, "boot") && nrow(pred.data) > 200)
    cat("\n[13] NOTE: ci_method = 'boot' with", nrow(pred.data),
        "prediction points -- this may be slow (", nboots, "refits x", nrow(pred.data), "predictions each).\n")

  ci.plot <- compute_ci(pred.data)
  cat("\n[13] ci.plot (back-transformed) head:\n"); print(head(ci.plot))
  if (grouping) {
    cat("[13] fit range per group (should NOT be a single repeated value):\n")
    group_col <- if (!isFALSE(factor_var)) factor_var else slice_final
    for (v in group_vals) {
      sel <- pred.data[[group_col]] == v
      cat("  group =", v, "| fit range:", range(ci.plot[sel, "fit"]), "\n")
    }
  } else {
    cat("[13] fit range (no grouping):", range(ci.plot[, "fit"]), "\n")
  }

  pa <- compute_point_aesthetics(xdata, color_var, size_var, size_factor)

  set_plot_mar(c(4, 4, 0.5, right_mar))

  plot(NA, xlim = range(xdata[[xvar]], na.rm = TRUE),
       ylim = range(c(y_plot, ci.plot[, "lwr"], ci.plot[, "upr"])),
       xlab = xlabel, ylab = ylabel, bty = "l",
       xaxt = if (xvar_scaled || xvar_log) "n" else "s",
       yaxt = if (yvar_log) "n" else "s")
  cat("\n[15] Empty plot created.\n")

  if (!grouping) {
    cat("[16] Drawing single CI ribbon + line + points.\n")
    polygon(x = c(pred.data[[xvar]], rev(pred.data[[xvar]])),
            y = c(ci.plot[, "lwr"], rev(ci.plot[, "upr"])),
            col = adjustcolor("steelblue", alpha.f = 0.2), border = NA)
    lines(pred.data[[xvar]], ci.plot[, "fit"], col = "steelblue4", lwd = 2)
    draw_points(xdata[[xvar]], y_plot, "steelblue4", point_style, pa$point_cols, pa$point_cex)

  } else if (!isFALSE(factor_var)) {
    pchs <- 1:n_levels; ltys <- 2:(n_levels + 1)
    cols <- c("steelblue4", "firebrick4", "forestgreen", "darkorange3", "purple4", "goldenrod3")[1:n_levels]
    cat("[16] Drawing factor-grouped CI ribbons + lines + points.\n")
    for (i in seq_along(group_vals)) {
      lev <- group_vals[i]; sel <- pred.data[[factor_var]] == lev
      polygon(x = c(x.seq, rev(x.seq)), y = c(ci.plot[sel, "lwr"], rev(ci.plot[sel, "upr"])),
              col = adjustcolor(cols[i], alpha.f = 0.2), border = NA)
      lines(x.seq, ci.plot[sel, "fit"], lwd = 2, lty = ltys[i], col = cols[i])
      sel_data <- xdata[[factor_var]] == lev
      draw_points(xdata[[xvar]], y_plot, cols[i], point_style, pa$point_cols, pa$point_cex,
                  sel = sel_data, pch_shape = pchs[i])
    }
    draw_legend(legend = group_vals, pch = pchs, lty = ltys, col = cols, bty = "n",
                title = factor_label, title.font = 2)

  } else {
    ltys <- c(2, 1, 4); cols <- c("firebrick4", "steelblue4", "forestgreen")
    labels <- c("-1 SD", "Mean", "+1 SD")
    cat("[16] Drawing slice-grouped CI ribbons + lines (single point cloud).\n")
    for (i in seq_along(group_vals)) {
      val <- group_vals[i]; sel <- pred.data[[slice_final]] == val
      polygon(x = c(x.seq, rev(x.seq)), y = c(ci.plot[sel, "lwr"], rev(ci.plot[sel, "upr"])),
              col = adjustcolor(cols[i], alpha.f = 0.2), border = NA)
      lines(x.seq, ci.plot[sel, "fit"], lwd = 2, lty = ltys[i], col = cols[i])
    }
    draw_points(xdata[[xvar]], y_plot, cols[2], point_style, pa$point_cols, pa$point_cex)
    draw_legend(legend = labels, lty = ltys, col = cols, bty = "n",
                title = slice_label, title.font = 2)
  }

  draw_aes_legends(pa, color_var, size_var, extra_var_label, size_factor, draw_legend)

  if (xvar_scaled || xvar_log) {
    pre_vals <- xdata[[xvar_prescale]]; mean_val <- mean(pre_vals, na.rm = TRUE); sd_val <- sd(pre_vals, na.rm = TRUE)
    raw_vals <- xdata[[xvar_raw]]
    if (xvar_log) {
      lo <- floor(log10(min(raw_vals, na.rm = TRUE))); hi <- ceiling(log10(max(raw_vals, na.rm = TRUE)))
      candidates <- as.vector(outer(c(1, 2, 5), 10^(lo:hi)))
      nice_raw <- sort(candidates[candidates >= min(raw_vals, na.rm = TRUE) & candidates <= max(raw_vals, na.rm = TRUE)])
    } else nice_raw <- pretty(raw_vals)
    transformed_nice <- if (xvar_log) log(nice_raw) else nice_raw
    x.at <- if (xvar_scaled) (transformed_nice - mean_val) / sd_val else transformed_nice
    cat("\n[18] x-axis ticks:\n"); print(nice_raw)
    axis(side = 1, at = x.at, labels = nice_raw, las = 1)
  } else cat("\n[18] No x-axis back-transform needed.\n")

  if (yvar_log) {
    axis(side = 2, at = log(compute_log_ticks(y_plot, log_ybase, "19")), las = 1,
         labels = compute_log_ticks(y_plot, log_ybase))
  } else cat("\n[19] No y-axis back-transform needed.\n")

  cat("\n========== plot_glmm END ==========\n\n")
  invisible(plot.res)
}
