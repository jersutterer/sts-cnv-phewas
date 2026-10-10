################################################################################
# functions.R
# Helper functions shared by the PheWAS, LabWAS, and sex hormone scripts.
################################################################################

# Convert p-values exported as text (e.g. "<2e-16") to numeric.
# Values reported as "<x" are set to x.
clean_pvalue <- function(p) {
  as.numeric(gsub("<", "", as.character(p)))
}

# Read one sheet per sex x CNV stratum into a named list of data frames.
read_strata <- function(file, analysis) {
  set_names(groups) %>%
    map(~ read_excel(file, sheet = sheet_name(.x, analysis)) %>%
          mutate(group = group_labels[[.x]]))
}

# Apply multiple-testing correction within one stratum.
# Bonferroni threshold = alpha / number of tests passing QC in that stratum;
# FDR by Benjamini-Hochberg.
add_multiple_testing <- function(df, alpha = params$alpha, q_threshold = params$fdr_q) {
  bonf <- alpha / nrow(df)
  df %>%
    arrange(p) %>%
    mutate(
      fdr_q          = p.adjust(p, method = "BH"),
      sig_bonferroni = p < bonf,
      sig_fdr        = fdr_q < q_threshold,
      row_id         = row_number()
    )
}

# Bonferroni threshold for a stratum (every row of a QC-filtered stratum is a test)
bonf_threshold <- function(df, alpha = params$alpha) alpha / nrow(df)

# Print a one-line summary per stratum (counts reported in Results / Table S1).
report_strata <- function(data_list, unit = "tests") {
  imap(data_list, function(df, nm) {
    cat(sprintf("%-20s %5d %s | Bonferroni p < %.2e: %3d | FDR q < 0.05: %3d\n",
                nm, nrow(df), unit, bonf_threshold(df),
                sum(df$sig_bonferroni, na.rm = TRUE),
                sum(df$sig_fdr, na.rm = TRUE)))
  })
  invisible(NULL)
}

# Manhattan plot used for Figures 1, 2, S3 and S4.
#   category_col : column holding the phenotype / laboratory category
#   effect_col   : column whose sign gives the effect direction (log OR or beta)
#   label_data   : rows to label (already containing a `label_text` column)
#   shape_labels : legend text for upward / downward triangles
#   title        : NULL for journal figures (eBioMedicine: no titles in artwork)
manhattan_plot <- function(data, category_col, effect_col, label_data = NULL,
                           shape_labels, x_title, title = NULL, caption = NULL) {

  plot_data <- data %>%
    filter(!is.na(p), is.finite(p), p > 0) %>%
    mutate(
      log10p   = -log10(p),
      category = tools::toTitleCase(tolower(as.character(.data[[category_col]]))),
      effect   = .data[[effect_col]]
    )

  bonf <- bonf_threshold(data)
  fdr_line <- if (any(plot_data$sig_fdr, na.rm = TRUE)) {
    max(plot_data$p[plot_data$sig_fdr], na.rm = TRUE)
  } else NA

  # Order categories so those with Bonferroni-significant hits are spread
  # across the x-axis rather than clustered (improves label legibility).
  cat_sig <- plot_data %>%
    group_by(category) %>%
    summarise(has_sig = any(sig_bonferroni, na.rm = TRUE), .groups = "drop")
  sig_cats    <- cat_sig$category[cat_sig$has_sig]
  nonsig_cats <- cat_sig$category[!cat_sig$has_sig]

  if (length(sig_cats) > 0 && length(nonsig_cats) > 0) {
    ratio <- length(nonsig_cats) / length(sig_cats)
    cat_order <- character(0)
    for (i in seq_along(sig_cats)) {
      cat_order <- c(cat_order, sig_cats[i])
      from <- ceiling((i - 1) * ratio) + 1
      to   <- min(ceiling(i * ratio), length(nonsig_cats))
      if (from <= to) cat_order <- c(cat_order, nonsig_cats[from:to])
    }
    cat_order <- c(cat_order, setdiff(nonsig_cats, cat_order))
  } else {
    cat_order <- c(sig_cats, nonsig_cats)
  }
  plot_data$category_factor <- factor(plot_data$category, levels = cat_order)

  # Random order within category (fixed seed for reproducibility)
  set.seed(123)
  plot_data <- plot_data %>%
    group_by(category_factor) %>%
    mutate(random_order = sample(row_number())) %>%
    arrange(category_factor, random_order) %>%
    ungroup()

  # x positions with a gap between categories
  cat_info <- plot_data %>%
    count(category_factor, name = "n_points") %>%
    mutate(
      spacing         = ceiling(n_points * 0.2),
      cum_points      = cumsum(n_points + spacing),
      category_start  = c(1, head(cum_points, -1) + 1),
      category_center = category_start + (n_points - 1) / 2
    )

  plot_data <- plot_data %>%
    left_join(cat_info, by = "category_factor") %>%
    group_by(category_factor) %>%
    mutate(x_pos = category_start + row_number() - 1) %>%
    ungroup() %>%
    mutate(
      color_group = rep(c("blue", "gray"), length.out = nrow(cat_info))[as.integer(category_factor)],
      shape_type  = ifelse(effect > 0, "up", "down")
    )

  p <- ggplot(plot_data, aes(x = x_pos, y = log10p)) +
    geom_point(aes(color = color_group, fill = color_group, shape = shape_type),
               size = 2.5, alpha = 0.7) +
    scale_color_manual(values = c(blue = "#4472C4", gray = "#7F7F7F")) +
    scale_fill_manual(values  = c(blue = "#4472C4", gray = "#7F7F7F")) +
    scale_shape_manual(values = c(up = 24, down = 25), labels = shape_labels,
                       name = "Effect Direction") +
    geom_hline(yintercept = -log10(bonf), linetype = "solid",
               color = "red", linewidth = 0.8) +
    scale_x_continuous(breaks = cat_info$category_center, labels = cat_order,
                       expand = c(0.01, 0.01)) +
    labs(title = title, caption = caption, x = x_title,
         y = expression(-log[10](p))) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
      axis.text.y  = element_text(size = 12),
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text  = element_text(size = 10),
      panel.grid   = element_blank(),
      plot.margin  = margin(10, 20, 30, 20, unit = "pt")
    ) +
    guides(color = "none", fill = "none")

  if (!is.na(fdr_line)) {
    p <- p + geom_hline(yintercept = -log10(fdr_line), linetype = "dashed",
                        color = "blue", linewidth = 0.6)
  }

  if (!is.null(label_data) && nrow(label_data) > 0) {
    label_data <- semi_join(plot_data, label_data, by = "row_id") %>%
      left_join(select(label_data, row_id, label_text), by = "row_id")
    p <- p + geom_text_repel(data = label_data, aes(label = label_text),
                             size = 3.5, max.overlaps = 40, segment.size = 0.3,
                             segment.color = "grey50", show.legend = FALSE)
  }
  p
}

# Save a figure as 300-dpi PNG plus vector PDF.
save_figure <- function(plot, name, width, height) {
  ggsave(out_path(paste0(name, ".png")), plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(out_path(paste0(name, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf)
  cat("Saved:", name, "(.png, .pdf)\n")
}

# Number formatting used in tables
fmt_p <- function(p, digits = 3) {
  ifelse(p < 0.001, sprintf("%.2e", p), formatC(p, format = "f", digits = digits))
}
