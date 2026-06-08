#' Plot Variable Relative Importance from ABT Model
#'
#' Generates a publication-ready horizontal bar plot of variable relative importance
#' derived from an Adaptive Boosting Tree (ABT) model. The plot is sorted by influence
#' (most important variables at the top), supports custom axis labels, and automatically
#' exports to both raster (PNG) and vector (PDF) formats.
#'
#' @param fit_softcorals_gbm An ABT model object returned by \code{\link{ABT}}. The
#'   function extracts the \code{var.importance} component, which must be a
#'   \code{data.frame} with columns \code{rel.inf} and \code{var}.
#' @param dict_label Optional named \code{list} or \code{environment} providing
#'   human-readable labels for variables. Names must match the values in the
#'   \code{var} column of the importance table. Variables not found in
#'   \code{dict_label} retain their original names with a warning. Default is
#'   \code{NULL} (no relabeling).
#' @param out_plot Character string giving the output directory path. Default is the
#'   current working directory (\code{"."}). Non-existent directories are created
#'   recursively via \code{\link{dir.create}}.
#' @param plot_title Character. Plot title displayed at the top center. Default is
#'   \code{"Relative Influence of Variables"}.
#' @param x_label Character. X-axis label describing the importance metric. Default
#'   is \code{"Relative Influence"}.
#' @param y_label Character. Y-axis label describing the variables shown. Default is
#'   \code{"Variable"}.
#' @param width Numeric. Plot width in inches. Default is \code{10}.
#' @param height Numeric. Plot height in inches. Default is \code{8}.
#' @param dpi Numeric. Resolution of the PNG output in dots per inch. Default is
#'   \code{300}, suitable for print-quality figures. Higher values (e.g., 600) are
#'   recommended for journal submissions.
#'
#' @return Invisibly returns the \code{ggplot} object, allowing further post-hoc
#'   customization (e.g., adding vertical reference lines or annotations). The
#'   plot is also printed to the active graphics device and saved to disk as:
#'   \itemize{
#'     \item \code{ABT.png} — raster image (resolution controlled by \code{dpi}).
#'     \item \code{ABT.pdf} — vector graphic (resolution-independent).
#'   }
#'
#' @details
#' \strong{Plot layout}
#' The function produces a horizontal bar chart with the following styling choices:
#' \itemize{
#'   \item Variables are ordered by decreasing \code{rel.inf} so that the most
#'     influential predictor appears at the top of the Y axis.
#'   \item Each bar is filled with a discrete color mapped to the variable name;
#'     the legend is suppressed to reduce clutter.
#'   \item The classic theme (\code{theme_classic}) is used with bold axis text
#'     and titles for readability in presentations and manuscripts.
#'   \item A thin black border is drawn around the plotting panel.
#' }
#'
#' \strong{Automatic axis scaling}
#' The X-axis breaks are computed heuristically to yield approximately six evenly
#' spaced, human-friendly ticks (multiples of 1, 2, 5, or 10 in the appropriate
#' order of magnitude). If all relative influences are zero or negative, the axis
#' is pinned to \code{[0, 1]}.
#'
#' \strong{Label mapping}
#' When \code{dict_label} is supplied, the function uses \code{scale_y_discrete}
#' with a named vector (breaks = original names, labels = custom names). This
#' preserves the internal factor ordering while displaying friendly text, and
#' avoids the pitfalls of direct string replacement in the data frame.
#'
#' @section File output:
#' Both files are written silently (no confirmation message). If the directory is
#' not writable, \code{ggsave} will raise an error. Existing files with the same
#' name are overwritten without warning.
#'
#' @examples
#' # Example usage:
#' data(fit_softcorals_gbm)
#' abt_plot(fit_softcorals_gbm)		
#'
#' @export
#' 

abt_plot <- function(fit_softcorals_gbm, dict_label = NULL, out_plot = '.', 
                     plot_title = "Relative Influence of Variables", 
                     x_label = "Relative Influence", y_label = "Variable",
                     width = 10, height = 8, dpi = 300) 
{
  abt_summary <- fit_softcorals_gbm$var.importance
  
  # 检查必需列
  if (!is.data.frame(abt_summary) || !all(c("rel.inf", "var") %in% colnames(abt_summary))) {
    stop("var.importance must be a data.frame with columns 'rel.inf' and 'var'")
  }
  
  # 检查 dict_label 类型（支持 list 和 environment）
  if (!is.null(dict_label) && !is.list(dict_label) && !is.environment(dict_label)) {
    stop("dict_label must be a list, environment, or NULL")
  }
  
  # 检查输出目录
  if (!dir.exists(out_plot)) {
    dir.create(out_plot, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(out_plot)) {
    stop("out_plot must be a valid directory path")
  }
  
  # 按相对影响排序（降序，让最重要的在上方）
  abt_summary <- abt_summary[order(abt_summary$rel.inf, decreasing = TRUE), ]
  
  # 将 var 转为因子，固定顺序
  abt_summary$var <- factor(abt_summary$var, levels = rev(as.character(abt_summary$var)))
  
  # 获取变量名
  var_names <- as.character(abt_summary$var)
  
  # 应用自定义标签
  if (!is.null(dict_label)) {
    custom_labels <- sapply(var_names, function(v) {
      if (exists(v, envir = if (is.environment(dict_label)) dict_label else as.environment(dict_label))) {
        dict_label[[v]]
      } else {
        warning("Label '", v, "' not found in dict_label. Using original label.")
        v
      }
    }, USE.NAMES = FALSE)
  } else {
    custom_labels <- var_names
  }
  
  # 构建标签映射（确保 breaks 和 labels 一一对应）
  label_map <- setNames(custom_labels, var_names)
  
  # 智能计算 x 轴刻度间隔
  max_rel_inf <- max(abt_summary$rel.inf, 0)
  if (max_rel_inf <= 0) {
    x_breaks <- c(0)
    x_limits <- c(0, 1)
  } else {
    # 自动选择合适的刻度间隔
    rough_interval <- max_rel_inf / 6
    magnitude <- 10^floor(log10(rough_interval))
    normalized <- rough_interval / magnitude
    if (normalized <= 1) {
      interval <- magnitude
    } else if (normalized <= 2) {
      interval <- 2 * magnitude
    } else if (normalized <= 5) {
      interval <- 5 * magnitude
    } else {
      interval <- 10 * magnitude
    }
    x_breaks <- seq(0, ceiling(max_rel_inf / interval) * interval, by = interval)
    x_limits <- c(0, max(x_breaks) * 1.05)  # 5% 边距
  }
  
  # 创建条形图
  p <- ggplot2::ggplot(abt_summary, ggplot2::aes(x = rel.inf, y = var, fill = var)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7, colour = NA) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(color = 'black', fill = NA),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 20), 
      axis.title = ggplot2::element_text(face = "bold", size = 16),
      axis.text = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "bold", size = 12),
      plot.margin = ggplot2::unit(c(1, 2, 1, 1), "lines")
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0), 
      breaks = x_breaks, 
      limits = x_limits
    ) +
    ggplot2::scale_y_discrete(
      labels = label_map
    ) +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = y_label
    )
  
  # 打印图形
  print(p)
  
  # 保存图形
  out_plot_png <- file.path(out_plot, 'ABT.png')
  out_plot_pdf <- file.path(out_plot, 'ABT.pdf')
  
  ggplot2::ggsave(out_plot_png, plot = p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(out_plot_pdf, plot = p, width = width, height = height)
  
  invisible(p)
}