#' Prepare Microbiome Data for ABT Analysis
#'
#' Reads environmental and OTU/ASV data from files, processes them, and prepares
#' a pairwise difference data frame for Adaptive Boosting Tree (ABT) analysis.
#' Supports both pre-computed distance matrices and raw OTU/ASV count tables with
#' on-the-fly distance computation via \code{\link[vegan]{vegdist}}.
#'
#' @param env_data A \code{data.frame} or \code{matrix} of environmental variables.
#'   Rows must correspond to samples (in the same order as columns of \code{otu_data}),
#'   columns to environmental predictors (e.g., temperature, pH, salinity).
#' @param otu_data A \code{data.frame} or \code{matrix} of microbial community data.
#'   When \code{datatype = "otu/asv"}, rows are taxa (OTUs/ASVs) and columns are samples.
#'   When \code{datatype = "matrix"}, a square symmetric distance matrix with rows and
#'   columns both representing samples.
#' @param datatype Character string specifying the input format. One of:
#'   \describe{
#'     \item{"otu/asv"}{(default) Raw abundance table; a distance matrix will be computed
#'       using the method specified in \code{method}.}
#'     \item{"matrix"}{Pre-computed distance matrix; \code{method} is ignored.}
#'   }
#' @param method Character string specifying the ecological distance metric. Required
#'   when \code{datatype = "otu/asv"}. See \code{\link[vegan]{vegdist}} for the full list.
#'   Commonly used methods in microbial ecology include:
#'   \describe{
#'     \item{"bray"}{Bray–Curtis dissimilarity (recommended for abundance data; default
#'       in many pipelines).}
#'     \item{"jaccard"}{Jaccard distance (presence/absence only; ignores abundance).}
#'     \item{"euclidean"}{Euclidean distance (suitable for Hellinger- or CLR-transformed
#'       data).}
#'     \item{"manhattan"}{Manhattan (L1) distance (robust to outliers).}
#'     \item{"canberra"}{Canberra distance (weights rare taxa heavily).}
#'     \item{"horn"}{Horn–Morisita index (accounts for abundance and sample size).}
#'     \item{"gower"}{Gower distance (handles mixed variable types).}
#'   }
#' @param include_pairs Logical. If \code{TRUE} (default), the returned data frame
#'   includes identifier columns \code{Sample_i} and \code{Sample_j} for each pair.
#'   Set to \code{FALSE} to drop them and save memory when sample identity is not
#'   needed downstream.
#'
#' @return A \code{data.frame} with one row per sample pair (N\eqn{^2}{^2} rows for N samples).
#'   Columns contain:
#'   \describe{
#'     \item{Sample_i, Sample_j}{Character identifiers of the two samples forming the pair
#'       (present only if \code{include_pairs = TRUE}).}
#'     \item{env_diff}{Difference of each environmental variable (\code{env_j} - \code{env_i}),
#'       preserving column names from \code{env_data}.}
#'     \item{otu_value}{Dissimilarity/distance value between \code{Sample_i} and \code{Sample_j}
#'       extracted from the computed or supplied distance matrix.}
#'   }
#'   The object carries the following \code{\link{attributes}}:
#'   \itemize{
#'     \item \code{datatype} — input data type ("otu/asv" or "matrix").
#'     \item \code{method} — distance method used (\code{NULL} for \code{datatype = "matrix"}).
#'     \item \code{n_samples} — number of samples (N).
#'     \item \code{n_pairs} — total number of rows (N\eqn{^2}{^2}).
#'   }
#'
#' @details
#' This function transforms a traditional “samples × variables” data structure into a
#' pairwise-difference framework required by ABT models. The key idea is that each row
#' represents the difference between two samples, allowing the model to learn how
#' environmental gradients drive community turnover.
#'
#' \strong{Data alignment requirements}
#' \itemize{
#'   \item \code{nrow(env_data)} must equal \code{ncol(otu_data)} (same number of samples).
#'   \item When \code{datatype = "matrix"}, \code{otu_data} must be square
#'     (\code{nrow == ncol}) and ideally symmetric.
#'   \item If both \code{env_data} has row names and \code{otu_data} has column names,
#'     they are compared; a warning is issued if they do not match, but execution
#'     continues.
#' }
#'
#' \strong{Distance computation (\code{datatype = "otu/asv"})}
#' The OTU/ASV table is internally transposed so that samples become rows and taxa
#' become columns before calling \code{vegan::vegdist}. This matches the conventional
#' input orientation of the vegan package.
#'
#' \strong{Memory note}
#' For N samples, the output has N\eqn{^2}{^2} rows. With 100 samples this is 10,000 rows;
#' with 1,000 samples it is 1,000,000 rows. Consider subsetting or block-wise processing
#' for very large datasets.
#'
#' @examples
#' # Example usage:
#' data(env)
#' data(otu)
#' datam <- mirco_data(env, otu, header = TRUE, datatype = 'otu/asv', method = "bray")
#' @export

mirco_data <- function(env_data = NULL, otu_data = NULL, tree = NULL, 
                       datatype = NULL, method = NULL, include_pairs = TRUE) 
{
  # 1. 匹配 datatype 参数
  datatype <- match.arg(datatype, choices = c("otu/asv", "tre/nwk"))
  
  # 2. 输入验证
  if (is.null(env_data)) {
    stop("env_data must be provided")
  }
  if (!is.matrix(env_data) && !is.data.frame(env_data)) {
    stop("env_data must be a matrix or data.frame")
  }
  
  if (is.null(otu_data)) {
    stop("otu_data must be provided")
  }
  if (!is.matrix(otu_data) && !is.data.frame(otu_data)) {
    stop("otu_data must be a matrix or data.frame")
  }
  
  if (is.null(method) || !method %in% c("bray", "jaccard", "weighted unifrac", "unweighted unifrac")) {
    stop("method must be one of: 'bray', 'jaccard', 'weighted unifrac', 'unweighted unifrac'")
  }
  if (!is.character(method) || length(method) != 1) {
    stop("method must be a single character string")
  }
  
  # 3. 确保数据框格式
  env_data <- as.data.frame(env_data)
  otu_data <- as.data.frame(otu_data)
  
  # 4. 检查样本数量：env 行数 = otu 列数（样本为列）
  if (nrow(env_data) != ncol(otu_data)) {
    stop(sprintf("Sample count mismatch: env_data has %d rows, otu_data has %d columns.", 
                 nrow(env_data), ncol(otu_data)))
  }
  
  # 5. 检查样本名是否匹配
  if (!is.null(rownames(env_data)) && !is.null(colnames(otu_data))) {
    if (!all(rownames(env_data) == colnames(otu_data))) {
      warning("Sample names in env_data (rows) and otu_data (columns) do not match.")
    }
  }
  
  # 6. 转置 OTU 数据：样本为行，OTU 为列
  otu_t <- t(otu_data)
  
  # 7. 计算距离矩阵
  if (datatype == 'otu/asv') {
    if (!method %in% c("bray", "jaccard")) {
      stop("For 'otu/asv' datatype, method must be 'bray' or 'jaccard'")
    }
	if (!is.null(tree)) {
      warning("Phylogenetic tree is not used in the current mode")
    }
    otu_dist <- vegan::vegdist(otu_t, method = method)
    
  } else if (datatype == 'tre/nwk') {
    if (!method %in% c("weighted unifrac", "unweighted unifrac")) {
      stop("For 'tre/nwk' datatype, method must be 'weighted unifrac' or 'unweighted unifrac'")
    }
    if (is.null(tree)) {
      stop("tree file path must be provided for 'tre/nwk' datatype")
    }
    if (!is.character(tree) || length(tree) != 1) {
      stop("tree must be a single character string (file path)")
    }
    if (!file.exists(tree)) {
      stop("tree file does not exist: ", tree)
    }
    
    library(ape)
    tree_data <- read.tree(tree)
    
    # 检查树和OTU表的一致性
    common_otus <- intersect(colnames(otu_t), tree_data$tip.label)
    if (length(common_otus) == 0) {
      stop("No matching OTU names between otu_data and tree")
    }
    if (length(common_otus) < ncol(otu_t)) {
      warning(sprintf("Tree is missing %d OTUs from the table. Using %d common OTUs.", 
                      ncol(otu_t) - length(common_otus), length(common_otus)))
      otu_t <- otu_t[, common_otus, drop = FALSE]
      tree_data <- ape::keep.tip(tree_data, common_otus)
    }
    
    # 计算 UniFrac
    unifrac <- GUniFrac::GUniFrac(otu_t, tree_data)$unifracs
    
    if (method == "weighted unifrac") {
      otu_dist <- as.dist(unifrac[, , 'd_1'])
    } else if (method == "unweighted unifrac") {
      otu_dist <- as.dist(unifrac[, , 'd_UW'])
    }
  }
  
  # 8. 转为矩阵
  otu_matrix <- as.matrix(otu_dist)
  
  # 9. 获取样本名
  sample_names <- rownames(otu_matrix)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(nrow(otu_matrix)))
    rownames(otu_matrix) <- colnames(otu_matrix) <- sample_names
  }
  
  # 10. 距离矩阵转为长格式（只保留 i < j，避免重复）
  n <- nrow(otu_matrix)
  pair_idx <- expand.grid(i = seq_len(n), j = seq_len(n))
  pair_idx <- pair_idx[pair_idx$i < pair_idx$j, , drop = FALSE]  # 只保留上三角
  
  otu_long <- data.frame(
    Sample_i = sample_names[pair_idx$i],
    Sample_j = sample_names[pair_idx$j],
    otu_value = otu_matrix[cbind(pair_idx$i, pair_idx$j)],
    stringsAsFactors = FALSE
  )
  
  # 11. 计算环境变量差异矩阵（只计算 i < j 的对）
  n_pairs <- nrow(pair_idx)
  p_env <- ncol(env_data)
  env_names <- colnames(env_data)
  
  env_i <- env_data[pair_idx$i, , drop = FALSE]
  env_j <- env_data[pair_idx$j, , drop = FALSE]
  env_diff <- as.data.frame(env_i - env_j)
  colnames(env_diff) <- env_names
  
  # 12. 组合结果
  if (include_pairs) {
    ABT_data <- cbind(
      Sample_i = otu_long$Sample_i,
      Sample_j = otu_long$Sample_j,
      env_diff,
      otu_value = otu_long$otu_value
    )
  } else {
    ABT_data <- cbind(env_diff, otu_value = otu_long$otu_value)
  }
  
  # 13. 添加属性
  attr(ABT_data, "datatype") <- datatype
  attr(ABT_data, "method") <- method
  attr(ABT_data, "n_samples") <- n
  attr(ABT_data, "n_pairs") <- n_pairs
  
  return(ABT_data)
}

