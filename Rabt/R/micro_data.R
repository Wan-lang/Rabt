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

micro_data <- function(env_data, otu_data, datatype = c('otu/asv', 'matrix'), 
                       method = NULL, include_pairs = TRUE) 
{
  # 匹配 datatype 参数
  datatype <- match.arg(datatype)
  
  # 输入验证
  if (!is.matrix(env_data) && !is.data.frame(env_data)) {
    stop("env_data must be a matrix or data.frame")
  }
  if (!is.matrix(otu_data) && !is.data.frame(otu_data)) {
    stop("otu_data must be a matrix or data.frame")
  }
  
  # 确保 env_data 是数据框
  env_data <- as.data.frame(env_data)
  otu_data <- as.data.frame(otu_data)
  
  # 检查样本数量是否匹配
  if (nrow(env_data) != ncol(otu_data)) {
    stop(sprintf("Sample count mismatch: env_data has %d rows, otu_data has %d columns. 
                  For OTU data, samples should be columns.", 
                 nrow(env_data), ncol(otu_data)))
  }
  
  # 检查样本名是否匹配
  if (!is.null(rownames(env_data)) && !is.null(colnames(otu_data))) {
    if (!all(rownames(env_data) == colnames(otu_data))) {
      warning("Sample names in env_data (rows) and otu_data (columns) do not match. 
               Ensure samples are aligned.")
    }
  }
  
  # 处理 OTU/ASV 数据
  if (datatype == 'otu/asv') {
    if (is.null(method)) {
      stop("Please provide a distance matrix method for 'otu/asv' datatype.")
    }
    if (!is.character(method) || length(method) != 1) {
      stop("method must be a single character string")
    }
    
    # 转置 OTU 数据：样本为行，OTU 为列
    otu_t <- t(otu_data)
    
    # 计算距离矩阵
    otu_dist <- tryCatch({
      vegan::vegdist(otu_t, method = method)
    }, error = function(e) {
      stop("Error computing distance matrix: ", e$message, 
           "\nSupported methods include: manhattan, euclidean, canberra, bray, 
            kulczynski, jaccard, gower, altGower, morisita, horn, mountford, 
            raup, binomial, chao, cao, mahalanobis, and more. 
            See ?vegan::vegdist for full list.")
    })
    
    otu_matrix <- as.matrix(otu_dist)
    
  } else {  # datatype == 'matrix'
    # 假设 otu_data 已经是距离矩阵
    if (nrow(otu_data) != ncol(otu_data)) {
      stop("For datatype = 'matrix', otu_data must be a square distance matrix")
    }
    if (!all(rownames(otu_data) == colnames(otu_data))) {
      warning("Row and column names of otu_data distance matrix do not match")
    }
    otu_matrix <- as.matrix(otu_data)
  }
  
  # 获取样本名
  sample_names <- rownames(otu_matrix)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(nrow(otu_matrix)))
  }
  
  # 将距离矩阵转为长格式（使用基础 R，避免 reshape2 依赖）
  n <- nrow(otu_matrix)
  pair_idx <- expand.grid(i = seq_len(n), j = seq_len(n))
  
  otu_long <- data.frame(
    Sample_i = sample_names[pair_idx$i],
    Sample_j = sample_names[pair_idx$j],
    value = otu_matrix[cbind(pair_idx$i, pair_idx$j)],
    stringsAsFactors = FALSE
  )
  
  # 计算环境变量差异矩阵（向量化，避免循环）
  n_env <- nrow(env_data)
  p_env <- ncol(env_data)
  env_names <- colnames(env_data)
  
  # 使用 rep() 和矩阵运算替代 outer() 循环
  env_i <- env_data[rep(seq_len(n_env), each = n_env), , drop = FALSE]
  env_j <- env_data[rep(seq_len(n_env), times = n_env), , drop = FALSE]
  env_diff <- as.data.frame(env_i - env_j)
  colnames(env_diff) <- env_names
  
  # 组合结果
  if (include_pairs) {
    ABT_data <- cbind(
      Sample_i = otu_long$Sample_i,
      Sample_j = otu_long$Sample_j,
      env_diff,
      otu_value = otu_long$value
    )
  } else {
    ABT_data <- cbind(env_diff, otu_value = otu_long$value)
  }
  
  # 添加属性信息
  attr(ABT_data, "datatype") <- datatype
  attr(ABT_data, "method") <- method
  attr(ABT_data, "n_samples") <- n_env
  attr(ABT_data, "n_pairs") <- nrow(ABT_data)
  
  return(ABT_data)
}




