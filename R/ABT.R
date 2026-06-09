#' Adaptive Boosting Trees (ABT) Model Fitting
#'
#' Adaptive boosting tree model based on Gradient Boosting Machine (GBM) with
#' cross-validation for optimal tree number selection. Suitable for ecology,
#' microbiome, and environmental-community relationship analysis.
#'
#' @param formula Model formula specifying response and predictors, e.g., y ~ x1 + x2
#' @param data Data frame or list containing model data
#' @param distribution Error distribution type. Options: "gaussian" (regression),
#'   "bernoulli" (binary classification), "poisson" (count data), "adaboost",
#'   "laplace" (robust regression), "coxph" (survival analysis)
#' @param weights Optional observation weights vector for handling imbalanced samples
#'   or differential importance
#' @param offset Optional linear predictor offset term, used for Poisson or Cox models
#' @param var.monotone Vector of monotonicity constraints: -1 (decreasing), 0 (none),
#'   1 (increasing)
#' @param n.trees Total number of boosting trees, default 200. More trees may improve
#'   accuracy but increase overfitting risk
#' @param interaction.depth Tree interaction depth, default 3. Controls interaction
#'   complexity; 1 means tree stumps (no interactions)
#' @param n.minobsinnode Minimum observations in terminal nodes. Default auto-computed as
#'   max(2, ceiling(log10(nrow(data)) + 1))
#' @param shrinkage Learning rate / shrinkage factor, default 0.05. Smaller values
#'   require more trees but generally generalize better
#' @param bag.fraction Fraction of data randomly sampled for each tree, default 0.5
#'   (stochastic gradient boosting). Set to 1 to use all data
#' @param keep.data Whether to retain training data, default FALSE to save memory
#' @param verbose Whether to output detailed training information, default FALSE
#' @param stratify Whether to stratify response variable during CV sampling, default TRUE
#'   (recommended for imbalanced data)
#' @param monitor Whether to print model summary and monitoring information, default TRUE
#' @param use Model selection strategy: "best" uses optimal tree count, "all" uses all trees
#' @param seed Random seed for reproducibility. Default 0 means no seed set
#' @param na.omit.y Whether to remove missing values in response variable y, default TRUE
#' @param cv.folds Number of cross-validation folds, default 5. Must be integer > 1;
#'   set to sample count for leave-one-out CV
#' @param train.fraction Fraction of data used for training, default 1 (use all data)
#' @param class.stratify.cv Stratification strategy for classification CV. NULL for
#'   auto-inference
#' @param n.cores Number of cores for parallel computation. NULL for auto-detection
#'
#' @return Returns an ABT-class object containing:
#' \itemize{
#'   \item gbm.model - Final GBM model trained on full data
#'   \item cv.predictions - Cross-validation predictions (each sample predicted in its
#'     validation fold)
#'   \item cv.folds - Cross-validation fold assignments for each sample
#'   \item best.n.trees - Optimal tree count determined by cross-validation
#'   \item best.cv.error - Cross-validation error at optimal tree count
#'   \item cv.error - Average cross-validation error vector across tree counts
#'   \item cv.error.matrix - Validation error matrix (folds x tree counts)
#'   \item shrinkage, interaction.depth, n.minobsinnode, bag.fraction - Model parameters
#'   \item distribution - Distribution type used
#'   \item var.importance - Variable relative importance data frame
#'   \item formula, data - Original formula and data
#' }
#'
#' @details
#' The ABT function automatically determines optimal iteration count (tree number) via
#' k-fold cross-validation to prevent overfitting. Stratified sampling (stratify=TRUE)
#' ensures consistent response distribution across folds, particularly important for
#' classification tasks or skewed response distributions.
#'
#' Supported distributions:
#' \itemize{
#'   \item gaussian: Continuous response (regression)
#'   \item bernoulli: Binary classification (0/1 response)
#'   \item poisson: Count data (non-negative integers)
#'   \item adaboost: AdaBoost exponential loss classification
#'   \item laplace: Robust regression (insensitive to outliers)
#'   \item coxph: Survival analysis (right-censored data)
#' }
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' data(softcorals)
#' fit_softcorals_gbm <- ABT(formula = Richness~Across+Along+Visibility+Slope+Flow+Wave+Sediment,
#'                           data = softcorals, n.trees = 100, cv.folds = 5)
#'
#' @export

ABT <- function(formula = formula(data), data = list(), distribution = "gaussian", 
                weights, offset = NULL, var.monotone = NULL, n.trees = 200, 
                interaction.depth = 3, n.minobsinnode = NULL,
                shrinkage = 0.05, bag.fraction = 0.5, keep.data = FALSE, 
                verbose = FALSE, stratify = TRUE, monitor = TRUE, 
                use = c("best", "all")[1], seed = 0, na.omit.y = TRUE, 
                cv.folds = 5, train.fraction = 1, class.stratify.cv = NULL) 
{
  # 检查 cv.folds 是否为大于1的正整数
  if (!is.numeric(cv.folds) || cv.folds <= 1 || cv.folds != as.integer(cv.folds)) {
    stop("cv.folds must be an integer greater than 1.")
  }
  
  # 设置 n.minobsinnode 默认值
  if (is.null(n.minobsinnode)) {
    n.minobsinnode <- max(2, ceiling(log10(max(nrow(data), 10)) + 1))
  }
  
  # 提取模型框架
  m <- match.call(expand.dots = FALSE)
  params <- c("", "formula", "data", "weights", "offset")
  m <- m[names(m) %in% params]
  m[[1]] <- as.name("model.frame")
  m$na.action <- na.pass
  m <- eval(m, parent.frame())
  
  # 提取 terms 和属性
  Terms <- attr(m, "terms")
  a <- attributes(Terms)
  
  # 提取响应变量和预测变量
  y <- model.response(m)
  x <- model.frame(delete.response(Terms), data, na.action = na.pass)
  
  # 提取权重
  w <- model.extract(m, weights)
  if (is.null(w)) w <- rep(1, length(y))
  
  # 处理响应变量中的缺失值
  if (na.omit.y && any(nays <- is.na(y))) {
    y <- y[!nays]
    x <- x[!nays, , drop = FALSE]
    w <- w[!nays]
    if (!is.null(offset)) offset <- offset[!nays]
  }
  
  # 检查 x 中是否存在缺失值
  if (any(is.na(x))) {
    warning("Missing values found in predictors. Consider imputation or use na.omit.y = FALSE with complete data.")
  }
  
  # 检查 offset 项
  if (!is.null(model.extract(m, offset))) {
    stop("In gbm the offset term needs to be specified using the offset parameter.")
  }
  
  # 提取变量名和响应名
  var.names <- a$term.labels
  response.name <- dimnames(attr(terms(formula), "factors"))[[1]][1]
  
  # 初始化
  nrows <- nrow(x)
  preds <- rep(NA, length = nrows)
  
  # 定义支持的分布
  distributions <- c("bernoulli", "gaussian", "poisson", "adaboost", "laplace", "coxph")
  distribution <- distributions[pmatch(distribution, distributions)]
  if (is.na(distribution)) {
    stop("Unsupported distribution. Choose from: ", paste(distributions, collapse = ", "))
  }
  
  # 定义 use 参数
  btype <- c("best", "all")
  use <- btype[pmatch(use, btype)]
  if (is.na(use)) {
    stop("use must be either 'best' or 'all'")
  }
  
  # 初始化分布列表
  if (is.character(distribution)) {
    distribution <- list(name = distribution)
  }
  
  class.stratify.cv <- gbm::getStratify(class.stratify.cv, d = distribution)
  
  # 设置随机种子
  if (seed > 0) {
    set.seed(seed)
  }
  
  # 创建交叉验证分组
  if (length(cv.folds) == nrows) {
    nsamps <- cv.folds
    cv.folds <- max(nsamps)
  } else {
    if (!stratify) {
      nsamps <- sample(rep(1:cv.folds, length = nrows))
    } else {
      nsamps <- rep(sample(1:cv.folds), ceiling(nrows / cv.folds))[1:nrows]
      nsamps <- nsamps[order(y)]
    }
  }
  
  # 存储每个 fold 的验证误差（用于确定最优树数量）
  valid.err.matrix <- matrix(NA, nrow = n.trees, ncol = cv.folds)
  
  # 执行交叉验证
  for (i in 1:cv.folds) {
    nTrain <- length(trainees <- which(nsamps != i))
    non.trainees <- which(nsamps == i)
    
    # 分割数据
    temp.x <- x[c(trainees, non.trainees), , drop = FALSE]
    temp.y <- y[c(trainees, non.trainees)]
    temp.w <- w[c(trainees, non.trainees)]
    temp.offset <- if (!is.null(offset)) offset[c(trainees, non.trainees)] else NULL
    
    # 拟合 GBM 模型
    gbm.obj <- gbm::gbm.fit(temp.x, temp.y, offset = temp.offset,
                            distribution = distribution, w = temp.w, 
                            var.monotone = var.monotone,
                            n.trees = n.trees, interaction.depth = interaction.depth,
                            n.minobsinnode = n.minobsinnode, shrinkage = shrinkage,
                            bag.fraction = bag.fraction, nTrain = nTrain,
                            keep.data = keep.data, verbose = verbose, 
                            var.names = var.names,
                            response.name = response.name)
    
    # 存储验证误差
    valid.err.matrix[, i] <- gbm.obj$valid.error[1:n.trees]
    
    # 对当前 fold 的验证集进行预测
    preds[non.trainees] <- gbm::predict.gbm(gbm.obj, 
                                            newdata = x[non.trainees, , drop = FALSE],
                                            n.trees = n.trees)
  }
  
  # 计算平均验证误差
  valid.err <- rowMeans(valid.err.matrix, na.rm = TRUE)
  
  # 确定最优树数量
  best <- which.min(valid.err)
  best.error <- valid.err[best]
  
  # 使用全量数据训练最终模型（使用最优树数量）
  final.gbm <- gbm::gbm.fit(x, y, offset = offset,
                            distribution = distribution, w = w, 
                            var.monotone = var.monotone,
                            n.trees = n.trees, interaction.depth = interaction.depth,
                            n.minobsinnode = n.minobsinnode, shrinkage = shrinkage,
                            bag.fraction = bag.fraction, nTrain = nrows,
                            keep.data = keep.data, verbose = verbose, 
                            var.names = var.names,
                            response.name = response.name)
  
  # 监控输出
  if (monitor) {
    cat("\n========================================\n")
    cat("ABT Model Summary\n")
    cat("========================================\n")
    cat("N trees =", n.trees, "\n")
    cat("Depth =", interaction.depth, "\n")
    cat("Minimum node size =", n.minobsinnode, "\n")
    cat("Distribution =", distribution$name, "\n")
    cat("Shrinkage =", signif(shrinkage, 3), "\n")
    cat("Bag Fraction =", bag.fraction, "\n")
    cat("CV Folds =", cv.folds, "\n")
    cat("----------------------------------------\n")
    cat("Best size =", best, "\n")
    cat("Best CV error =", signif(best.error, 5), "\n")
    cat("========================================\n\n")
    
    # 调用修正后的 monitor 函数
    monitor.gbm(final.gbm, cv.error = valid.err, best = best)
  }
  
  # 计算变量重要性
  var.importance <- gbm::summary.gbm(final.gbm, plotit = FALSE)
  
  # 构建返回对象
  result <- list(
    gbm.model = final.gbm,
    cv.predictions = preds,
    cv.folds = nsamps,
    best.n.trees = best,
    best.cv.error = best.error,
    cv.error = valid.err,
    cv.error.matrix = valid.err.matrix,
    shrinkage = shrinkage,
    interaction.depth = interaction.depth,
    n.minobsinnode = n.minobsinnode,
    bag.fraction = bag.fraction,
    distribution = distribution$name,
    var.importance = var.importance,
    formula = formula,
    data = data
  )
  
  class(result) <- "ABT"
  
  return(result)
}

monitor.gbm <- function(obj, cv.error = NULL, best = NULL) 
{
  nt <- 1:obj$n.trees
  
  # 确定监控间隔（对数尺度）
  if (length(nt) <= 20) {
    ind <- nt
  } else {
    ind <- unique(c(1, round(seq(1, length(nt), length.out = 20)), max(nt)))
  }
  
  # 提取训练误差
  train.err <- obj$train.error[ind]
  
  # 创建数据框
  if (!is.null(cv.error) && length(cv.error) == obj$n.trees) {
    val.err <- cv.error[ind]
    df <- data.frame(
      Iteration = ind,
      TrainError = signif(train.err, 5),
      ValidError = signif(val.err, 5)
    )
    cat("Training & Cross-Validation Error:\n")
  } else if (!is.null(obj$valid.error) && length(obj$valid.error) > 0) {
    val.err <- obj$valid.error[ind]
    df <- data.frame(
      Iteration = ind,
      TrainError = signif(train.err, 5),
      ValidError = signif(val.err, 5)
    )
    cat("Training & Validation Error:\n")
  } else {
    df <- data.frame(
      Iteration = ind,
      TrainError = signif(train.err, 5)
    )
    cat("Training Error:\n")
  }
  
  print(df, row.names = FALSE)
  
  # 标记最优迭代
  if (!is.null(best)) {
    cat("\nOptimal iteration:", best, "\n")
  }
  
  cat("\n")
}
