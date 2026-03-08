#############################################################引物对##########################


#' Generate Simulated Primer Pair Data
#'
#' Creates synthetic primer pair data for testing multiplex PCR ILP models.
#' Each primer pair consists of a 5' and 3' primer for a specific locus.
#'
#' @param n_pairs Number of primer pairs to generate (default: 2000)
#' @param n_loci Number of distinct loci (default: 24)
#' @param num_channels Number of fluorescence channels (default: 3)
#' @param size_gap Minimum size gap between amplicons (default: 10)
#' @param seed Random seed for reproducibility (default: 1123)
#'
#' @return A list containing:
#'   \item{primers_5}{Character vector of 5' primer sequences}
#'   \item{primers_3}{Character vector of 3' primer sequences}
#'   \item{ids_5}{Character vector of 5' primer IDs}
#'   \item{ids_3}{Character vector of 3' primer IDs}
#'   \item{str_loci_all}{Character vector of locus assignments}
#'   \item{dis_5}{Numeric vector of 5' primer distances}
#'   \item{dis_3}{Numeric vector of 3' primer distances}
#'   \item{temp_5}{Numeric vector of 5' primer annealing temperatures}
#'   \item{temp_3}{Numeric vector of 3' primer annealing temperatures}
#'   \item{weights}{Numeric vector of primer pair weights}
#'   \item{amplicon_size_low}{Numeric vector of minimum amplicon sizes}
#'   \item{amplicon_size_high}{Numeric vector of maximum amplicon sizes}
#'   \item{compatible_pairs}{Character vector of compatible primer pair combinations}
#'   \item{num_channels}{Number of channels}
#'   \item{f_color_lenght_max}{Placeholder for max channel length}
#'   \item{size_gap}{Size gap parameter}
#'
#' @examples
#' \dontrun{
#' primer_data <- generate_primer_pair_data(n_pairs = 100, n_loci = 10)
#' }
#'
#' @export
generate_primer_pair_data <- function(n_pairs = 2000, 
                                 n_loci = 24, 
                                 num_channels = 3,
                                 size_gap = 10,
                                 seed = 1123) {
  
  # 设置随机种子
  set.seed(seed)
  
  # 4. 生成随机引物序列（8bp）
  generate_random_primer <- function(n) {
    bases <- c("A", "T", "C", "G")
    sapply(1:n, function(x) paste(sample(bases, 8, replace = TRUE), collapse = ""))
  }
  
  primers_5 <- generate_random_primer(n_pairs)
  primers_3 <- generate_random_primer(n_pairs)
  # 5. ID和位点分配
  ids <- paste0("P", 1:n_pairs)  # 引物对ID
  ids_5 <- paste0("P", 1:n_pairs, "_5")
  ids_3 <- paste0("P", 1:n_pairs, "_3")
  # 6. 基本属性
  dis_5 <- rep(1.0, n_pairs)
  dis_3 <- rep(1.0, n_pairs)
  temp_5 <- rep(60.0, n_pairs)
  temp_3 <- rep(60.0, n_pairs)
  
  
  
  # 确保每个位点至少有1个引物对
  primer_per_locus <- ceiling(n_pairs / n_loci)
  
  str_loci_all <- character(n_pairs)
  for (i in 1:n_loci) {
    start_idx <- (i - 1) * primer_per_locus + 1
    end_idx <- min(i * primer_per_locus, n_pairs)
    
    if (start_idx <= n_pairs) {
      str_loci_all[start_idx:end_idx] <- paste0("LOCUS", i)
    }
  }
  

  # 7. 扩增子大小 - 每个位点内使用相同的扩增子范围
  amplicon_size_low <- integer(n_pairs)
  amplicon_size_high <- integer(n_pairs)
  
  # 为每个位点生成一个基础扩增子范围
  for (i in 1:n_loci) {
    # 找出当前位点的引物对索引
    locus_indices <- which(str_loci_all == paste0("LOCUS", i))
    
    if (length(locus_indices) > 0) {
      # 为这个位点生成一个基础扩增子大小
      base_low <- sample(50:200, length(locus_indices), replace = TRUE)
      base_high <- base_low + sample(30:150, 1, replace = TRUE)
      
      # 将该位点内所有引物对设为相同的扩增子范围
      amplicon_size_low[locus_indices] <- base_low
      amplicon_size_high[locus_indices] <- base_high
    }
  }
  weights <- (10000 + floor(10000 / amplicon_size_low)) / 100
  
  # ============ 第一部分：生成兼容对关系 ============
  cat("\n=== 生成兼容对关系 ===\n")
  
  # 生成ids_merge_raw
  ids_merge_raw <- paste0(str_loci_all, "$P", 1:n_pairs, "_5&", str_loci_all, "$P", 1:n_pairs, "_3")
  
  compatible_pairs <- character(0)
  loci <- unique(str_loci_all)
  n_loci_unique <- length(loci)
  
  # 为每个位点生成引物ID列表
  loci_primer_ids <- lapply(loci, function(loc) {
    idx <- which(str_loci_all == loc)
    ids_merge_raw[idx]
  })
  
  # 计算总迭代次数
  total_iterations <- (n_loci_unique - 1) * n_loci_unique / 2
  
  # 创建进度条
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3, width = 50, char = "=")
  iter_counter <- 0
  
  # 生成跨位点的兼容对
  compatible_pairs_list <- vector("list", total_iterations * 2)
  counter <- 1
  
  for (i in 1:(n_loci_unique - 1)) {
    for (j in (i + 1):n_loci_unique) {
      iter_counter <- iter_counter + 1
      setTxtProgressBar(pb, iter_counter)
      
      n_from_i <- min(floor(length(loci_primer_ids[[i]]) / 1.3), length(loci_primer_ids[[i]]))
      n_from_j <- min(floor(length(loci_primer_ids[[j]]) / 1.3), length(loci_primer_ids[[j]]))
      
      from_i <- loci_primer_ids[[i]][sample(1:length(loci_primer_ids[[i]]), n_from_i)]
      from_j <- loci_primer_ids[[j]][sample(1:length(loci_primer_ids[[j]]), n_from_j)]
      
      # 生成所有组合
      combinations <- expand.grid(a = from_i, b = from_j)
      
      # 添加正反两种组合
      compatible_pairs_list[[counter]] <- paste0(combinations$a, "#", combinations$b)
      counter <- counter + 1
      
      compatible_pairs_list[[counter]] <- paste0(combinations$b, "#", combinations$a)
      counter <- counter + 1
    }
  }
  close(pb)
  
  compatible_pairs <- unlist(compatible_pairs_list)
  
  cat("  总兼容对数量:", length(compatible_pairs), "\n")
  cat("  前3个兼容对示例:\n")
  if (length(compatible_pairs) > 0) {
    for (i in 1:min(3, length(compatible_pairs))) {
      cat("    ", compatible_pairs[i], "\n")
    }
  }
  
  # 将所有数据组合成一个list
  primer_data <- list(
    primers_5 = primers_5,
    primers_3 = primers_3,
    ids_5 = ids_5,
    ids_3 = ids_3,
    str_loci_all = str_loci_all,
    dis_5 = dis_5,
    dis_3 = dis_3,
    temp_5 = temp_5,
    temp_3 = temp_3,
    weights = weights,
    amplicon_size_low = amplicon_size_low,
    amplicon_size_high = amplicon_size_high,
    compatible_pairs = compatible_pairs,
    num_channels = num_channels,
    f_color_lenght_max = rep(0, num_channels),
    size_gap = size_gap
  )
  
  return(primer_data)
}



# 函数2：ILP生成

#' Generate ILP  for Primer Pair Data
#'
#' Executes the ILP model generation for primer pair data created by
#' \code{generate_primer_pair_data}.
#'
#' @param primer_data List returned by \code{generate_primer_pair_data}
#' @param lower_bound Optional lower bound for objective function (default: 0)
#'
#' @return Result from \code{generate_correct_multiplex_ilp_locus_C_fixed_binary}
#'
#' @examples
#' \dontrun{
#' primer_data <- generate_primer_pair_data(n_pairs = 100, n_loci = 10)
#' ilp_result <- gen_ilp_primer_pair(primer_data)
#' }
#'
#' @export
gen_ilp_primer_pair <- function(primer_data, lower_bound = 0) {
  
  cat("\n=== 运行ILP函数 ===\n")
  
  # 从list中提取参数并运行ILP函数
  ilp_result <- generate_correct_multiplex_ilp_locus_C_fixed_binary(
    primers_5 = primer_data$primers_5,
    primers_3 = primer_data$primers_3,
    ids_5 = primer_data$ids_5,
    ids_3 = primer_data$ids_3,
    str_loci_all = primer_data$str_loci_all,
    dis_5 = primer_data$dis_5,
    dis_3 = primer_data$dis_3,
    temp_5 = primer_data$temp_5,
    temp_3 = primer_data$temp_3,
    weights = primer_data$weights,
    amplicon_size_low = primer_data$amplicon_size_low,
    amplicon_size_high = primer_data$amplicon_size_high,
    compatible_pairs = primer_data$compatible_pairs,
    num_channels_v = primer_data$num_channels,
    size_gap_v = primer_data$size_gap,
    lower_bound = lower_bound
  )
  
  return(ilp_result)
}



#############################################################单独引物##########################
# 函数1：生成单个引物数据

#' Generate Simulated Single Primer Data
#'
#' Creates synthetic single primer data for testing multiplex PCR ILP models.
#' Each primer is either 5' or 3' end and belongs to a specific locus.
#'
#' @param n_primers Number of primers to generate (default: 100)
#' @param n_loci Number of distinct loci (default: 6)
#' @param num_channels Number of fluorescence channels (default: 2)
#' @param size_gap Minimum size gap between amplicons (default: 10)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A list containing:
#'   \item{primers}{Character vector of primer sequences}
#'   \item{primer_ids}{Character vector of primer IDs}
#'   \item{primer53}{Character vector indicating "5" or "3" end}
#'   \item{primer_temp}{Numeric vector of annealing temperatures}
#'   \item{locus_ids}{Character vector of locus assignments}
#'   \item{distances}{Numeric vector of primer distances}
#'   \item{copy_length_max}{Numeric vector of maximum copy lengths}
#'   \item{copy_length_min}{Numeric vector of minimum copy lengths}
#'   \item{repeat_distances}{Numeric vector of repeat region distances}
#'   \item{compatible_pairs}{Character vector of compatible primer combinations}
#'   \item{num_channels}{Number of channels}
#'   \item{size_gap}{Size gap parameter}
#'
#' @examples
#' \dontrun{
#' primer_data <- generate_single_primer_data(n_primers = 50, n_loci = 5)
#' }
#'
#' @export
generate_single_primer_data <- function(n_primers = 100,
                                        n_loci = 6,
                                        num_channels = 2,
                                        size_gap = 10,
                                        seed = 123) {
  
  cat("\n=== 参数设置 ===\n")
  set.seed(seed)
  
  # 3. 生成随机引物数据
  cat("\n=== 生成随机引物数据 ===\n")
  
  # 生成随机引物序列（8bp）
  generate_random_primer <- function(n) {
    bases <- c("A", "T", "C", "G")
    sapply(1:n, function(x) paste(sample(bases, 8, replace = TRUE), collapse = ""))
  }
  
  primers <- generate_random_primer(n_primers)
  primer_ids <- paste0("P", 1:n_primers)
  primer_temp <- rep(52, length(primers))
  distances <- sample(10:100, n_primers, replace = TRUE)
  
  
  # 分配5'和3'端（确保每个位点都有5'和3'引物）
  cat("分配5'和3'端引物...\n")
  primer53 <- character(n_primers)
  
  # 先分配位点
  primer_per_locus <- ceiling(n_primers / n_loci)
  locus_ids <- character(n_primers)
  loci <- paste0("LOCUS", 1:n_loci)
  
  for (i in 1:n_loci) {
    start_idx <- (i - 1) * primer_per_locus + 1
    end_idx <- min(i * primer_per_locus, n_primers)
    
    if (start_idx <= n_primers) {
      locus_ids[start_idx:end_idx] <- loci[i]
    }
  }
  
  
  # 为每个位点分配5'和3'引物
  for (locus in unique(locus_ids)) {
    idx <- which(locus_ids == locus)
    n_in_locus <- length(idx)
    
    # 确保至少有一个5'和一个3'引物
    if (n_in_locus >= 2) {
      # 随机分配约一半为5'，一半为3'
      n_5 <- max(1, round(n_in_locus / 2))
      n_3 <- n_in_locus - n_5
      
      primer53[idx[1:n_5]] <- "5"
      primer53[idx[(n_5+1):n_in_locus]] <- "3"
    } else {
      # 如果只有一个引物，随机分配为5'或3'
      primer53[idx] <- sample(c("5", "3"), 1)
    }
  }
  
  cat("引物端分布:\n")
  cat("  5'端引物:", sum(primer53 == "5"), "\n")
  cat("  3'端引物:", sum(primer53 == "3"), "\n")
  
  # 4. 生成其他参数
  cat("\n=== 生成其他参数 ===\n")
  
  # 实际计算唯一的位点数
  actual_n_loci <- length(unique(locus_ids))
  cat("实际位点数量:", actual_n_loci, "\n")
  
  # 位点长度范围 - 这是位点的属性，同一locus应该相同
  locus_properties <- data.frame(
    locus = unique(locus_ids),
    copy_min = sample(20:50, actual_n_loci, replace = TRUE),
    copy_max = sample(30:40, actual_n_loci, replace = TRUE),  # 确保max > min
    repeat_dist = sample(10:40, actual_n_loci, replace = TRUE)
  )
  
  # 确保copy_max > copy_min
  locus_properties$copy_max <- (locus_properties$copy_min + locus_properties$repeat_dist)
  
  # 现在为每个引物分配其所在locus的属性
  copy_length_min <- numeric(n_primers)
  copy_length_max <- numeric(n_primers)
  repeat_distances <- numeric(n_primers)
  
  for (i in 1:n_primers) {
    locus <- locus_ids[i]
    props <- locus_properties[locus_properties$locus == locus, ]
    
    copy_length_min[i] <- props$copy_min
    copy_length_max[i] <- props$copy_max
    repeat_distances[i] <- props$repeat_dist
  }
  
  # 5. 生成兼容对关系
  cat("\n=== 生成兼容对关系 ===\n")
  
  compatible_pairs <- character(0)
  unique_loci <- unique(locus_ids)
  n_unique_loci <- length(unique_loci)
  
  # 为每个位点生成引物ID列表
  loci_primer_ids <- lapply(unique_loci, function(loc) {
    primer_ids[which(locus_ids == loc)]  # 使用引物ID
  })
  
  # 计算总迭代次数
  total_iterations <- (n_unique_loci - 1) * n_unique_loci / 2
  
  # 创建进度条
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3, width = 50, char = "=")
  iter_counter <- 0
  
  # 生成跨位点的兼容对
  compatible_pairs_list <- list()
  counter <- 1
  
  for (i in 1:(n_unique_loci-1)) {
    for (j in (i+1):n_unique_loci) {
      iter_counter <- iter_counter + 1
      setTxtProgressBar(pb, iter_counter)
      
      n_from_i <- min(floor(length(loci_primer_ids[[i]])/1.1), length(loci_primer_ids[[i]]))
      n_from_j <- min(floor(length(loci_primer_ids[[j]])/1.1), length(loci_primer_ids[[j]]))
      
      from_i <- sample(loci_primer_ids[[i]], n_from_i)
      from_j <- sample(loci_primer_ids[[j]], n_from_j)
      
      # 生成所有组合
      combinations <- expand.grid(a = from_i, b = from_j)
      
      # 添加正反两种组合
      pairs <- c(paste0(combinations$a, "&", combinations$b),
                 paste0(combinations$b, "&", combinations$a))
      
      compatible_pairs_list[[counter]] <- pairs
      counter <- counter + 1
    }
  }
  close(pb)
  
  compatible_pairs <- unlist(compatible_pairs_list)
  
  cat("总兼容对数量:", length(compatible_pairs), "\n")
  cat("前3个兼容对示例:\n")
  if (length(compatible_pairs) > 0) {
    for (i in 1:min(3, length(compatible_pairs))) {
      cat("  ", compatible_pairs[i], "\n")
    }
  }
  
  # 将所有数据组合成一个list
  primer_data <- list(
    primers = primers,
    primer_ids = primer_ids,
    primer53 = primer53,
    primer_temp = primer_temp,
    locus_ids = locus_ids,
    distances = distances,
    copy_length_max = copy_length_max,
    copy_length_min = copy_length_min,
    repeat_distances = repeat_distances,
    compatible_pairs = compatible_pairs,
    num_channels = num_channels,
    size_gap = size_gap
  )
  
  return(primer_data)
}

# 函数2：运行单个引物ILP分析

#' Generate ILP for Single Primer Data
#'
#' Executes the ILP model generation for single primer data created by
#' \code{generate_single_primer_data}.
#'
#' @param primer_data List returned by \code{generate_single_primer_data}
#'
#' @return Result from \code{generate_multiplex_ilp_single_primers_direct}
#'
#' @examples
#' \dontrun{
#' primer_data <- generate_single_primer_data(n_primers = 50, n_loci = 5)
#' ilp_result <- gen_ilp_single_primer(primer_data)
#' }
#'
#' @export
gen_ilp_single_primer <- function(primer_data) {
  
  cat("\n=== 运行单个引物ILP模型 ===\n")
  
  # 从list中提取参数并运行ILP函数
  result <- generate_multiplex_ilp_single_primers_direct(
    primers = primer_data$primers,
    primer_ids = primer_data$primer_ids,
    primer53 = primer_data$primer53,
    primer_temp = primer_data$primer_temp,
    locus_ids = primer_data$locus_ids,
    distances = primer_data$distances,
    copy_length_max = primer_data$copy_length_max,
    copy_length_min = primer_data$copy_length_min,
    repeat_distances = primer_data$repeat_distances,
    compatible_pairs = primer_data$compatible_pairs,
    num_channels_v = primer_data$num_channels,
    size_gap_v = primer_data$size_gap
  )
  
  return(result)
}


