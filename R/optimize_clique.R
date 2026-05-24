#' Optimize clique constraints in LP file
#'
#' Replaces pairwise constraints (x_i + x_j <= 1) with optimized clique constraints.
#'
#' @param in_file Input LP file path. Default "ilp_problem.lp".
#' @param out_file Output LP file path. Default "ilp_problem_opti.lp".
#'
#' @return Invisibly returns number of new constraints added.
#'
#' @export
optimize_clique <- function(in_file = "ilp_problem.lp", out_file = "ilp_problem_opti.lp") {
  lp <- read_lines(in_file)
  
  #pattern <- "^\\s*[A-Za-z0-9_]+\\s*\\+\\s*[A-Za-z0-9_]+\\s*<=\\s*1\\s*$"
  #idx <- which(str_detect(lp, pattern))
  
  idx <- which(grepl("^\\s*[A-Za-z0-9_]+\\s*\\+\\s*[A-Za-z0-9_]+\\s*<=\\s*1\\s*$", lp, perl = TRUE))
  
  if (length(idx) == 0) {
    write_lines(lp, out_file)
    return(invisible(0))
  }
  
  # 只对输入convert_constraints的内容排序
  cat(sprintf("[%s] grouping constraints...\n", length(idx) ) );
  
  #lp_Y_Y <- (lp[idx])
  
  #nums <- do.call(rbind, strsplit( (lp[idx]), "\\D+", perl = TRUE) );
  #dt <- data.table(idx = seq_along( (lp[idx]) ), n1 = as.numeric(nums[,2]), n2 = as.numeric(nums[,3]) );rm(nums)
  parts <- tstrsplit( (lp[idx]), "\\D+", perl = TRUE)
  dt <- data.table(idx = seq_along( (lp[idx]) ),
                   n1 = as.numeric(parts[[2]]),
                   n2 = as.numeric(parts[[3]]))

  # numbers <- stri_extract_all_regex(lp[idx], "\\d+", simplify = TRUE)
  # 
  # dt <- data.table(
  #   idx = seq_along(lp[idx]),
  #   n1 = as.numeric(numbers[, 1]),  # 第一个数字 (100)
  #   n2 = as.numeric(numbers[, 2])   # 第二个数字 (10000)
  # )
  
  
  
  #batch_size <- 10000
  #分组
  small_bin <- 100000
  
  dt[, `:=`(s1 = n1 %/% small_bin, s2 = n2 %/% small_bin)]
  dt[, `:=`(b1 = s1 %/% 2, b2 = s2 %/% 2)]
  
  # 使用整数键而不是字符串
  dt[, group_key := b1]
  dt[b1 != b2, group_key := -(s1 * 1000000L + s2)]
  
  result <- split(dt$idx, dt$group_key)
  
  
  
  all_new_cons <- c()
  cat(sprintf("[%s] total batch num\n", length(result) ) );
  write_lines(lp[1:(idx[1]-1)], out_file)
  
  
  for (i in 1:length(result) ) {
    cat(sprintf("batch [%s] | ", i ) );
    
    all_new_cons <- c( convert_constraints( (lp[idx])[result[[i]] ] ) );
    write_lines(all_new_cons, out_file, append=TRUE)
    
  }
  
  # 使用原始idx，保持文件结构不变
  #new_lp <- c(lp[1:(idx[1]-1)], all_new_cons, lp[setdiff(1:length(lp), c(idx, 1:(idx[1]-1)))])
  
  write_lines(lp[setdiff(1:length(lp), c(idx, 1:(idx[1]-1)))], out_file, append=TRUE)
  invisible(length(all_new_cons))
}

# 使用
# optimize_clique()  # 使用默认值
# optimize_clique("clique_fast.lp", "clique_optimized3.lp")

