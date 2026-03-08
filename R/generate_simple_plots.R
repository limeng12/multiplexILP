
#' Generate Simple Visualization Plots for ILP Results
#'
#' Creates basic plots to visualize the selected primer pairs and their
#' amplicon size distribution across channels.
#'
#' @param results Data frame of selected primer pairs (from parser functions)
#' @param output_prefix Prefix for output plot files
#'
#' @return A list containing:
#'   \item{amplicon_plot}{ggplot2 object of amplicon size plot}
#'
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' parsed_result <- parse_scip_locus_C_kt_solution()
#' plots <- generate_simple_plots(parsed_result$selected_primers, "my_result")
#' }
#'
#' @export
generate_simple_plots <- function(results, output_prefix) {
  cat("\n=== GENERATING SIMPLE VISUALIZATIONS ===\n")
  
  if (nrow(results) == 0) {
    cat("No results to visualize\n")
    return(NULL)
  }
  
  library(ggplot2)
  library(dplyr)
  
  # 1. 扩增子大小图（最核心的图）
  p1 <- ggplot(results, 
               aes(y = factor(channel), 
                   x = amplicon_mid, 
                   color = locus, 
                   fill = locus)) +
    
    # 用矩形表示范围（横向）- 调整高度更低
    geom_rect(aes(xmin = amplicon_low, 
                  xmax = amplicon_high,
                  ymin = as.numeric(factor(channel)) - 0.15,  # 高度降低：0.3→0.15
                  ymax = as.numeric(factor(channel)) + 0.15),
              alpha = 0.4,  # 透明度降低：0.2→0.4（颜色更深）
              color = NA) +
    
    # 添加中点标记
    geom_point(size = 4) +
    
    # 添加标签
    geom_text(aes(label = paste0(locus, "\n", amplicon_low, "-", amplicon_high, "bp")),
              hjust = -0.1,
              size = 3) +
    
    # 调整坐标轴和标签
    labs(
      title = "Multiplex PCR Design - Amplicon Sizes",
      x = "Amplicon Size (bp)",
      y = "Channel"
      # 移除了color和fill图例标题
    ) +
    
    # 横向主题
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"  # 完全移除图例
    ) +
    
    # 调整标签位置，防止重叠
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  #print(p1)
  #print(p1_horizontal)
  
  ggsave(paste0(output_prefix, "_amplicon_sizes.pdf"), p1, width = 8, height = 2*max(results$channel) );
  cat("Saved: ", output_prefix, "_amplicon_sizes.pdf\n", sep = "")
  
  # 4. 一个汇总图（仪表板风格）
  # create_summary_plot(results, output_prefix)
  
  cat("\nVisualizations generated successfully!\n")
  return(list(amplicon_plot = p1))
}
