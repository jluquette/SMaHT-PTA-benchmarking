library(Cairo)
library(GenomicRanges)
library(data.table)
library(tidyverse)
library(IRdisplay)
library(dplyr)
library(patchwork)


add_stuff_to_seg_merged <- function(seg_merged_data){
    seg_merged_data$CN_total <- NA
    
    seg_merged_data$CN <- as.character(seg_merged_data$CN)
    
    # Split CN into two numeric components
    cn_split <- str_split(seg_merged_data$CN, "\\|")
    
    # Extract alleles as numeric values
    seg_merged_data$CN_A <- as.numeric(sapply(cn_split, `[`, 1))
    seg_merged_data$CN_B <- as.numeric(sapply(cn_split, `[`, 2))
    
    # Compute total copy number
    sex_chr_mask <- seg_merged_data$chrom %in% c("X", "Y")
    seg_merged_data$CN_A[sex_chr_mask & is.na(seg_merged_data$CN_A)] <- 0
    seg_merged_data$CN_B[sex_chr_mask & is.na(seg_merged_data$CN_B)] <- 0
    
    seg_merged_data$CN_total <- seg_merged_data$CN_A + seg_merged_data$CN_B
     
    #seg_merged_data$CN_total[seg_merged_data$chrom %in% c("X", "Y")] <- as.integer(seg_merged_data$CN[seg_merged_data$chrom %in% c("X", "Y")])
    seg_merged_data$is_LOH <- (seg_merged_data$chrom %in% c(1:22)) & (seg_merged_data$CN_A == 0 | seg_merged_data$CN_B == 0)
    seg_merged_data$is_CNV <- (seg_merged_data$chrom %in% c(1:22)) & (seg_merged_data$CN != "1|1")
    seg_merged_data$length <- seg_merged_data$end - seg_merged_data$start
    seg_merged_data$is_gain <- (seg_merged_data$chrom %in% c(1:22)) & seg_merged_data$CN_total>2
    seg_merged_data$is_loss <- (seg_merged_data$chrom %in% c(1:22)) & seg_merged_data$CN_total<2

    seg_merged_data
    }



plot_cell2_supplementary <- function(cell_name, title, cell_data, seg_merged_data, cnv_presence = NULL,
                       output_dir = NULL, print_plot = TRUE, chrom_list = c(1:22, "X", "Y")) {
  category_colors <- c(
    "Aneuploidy" = "#F0E442",  # yellow
    "CN-LOH"     = "#009E73",  # green
    "Deletion"   = "#56B4E9",  # blue
    "Duplication"= "#D55E00",   # red-orange
    "Failed QC" = "#000000" # black
  )

  max_y <- 4
  max_y_offset <- 0.25
    
  cell_df <- cell_data[cell_data$cell_name == cell_name & cell_data$CHROM %in% chrom_list, ]
  seg_df <- seg_merged_data[seg_merged_data$cell_name == cell_name & seg_merged_data$chrom %in% chrom_list, ]


  
    
  chrom_order <- as.character(chrom_list)
  
  chrom_bin_overrides <- tibble(CHROM = chrom_order) %>%
    mutate(force_n_bins = case_when(
      CHROM == "X" ~ 200,
      CHROM == "Y" ~ 100,
      TRUE ~ NA_real_
    ))
  
  cell_df$CHROM <- factor(cell_df$CHROM, levels = chrom_order)
  seg_df$chrom <- factor(seg_df$chrom, levels = chrom_order)
  
  chrom_sizes <- cell_df %>%
    count(CHROM, name = "actual_n_bins") %>%
    full_join(chrom_bin_overrides, by = "CHROM") %>%
    mutate(n_bins = ifelse(!is.na(force_n_bins), force_n_bins, actual_n_bins)) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(n_bins), default = 0))
  
  cell_df <- cell_df %>%
    left_join(chrom_sizes, by = "CHROM") %>%
    group_by(CHROM) %>%
    arrange(START) %>%
    mutate(
      bin_index = row_number(),
      bin_spacing = n_bins / n(),
      x = offset + (bin_index - 1) * bin_spacing
    ) %>%
    ungroup()
  
  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  x_breaks <- chrom_centers$center
  x_labels <- chrom_centers$CHROM
  
  bin_lookup_start <- cell_df %>% select(CHROM, START, x)
  bin_lookup_end <- cell_df %>% select(CHROM, END, x) %>% distinct()
  
  seg_df <- seg_df %>%
    left_join(bin_lookup_start, by = c("chrom" = "CHROM", "start" = "START")) %>%
    rename(x_start = x) %>%
    left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
    rename(x_end = x) %>%
    filter(!is.na(x_start), !is.na(x_end))
  
  # === CN LINE PLOTTING ===
  # Create stepwise line segments for CN
  cn_path_df <- seg_df %>%
    mutate(segment_id = row_number()) %>%
    pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
    arrange(segment_id, x) %>%
    mutate(y = CN_total)
  
  # === VERTICAL JUMPS ONLY IF CN CHANGES OR CHROM CHANGES ===
  jump_df <- seg_df %>%
    arrange(chrom, start) %>%
    mutate(
      next_CN = lead(CN_total),
      next_chrom = lead(chrom),
      next_x = lead(x_start)
    ) %>%
    filter(!is.na(next_CN)) %>%
    filter(
      chrom != next_chrom |
        abs(CN_total - next_CN) > 1e-3  # ignore tiny differences
    ) %>%
    transmute(
      x = x_end,
      xend = x,
      y = pmin(CN_total, max_y),
      yend = pmin(next_CN, max_y)
    )
  
  chrom_backgrounds <- chrom_sizes %>%
    mutate(
      xmin = offset,
      xmax = offset + n_bins,
      fill = as.factor(row_number() %% 2)
    )

    if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
      x_min <- min(cell_df$x, na.rm = TRUE)
      x_max <- max(cell_df$x, na.rm = TRUE)
      
      n_squares <- nrow(cnv_presence)
      
      # Let's reserve about 5% of the x-axis width for all squares combined,
      # then divide that by number of squares to get spacing.
      total_space <- (x_max - x_min) * 0.1
      
      x_spacing <- total_space / max(n_squares - 1, 1)  # avoid division by zero
      
      cnv_presence <- cnv_presence %>%
        mutate(
          # Place squares starting from right (x_max - margin) going left
          x = x_max - (row_number() - 0.75) * x_spacing,
          y = max_y + 0.5
        )
    }
  cell_df <- cell_df %>% mutate(above_max_y=CN_signal>max_y)
  cell_df$CN_signal <- pmin(cell_df$CN_signal, max_y)
  # === CN PLOT ===
  p1 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = CN_signal, shape=factor(above_max_y)),
               color="#999933", size = 0.1, alpha = 0.9, stroke=0.1) +
    scale_shape_manual(values=c(16, 2), guide="none") + #(16, 2)
    geom_line(data = cn_path_df, aes(x = x, y = y, group = segment_id),
              color = "black", size = 0.6) +
    geom_segment(data = jump_df, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black") +
    {
    if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
      geom_point(data = cnv_presence, aes(x = x, y = y, color = category),
                 shape = 15, size = 3)  # bigger squares here
    } else {
      NULL
    }
      } +
      scale_color_manual(values = category_colors) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
      coord_cartesian(ylim = c(0, max_y + max_y_offset), clip = "off") +  # extend y limit to 7 here
      scale_y_continuous(breaks = c(0, 2, 4)) +
      labs(x = NULL, y = "CN") +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_text(size=10, color="black"),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 25, r = 5, b = 0, l = 5),  # add top margin to fit squares and title
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
        axis.line = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.position = "none"
      ) +
      ggtitle(title)
  
      # === BAF PLOT ===
      p2 <- ggplot() +
        geom_rect(data = chrom_backgrounds,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                  alpha = 0.5, inherit.aes = FALSE) +
        scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
        geom_point(data = cell_df, aes(x = x, y = pBAF), color = "#3399FF", size = 0.1/4, alpha = 0.9, stroke=0.1/4) +    
        scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(x = "Chromosome", y = "BAF") +
        theme_bw(base_size = 10) +
        theme(
          panel.grid = element_blank(),
          axis.title.y = element_text(size=10, color="black"),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
          axis.line = element_blank()
        )
  
      plot_output <- p1 / plot_spacer() / p2 + plot_layout(heights = c(1, 0.1, 1))
  
  if (!is.null(output_dir)){
      ggsave(
        file.path(output_dir, paste0(cell_name, "_track.png")),
        plot_output,
        width = 10, height = 4, dpi = 300
      )
      ggsave(
        file.path(output_dir, paste0(cell_name, "_track.pdf")),
        plot_output,
        width = 10, height = 4, device = cairo_pdf
      )
      
  }
  if (print_plot) {
    print(plot_output)
  }
  return(list(rdr=p1, baf=p2))
}


plot_cell2_supplementary2 <- function(cell_name, title, cell_data, seg_merged_data, cnv_presence = NULL,
                       output_dir = NULL, print_plot = TRUE, chrom_list = c(1:22, "X", "Y")) {
  category_colors <- c(
    "Aneuploidy" = "#F0E442",  # yellow
    "CN-LOH"     = "#009E73",  # green
    "Deletion"   = "#56B4E9",  # blue
    "Duplication"= "#D55E00",   # red-orange
    "Failed QC" = "#000000" # black
  )

  max_y <- 4
  max_y_offset <- 0.25
    
  cell_df <- cell_data[cell_data$cell_name == cell_name & cell_data$CHROM %in% chrom_list, ]
  seg_df <- seg_merged_data[seg_merged_data$cell_name == cell_name & seg_merged_data$chrom %in% chrom_list, ]

  chrom_order <- as.character(chrom_list)
  
  chrom_bin_overrides <- tibble(CHROM = chrom_order) %>%
    mutate(force_n_bins = case_when(
      CHROM == "X" ~ 200,
      CHROM == "Y" ~ 100,
      TRUE ~ NA_real_
    ))
  
  cell_df$CHROM <- factor(cell_df$CHROM, levels = chrom_order)
  seg_df$chrom <- factor(seg_df$chrom, levels = chrom_order)
  
  chrom_sizes <- cell_df %>%
    count(CHROM, name = "actual_n_bins") %>%
    full_join(chrom_bin_overrides, by = "CHROM") %>%
    mutate(n_bins = ifelse(!is.na(force_n_bins), force_n_bins, actual_n_bins)) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(n_bins), default = 0))
  
  cell_df <- cell_df %>%
    left_join(chrom_sizes, by = "CHROM") %>%
    group_by(CHROM) %>%
    arrange(START) %>%
    mutate(
      bin_index = row_number(),
      bin_spacing = n_bins / n(),
      x = offset + (bin_index - 1) * bin_spacing
    ) %>%
    ungroup()
  
  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  x_breaks <- chrom_centers$center
  x_labels <- chrom_centers$CHROM
  
  bin_lookup_start <- cell_df %>% select(CHROM, START, x)
  bin_lookup_end <- cell_df %>% select(CHROM, END, x) %>% distinct()
  
  seg_df <- seg_df %>%
    left_join(bin_lookup_start, by = c("chrom" = "CHROM", "start" = "START")) %>%
    rename(x_start = x) %>%
    left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
    rename(x_end = x) %>%
    filter(!is.na(x_start), !is.na(x_end))
  
  # === CN LINE PLOTTING ===
  # Create stepwise line segments for CN
  cn_path_df <- seg_df %>%
    mutate(segment_id = row_number()) %>%
    pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
    arrange(segment_id, x) %>%
    mutate(y = CN_total)
  
  # === VERTICAL JUMPS ONLY IF CN CHANGES OR CHROM CHANGES ===
  jump_df <- seg_df %>%
    arrange(chrom, start) %>%
    mutate(
      next_CN = lead(CN_total),
      next_chrom = lead(chrom),
      next_x = lead(x_start)
    ) %>%
    filter(!is.na(next_CN)) %>%
    filter(
      chrom != next_chrom |
        abs(CN_total - next_CN) > 1e-3  # ignore tiny differences
    ) %>%
    transmute(
      x = x_end,
      xend = x,
      y = pmin(CN_total, max_y),
      yend = pmin(next_CN, max_y)
    )
  
  chrom_backgrounds <- chrom_sizes %>%
    mutate(
      xmin = offset,
      xmax = offset + n_bins,
      fill = as.factor(row_number() %% 2)
    )

  # Get x-axis limits for proper positioning
  x_min <- min(cell_df$x, na.rm = TRUE)
  x_max <- max(cell_df$x, na.rm = TRUE)

  if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
    n_squares <- nrow(cnv_presence)
    
    # Let's reserve about 5% of the x-axis width for all squares combined,
    # then divide that by number of squares to get spacing.
    total_space <- (x_max - x_min) * 0.1
    
    x_spacing <- total_space / max(n_squares - 1, 1)  # avoid division by zero
    
    cnv_presence <- cnv_presence %>%
      mutate(
        # Place squares starting from right (x_max - margin) going left
        x = x_max - (row_number() - 0.75) * x_spacing,
        y = max_y + 1.5
      )
  }
  
  cell_df <- cell_df %>% mutate(above_max_y=CN_signal>max_y)
  cell_df$CN_signal <- pmin(cell_df$CN_signal, max_y)
  
  # === CN PLOT ===
  p1 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = CN_signal, shape=factor(above_max_y)),
               color="#999933", size = 0.1, alpha = 0.9, stroke=0.1) +
    scale_shape_manual(values=c(16, 2), guide="none") + #(16, 2)
    geom_line(data = cn_path_df, aes(x = x, y = y, group = segment_id),
              color = "black", size = 0.6) +
    geom_segment(data = jump_df, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black") +
    {
      if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
        geom_point(data = cnv_presence, aes(x = x, y = y, color = category),
                   shape = 15, size = 3)  # bigger squares here
      } else {
        NULL
      }
    } +
    # Add title as text annotation at the left side
    annotate("text", x = x_min, y = max_y + 1.5, label = title, 
             hjust = 0, vjust = 0.5, size = 4.2) +
    scale_color_manual(values = category_colors) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, max_y + max_y_offset), xlim = c(x_min, x_max), clip = "off") +
    scale_y_continuous(breaks = c(0, 2, 4)) +
    labs(x = NULL, y = "CN") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size=10, color="black"),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 25, r = 5, b = 0, l = 5),  # add top margin to fit title and squares
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank(),
      plot.title = element_blank(),  # Remove default title since we're using annotate
      legend.position = "none"
    )
  
  # === BAF PLOT ===
  p2 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = pBAF), color = "#3399FF", size = 0.1/4, alpha = 0.9, stroke=0.1/4) +    
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Chromosome", y = "BAF") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size=10, color="black"),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank()
    )

  plot_output <- p1 / plot_spacer() / p2 + plot_layout(heights = c(1, 0.1, 1))

  if (!is.null(output_dir)){
    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.png")),
      plot_output,
      width = 10, height = 4, dpi = 300
    )
    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.pdf")),
      plot_output,
      width = 10, height = 4, device = cairo_pdf
    )
  }
  
  if (print_plot) {
    print(plot_output)
  }
  
  return(list(rdr=p1, baf=p2))
}



plot_cell2_supplementary2_qc_fail <- function(cell_name, title, cell_data, seg_merged_data, cnv_presence = NULL,
                       output_dir = NULL, print_plot = TRUE, chrom_list = c(1:22, "X", "Y")) {
  category_colors <- c(
    "Aneuploidy" = "#F0E442",  # yellow
    "CN-LOH"     = "#009E73",  # green
    "Deletion"   = "#56B4E9",  # blue
    "Duplication"= "#D55E00",   # red-orange
    "Failed QC" = "#000000" # black
  )

  max_y <- 4
  max_y_offset <- 0.25
    
  cell_df <- cell_data[cell_data$cell_name == cell_name & cell_data$CHROM %in% chrom_list, ]
  seg_df <- seg_merged_data[seg_merged_data$cell_name == cell_name & seg_merged_data$chrom %in% chrom_list, ]


  chrom_order <- as.character(chrom_list)
  
  chrom_bin_overrides <- tibble(CHROM = chrom_order) %>%
    mutate(force_n_bins = case_when(
      CHROM == "X" ~ 200,
      CHROM == "Y" ~ 100,
      TRUE ~ NA_real_
    ))
  
  cell_df$CHROM <- factor(cell_df$CHROM, levels = chrom_order)
  seg_df$chrom <- factor(seg_df$chrom, levels = chrom_order)
  
  chrom_sizes <- cell_df %>%
    count(CHROM, name = "actual_n_bins") %>%
    full_join(chrom_bin_overrides, by = "CHROM") %>%
    mutate(n_bins = ifelse(!is.na(force_n_bins), force_n_bins, actual_n_bins)) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(n_bins), default = 0))
  
  cell_df <- cell_df %>%
    left_join(chrom_sizes, by = "CHROM") %>%
    group_by(CHROM) %>%
    arrange(START) %>%
    mutate(
      bin_index = row_number(),
      bin_spacing = n_bins / n(),
      x = offset + (bin_index - 1) * bin_spacing
    ) %>%
    ungroup()

  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  x_breaks <- chrom_centers$center
  x_labels <- chrom_centers$CHROM
  
  bin_lookup_start <- cell_df %>% select(CHROM, START, x)
  bin_lookup_end <- cell_df %>% select(CHROM, END, x) %>% distinct()
  
  seg_df <- seg_df %>%
    left_join(bin_lookup_start, by = c("chrom" = "CHROM", "start" = "START")) %>%
    rename(x_start = x) %>%
    left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
    rename(x_end = x) %>%
    filter(!is.na(x_start), !is.na(x_end))
  
  # === CN LINE PLOTTING ===
  # Create stepwise line segments for CN
  cn_path_df <- seg_df %>%
    mutate(segment_id = row_number()) %>%
    pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
    arrange(segment_id, x) %>%
    mutate(y = CN_total)

  cn_path_df <- cn_path_df %>% mutate(dash_line=CN_total>max_y, y=pmin(y, max_y))
  
  # === VERTICAL JUMPS ONLY IF CN CHANGES OR CHROM CHANGES ===
  jump_df <- seg_df %>%
    arrange(chrom, start) %>%
    mutate(
      next_CN = lead(CN_total),
      next_chrom = lead(chrom),
      next_x = lead(x_start)
    ) %>%
    filter(!is.na(next_CN)) %>%
    filter(
      chrom != next_chrom |
        abs(CN_total - next_CN) > 1e-3  # ignore tiny differences
    ) %>%
    transmute(
      x = x_end,
      xend = x,
      y = pmin(CN_total, max_y),
      yend = pmin(next_CN, max_y)
    )
  
  chrom_backgrounds <- chrom_sizes %>%
    mutate(
      xmin = offset,
      xmax = offset + n_bins,
      fill = as.factor(row_number() %% 2)
    )

  # Get x-axis limits for proper positioning
  x_min <- min(cell_df$x, na.rm = TRUE)
  x_max <- max(cell_df$x, na.rm = TRUE)

  if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
    n_squares <- nrow(cnv_presence)
    
    # Let's reserve about 5% of the x-axis width for all squares combined,
    # then divide that by number of squares to get spacing.
    total_space <- (x_max - x_min) * 0.1
    
    x_spacing <- total_space / max(n_squares - 1, 1)  # avoid division by zero
    
    cnv_presence <- cnv_presence %>%
      mutate(
        # Place squares starting from right (x_max - margin) going left
        x = x_max - (row_number() - 0.75) * x_spacing,
        y = max_y + 1.5
      )
  }
  
  cell_df <- cell_df %>% mutate(above_max_y=CN_signal>max_y)
  cell_df$CN_signal <- pmin(cell_df$CN_signal, max_y)

  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  # Only label even-numbered chromosomes and sex chromosomes
  x_breaks <- chrom_centers$center[chrom_centers$CHROM %in% c(as.character(seq(2,22,2)), "X", "Y")]
  x_labels <- chrom_centers$CHROM[chrom_centers$CHROM %in% c(as.character(seq(2,22,2)), "X", "Y")]

    
  # === CN PLOT ===
  p1 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = CN_signal, shape=factor(above_max_y)),
               color="black", size = 0.75, alpha = 0.9, stroke=0.1) +
    scale_shape_manual(values=c(16, 2), guide="none") +
    {
      if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
        geom_point(data = cnv_presence, aes(x = x, y = y, color = category),
                   shape = 15, size = 3)
      } else {
        NULL
      }
    } +
    ggtitle(title) +
    scale_color_manual(values = category_colors) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, max_y + max_y_offset), xlim = c(x_min, x_max), clip = "off") +
    scale_y_continuous(breaks = c(0, 2, 4)) +
    labs(x = NULL, y = "CN") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size=35, color="black"),
      axis.text.y = element_text(size=30, color="black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=30, color="black"),
      plot.margin = margin(t = 25, r = 5, b = 0, l = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank(),
      plot.title = element_text(size=40, color="black", hjust=0.5),
      legend.position = "none"
    )
  
  # === BAF PLOT ===
  p2 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = pBAF), color = "#3399FF", size = 0.1/4, alpha = 0.9, stroke=0.1/4) +    
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Chromosome", y = "BAF") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size=10, color="black"),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank()
    )

  #plot_output <- p1 / plot_spacer() / p2 + plot_layout(heights = c(1, 0.1, 1))
  plot_output <- p1
    
  if (!is.null(output_dir)){
    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.png")),
      plot_output,
      width = 16, height = 4, dpi = 300
    )
    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.pdf")),
      plot_output,
      width = 16, height = 4, device = cairo_pdf
    )
  }
  
  if (print_plot) {
    print(plot_output)
  }
  
  return(list(rdr=p1, baf=p2))
}



plot_cell2_supplementary2_qc_fail_scan2 <- function(cell_name, title, cell_data, seg_merged_data, cnv_presence = NULL,
                       output_dir = NULL, print_plot = TRUE, chrom_list = c(1:22, "X", "Y")) {
  category_colors <- c(
    "Aneuploidy" = "#F0E442",  # yellow
    "CN-LOH"     = "#009E73",  # green
    "Deletion"   = "#56B4E9",  # blue
    "Duplication"= "#D55E00",   # red-orange
    "Failed QC" = "#000000" # black
  )

    
  cell_df <- cell_data[cell_data$cell_name == cell_name & cell_data$CHROM %in% chrom_list, ]
  seg_df <- seg_merged_data[seg_merged_data$cell_name == cell_name & seg_merged_data$chrom %in% chrom_list, ]

  max_y <- quantile(cell_df$CN_signal, 0.99)[["99%"]] * 1.1
  max_y_offset <- 0.25

  chrom_order <- as.character(chrom_list)
  
  chrom_bin_overrides <- tibble(CHROM = chrom_order) %>%
    mutate(force_n_bins = case_when(
      CHROM == "X" ~ 400,
      CHROM == "Y" ~ 200,
      TRUE ~ NA_real_
    ))
  
  cell_df$CHROM <- factor(cell_df$CHROM, levels = chrom_order)
  seg_df$chrom <- factor(seg_df$chrom, levels = chrom_order)
  
  chrom_sizes <- cell_df %>%
    count(CHROM, name = "actual_n_bins") %>%
    full_join(chrom_bin_overrides, by = "CHROM") %>%
    mutate(n_bins = ifelse(!is.na(force_n_bins), force_n_bins, actual_n_bins)) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(n_bins), default = 0))
  print(chrom_sizes)
    
  cell_df <- cell_df %>%
    left_join(chrom_sizes, by = "CHROM") %>%
    group_by(CHROM) %>%
    arrange(START) %>%
    mutate(
      bin_index = row_number(),
      bin_spacing = n_bins / n(),
      x = offset + (bin_index - 1) * bin_spacing
    ) %>%
    ungroup()

  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  x_breaks <- chrom_centers$center
  x_labels <- chrom_centers$CHROM
  
  bin_lookup_start <- cell_df %>% select(CHROM, START, x)
  bin_lookup_end <- cell_df %>% select(CHROM, END, x) %>% distinct()
  
  seg_df <- seg_df %>%
    left_join(bin_lookup_start, by = c("chrom" = "CHROM", "start" = "START")) %>%
    rename(x_start = x) %>%
    left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
    rename(x_end = x) %>%
    filter(!is.na(x_start), !is.na(x_end))
  
  # === CN LINE PLOTTING ===
  # Create stepwise line segments for CN
  cn_path_df <- seg_df %>%
    mutate(segment_id = row_number()) %>%
    pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
    arrange(segment_id, x) %>%
    mutate(y = CN_total)

  cn_path_df <- cn_path_df %>% mutate(dash_line=CN_total>max_y, y=pmin(y, max_y))
  
  # === VERTICAL JUMPS ONLY IF CN CHANGES OR CHROM CHANGES ===
  jump_df <- seg_df %>%
    arrange(chrom, start) %>%
    mutate(
      next_CN = lead(CN_total),
      next_chrom = lead(chrom),
      next_x = lead(x_start)
    ) %>%
    filter(!is.na(next_CN)) %>%
    filter(
      chrom != next_chrom |
        abs(CN_total - next_CN) > 1e-3  # ignore tiny differences
    ) %>%
    transmute(
      x = x_end,
      xend = x,
      y = pmin(CN_total, max_y),
      yend = pmin(next_CN, max_y)
    )
  
  chrom_backgrounds <- chrom_sizes %>%
    mutate(
      xmin = offset,
      xmax = offset + n_bins,
      fill = as.factor(row_number() %% 2)
    )

  # Get x-axis limits for proper positioning
  x_min <- min(cell_df$x, na.rm = TRUE)
  x_max <- max(cell_df$x, na.rm = TRUE)

  if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
    n_squares <- nrow(cnv_presence)
    
    # Let's reserve about 5% of the x-axis width for all squares combined,
    # then divide that by number of squares to get spacing.
    total_space <- (x_max - x_min) * 0.1
    
    x_spacing <- total_space / max(n_squares - 1, 1)  # avoid division by zero
    
    cnv_presence <- cnv_presence %>%
      mutate(
        # Place squares starting from right (x_max - margin) going left
        x = x_max - (row_number() - 0.75) * x_spacing,
        y = max_y + 1.5
      )
  }
  
  cell_df <- cell_df %>% mutate(above_max_y=CN_signal>max_y)
  cell_df$CN_signal <- pmin(cell_df$CN_signal, max_y)

  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  # Only label even-numbered chromosomes and sex chromosomes
  x_breaks <- chrom_centers$center[chrom_centers$CHROM %in% c(as.character(seq(2,22,2)),  "Y")]
  x_labels <- chrom_centers$CHROM[chrom_centers$CHROM %in% c(as.character(seq(2,22,2)), "Y")]

    
  # === CN PLOT ===
  p1 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = CN_signal, shape=factor(above_max_y)),
               color="black", size = 0.75, alpha = 0.9, stroke=0.1) +
    scale_shape_manual(values=c(16, 2), guide="none") +
    {
      if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
        geom_point(data = cnv_presence, aes(x = x, y = y, color = category),
                   shape = 15, size = 3)
      } else {
        NULL
      }
    } +
    ggtitle(title) +
    scale_color_manual(values = category_colors) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, max_y + max_y_offset), xlim = c(x_min, x_max), clip = "off") +
    labs(x = NULL, y = "Read count") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size=35, color="black"),
      axis.text.y = element_text(size=25, color="black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=25, color="black"),
      plot.margin = margin(t = 25, r = 5, b = 0, l = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank(),
      plot.title = element_text(size=40, color="black", hjust=0.5),
      legend.position = "none"
    )
  
  # === BAF PLOT ===
  p2 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = pBAF), color = "#3399FF", size = 0.1/4, alpha = 0.9, stroke=0.1/4) +    
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Chromosome", y = "BAF") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size=10, color="black"),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank()
    )

  #plot_output <- p1 / plot_spacer() / p2 + plot_layout(heights = c(1, 0.1, 1))
  plot_output <- p1
    
  if (!is.null(output_dir)){
    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.png")),
      plot_output,
      width = 16, height = 4, dpi = 300
    )
    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.pdf")),
      plot_output,
      width = 16, height = 4, device = cairo_pdf
    )
  }
  
  if (print_plot) {
    print(plot_output)
  }
  
  return(list(rdr=p1, baf=p2))
}


plot_phased_hets <- function(phased_hets, title, cyto_file = "refs/hg38_cytoBand.txt",
                             chrom_list = c(1:22, "X", "Y")) {
  library(dplyr)
  library(ggplot2)
  library(readr)

  # Read chromosome lengths from cytoBand.txt
  cyto <- read_tsv(cyto_file,
                   col_names = c("chrom", "start", "end", "name", "gieStain"))

  chrom_lengths <- cyto %>%
    group_by(chrom) %>%
    summarise(LENGTH = max(end), .groups = "drop") %>%
    mutate(CHROM = gsub("^chr", "", chrom)) %>%
    filter(CHROM %in% chrom_list) %>%
    arrange(factor(CHROM, levels = chrom_list)) %>%
    select(CHROM, LENGTH)

  # Compute cumulative offsets
  chrom_lengths <- chrom_lengths %>%
    mutate(offset = lag(cumsum(LENGTH), default = 0),
           midpoint = offset + LENGTH / 2)

  # Add adjusted x coordinate
  phased_hets <- phased_hets %>%
    inner_join(chrom_lengths, by = c("CHROM")) %>%
    mutate(x_coord = POS + offset)

  # Background rectangles for alternating chromosomes
  chrom_backgrounds <- chrom_lengths %>%
    mutate(fill = rep(c("white", "#D0D0D0"), length.out = n()))

  # Plot
  p <- ggplot() +
    # backgrounds
    geom_rect(data = chrom_backgrounds,
              aes(xmin = offset,
                  xmax = offset + LENGTH,
                  ymin = -Inf, ymax = Inf,
                  fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_identity() +

    # SNP points
    geom_point(data = phased_hets,
               aes(x = x_coord, y = TOTAL),
               size = 0.4, alpha = 0.6) +

    # x-axis labels: only even numbered chromosomes
    scale_x_continuous(
      breaks = chrom_lengths$midpoint[chrom_lengths$CHROM %in% as.character(seq(2,22,2))],
      labels = chrom_lengths$CHROM[chrom_lengths$CHROM %in% as.character(seq(2,22,2))],
      expand = c(0, 0)
    ) +
    ggtitle(title) +
    labs(x = "Chromosome", y = "TOTAL") +

    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold")
    )

  return(p)
}

plot_binned_phased_hets <- function(phased_hets, title, cyto_file = "refs/hg38_cytoBand.txt",
                                    chrom_list = c(1:22, "X", "Y"),
                                    binsize = 50000,
                                    yvar = c("TOTAL", "pBAF", "RDR_MEAN")) {
  library(dplyr)
  library(ggplot2)
  library(readr)

  yvar <- match.arg(yvar)

  # read cytoband to get chromosome lengths
  cyto <- read_tsv(cyto_file,
                   col_names = c("chrom", "start", "end", "name", "gieStain"))
  chrom_lengths <- cyto %>%
    group_by(chrom) %>%
    summarise(LENGTH = max(end), .groups = "drop") %>%
    mutate(CHROM = gsub("^chr", "", chrom)) %>%
    filter(CHROM %in% chrom_list) %>%
    arrange(factor(CHROM, levels = chrom_list)) %>%
    select(CHROM, LENGTH) %>%
    mutate(offset = lag(cumsum(LENGTH), default = 0),
           midpoint = offset + LENGTH/2)

  # bin phased hets
  phased_hets_binned <- phased_hets %>%
    filter(CHROM %in% chrom_list) %>%
    mutate(bin_index = floor(POS / binsize),
           bin_start = bin_index * binsize,
           bin_end = bin_start + binsize,
           bin_mid = (bin_start + bin_end) / 2) %>%
    group_by(CHROM, bin_index, bin_start, bin_end, bin_mid) %>%
    summarise(VALUE = mean(.data[[yvar]], na.rm = TRUE),
              .groups = "drop") %>%
    inner_join(chrom_lengths, by = c("CHROM")) %>%
    mutate(x_coord = bin_mid + offset)

  # backgrounds
  chrom_backgrounds <- chrom_lengths %>%
    mutate(fill = rep(c("white", "#D0D0D0"), length.out = n()))

  # plot
  ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = offset, xmax = offset + LENGTH,
                  ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_identity() +
    geom_point(data = phased_hets_binned,
               aes(x = x_coord, y = VALUE),
               size = 0.1, alpha = 0.8) +
    scale_x_continuous(
      breaks = chrom_lengths$midpoint[chrom_lengths$CHROM %in% as.character(seq(2,22,2))],
      labels = chrom_lengths$CHROM[chrom_lengths$CHROM %in% as.character(seq(2,22,2))],
      expand = c(0,0)
    ) +
    ggtitle(title) +
    labs(x = "Chromosome", y = yvar) +
    theme_classic(base_size = 16)
}



plot_cell2_truncated <- function(cell_name, title, cell_data, seg_merged_data, cnv_presence = NULL,
                       output_dir = NULL, print_plot = TRUE, chrom_list = c(1:22, "X", "Y")) {
  category_colors <- c(
    "Aneuploidy" = "#F0E442",  # yellow
    "CN-LOH"     = "#009E73",  # green
    "Deletion"   = "#56B4E9",  # blue
    "Duplication"= "#D55E00"   # red-orange
  )

  cell_type_colors <- c(
  "colon" = "#984ea3",  # Purple
  "lung"  = "#e31a1c"   # Red
)
  # Okabe-Ito palette (colorblind-safe colors)
  site_linetypes <- c("BCH-Broad" = "solid", 
                      "Mayo-WashU" = "twodash", 
                      "NYGC" = "blank",
                      "Yale-BCM" = "dotted", 
                      "Yonsei" = "longdash")

    
  max_y <- 4
  max_y_offset <- 0.25
    
  cell_df <- cell_data[cell_data$cell_name == cell_name & cell_data$CHROM %in% chrom_list, ]
  seg_df <- seg_merged_data[seg_merged_data$cell_name == cell_name & seg_merged_data$chrom %in% chrom_list, ]

  cell_type <- cell_df$cell_type[1]
  point_color <- cell_type_colors[[cell_type]]

  site_name <- cell_df$site_name[1]
  panel_border_linetype=site_linetypes[[site_name]]
    
  chrom_order <- as.character(chrom_list)
  
  chrom_bin_overrides <- tibble(CHROM = chrom_order) %>%
    mutate(force_n_bins = case_when(
      CHROM == "X" ~ 200,
      CHROM == "Y" ~ 100,
      TRUE ~ NA_real_
    ))
  
  cell_df$CHROM <- factor(cell_df$CHROM, levels = chrom_order)
  seg_df$chrom <- factor(seg_df$chrom, levels = chrom_order)
  
  chrom_sizes <- cell_df %>%
    count(CHROM, name = "actual_n_bins") %>%
    full_join(chrom_bin_overrides, by = "CHROM") %>%
    mutate(n_bins = ifelse(!is.na(force_n_bins), force_n_bins, actual_n_bins)) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(n_bins), default = 0))
  
  cell_df <- cell_df %>%
    left_join(chrom_sizes, by = "CHROM") %>%
    group_by(CHROM) %>%
    arrange(START) %>%
    mutate(
      bin_index = row_number(),
      bin_spacing = n_bins / n(),
      x = offset + (bin_index - 1) * bin_spacing
    ) %>%
    ungroup()
  
  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  x_breaks <- chrom_centers$center
  x_labels <- chrom_centers$CHROM
  
  bin_lookup_start <- cell_df %>% select(CHROM, START, x)
  bin_lookup_end <- cell_df %>% select(CHROM, END, x) %>% distinct()
  
  seg_df <- seg_df %>%
    left_join(bin_lookup_start, by = c("chrom" = "CHROM", "start" = "START")) %>%
    rename(x_start = x) %>%
    left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
    rename(x_end = x) %>%
    filter(!is.na(x_start), !is.na(x_end))
  
  # === CN LINE PLOTTING ===
  # Create stepwise line segments for CN
  cn_path_df <- seg_df %>%
    mutate(segment_id = row_number()) %>%
    pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
    arrange(segment_id, x) %>%
    mutate(y = CN_total)
  
  # === VERTICAL JUMPS ONLY IF CN CHANGES OR CHROM CHANGES ===
  jump_df <- seg_df %>%
    arrange(chrom, start) %>%
    mutate(
      next_CN = lead(CN_total),
      next_chrom = lead(chrom),
      next_x = lead(x_start)
    ) %>%
    filter(!is.na(next_CN)) %>%
    filter(
      chrom != next_chrom |
        abs(CN_total - next_CN) > 1e-3  # ignore tiny differences
    ) %>%
    transmute(
      x = x_end,
      xend = x,
      y = pmin(CN_total, max_y),
      yend = pmin(next_CN, max_y)
    )
  
  chrom_backgrounds <- chrom_sizes %>%
    mutate(
      xmin = offset,
      xmax = offset + n_bins,
      fill = as.factor(row_number() %% 2)
    )

    if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
      x_min <- min(cell_df$x, na.rm = TRUE)
      x_max <- max(cell_df$x, na.rm = TRUE)
      
      n_squares <- nrow(cnv_presence)
      
      # Let's reserve about 5% of the x-axis width for all squares combined,
      # then divide that by number of squares to get spacing.
      total_space <- (x_max - x_min) * 0.1
      
      x_spacing <- total_space / max(n_squares - 1, 1)  # avoid division by zero
      
      cnv_presence <- cnv_presence %>%
        mutate(
          # Place squares starting from right (x_max - margin) going left
          x = x_min + (row_number() - 0.75) * x_spacing,
          y = max_y 
        )
    }
  cell_df <- cell_df %>% mutate(above_max_y=CN_signal>max_y)
  cell_df$CN_signal <- pmin(cell_df$CN_signal, max_y)
  # === CN PLOT ===
  p1 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = CN_signal, shape=factor(above_max_y)),
               color=point_color, size = 0.3, alpha = 0.9) + #color = "#999933"
    scale_shape_manual(values=c(16, 2), guide="none") + #(16, 2)
    geom_line(data = cn_path_df, aes(x = x, y = y, group = segment_id),
              color = "black", size = 0.6) +
    geom_segment(data = jump_df, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black") +
    {
    if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
      geom_point(data = cnv_presence, aes(x = x, y = y, color = category),
                 shape = 15, size = 4)  # bigger squares here
    } else {
      NULL
    }
      } +
      scale_color_manual(values = category_colors) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
      coord_cartesian(ylim = c(0, max_y + max_y_offset), clip = "off") +  # extend y limit to 7 here
      scale_y_continuous(breaks = c(0, 2, 4)) +
      labs(x = NULL, y = "Copy number") +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 25, r = 5, b = 0, l = 5),  # add top margin to fit squares and title
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5, linetype=panel_border_linetype),
        axis.line = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none"
      ) +
      ggtitle(title)
  
 
  
  plot_output <- p1

  if (!is.null(output_dir)){
      ggsave(
        file.path(output_dir, paste0(cell_name, "_track.png")),
        plot_output,
        width = 10, height = 4, dpi = 300
      )
  }
  if (print_plot) {
    print(plot_output)
  }
  return(plot_output)
}


## BE CAREFUL with GenomicRanges, it will override some dplyr functions here and break function
## load GenomicRanges first
plot_cell2 <- function(cell_name, title, cell_data, seg_merged_data, cnv_presence = NULL,
                       output_dir = NULL, print_plot = TRUE, chrom_list = c(1:22, "X", "Y"), 
                       color_by_CN=FALSE, alt_seg_df=NULL) {
  category_colors <- c(
    "Aneuploidy" = "#F0E442",  # yellow
    "CN-LOH"     = "#009E73",  # green
    "Deletion"   = "#56B4E9",  # blue
    "Duplication"= "#D55E00"   # red-orange
  )
  
  cell_df <- cell_data[cell_data$cell_name == cell_name & cell_data$CHROM %in% chrom_list, ]
  seg_df <- seg_merged_data[seg_merged_data$cell_name == cell_name & seg_merged_data$chrom %in% chrom_list, ]
  
  chrom_order <- as.character(chrom_list)

  chrom_bin_overrides <- tibble(CHROM = chrom_order) %>%
    mutate(force_n_bins = case_when(
      CHROM == "X" ~ 200,
      CHROM == "Y" ~ 100,
      TRUE ~ NA_real_
    ))
  
  cell_df$CHROM <- factor(cell_df$CHROM, levels = chrom_order)
  seg_df$chrom <- factor(seg_df$chrom, levels = chrom_order)
  
  chrom_sizes <- cell_df %>%
    count(CHROM, name = "actual_n_bins") %>%
    full_join(chrom_bin_overrides, by = "CHROM") %>%
    mutate(n_bins = ifelse(!is.na(force_n_bins), force_n_bins, actual_n_bins)) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(n_bins), default = 0))
  
  cell_df <- cell_df %>%
    left_join(chrom_sizes, by = "CHROM") %>%
    group_by(CHROM) %>%
    arrange(START) %>%
    mutate(
      bin_index = row_number(),
      bin_spacing = n_bins / n(),
      x = offset + (bin_index - 1) * bin_spacing
    ) %>%
    ungroup()

    
  chrom_centers <- cell_df %>%
    group_by(CHROM) %>%
    summarise(center = mean(x)) %>%
    arrange(CHROM)
  
  x_breaks <- chrom_centers$center
  x_labels <- chrom_centers$CHROM
  
  bin_lookup_start <- cell_df %>% select(CHROM, START, x)
  bin_lookup_end <- cell_df %>% select(CHROM, END, x) %>% distinct()
  seg_df <- seg_df %>%
    left_join(bin_lookup_start, by = c("chrom" = "CHROM", "start" = "START")) %>%
    rename(x_start = x) %>%
    left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
    rename(x_end = x) %>%
    filter(!is.na(x_start), !is.na(x_end))
  
  # === CN LINE PLOTTING ===
  # Create stepwise line segments for CN
  cn_path_df <- seg_df %>%
    mutate(segment_id = row_number()) %>%
    pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
    arrange(segment_id, x) %>%
    mutate(y = CN_total)
  
  # === VERTICAL JUMPS ONLY IF CN CHANGES OR CHROM CHANGES ===
  jump_df <- seg_df %>%
    arrange(chrom, start) %>%
    mutate(
      next_CN = lead(CN_total),
      next_chrom = lead(chrom),
      next_x = lead(x_start)
    ) %>%
    filter(!is.na(next_CN)) %>%
    filter(
      chrom != next_chrom |
        abs(CN_total - next_CN) > 1e-3  # ignore tiny differences
    ) %>%
    transmute(
      x = x_end,
      xend = x,
      y = CN_total,
      yend = next_CN
    )
  
  chrom_backgrounds <- chrom_sizes %>%
    mutate(
      xmin = offset,
      xmax = offset + n_bins,
      fill = as.factor(row_number() %% 2)
    )
  
  # Compute max x and max y for positioning the CNV squares
    max_x <- max(cell_df$x, na.rm = TRUE)
    
    # Increase y limit for squares to be above title
    max_y <- 6
    
    if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
      x_min <- min(cell_df$x, na.rm = TRUE)
      x_max <- max(cell_df$x, na.rm = TRUE)
      
      n_squares <- nrow(cnv_presence)
      
      # Let's reserve about 5% of the x-axis width for all squares combined,
      # then divide that by number of squares to get spacing.
      total_space <- (x_max - x_min) * 0.05  # 5% of full width
      
      x_spacing <- 2 * total_space / max(n_squares - 1, 1)  # avoid division by zero
      
      cnv_presence <- cnv_presence %>%
        mutate(
          # Place squares starting from right (x_max - margin) going left
          x = x_max - 0.05 * (x_max - x_min) - (row_number() - 1) * x_spacing,
          y = max_y + 2
        )
    }

  # === If alt_seg_df provided, map to x coords ===
  alt_path_df <- NULL
  if (!is.null(alt_seg_df)) {
      alt_seg_df <- alt_seg_df %>%
        filter(cell_name == cell_name, chrom %in% chrom_list)
    
      # Use bin lookups for mapping genomic coords to x
      # bin_lookup_start has CHROM, START, x
      # bin_lookup_end   has CHROM, END, x
      bin_starts <- bin_lookup_start %>% arrange(START)
      bin_ends   <- bin_lookup_end   %>% arrange(END)
    
      map_to_x <- function(chrom, pos, lookup) {
        # filter lookup for chrom
        chr_bins <- lookup %>% filter(CHROM == chrom)
        if (nrow(chr_bins) == 0) return(NA_real_)
        # find interval (binary search)
        idx <- findInterval(pos, chr_bins[[2]])  # second column is START or END
        # clamp
        idx <- max(1, min(idx, nrow(chr_bins)))
        return(chr_bins$x[idx])
      }
    
      alt_seg_df <- alt_seg_df %>%
        rowwise() %>%
        mutate(
          x_start = map_to_x(chrom, start, bin_lookup_start),
          x_end   = map_to_x(chrom, end,   bin_lookup_end)
        ) %>%
        ungroup() %>%
        filter(!is.na(x_start), !is.na(x_end))
    
      alt_path_df <- alt_seg_df %>%
        mutate(segment_id = row_number()) %>%
        tidyr::pivot_longer(cols = c(x_start, x_end), names_to = "point", values_to = "x") %>%
        arrange(segment_id, x) %>%
        mutate(y = CN_total)
    }

  # === CN PLOT ===
  p1 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    geom_point(data = cell_df, aes(x = x, y = CN_signal),
               color = "#999933", size = 0.3, alpha = 0.9) +
    geom_line(data = cn_path_df, aes(x = x, y = y, group = segment_id),
              color = "black", size = 0.6) +
    geom_segment(data = jump_df, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black") +
    {
    if (!is.null(cnv_presence) && nrow(cnv_presence) > 0) {
      geom_point(data = cnv_presence, aes(x = x, y = y, color = category),
                 shape = 15, size = 8)  # bigger squares here
    } else {
      NULL
    }
      } +
    {
      if (!is.null(alt_path_df)) {
        geom_line(data = alt_path_df,
                  aes(x = x, y = y, group = segment_id),
                  color = "blue", size = 0.8, linetype = "dashed")
      } else NULL
    } +
      scale_color_manual(values = category_colors) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
      coord_cartesian(ylim = c(0, max_y + 1), clip = "off") +  # extend y limit to 7 here
      scale_y_continuous(breaks = 0:6) +
      labs(x = NULL, y = "Copy number") +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 25, r = 5, b = 0, l = 5),  # add top margin to fit squares and title
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.position = "none"
      ) +
      ggtitle(title, size)
  
  # === BAF PLOT ===
  p2 <- ggplot() +
    geom_rect(data = chrom_backgrounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    {
    if (!color_by_CN){
        geom_point(data = cell_df, aes(x = x, y = pBAF),
               color = "#3399FF", size = 0.3, alpha = 0.9)
        } else{
         geom_point(data = cell_df, aes(x = x, y = pBAF, color=factor(CN)),
                                        size = 0.3, alpha = 0.9)
        }
        } +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Chromosome", y = "B-allele freq") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank()
    )
  
  plot_output <- p1 / plot_spacer() / p2 + plot_layout(heights = c(1, 0.1, 1))

  if (!is.null(output_dir)){
      ggsave(
        file.path(output_dir, paste0(cell_name, "_track.png")),
        plot_output,
        width = 10, height = 4, dpi = 300
      )
  }
  if (print_plot) {
    print(plot_output)
  }
  return(plot_output)
}




plot_cell <- function(cell_name, title, cell_data, seg_merged_data, output_dir="cell_plots"){
    y_axis_title_fontsize = 16
    x_axis_title_fontsize = 16
    y_axis_label_fontsize = 12
    x_axis_label_fontsize = 10
    title_fontsize=20
    cell_df <- cell_data[cell_data$cell_name==cell_name & cell_data$CHROM %in% as.character(c(1:22)),]
    seg_df <- seg_merged_data[seg_merged_data$cell_name==cell_name & seg_merged_data$chrom %in% as.character(c(1:22)),]

    # Convert chromosomes to factor with natural ordering
    cell_df$CHROM <- factor(cell_df$CHROM, levels = as.character(1:22))
    seg_df$chrom <- factor(seg_df$chrom, levels = as.character(1:22))


    # Step 1: count bins per chromosome
    chrom_sizes <- cell_df %>%
      count(CHROM, name = "n_bins") %>%
      arrange(as.numeric(as.character(CHROM))) %>%
      mutate(offset = lag(cumsum(n_bins), default = 0))  # cumulative offset
    
    # Step 2: join offsets into main df and compute x
    cell_df <- cell_df %>%
      left_join(chrom_sizes, by = "CHROM") %>%
      group_by(CHROM) %>%
      arrange(START) %>%
      mutate(bin_index = row_number(),
             x = bin_index + offset) %>%
      ungroup()

    
    chrom_centers <- cell_df %>%
      group_by(CHROM) %>%
      summarise(center = mean(x)) %>%
      arrange(as.numeric(as.character(CHROM)))



    # Create a mapping of START positions to x-coordinates
    bin_lookup <- cell_df %>%
      select(CHROM, START, x)

    seg_df <- seg_df %>%
      left_join(bin_lookup, by = c("chrom" = "CHROM", "start" = "START")) %>%
      dplyr::rename(x_start = x)
    
    # Do the same for segment end
    bin_lookup_end <- cell_df %>%
      select(CHROM, END, x) %>%
      distinct()
    
    seg_df <- seg_df %>%
      left_join(bin_lookup_end, by = c("chrom" = "CHROM", "end" = "END")) %>%
      dplyr::rename(x_end = x)

   # Sort segments by chromosome and start position
    seg_df <- seg_df %>%
      arrange(as.numeric(as.character(chrom)), x_start)
    
    # Identify jumps between segments
    jump_df <- seg_df %>%
      mutate(next_CN = lead(CN_total),
             next_chrom = lead(chrom),
             next_x = lead(x_start)) %>%
      select(x = x_end, nextCN = next_CN, CN = CN_total)
 
    # Vertical chromosome boundaries
    vlines <- chrom_sizes$offset
    vlines <- vlines[2:length(vlines)] 
    
    # Chromosome tick labels
    x_breaks <- chrom_centers$center
    x_labels <- chrom_centers$CHROM

    # Prepare alternating background rectangles
    chrom_backgrounds <- chrom_sizes %>%
      mutate(xmin = offset,
             xmax = offset + n_bins,
             fill = as.factor((as.numeric(CHROM) %% 2)))  # alternating fill

    
    # CN plot (top)
   p1 <- ggplot() +
      # Alternating chromosome backgrounds, now darker
      geom_rect(data = chrom_backgrounds,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                alpha = 0.5, inherit.aes = FALSE) +
      scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +  # Darker gray
    
      # CN scatter + segments
      geom_point(data = cell_df, aes(x = x, y = CN_signal),
                 color = "#999933", size = 0.3, alpha = 0.9) +
      geom_segment(data = seg_df, aes(x = x_start, xend = x_end, y = CN_total, yend = CN_total),
                   color = "black", size = 0.6) +
    
      # Axis
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
      coord_cartesian(ylim = c(0, 10)) +
      scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
      geom_segment(data = jump_df, aes(x = x, xend=x, y=CN, yend=nextCN)) +
      labs(x = NULL, y = "CN") +
      ggtitle(title) +
      # Theme with full box
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_text(color = "black", size = y_axis_title_fontsize),
        axis.text.y = element_text(color = "black", size = y_axis_label_fontsize),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # Full box
        axis.line = element_blank(),
        plot.title = element_text(color = "black", size = title_fontsize, hjust = 0.5)
      )


    
    # BAF plot (bottom)
    p2 <- ggplot() +
      # Alternating chromosome backgrounds, same as above
      geom_rect(data = chrom_backgrounds,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                alpha = 0.5, inherit.aes = FALSE) +
      scale_fill_manual(values = c("white", "#D0D0D0"), guide = "none") +
    
      # BAF points
      geom_point(data = cell_df, aes(x = x, y = pBAF),
                 color = "#3399FF", size = 0.3, alpha = 0.9) +
    
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      labs(x = "Chromosome", y = "BAF") +
    
      # Theme with full box
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_text(color = "black", size = y_axis_title_fontsize),
        axis.text.y = element_text(color = "black", size = y_axis_label_fontsize),
        axis.title.x = element_text(color = "black", size = x_axis_title_fontsize),
        axis.text.x = element_text(color = "black", size = x_axis_label_fontsize),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # Full box
        axis.line = element_blank()
      )



    # Combine
    plot_output <- p1 / plot_spacer() / p2 + 
     plot_layout(heights = c(1, 0.1, 1))

    ggsave(
      file.path(output_dir, paste0(cell_name, "_track.png")),
      plot_output,
      width = 10, height = 4, dpi = 300
    )


    print(plot_output, width = 10, height = 4)

    }


compute_mapd <- function(rdr_vector) {
  log2_rdr <- log2(rdr_vector)
  diffs <- diff(log2_rdr)
  median(abs(diffs), na.rm = TRUE)
}

get_mapd_per_cell <- function(cell_data){
    # Assuming bins are ordered by CHROM and START
    mapd_per_cell <- cell_data %>%
      group_by(cell_name) %>%
      arrange(cell_name, CHROM, START) %>%  # Ensure bins are in genomic order
      summarize(MAPD = compute_mapd(RDR), .groups = "drop")
    mapd_per_cell
}

get_site_table <- function(batch_names){
    site_dictionary <- c("bcm"="Yale-BCM", "washu"="Mayo-WashU", "broad"="BCH-Broad", "yonsei"="Yonsei", "nygc"="NYGC")
    site_table <- data.table()  # initialize as data.table
    
    for (batch_name in batch_names) {
        message(paste("Reading site info for", batch_name))
        site_table_batch <- fread(file.path(prefix, batch_name, "site_metadata.csv"), header = TRUE)  # fread returns data.table
        site_table <- rbind(site_table, site_table_batch, use.names = TRUE, fill = TRUE)
    }
    
    site_table[, site_name := site_dictionary[site]]  # data.table style column creation
    site_table
}    

expand_df_with_cell_metadata <- function(cell_data){
    site_df <- as.data.frame(get_site_table(unique(cell_data$batch)))
    mapd_df <- get_mapd_per_cell(cell_data)
    
    cell_data <- merge(cell_data, site_df, by="cell_name")
    cell_data <- merge(cell_data, mapd_df, by="cell_name")
    return(cell_data)
    }

get_phased_hets_data_for_cell <- function(prefix, batch_name, binsize, lambda_str, sample, cell_name, chromosome="", region=c(start=0, end=Inf)){
    base_dir <- file.path(prefix, batch_name, "hiscanner_dsa_multisample-false", paste0("binsize_", binsize, "_lambda_", lambda_str))
    phased_hets_fname <- file.path(base_dir, sample, "hiscanner_output", "phased_hets", paste0(cell_name, ".hetsnp.txt"))
    phased_hets_df <- fread(phased_hets_fname)
    if (chromosome == ""){
        return (phased_hets_df)
        }
    return (phased_hets_df[phased_hets_df$CHROM == chromosome & phased_hets_df$POS >= region[["start"]] & phased_hets_df$POS <= region[["end"]]])
    }

get_phased_hets_data_for_cell_chr_list <- function(prefix, batch_name, binsize, lambda_str, sample, cell_name, chromosomes=NULL){
    base_dir <- file.path(prefix, batch_name, "hiscanner_dsa_multisample-false", paste0("binsize_", binsize, "_lambda_", lambda_str))
    phased_hets_fname <- file.path(base_dir, sample, "hiscanner_output", "phased_hets", paste0(cell_name, ".hetsnp.txt"))
    phased_hets_df <- fread(phased_hets_fname)
    if (is.null(chromosomes)){
        return (phased_hets_df)
        }
    return (phased_hets_df[phased_hets_df$CHROM %in% chromosomes])
    }


get_cell_data_for_cell <- function(prefix, batch_name, binsize, lambda_str, lambda, sample, cell_name, file_suffix=""){
    base_dir <- file.path(prefix, batch_name, "hiscanner_multisample-false", paste0("binsize_", binsize, "_lambda_", lambda_str))
    
    cell_fname <- file.path(base_dir, sample, "hiscanner_output", paste0("final_calls_lambda", lambda), paste0(cell_name, file_suffix, ".txt"))
    cell_data <- fread(cell_fname)
    cell_data$sample <- sample
    cell_data$cell_name <- cell_name
    cell_data$batch_name <- batch_name
    
    return(cell_data)
    }

    
get_cell_data <- function(prefix, batch_names, binsize, lambda_str, lambda, file_suffix="", hiscanner_dir_name="hiscanner_multisample-false"){
    # ----------------------------
    # Construct base directory and find *_input_table.txt files
    # ----------------------------
    all_data <- data.table()
    for (batch_name in batch_names){
        base_dir <- file.path(prefix, batch_name, hiscanner_dir_name, paste0("binsize_", binsize, "_lambda_", lambda_str))
        message(base_dir)
        
        input_table_files <- list.files(
          base_dir,
          pattern = "_input_table.txt$",
          full.names = TRUE,
          recursive = TRUE
        )
        
        # Restrict to final_calls_lambda{lambda}
        input_table_files <- input_table_files[grepl(paste0("final_calls_lambda", lambda), input_table_files)]
        message("Found ", length(input_table_files), " input_table.txt files")
        
        # ----------------------------
        # Extract matching {cell_name}.txt files and metadata
        # ----------------------------
        extract_data_from_cell_txt <- function(input_path) {
          # Get metadata from folder structure
          parts <- strsplit(input_path, "/")[[1]]
          sample <- parts[length(parts) - 3]
          
          # Get cell name from input file
          file_basename <- basename(input_path)
          prefix <- sub("_input_table.txt$", "", file_basename)
          cell_name <- prefix
        
          # Construct path to {cell_name}.txt file
          cell_txt_path <- file.path(dirname(input_path), paste0(prefix, file_suffix, ".txt"))
        
          # Skip if the file doesn’t exist
          if (!file.exists(cell_txt_path)) {
            warning("Missing .txt file for: ", cell_name)
            return(NULL)
          }
        
          # Read only the {cell_name}.txt file
          dt <- fread(cell_txt_path)
          dt[, `:=`(
            sample = sample,
            cell_name = cell_name
          )]
        
          return(dt)
        }
        
        # ----------------------------
        # Load all {cell_name}.txt files
        # ----------------------------
        message(batch_name)
        message("Reading cell_name.txt files...")
        all_data_batch <- rbindlist(lapply(input_table_files, extract_data_from_cell_txt), fill = TRUE)
        all_data_batch$batch_name <- batch_name
        all_data_batch$cell_type <- strsplit(batch_name, "_")[[1]][2]
        message("Loaded data for ", nrow(all_data_batch), " rows")
        all_data <- rbind(all_data, all_data_batch)
    }
    return(as.data.frame(all_data))
}



display_cell_image <- function(prefix, batch_name, binsize, lambda_str, lambda, sample, cell_name){
    base_dir <- file.path(prefix, batch_name, "hiscanner_multisample-false", paste0("binsize_", binsize, "_lambda_", lambda_str))
    
    track_path <- file.path(base_dir, sample, "hiscanner_output", paste0("final_calls_lambda", lambda), paste0(cell_name, "_track.png"))
    if (file.exists(track_path)) {
        #cat("\n####", cell_id, "- Track Image\n")
        text_html <- paste0(
        "<div style='font-size:18px; font-weight:bold;'>",
        cell_name,
        "</div>"
          )
        display_html(text_html)
        #display_text(paste(FANS_smple_ID, paste0("(MAPD: ", mapd, ";"), paste0("coverage: ", coverage, ")")))
        display_png(file = track_path, width=800, height=800)
      }
    }