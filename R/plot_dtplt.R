#' @title plot_dtplt
#'
#' @description Generates multi-feature plot(s).
#'
#' @param dtplt_data A data frame with the following columns: Gene, Grp, Split, Avg, Avg.SC, Pct, Pct.SC
#' @param split A logical value indicating the data is split by gene
#'
#' @return A ggplot object or a list of ggplot objects
#'
#' @export
#'
plot_dtplt <- function(dtplt_data, split = FALSE) {
  if (split) {
    plots <- list()
    for (gene in unique(dtplt_data$Gene)) {
      plots[[gene]] <- ggplot(dt_data_split[ which(dt_data_split$Gene==gene), ], aes(x = Split, y = Grp)) +
        geom_point(aes(size = Pct.SC, colour = Avg.SC)) + coord_flip() +
        labs(title = str_to_title(gene)) +
        scale_color_gradient(name = "Avg Exp", low = "lightgray", high = "blue"#, labels = NULL#, breaks = c(-0.4, 2.2), labels = c('lo', 'hi')
        ) +
        guides(colour = guide_colorbar(order=1, ticks.colour = NA, barheight = 1, barwidth = 4.5,
                                       title.position="top", title.hjust = 0.5),
               size = guide_legend(order=2, reverse=F, title.position="top", title.hjust = 0.5,
                                   label.position = "bottom", label.vjust = unit(-0.3, "points"))) +
        scale_size_continuous(name = "% Exp"#, labels = NULL#, breaks = c(0, 0.75, 1.25, 2), labels = c('lo', '', '', 'hi')
        ) +
        theme(plot.title = element_text(face = 'bold', size = 18, hjust = 0.5),
              panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
              panel.background = element_rect(fill = 'white'),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 0),
              axis.text.y = element_text(size = 12, colour = "black"), #, face = "italic"
              axis.title = element_blank(),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 11.5),
              legend.key = element_rect(fill = NA),
              legend.spacing = unit(0.5, "points"),
              legend.box.spacing = unit(0.75, "points"),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.spacing.x = unit(0.1, 'points')
        ) + scale_x_discrete(position = "bottom") +
        scale_y_discrete(position = "right")
    }
    return (plots)
  } else {
    plot <- ggplot(dt_data, aes(x = Gene, y = Grp)) +
      geom_point(aes(size = Pct.SC, colour = Avg.SC)) + coord_flip() +
      #labs(title = str_to_title(gene)) +
      scale_color_gradient(name = "Avg Exp", low = "lightgray", high = "blue"#, labels = NULL#, breaks = c(-0.4, 2.2), labels = c('lo', 'hi')
      ) +
      guides(colour = guide_colorbar(order=1, ticks.colour = NA, barheight = 1, barwidth = 4.5,
                                     title.position="top", title.hjust = 0.5),
             size = guide_legend(order=2, reverse=F, title.position="top", title.hjust = 0.5,
                                 label.position = "bottom")) + #, label.vjust = unit(-0.3, "points"))) +
      scale_size_continuous(name = "% Exp"#, labels = NULL#, breaks = c(-1, -0.33, 0.33, 1), labels = c('lo', '', '', 'hi')
      ) +
      theme(plot.title = element_text(face = 'bold', size = 18, hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
            panel.background = element_rect(fill = 'white'),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            #axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 0),
            axis.text.x = element_text(size = 12, colour = "black", angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 12, colour = "black"), #, face = "italic"
            axis.title = element_blank(),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11.5),
            legend.key = element_rect(fill = NA),
            legend.spacing = unit(0.5, "points"),
            legend.box.spacing = unit(0.75, "points"),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.spacing.x = unit(0.1, 'points')
      ) + scale_x_discrete(position = "bottom") +
      scale_y_discrete(position = "right") #labels=function(x){sub("\\s", "\n", x)}
  return (plot)
  }
}
