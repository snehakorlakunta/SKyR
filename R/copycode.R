#' @title copycode
#'
#' @description Adds function code to clipboard to allow pasting with formatting.
#'
#' @return Adds selected code to clipboard.
#'
#' @export
#'
copycode <- function() {
  select <- menu(c(
    "Assign cell cluster labels",
    "Get cell counts",
    "Get proportion plot"),
    title = "Which code do you want to copy?")
  if (select == 0) {stop()}

  # Functions #####
  functions = list()
  fxn_ct <- 0

  # formatting for new function:
  #fxn_ct <- fxn_ct + 1
  #functions[fxn_ct] <-
  #  ""

  ## Assign cell cluster labels ################################################
  fxn_ct <- fxn_ct + 1
  functions[fxn_ct] <-
    "cluster_labels <- c(\"\",
                    \"\",
                    \"\")

names(cluster_labels) <- levels(seurat_object)\nseurat_object <- RenameIdents(seurat_object, cluster_labels)\nDimPlot(seurat_object, reduction = \"umap.harmony\", label = TRUE) + NoLegend()

seurat_object$CellType <- Idents(seurat_object)\nseurat_object$CellType <- factor(seurat_object$CellType)
    "

  ## Get cell counts ###########################################################
  fxn_ct <- fxn_ct + 1
  functions[fxn_ct] <-
    "library(data.table)\n\nmd <- seurat_object@meta.data %>% as.data.table\ncounts <- md[, .N, by = c(\"orig.ident\", \"CellType\")]\n\n# with additional casting after the counting\n#md[, .N, by = c(\"orig.ident\", \"CellType\")] %>% dcast(., orig.ident ~ CellType, value.var = \"N\")\n\n# to order cells\ncounts <- counts[order(counts[[order]], decreasing = decr_order), ]"

  ## Get proportion plot #######################################################
  fxn_ct <- fxn_ct + 1
  functions[fxn_ct] <-
    "seurat_object@meta.data %>%
      #filter(Condition %in% c(*CONDITIONS*)) %>%
      #filter(CellType %in% c(*CELLTYPES*)) %>%
      group_by(Condition, CellType) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(Condition) %>%
      mutate(total = sum(count)) %>%
      ungroup() %>%
      mutate(percentage = ((count / total)*100)) %>%
      #ggplot(aes(x=Condition, y=percentage, fill=CellType)) + #, label=count)) +
      ggplot(aes(x=Condition, y=percentage, fill=CellType, label=sprintf(\"%0.2f\", round(percentage, digits = 2)))) +
      geom_bar(stat=\"identity\", width=0.9, colour=\"white\") + #NoLegend() +
      geom_text(size = 3, position = position_stack(vjust = 0.5)) +
      scale_y_continuous(limits = c(-5,105), expand = c(0, 0)) +
      scale_fill_discrete(name = \"Cell Type\")+#, limits = c()) +
      labs(y = \"Percentage of cell types\") +
      theme(plot.title = element_text(size = 10, face = \"bold\"),
            axis.ticks.x = element_blank(),
            axis.line = element_line(color = \"black\"),
            panel.background = element_rect(fill = \"white\"),
            panel.grid = element_blank())
      "

  # Write to clipboard #####
  utils::writeClipboard(toString(functions[select]))
}

