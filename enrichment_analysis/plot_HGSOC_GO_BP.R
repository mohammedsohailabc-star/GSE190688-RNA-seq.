library(ggplot2)

# Input file
input_file <- "~/GSE190688/enrichment_analysis/HGSOC_up_GO_BP.csv"

# Output folder
output_dir <- "~/GSE190688/enrichment_analysis/GO_plots"

# Read GO results
go <- read.csv(input_file)

# Convert GeneRatio from fraction to decimal
go$GeneRatio_num <- sapply(strsplit(go$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

# Sort by adjusted p-value
go <- go[order(go$p.adjust), ]

# Function to create a dot plot
make_plot <- function(data, title, filename) {

  data$Description <- factor(
    data$Description,
    levels = rev(data$Description)
  )

  p <- ggplot(data, aes(
    x = GeneRatio_num,
    y = Description,
    size = Count,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = title,
      x = "Gene Ratio",
      y = "Biological Process",
      size = "Gene Count",
      color = "Adjusted p-value"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    plot = p, width = 12, height = 9, dpi = 300
  )

  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot = p, width = 12, height = 9
  )
}

# Plot top 20 terms
top20 <- head(go, 20)

make_plot(
  top20,
  "HGSOC-upregulated: Top 20 GO Biological Processes",
  "HGSOC_GO_BP_Top20"
)

# Plot all enriched terms
make_plot(
  go,
  "HGSOC-upregulated: All GO Biological Processes",
  "HGSOC_GO_BP_All"
)

cat("Both GO plots saved successfully!\n")
