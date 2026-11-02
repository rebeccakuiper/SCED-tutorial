################################################################################
# Plot of the person-specific GORICA weights (run after the ENET_Pk scripts)
#
# Works with restriktor 0.6-x (CRAN) and 0.7+: weights are selected by column
# name ('gorica.weights'), not by position.
################################################################################

library(ggplot2)

# Names of the persons, in the order shown on the x-axis
# (objects output_P1, output_P3, ... must exist in the workspace)
P_names <- c("P1", "P3", "P4", "P5", "P6", "P7", "P9", "P10", "P11")

# Optional: order persons by the support for one hypothesis, e.g.
# order_by <- "unconstrained"
# (this gives P9, P3, P11, P1, P4, P5, P6, P10, P7 for the paper's results).
order_by <- NULL

# Label for y-axis
y_label <- "GORICA weights"


## Collect the weights -------------------------------------------------------------
outputs <- mget(paste0("output_", P_names), envir = globalenv(), ifnotfound = list(NULL))
missing_P <- P_names[vapply(outputs, is.null, logical(1))]
if (length(missing_P) > 0) {
  stop("Run the analysis scripts first; not found: ",
       paste0("output_", missing_P, collapse = ", "))
}

plot_data <- do.call(rbind, Map(function(out, p) {
  res <- as.data.frame(out$result)
  data.frame(Person = p, Hypothesis = res$model, weight = res[["gorica.weights"]])
}, outputs, P_names))
rownames(plot_data) <- NULL

# Same set of hypotheses for everyone?
namesH <- unique(plot_data$Hypothesis)
if (any(table(plot_data$Person) != length(namesH))) {
  warning("Not all persons have the same set of hypotheses.")
}

if (!is.null(order_by)) {
  w <- plot_data[plot_data$Hypothesis == order_by, ]
  P_names <- w$Person[order(w$weight)]
}
plot_data$Person     <- factor(plot_data$Person, levels = P_names)
plot_data$Hypothesis <- factor(plot_data$Hypothesis, levels = namesH)


## Plot -----------------------------------------------------------------------------
angle_x <- 30
vjust_x <- 1 # 0: bottom-aligned to 1: top-aligned
hjust_x <- 1 # right-aligned

p_weights <- ggplot(plot_data, aes(x = Person, y = weight,
                                   color = Hypothesis, shape = Hypothesis)) +
  geom_point(size = 3) +   # colour AND shape, so it also works in black-and-white
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
  scale_x_discrete(expand = c(0, 0.05)) +
  labs(x = "Person",
       y = y_label,
       title = "Person-specific GORICA weights",
       subtitle = "for all trajectory hypotheses") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(unit(0, "cm")),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.text = element_text(size = 12),
    axis.text.x  = element_text(size = 12, angle = angle_x, vjust = vjust_x, hjust = hjust_x),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14, vjust = -3),
    axis.title.y = element_text(size = 14, vjust = 5),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )
p_weights

# To save (e.g., for the paper):
# ggsave("Figure5_GORICA_weights.png", p_weights, width = 8, height = 5, dpi = 300)
