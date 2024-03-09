# plot the simulation results

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(latex2exp)

###### Boxplots for estimators #######

tau_mc_1<-read.csv('tau_mc_1.csv')
tau_mc_0<-read.csv('tau_mc_0.csv')

df <- data.frame(
  strategy = c(rep('A',400),rep('B',400)),
  stage = rep(c('Stage 1','Stage 2','Stage 3','Stage 4'), times = 200),
  value = c(c(t(tau_mc_1)),c(t(tau_mc_0)))
)

# Define custom colors
custom_colors <- c("A" = "steelblue", "B" = "darkseagreen")

# Create a boxplot with facet_wrap and custom colors
p1<-ggplot(df, aes(x = strategy, y = value, fill = strategy)) +
  geom_boxplot() +
  facet_wrap(~ stage, ncol = 4) +
  labs(title = "Boxplots for Design-based Estimators", x = "Stage", y = "Treatment Effect Estimates") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(10, 10),
    axis.text = element_text(size = rel(1.1)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1)),
    plot.title = element_text(size = rel(1.25), hjust = 0.5),
    plot.margin = margin(0.5, 0.2, 0.5, 0.2, "cm")
  ) +
  theme_bw() +
  scale_fill_manual(values = custom_colors)

ggsave("comparison_1.pdf", p1, width = 4.5, height = 4, units = "in")

##### Sample size v.s. Design Strategy ########

df_max<-read.csv('df_max.csv')

p_max <- df_max %>%
  ggplot(aes(
    x = m,
    y = bias,
    color = Design,
    linetype = Design,
    shape = Design
  )) +
  geom_line(position = position_dodge(width = 0.1), size = 1) +
  geom_ribbon(aes(
    ymin = lower_ci,
    ymax = upper_ci,
    color = Design,
    fill = Design
  ),
  alpha = 0.3) +
  geom_point(size = 4, position = position_dodge(width = 0.1)) +
  scale_color_manual(values = c("steelblue", "darkseagreen", "tan1")) +
  scale_fill_manual(values = c("steelblue", "darkseagreen", "tan1")) +
  scale_shape_manual(values = c(15, 16, 17, 4, 7, 5)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.05, decimal.mark = '.')) +
  ylab('Propensity score / Subgroup Proportion') +
  xlab("Sample Size") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(10, 10),
    axis.text = element_text(size = rel(1.1)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1)),
    plot.title = element_text(size = rel(1.25), hjust = 0.5),
    plot.margin = margin(0.5, 0.2, 0.5, 0.2, "cm")
  ) +
  theme_bw() +
  ggtitle("Sample Size vs. Design Strategy") +
  theme(
    plot.title = element_text(size = rel(1.1), hjust = 0.5),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2))
  ) 
  
  # # Add horizontal lines with labels
  # + geom_hline(
  #   yintercept = c(p.star[1], e.star[1]),
  #   linetype = "dashed",
  #   color = "tan1",
  #   size = 1
  # )

ggsave("sample_size.pdf", p_max, width = 4.5, height = 4, units = "in")


##### Design strategies comparison between Problems 2 and 3 

design1_star<-read.csv("design1_star.csv", header = TRUE)
design1_star<-as.vector(design1_star$x)
design2_star<-read.csv("design2_star.csv", header = TRUE)
design2_star<-as.vector(design2_star$x)

df2 <- data.frame(
  strategy = c(rep('A',4),rep('C',4)),
  stage = rep(c('Stratum 1','Stratum 2','Stratum 3','Stratum 4'), times = 2),
  value = c(design1_star[5:8],design2_star[5:8])
)

custom_colors <- c("A" = "steelblue", "C" = "tan1")

p2<-ggplot(df2, aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~ stage, ncol = 4) +
  labs(title = "Treatment  Allocation Comparison", x = "Stratum", y = "Propensity Score") +
  theme_minimal() +
  theme_bw() +
  scale_fill_manual(values = custom_colors)+
  ylim(0, 0.6)

ggsave("treatment_allocation.pdf", p2, width = 4.5, height = 4, units = "in")


df3 <- data.frame(
  strategy = c(rep('A',4),rep('C',4)),
  stage = rep(c('Stratum 1','Stratum 2','Stratum 3','Stratum 4'), times = 2),
  value = c(design1_star[1:4],design2_star[1:4])
)

p3<-ggplot(df3, aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~ stage, ncol = 4) +
  labs(title = "Subgroup Proportion Comparison", x = "Stratum", y = "Subgroup Proportion") +
  theme_minimal() +
  theme_bw() +
  scale_fill_manual(values = custom_colors)+
  ylim(0, 0.6)

ggsave("subgroup_proportion.pdf", p3, width = 4.5, height = 4, units = "in")

