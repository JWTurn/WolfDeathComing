cc.betas <- readRDS('data/derived-data/full_cc_betas.Rds')


cc.betas[, cor(case, control), .(term)]

cc.betas.long <- melt(cc.betas[,.(wolfID, term, control, case)])

cc.betas.long[,ICC::ICCest(wolfID, value), .(term)]

b.cor <- ggplot(cc.betas, aes(control, (case), fill = term)) +
  geom_smooth(aes(fill = term),# notch = TRUE, notchwidth = 0.7,
               outlier.color = NA, lwd = 0.6,
               alpha = 0.25) +
  geom_jitter(aes(color = term),
              position = position_jitterdodge(.35),
              size = 2, alpha = 0.4) +
  #ggtitle('Interaction with community identity') +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme(#legend.position = 'none',
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title=element_text(size = 16, hjust=0),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="black", size = 1, fill = "white"),
    strip.text = element_text(size = 14)) +
  xlab('control') +
  ylab('case') +
  ggtitle("Correlation") #+
  #scale_fill_manual(values = cbPalette) +
  #scale_color_manual(values = cbPalette)# + ylim(-2,2)
