source("packages.R")

load("cycle.thresholds.RData")

platesPlot <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ txt.path, labeller=function(var,val){
    gsub("/", "\n", val)
  })+
  geom_tile(aes(x=letter, y=number.fac, fill=Ct),
            data=cycle.thresholds)

pdf("figure-cycle-thresholds.pdf")
print(platesPlot)
dev.off()


