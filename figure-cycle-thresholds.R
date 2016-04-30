source("packages.R")

load("cycle.thresholds.RData")

TFAC <- function(Treatment){
  factor(
    ifelse(Treatment=="no treatment", "none", Treatment),
    c("none", "HT-DNA", "LPS", "Pam3CSK4"))
  factor(Treatment, c("no treatment","HT-DNA", "LPS", "Pam3CSK4"))
}
platesPlot <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ txt.path, labeller=function(var,val){
    gsub("/", "\n", val)
  })+
  scale_linetype_manual(values=c("1"="solid", "2"="dashed"))+
  geom_tile(aes(number.fac,
                factor(letter.fac, rev(levels(letter.fac))),
                fill=Ct,
                color=factor(`Biological replicate`),
                linetype=factor(`Technical replicate`)
                ),
            size=1,
            data=cycle.thresholds)+
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ txt.path, labeller=function(var,val){
    gsub("/", "\n", val)
  })+
  geom_tile(aes(number.fac,
                factor(letter.fac, rev(levels(letter.fac))),
                fill=Time,
                color=Treatment
                ),
            size=1,
            data=cycle.thresholds)+
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  scale_color_brewer(palette="Set1")

cycle.thresholds[, Time.fac := factor(Time, c(0, 0.5, 1, 2))]
cycle.thresholds[, Treatment.fac := TFAC(Treatment)]
replicate.means <- cycle.thresholds[, list(
  mean.relative.expression=mean(relative.expression),
  replicates=.N),
  by=.(Genotype, Time, Time.fac, Treatment.fac, Gene)]
KO <- replicate.means[Genotype=="KO",]
KO$WT <- replicate.means[Genotype=="WT",]$mean.relative.expression

dput(RColorBrewer::brewer.pal(Inf, "Blues"))
c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", 
  "#2171B5", "#08519C", "#08306B")

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(Gene ~ Treatment.fac, scales="free", space="free")+
  geom_tile(aes(Time.fac, Genotype, fill=mean.relative.expression),
            data=replicate.means)+
  scale_fill_gradient(low="#F7FBFF", high="#08306B")

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(Gene ~ Treatment.fac, scales="free", space="free")+
  geom_line(aes(Time.fac, relative.expression,
                color=Genotype,
                group=paste(`Technical replicate`, Genotype)),
            data=cycle.thresholds)+
  geom_point(aes(Time.fac, relative.expression,
                 color=Genotype),
             data=cycle.thresholds)

no.treatment <- cycle.thresholds[Treatment == "no treatment",]
treatments <- cycle.thresholds[Treatment != "no treatment",]
zero.treatments.list <- list(treatments)
for(Treatment in unique(treatments$Treatment)){
  no.treatment$Treatment <- Treatment
  zero.treatments.list[[Treatment]] <- no.treatment
}
zero.treatments <- do.call(rbind, zero.treatments.list)
zero.treatments[, Treatment.fac := TFAC(Treatment)]

genes.ordered <-
  KO[Treatment.fac=="no treatment",][order(mean.relative.expression-WT),]
KO[, Gene.fac := factor(Gene, genes.ordered$Gene)]
replicate.means[, Gene.fac := factor(Gene, genes.ordered$Gene)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ Treatment.fac, scales="free", space="free")+
  geom_tile(aes(Time.fac, Gene.fac, fill=mean.relative.expression-WT),
            data=KO)+
  scale_fill_gradient2()

viz <- list(
  title="H3F3B RNA mouse data",
  overview=ggplot()+
    ggtitle("Changes in RNA in KO relative to WT")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(. ~ Treatment.fac, scales="free")+
    geom_tile(aes(Time.fac, Gene.fac, fill=mean.relative.expression-WT,
                  clickSelects=Gene),
              data=KO)+
    xlab("hours after treatment")+
    ylab("Gene")+
    scale_fill_gradient2(),
  curves=ggplot()+
    ggtitle("Two technical replicates")+
    theme_bw()+
    xlab("hours after treatment")+
    ylab("relative expression (Ct difference from Gapdh)")+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(. ~ Treatment.fac, scales="free")+
    geom_line(aes(Time.fac, relative.expression,
                  color=Genotype,
                  showSelected=Gene,
                  group=paste(`Technical replicate`, Genotype)),
              data=zero.treatments),
  wtko=ggplot()+
    theme_bw()+
    xlab("hours after treatment")+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(height=700)+
    facet_grid(Gene.fac ~ Treatment.fac, scales="free")+
    geom_tile(aes(Time.fac, Genotype,
                  clickSelects=Gene,
                  fill=mean.relative.expression),
              data=replicate.means)+
    scale_fill_gradient(low="#F7FBFF", high="#08306B"),
  duration=list(Gene=2000))
animint2dir(viz, "figure-cycle-thresholds")

gg.curves <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(Gene ~ Treatment.fac, scales="free", space="free")+
  geom_line(aes(Time.fac, relative.expression,
                color=Genotype,
                group=paste(`Technical replicate`, Genotype)),
            data=zero.treatments)+
  ylab("relative expression (Ct difference from Gapdh)")+
  xlab("hours after treatment")

pdf("figure-cycle-thresholds-curves.pdf", h=8)
print(gg.curves)
dev.off()

pdf("figure-cycle-thresholds-diff.pdf")
print(viz$overview)
dev.off()

pdf("figure-cycle-thresholds.pdf", h=8)
print(viz$wtko)
dev.off()


