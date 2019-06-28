library(ggplot2)
library(grid)
library(gridExtra)

args <- commandArgs(TRUE)
in_file <- args[1]
out_file <- args[2]
coord_file <- args[3]

if (length(args) == 3) {
    axes11 = 'PC1'
    axes12 = 'PC2'
    axes21 = 'PC1'
    axes22 = 'PC3'
} else if (length(args) == 4) {
    axes11 = 'PC1'
    axes12 = 'PC2'
    axes2 = lapply(strsplit(args[4], ',')[[1]], function(x) paste0('PC', x))
    axes21 = axes2[[1]]
    axes22 = axes2[[2]]
} else {
    axes = strsplit(args[4:5], ',')
    axes1 = lapply(axes[[1]], function(x) paste0('PC', x))
    axes2 = lapply(axes[[2]], function(x) paste0('PC', x))
    axes11 = axes1[[1]]
    axes12 = axes1[[2]]
    axes21 = axes2[[1]]
    axes22 = axes2[[2]]
}




hgdp_data <- read.table(coord_file, header=T, sep='\t')
sample_data <- read.table(in_file, header=T, sep='\t', colClasses=c(rep('factor', 2), rep('numeric', 13)))
nsamples <- dim(sample_data)[[1]]

merged <- do.call(rbind, lapply(list(hgdp_data, sample_data), function (x) subset(x, select=intersect(names(hgdp_data),names(sample_data)))))


refhues <- hcl(h=seq(15, 375, length=8), l=75, c=100)[1:7]
samplehues <- hcl(h=seq(15, 375, length=nsamples+1), l=30, c=100)[1:nsamples]
myhues <- c(refhues, samplehues)
myshapes <- c(rep(16, 7), rep(4, nsamples))
mysizes <- c(rep(1.5, 7), rep(4, nsamples))

make_plot <- function(x, y)
    return(ggplot(aes_string(x=x, y=y, col='popID', shape='popID', size='popID'), data=merged)
           + geom_point()
           + theme(axis.ticks=element_blank(), axis.text=element_blank())
           + scale_color_manual(limits=levels(merged$popID), values=myhues, name='Population')
           + scale_shape_manual(limits=levels(merged$popID), values=myshapes, name='Population')
           + scale_size_manual(limits=levels(merged$popID), values=mysizes, name='Population')
           + guides(color=guide_legend(override.aes=list(size=4), title.position='top', title.hjust=0.5)))

p1 <- make_plot(axes11, axes12)
p2 <- make_plot(axes21, axes22)

png(file=out_file, width=1000, height=525)

g <- ggplotGrob(p1 + theme(legend.position='bottom'))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == 'guide-box')]]
lheight <- sum(legend$height)
args <- lapply(list(p1, p2), function(x) x + theme(legend.position='none'))
args$ncol <- 2
grid.arrange(do.call(arrangeGrob, args), legend, ncol=1, heights=unit.c(unit(1, 'npc') - lheight, lheight))

dev.off()
