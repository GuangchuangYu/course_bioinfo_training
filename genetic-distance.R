library(ggplot2)
library(ggtree)


set.seed(2019)
x = rtree(10)
x$tip.label = sub("t", "virus", x$tip.label)
x$edge.length = round(x$edge.length, 2)

p <- ggtree(x) + geom_tiplab(align = TRUE, size =5) 

n = c(2, 3, 14, 15)
d2 = subset(p$data, p$data$node %in% n)

d2$xend = p$data$x[d2$parent]
p2 <- p + 
    geom_segment(aes(x=x,xend=xend, y=y, yend=y), 
    data=d2, color='firebrick', size=2, lineend='round') +
    xlim(NA, 3.05) + 
    geom_label2(aes(x = branch, label=branch.length,
                    subset = (node != 11)
                    )
                )

ggsave(p2, file="genetic-distance.png",
    width = 8, height = 6)


 
