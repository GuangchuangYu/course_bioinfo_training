library("DiagrammeR")


grViz("digraph course {
rankdir = LR
node [shape = box, style=filled]
layout = dot
compound =true
#color = crimson

git [label='git和github']

subgraph clusterA{

label = 'Sequence'
style = dashed
rank = same
fasta [label='fasta文件解析']
genbank [label='genbank文件转fasta']
gb2 [label='批量下载genbank序列']
alignment [label='使用动态规划实验两序列比对']
msa [label='多重序列比对']
}

subgraph clusterB{

label = 'Database'
style = dashed
rank = same
sqlite [label='SQLite简介']
sqlite2 [label='SQLite实现课表查询']
}

subgraph clusterC{

label = 'R'
style = dashed
rank = same
vis [label='R语言数据可视化']
vis2 [label='R语言数据可视化']
rpkg [label='R包开发']
rpkg2 [label='R包开发']
}

subgraph clusterD{

label = 'NGS'
style = dashed
rank = same
rseq [label='R语言序列分析']
limma [label='RNA-seq数据分析差异基因']
gsea [label='GSEA分析']
}

subgraph clusterE{

label = 'Reproducible Research'
style = dashed
rank = same
make [label='用make搭建流程']
rmd [label='可重复性研究和Rmarkdown介绍']
}

subgraph clusterF{

label = 'Phylogeny'
style = dashed
rank = same
tree [label='如何阅读解析进化树']
treevis [label='构建进化树及可视化']
merse [label='MERS冠状病毒的来源']
}

sqlite -> sqlite2
genbank->fasta
fasta->alignment
fasta->msa
msa->treevis->tree->merse
rseq->limma->gsea
gb2-> genbank
make -> genbank [lhead=clusterA];
msa->rmd[ltail=clusterA]
msa->git[ltail=clusterA]
msa->rpkg[ltail=clusterA]
rpkg->rpkg2
vis->vis2
fasta->rseq
rpkg2->git
}") -> x

yyplot::gv2file(x,  file = 'diagram.pdf')
yyplot::gv2file(x, file = 'diagram.png' )
