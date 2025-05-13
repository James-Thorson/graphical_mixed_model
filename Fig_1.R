
root_dir = R'(C:\Users\James.Thorson\Desktop\Git\GGMM)'

#
Date = Sys.Date()
run_dir = file.path(root_dir, Date )
  dir.create(run_dir)

library(igraph)
library(Matrix)
library(ape)

plot_graph<-
function( graph ){           # rev(seq_len(ncol(A)))
  A = (as_adjacency_matrix(graph))   # [,rev(seq_len(nrow(A)))]
  image( as.matrix(A)[,rev(seq_len(ncol(A)))], x = seq_len(nrow(A)), y = seq_len(ncol(A)),
         col = c("white","grey"), xaxt = "n", yaxt = "n",
         xlab = "", ylab = "" )
  axis(3, at = axTicks(3), label = names(V(graph)), las = 2)
  axis(2, at = axTicks(2), label = rev(names(V(graph))), las = 2)
  box()
}

g1 <- graph_from_literal(
  "t" -+ "t+1",
  "t+1" -+ "t+2",
  "t+2" -+ "t+3"
)
coords1 = layout_(g1, with_sugiyama())
coords1$layout = cbind( 1:4, 4:1 )

g2 <- graph_from_literal(
  "x,y" -+ "x,y+1",
  "x,y" -+ "x+1,y",
  "x,y" -+ "x,y-1",
  "x,y" -+ "x-1,y"
)
coords2 = layout_(g2, with_sugiyama())
coords2$layout = cbind( x = c(0, 0, 1, 0, -1), y = c(0, 1, 0, -1, 0) )

# Tree
set.seed(101)
tree = rtree(4)
tree$tip.label = gsub( tree$tip.label, pattern = "t", replacement = "s")
tree$node.label = paste0("s", Ntip(tree) + seq_len(Nnode(tree)) )
# Rescale edge lengths
height = node.depth.edgelength(tree)[seq_len(Ntip(tree))]
extra_height = max(height) - height
match_tip_to_edge = match(seq_len(Ntip(tree)),tree$edge[,2])
tree$edge.length[match_tip_to_edge] = tree$edge.length[match_tip_to_edge] + extra_height
g3 = as.igraph(tree)
plot(tree)
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_positions <- data.frame(label = tree$tip.label,
                            x = coords$xx[1:Ntip(tree)],
                            y = coords$yy[1:Ntip(tree)])
node_positions <- data.frame(node = (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode),
                             x = coords$xx[(Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)],
                             y = coords$yy[(Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)])
names(node_positions)[1] = "label"
node_positions[,1] = tree$node.label
positions = rbind( tip_positions, node_positions)
coords3 = layout_(g3, with_sugiyama())
coords3$layout = as.matrix(positions[ match(names(V(g3)),positions$label),c("x","y")])

# Interactions
g4 <- graph_from_literal(
  "A" -+ "B",
  "A" -+ "C",
  "B" -+ "D",
  "C" -+ "D"
)
coords4 = layout_(g4, with_sugiyama())

png( file.path(run_dir,"graph_examples.png"), width = 4, height = 8, res = 200, units="in" )
  par( mfrow = c(4,2), mar=c(0,4,4,0), mgp = c(1.75, 0.25, 0), tck=-0.02, oma = c(1,0,0,0) )
  for( g_i in 1:4 ){
    graph = list( g1, g2, g3, g4 )[[g_i]]
    coords = list( coords1, coords2, coords3, coords4 )[[g_i]]

    plot( graph, layout = coords, vertex.color = "white", vertex.frame.color = NA,
              vertex.label = names(V(graph)), vertex.label.cex = 1.5, edge.label.cex = 1.5, xlim=c(-1,1.2) )
    mtext( side = 2, line = 1,
           text = c("Time lag\n", "Two dimensional\nspatial diffusion", "Phylogeny\n", "Interactions\namong variables")[g_i] )
    if(g_i==1) mtext( side = 3, text = "Graph", line = 1.5 )
    plot_graph(graph)
    if(g_i==1) mtext( side = 3, text = "Path matrix", line = 1.5 )
  }
dev.off()

