my_frame

p = ggplot(data = my_better_frame, mapping = aes(cromo_i, cromo_j, fill = avg))
p = p + geom_raster() + scale_fill_gradient(low = "grey", high = "green")
p

my_avg_matrix = reshape2::acast(my_better_frame, cromo_i~cromo_j, value.var = "avg")
my_dist_avg   = dist(my_avg_matrix)
my_clust_avg  = hclust(my_dist_avg)
?hclust
cutree(my_clust_avg, k = 3)
my_dist_avg

library(igraph)

g = graph_from_adjacency_matrix(adjmatrix = as.matrix(Matrix::forceSymmetric(as.matrix(my_dist_avg), 
                                                                   uplo = "U"), 
                                  ), 
                                weighted = TRUE, 
                                mode = "undirected"
                                )

g
plot(g)
E(g)$weight<-E(g)$weight/max(E(g)$weight)
plot(g)

E(g)$weight
write.graph(g, "gtest.gml", "gml")
write.table(x = get.data.frame(g, "edges"), file = "gtest.txt", row.names = FALSE, quote = FALSE)

cluster_infomap(g, e.weights = 1/E(g)$weight)
cluster_fast_greedy(g)
cluster_leading_eigen(g, weights = 1/E(g)$weight)
cluster_leading_eigen(g, weights = E(g)$weight)
