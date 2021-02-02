ValidatCluster <- function(x, y, plot=FALSE) {
  # Find the accuracy of our clustering of the projection by use of user input.
  #
  # Args:
  #   x: adjacency matrix of the projected bipartite graph.
  #   y: an array of user specified features for nodes of the considered projection. (For exmple, the bipartite graph is between cell lines and drugs,
  #      and the cell lines are used as projection. Then the tissues of the cell lines can be given to the function by the user to look at the clusters
  #      of the cell lines and compare them with the calculated clustering)
  #   Plot: a binary value indicating whether the result of two clusters should be plotted side by side or not. Default is False.
  #
  # Returns:
  #   the measured differences between two clustering clustering.
  
  return()
}