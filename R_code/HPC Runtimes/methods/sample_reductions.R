fulltie_samplereduction <- function(x) {
  return(unique(x))
}

greedy_samplereduction_BIV <- function(x, m) {
  # Initialization
  n <- nrow(x);
  b <- g <- c();
  
  # Calculation of initial weights
  for (i in 1:n) {
    b <- c(b, (sum(x[, 1] == x[i, 1]) > 1) + (sum(x[, 2] == x[i, 2]) > 1));
    g <- c(g, max(sum(x[, 1] == x[i, 1]), sum(x[, 2] == x[i, 2])));
  }
  w <- b / g;
  
  # Order of weights that are greater than zero
  index <- order(w, decreasing = TRUE);
  index <- index[w[index] > 0];
  
  # Select a representative per tie in the greedy selection by removing it from
  # the index set, decreasing in b. Check for full ties, then partial ties.
  I <- J <- c()
  for (i in index) { if (!i %in% union(I, J)) {
    
    B.0 <- B.1 <- B.2 <- c();
    i.0 <- i.1 <- i.2 <- NA;
    
    B.0 <- index[x[index, 1] == x[i, 1] & x[index, 2] == x[i, 2]];
    B.0 <- B.0[!B.0 %in% union(I, J)];
    
    B.1 <- index[x[index, 1] == x[i, 1]];
    B.1 <- B.1[!B.1 %in% union(union(I, J), B.0)];
    
    B.2 <- index[x[index, 2] == x[i, 2]];
    B.2 <- B.2[!B.2 %in% union(union(I, J), B.0)];
    
    if (!any(is.na(B.0)) && !is.na(B.0[1])) { i.0 <- sample(c(B.0, B.0), 1); }
    if (!any(is.na(B.1)) && !is.na(B.1[1])) { i.1 <- sample(c(B.1, B.1), 1); }
    if (!any(is.na(B.2)) && !is.na(B.2[1])) { i.2 <- sample(c(B.2, B.2), 1); }

    if (
      length(J) <= m & m <= length(J) + 
      length(B.0[B.0!=i.0]) + length(B.1[B.1!=i.1]) + length(B.2[B.2!=i.2])
    ) {
      if (!is.null(J) & length(J) > 0) {
        return(x[-J,])
        stop()
      } else {
        return(x)
        stop()
      }
    } else {
      I <- union(I, c(i.0, i.1, i.2)[!is.na(c(i.0, i.1, i.2))]);
      J <- union(J, union(B.0[B.0!=i.0], union(B.1[B.1!=i.1], B.2[B.2!=i.1])))
      J <- J[!is.na(J)]
    }
  }}
  if (!is.null(J)) {
    return(x[-J,])
  } else {
    return(x)
  }
}
