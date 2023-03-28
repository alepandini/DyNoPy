######################################################################
# functions
coev_heatmap <- function(coev_mat, main_title, rescale = T, col_scheme = "Viridis"){
  if (rescale){
    s_coev_mat <- coev_mat / mean(coev_mat)
    s_coev_mat[s_coev_mat < 1] = NA
  }else{
    s_coev_mat <- coev_mat
  }
  nres <- nrow(coev_mat)
  xseq <- seq(5, nres, 5)
  tseq <- xseq / nres
  opar <- par(no.readonly = TRUE)
  par(pty = 's')
  par(las = 2)
  image(
    s_coev_mat,
    xaxt = 'n',
    yaxt = 'n',
    col = hcl.colors(12, col_scheme, rev = T),
    xlab = 'residue position',
    ylab = 'residue position',
    main = main_title
  )
  axis(1, xseq, at = tseq)
  axis(2, xseq, at = tseq)
  par(opar)
}

coev_barplot <- function(coev_mat, main_title, rescale = T){
  if (rescale){
    s_coev_mat = coev_mat / mean(coev_mat)
    s_coev_mat[s_coev_mat < 1] = 0
  }else{
    s_coev_mat = coev_mat
  }
  nres <- nrow(coev_mat)
  xseq <- seq(5, nres, 5)
  cum_coev_vec <- apply(s_coev_mat, 1, median)
  opar <- par(no.readonly = TRUE)
  par(pty = 'm')
  par(las = 2)
  xpos <- barplot(
    cum_coev_vec,
    border = 'grey',
    col = 'pink',
    xlab = 'residue position',
    #ylim = c(0, ceiling(max(cum_coev_vec)/10)*10),
    main = main_title
  )
  axis(1, xseq, at = xpos[xseq])
  par(opar)
}

coev_boxplot <- function(coev_mat, main_title, rescale = T){
  if (rescale){
    s_coev_mat = coev_mat / mean(coev_mat)
    s_coev_mat[s_coev_mat < 1] = NA
  }else{
    s_coev_mat = coev_mat
  }
  nres <- nrow(coev_mat)
  xseq <- seq(5, nres, 5)
  opar <- par(no.readonly = TRUE)
  par(pty = 'm')
  par(las = 2)
  bxp = boxplot(
    s_coev_mat,
    xlab = 'residue position',
    xaxt = 'n',
    outcex=0.5,
    outpch=16,
    main = main_title
  )
  axis(1, xseq)
  par(opar)
}
