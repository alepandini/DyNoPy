######################################################################
# functions
dynor_mat_heatmap <- function(input_mat, main_title, rescale = T, col_scheme = "Viridis"){
  if (rescale){
    s_input_mat <- input_mat / mean(input_mat)
    s_input_mat[s_input_mat < 1] = NA
  }else{
    s_input_mat <- input_mat
  }
  nres <- nrow(input_mat)
  xseq <- seq(5, nres, 5)
  tseq <- xseq / nres
  opar <- par(no.readonly = TRUE)
  par(pty = 's')
  par(las = 2)
  image(
    s_input_mat,
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

dynor_mat_impplot <- function(input_mat, main_title, rescale = T){
  if (rescale){
    s_input_mat = input_mat / mean(input_mat)
    s_input_mat[s_input_mat < 1] = NA
  }else{
    s_input_mat = input_mat
  }
  nres <- nrow(input_mat)
  xseq <- seq(5, nres, 5)
  mean_coev_vec <- apply(s_input_mat, 1, mean, na.rm = TRUE)
  mean_of_mean <- mean(mean_coev_vec)
  opar <- par(no.readonly = TRUE)
  par(pty = 'm')
  par(las = 2)
  plot(
    mean_coev_vec,
    col = c('grey', 'red')[(mean_coev_vec > mean_of_mean)+1],
    xlab = 'residue position',
    ylab = 'mean',
    ylim = c(min(mean_coev_vec), max(mean_coev_vec)),
    type = 'h',
    xaxt = 'n',
    main = main_title
  )
  abline(h = mean_of_mean, col = rgb(0,0,1,0.4))
  points(
    mean_coev_vec,
    col = c('grey', 'red')[(mean_coev_vec > mean_of_mean)+1],
    type = 'h'
  )
  axis(1, xseq)
  par(opar)
}

dynor_mat_boxplot <- function(input_mat, main_title, rescale = T, outflag = TRUE){
  if (rescale){
    s_input_mat = input_mat / mean(input_mat)
    s_input_mat[s_input_mat < 1] = NA
  }else{
    s_input_mat = input_mat
  }
  nres <- nrow(input_mat)
  xseq <- seq(5, nres, 5)
  par(pty = 'm')
  par(las = 2)
  bxp = boxplot(
    s_input_mat,
    xlab = 'residue position',
    xaxt = 'n',
    outcex=0.5,
    outpch=16,
    outline = outflag,
    main = main_title
  )
  axis(1, xseq)
}

dynor_mat_panel_boxplot <- function(input_mat, main_title, rescale = T){
    par(mfrow = c(2,1))
    dynor_mat_boxplot(input_mat, main_title, rescale, outflag = FALSE)
    dynor_mat_boxplot(input_mat, main_title, rescale, outflag = TRUE)
    par(mfrow = c(1,1))
}
