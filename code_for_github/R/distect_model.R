# DiSTect autologistic model. Y = 1 for SMS and 0 for CTR.
# A positive coefficient is higher expression in SMS.
# The last two columns of x are coordinates and are used only for the neighbor test.
# The radius is 1 coordinate unit. On full-resolution pixels, neighboring 16 um bins
# are about 56 px apart, so that term stays off. Array row and column turn some of it on.

stan_program <- "
data {
  int<lower=0> N; int<lower=0> P; int<lower=0> S;
  matrix[N,P+2] x; int y[N]; vector[S] u_mu; int label[N];
}
parameters {
  vector[P] beta; real<lower=0,upper=8> eta; vector<lower=0>[P] beta_gamma;
  real<lower=0,upper=1> w; vector<lower=-1,upper=1>[S] U;
  real<lower=0,upper=1> sigma_square; real<lower=0,upper=1> rho;
}
model{
  vector[N] mu; matrix[N,N] prob_neigh; matrix[S,S] covar_matrix;
  rho ~ uniform(0,1); sigma_square ~ inv_gamma(5,50);
  for (i in 1:S){ for (j in 1:S){
    if (i == j){ covar_matrix[i,j]=sigma_square; } else { covar_matrix[i,j]=sigma_square*rho; }
  }}
  U ~ multi_normal(u_mu, covar_matrix);
  for(i in 1:P){
    w ~ uniform(0,1); beta_gamma[i] ~ inv_gamma(5,50);
    target += log_sum_exp(log(1-w)+normal_lpdf(beta[i]|0,0.000001*beta_gamma[i]), log(w)+normal_lpdf(beta[i]|0,beta_gamma[i]));
  }
  for(i in 1:N) {
    for(j in 1:N){
      eta ~ uniform(0,8);
      if (label[i]==label[j] && j != i && sqrt(square(x[i,P+1]-x[j,P+1])+square(x[i,P+2]-x[j,P+2]))<=1){
        prob_neigh[i,j]=eta*y[j];
      }else{ prob_neigh[i,j]=0; }
    }
    mu[i] = dot_product(x[i,1:P],beta)+sum(prob_neigh[i,])+U[label[i]];
  }
  y ~ bernoulli_logit(mu);
}
"

plot_sms_ctr <- function(dat, title, fname) {
  dat$gene <- factor(dat$gene, levels = dat$gene)
  show_n <- nrow(dat)
  p <- ggplot2::ggplot(dat, ggplot2::aes(gene, posterior_mean, color = direction)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = cri_2.5, ymax = cri_97.5), linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c(higher_in_SMS = "#D7191C", higher_in_CTR = "#2C7BB6")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      legend.position = "top"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Red / positive = higher in SMS. Blue / negative = higher in CTR. Bars = 95% CrI.",
      x = NULL, y = "Posterior mean", color = NULL
    )
  ggplot2::ggsave(fname, p, width = max(10, 0.28 * show_n + 4), height = 5.8, dpi = 140)
  invisible(p)
}

run_sms_ctr <- function(obj, cells, tag, max_per_rep, out, stan_mod) {
  set.seed(42)
  md <- obj[[]][cells, , drop = FALSE]
  sampled <- unlist(lapply(unique(as.character(md$spatial_replicate)), function(rep) {
    ids <- cells[as.character(md$spatial_replicate) == rep]
    sample(ids, min(length(ids), max_per_rep))
  }), use.names = FALSE)

  sub <- subset(obj, cells = sampled)
  sub <- Seurat::JoinLayers(sub)
  sub <- Seurat::FindVariableFeatures(sub, selection.method = "vst", nfeatures = 200, verbose = FALSE)
  hvg <- Seurat::VariableFeatures(sub)
  expr <- t(as.matrix(Seurat::GetAssayData(sub, layer = "data")[hvg, sampled, drop = FALSE]))
  matrix_x <- as.matrix(cbind(
    expr,
    x = obj$coord_x[sampled],
    y = obj$coord_y[sampled]
  ))
  y_vector <- as.integer(obj$condition[sampled] == "SMS")
  label_list <- as.integer(factor(obj$spatial_replicate[sampled]))
  ok <- complete.cases(matrix_x)
  matrix_x <- matrix_x[ok, , drop = FALSE]
  y_vector <- y_vector[ok]
  label_list <- label_list[ok]

  fit <- rstan::vb(
    stan_mod,
    data = list(
      N = nrow(matrix_x), P = ncol(matrix_x) - 2, x = matrix_x, y = y_vector,
      label = label_list, S = length(unique(label_list)),
      u_mu = rep(0, length(unique(label_list)))
    ),
    algorithm = "meanfield", seed = 128, tol_rel_obj = 1e-04, iter = 10000
  )
  beta <- rstan::summary(fit)$summary
  beta <- beta[grep("^beta\\[", rownames(beta)), , drop = FALSE]
  n_genes <- ncol(matrix_x) - 2
  coef_df <- data.frame(
    gene = colnames(matrix_x)[seq_len(n_genes)],
    posterior_mean = beta[seq_len(n_genes), "mean"],
    cri_2.5 = beta[seq_len(n_genes), "2.5%"],
    cri_97.5 = beta[seq_len(n_genes), "97.5%"],
    stringsAsFactors = FALSE
  )
  coef_df$excludes_zero <- (coef_df$cri_97.5 < 0) | (coef_df$cri_2.5 > 0)
  coef_df$direction <- ifelse(coef_df$posterior_mean > 0, "higher_in_SMS", "higher_in_CTR")
  coef_df <- coef_df[order(-abs(coef_df$posterior_mean)), ]
  utils::write.csv(coef_df, file.path(out, paste0(tag, "_coefficients.csv")), row.names = FALSE)
  selected <- coef_df[coef_df$excludes_zero, ]
  utils::write.csv(selected, file.path(out, paste0(tag, "_selected.csv")), row.names = FALSE)
  message(tag, ": ", nrow(selected), " / ", nrow(coef_df), " genes called")
  if (nrow(selected)) {
    plot_sms_ctr(
      utils::head(selected, 40),
      paste0("SMS vs CTR, ", tag),
      file.path(out, paste0(tag, "_cri.png"))
    )
  }
  invisible(coef_df)
}

fit_three <- function(obj, out) {
  if (!exists("distect_stan_mod", inherits = TRUE)) {
    distect_stan_mod <<- rstan::stan_model(model_code = stan_program)
  }
  all_cells <- colnames(obj)
  upper <- colnames(obj)[!is.na(obj$layer_class) & obj$layer_class == "upper"]
  deep <- colnames(obj)[!is.na(obj$layer_class) & obj$layer_class == "deep"]
  run_sms_ctr(obj, all_cells, "all_cortex", 500, out, distect_stan_mod)
  run_sms_ctr(obj, upper, "upper_L1_L23", 400, out, distect_stan_mod)
  run_sms_ctr(obj, deep, "deep_L5_L6", 400, out, distect_stan_mod)
}
