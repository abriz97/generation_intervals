suppressPackageStartupMessages({
    library(Rcpp)
    library(data.table)
    library(ggplot2)
    library(igraph)
})


relative_path <- "src/implement_Rcpp.R" 
here::i_am(relative_path)
gitdir <- here::here() 
source(file.path(gitdir, 'R/utils.R'))
source(file.path(gitdir, 'R/gibbs_gi_distribution.R'))

overleaf_figs <- "/extraspace/latex/overleaf/generation-intervals/figures"

# DEFINE TRANSMISSION NETWORK

edges_dt <- data.table(
    SOURCE = c("A", "B", "B"),
    RECIPIENT = c("B", "C", "D")
)

nodes_dt <- data.table(
    ID = c("A", "B", "C", "D"),
    m = c(1, 2, .5, 1),
    M = c(4, 3, 4, 5)
)

network <- list(edges=edges_dt, nodes=nodes_dt)
check_network(list(edges=edges_dt, nodes=nodes_dt), check_connected = TRUE)
network$nodes <- propagate_minimum_infection_date(network=network, verbose=TRUE)
network$nodes <- propagate_maximum_infection_date(network=network, verbose=TRUE)


if(0){
    # specify priors for infection dates 
    priors_u <- make_uniform_priors(network)
    priors_g <- make_normal_priors(network, t_s = nodes_dt$m, t_e = nodes_dt$M)

    # Inference step
    posterior <- estimate_dates_network_gibbs(network, iters=10000, priors=priors_u, verbose=TRUE)
    plot_posterior_doi(DT=posterior)
}

# load Rcpp functions in file.path(gitdir, "C/gibbs.cpp")


sourceCpp(file.path(gitdir, "C/gibbs.cpp"))

cppdata <- make_cpp_data_gibbs_network(network=network)
I <- c(2, 2.5, 3, 4)

posterior_cpp <- with(cppdata,
    GibbsSamplerUniform(initial_values = I,
        n_iter = 100000,
        source = source,
        recipients = recipients,
        min_infection_date = min_infection_date,
        max_infection_date = max_infection_date
    )
)
colnames(posterior_cpp) <- LETTERS[seq_len(ncol(posterior_cpp))]
colnames(posterior_cpp)
p1 <- plot_posterior_doi(posterior_cpp)
p1
# p2 <- plot_posterior_doi(posterior)
# ggpubr::ggarrange(p1, p2)


# now consider a case in which the generation interval is gamma distributed

SHAPE = 1.5; RATE = 2

# now, what do we have to do at every Gibbs Sampling iteration...
# 1. find [a, b] as before. Cool. this is easy.
# 2. find the maximizers of each function in the product... 
y_s <- 2.1; y_r <- c(3, 3.5)
a <- 2.5; b <- 3

posterior <- estimate_dates_network_gibbs_generationinterval(
    network,
    iters=100, verbose=TRUE) 
plot_posterior_doi(DT=posterior)

priors_u <- make_uniform_priors(network)


if(0){
    # visuallise
    delta <- .0001
    x_range <- seq(from=a, to=b, by=delta)
    y <- sapply(x_range, function(x){
        x <- densities$multiplier * x + densities$offset
        c(2, dgamma(x[densities$family == "gamma"], shape=SHAPE, rate=RATE))
    })
    pdfs <- data.table(t(as.data.table(y)))
    pdfs[, all := V1 * V2 * V3 * V4, ]
    pdfs <- pdfs[, lapply(.SD, proportions)]
    pdfs <- pdfs[, lapply(.SD, function(x) {x / delta}) ]
    pdfs[, x:=x_range]

    dim(y)
    integral_y <- sum(y) * delta
    sprintf("integral of y is %f", integral_y)

    g <- ggplot(data.frame(y=samples), aes(y)) +
    geom_histogram(aes(y= after_stat(density) ), fill="grey80", color="black")  +
    geom_line( data= pdfs, aes(x=x, y=all, color="all")) +
    geom_line( data= pdfs, aes(x=x, y=V1, color="prior"), linetype="dashed") +
    geom_line( data= pdfs, aes(x=x, y=V2, color="source"), linetype="dashed") +
    geom_line( data= pdfs, aes(x=x, y=V3, color="recipient"), linetype="dashed") +
    geom_line( data= pdfs, aes(x=x, y=V4, color="recipient"), linetype="dashed") +
    theme_minimal(base_size=12) + 
    labs( x="samples", y="density", title="ARSP Example: pdf vs samples")

    filename <- file.path(overleaf_figs, "arsp_example_gamma.pdf")
    ggsave(plot=g, filename=filename, width=6, height=4.5)
    system(sprintf("zathura %s &", filename))

}
