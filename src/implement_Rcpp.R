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
make_cpp_data_gibbs_network <- function(network){

    # need to first set all character ids to integers in a meaningful manner.
    src <- make_source_vector(network)
    rec <- make_recipients_list(network)
    char2int_dict <- setNames(seq_along(rec) - 1, names(rec))

    # names(rec) <- seq_along(rec)
    na2minus1 <- function(x) {x[is.na(x)] <- -1; x}
    rec <- lapply(rec, function(x) unname(char2int_dict[x]) )
    rec <- lapply(rec, na2minus1)
    src <- unname(char2int_dict[src])
    src <- na2minus1(src)

    list(
        initial_values = network$nodes$m,
        source = src,
        recipients = rec, 
        min_infection_date = network$nodes$m,
        max_infection_date = network$nodes$M
    )
}


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


# now consider a case in which the generation interval is exponentialy
SHAPE = 1.5; RATE = 2
MODE = max(0, (SHAPE - 1) / RATE)
{
    # visualise generation interval distribution
    .seq <- seq(0, 10, by=.1)
    plot(.seq, dgamma(x=.seq, shape=SHAPE, rate=RATE), type="l")
    abline(v=MODE, col="red")
}

# now, what do we have to do at every Gibbs Sampling iteration...
# 1. find [a, b] as before. Cool. this is easy.
# 2. find the maximizers of each function in the product... 
y_s <- 2.1; y_r <- 3
a <- 2.5; b <- 3

# compute information for factors
# f_sup_prior <- 1/(3 - 2.5); int_prior <- 1/f_sup_prior
densities <- list(
    prior = c(mode = (a+b)/2, I = b - a, sup = (b-a)^(-1)),
    source = get_maximisers_truncated_gamma(min= a - y_s , max= b - y_s ) ,
    recipients = get_maximisers_truncated_gamma(min= y_r - b, max= y_r - a)
)

# compute the normalising constant of the product of pdfs
kappa <- prod( sapply(densities, `[`, "I") )
# now compute the supemums of each factor, and divide by the maximum
sups <- sapply(densities, `[`, "sup") 
# identify idx with maximum sup:
n0_idx <- which.max(sups)
sups_const <- prod(sups) / max(sups)

if( which.max(sups) == 2){
    source <- TRUE
    recipient <- FALSE
}


# initalise the acceptance probability
log_c_lambda <- log(sups_const)

# This implies we should use the source distribution for our rejection sampling...
yy <- c()

# can we compute the acceptance probability? It is the inverse of the maximum sup_f
accepted = FALSE
while ( TRUE ){
    
    if ( source ){
        sample <- sample_truncated_distribution(fname="gamma", m=a - y_s, M = b - y_s, fargs=list(shape=SHAPE, rate=RATE))
        sample <- y_s + sample

    }else if ( recipient ){
        sample <- sample_truncated_distribution(fname="gamma", m=y_r - b, M = y_r - a, fargs=list(shape=SHAPE, rate=RATE))
        sample <- y_r - sample # TODO to check
    }

    # Once we have a sample, we can finally compute the acceptance probability:
    # we need to evaluate the pdf of each term in the product, and divide by the maximum 
    log_pdf_ratio_at_sample <- {
        logs <- c(
             - log( densities$prior["I"] ),
            dgamma(x = sample - y_s, shape=SHAPE, rate=RATE, log=TRUE) - log(densities$source["I"] ),
            dgamma(x = y_r - sample, shape=SHAPE, rate=RATE, log=TRUE) - log(densities$recipients["I"])
        )
        sum(logs) - logs[[n0_idx]] - log_c_lambda
    }
    yy <- c(exp(log_pdf_ratio_at_sample), yy)
    print(yy[1])

    if ( exp(log_pdf_ratio_at_sample) > runif(1)){
        accepted <- TRUE
    }
}

# But this could still be right! The acceptance probabilities do NOT need to look
# uniformly distributed!!!!!!
hist(yy)

