get_onedimensional_ranges <- function(y_src, y_rec, A = a, B = b) {
    rbind(
        data.frame(type = "prior", min = A, max = B, offset = 0, multiplier = 1),
        if(length(y_src)){
            data.frame(type = "source", min = A - y_src, max = B - y_src, offset = y_src, multiplier = 1)
        }else{data.frame()},
        if(length(y_rec)){
            data.frame(type = "recipient", min = y_rec - B, max = y_rec - A, offset = y_rec, multiplier = -1)
        }else{data.frame()}
    )
}

get_maximisers_uniform <- function(min, max) {
    data.frame(family = "unif", mode = (min + max) / 2, I = max - min, sup = (max - min)^(-1))
}

get_maximisers_truncated_gamma <- function(shape = SHAPE, rate = RATE, min, max) {
    # should I rewrite that in terms of fargs instead of shape and rate?

    # check inputs
    stopifnot(
        "min and max need to have same length" = length(min) == length(max),
        "shape and rate need to have same length" = length(shape) == length(rate)
    )
    if (length(shape) > 1 && length(shape) != length(min)) {
        stop("Dimensions for gamma parameters and min/max do not match")
    }

    # Want to find mode, integral and sup f
    mode <- pmax(0, (shape - 1) / rate)
    if (length(shape) == 1 & length(min) > 1) {
        mode <- rep(mode, length(min))
    }

    # check that mode is in [min, max]
    mode[mode < min] <- min[mode < min]
    mode[mode > max] <- max[mode > max]
    # compute intergal over domain, and supremum
    integral <- pgamma(q = max, shape = shape, rate = rate)
    integral <- integral - pgamma(q = min, shape = shape, rate = rate)
    sup_f <- dgamma(x = mode, shape = shape, rate = rate) / integral

    return(data.frame(family = "gamma", mode = mode, I = integral, sup = sup_f))
}

evaluate_family_pdf_at_x <- function(x, family, fargs) {
    do.call(paste0("d", family), c(list(x = x), fargs))
}

one_step_arsp_gammagi <- function(shape, rate, y_s, y_r, a, b) {

    if( a == b ) return(a)
    stopifnot( a < b)

    densities <- get_onedimensional_ranges(y_src = y_s, y_rec = y_r, A = a, B = b)

    stopifnot("First data.frame entry must concern prior" = densities$type[1] == "prior")
    info <- with(densities, {
        idx_prior <- type == "prior"
        rbind(
            get_maximisers_uniform(min = min[idx_prior], max = max[idx_prior]),
            get_maximisers_truncated_gamma(min = min[!idx_prior], max = max[!idx_prior])
        )
    })
    densities <- cbind(densities, info)

    # compute the normalising constant of the product of pdfs
    kappa <- prod(densities$I)

    # identify idx with maximum sup
    idx_n0 <- which.max(densities$sup)
    idx_prior <- which(rownames(densities) == "prior")

    ## Rejection sampling step
    # accepted = FALSE

    while (TRUE) {
        proposed_sample <- with(densities[idx_n0, ], {
            out <- sample_truncated_distribution(
                m = min,
                M = max,
                fname = "gamma",
                fargs = list(shape = SHAPE, rate = RATE)
            )
            out * multiplier + offset
        })

        stopifnot(proposed_sample %between% c(a,b))

        # Compute acceptance probability

        d_accept <- densities[-c(idx_n0), ]

        log_pdf_ratio_at_sample <- with(d_accept, {
            x <- offset + multiplier * proposed_sample

            # initialise with zeros (~ uniform priors)
            pdfs <- rep(0, length(x))
            idx_gamma <- which(family == "gamma")
            pdfs[idx_gamma] <- dgamma(x = x[idx_gamma], shape = SHAPE, rate = RATE, log = TRUE) - log(I[idx_gamma])
            sum(pdfs - log(sup))
        })


        if (exp(log_pdf_ratio_at_sample) > runif(1)) {
            # accepted <- TRUE
            return(proposed_sample)
        }
    }
}


estimate_dates_network_gibbs_generationinterval <- function(network, iters, y0 = NA_real_, verbose = FALSE, shape, rate) {
    # network: list of data.tables with edges and nodes
    # iters: number of iterations
    ### REMOVED PRIORS
    ## priors: list specifying priors for each individual.
    #       Each element of the list is a list with the following elements:
    #       - fname: suffix of the distribution to sample from (eg. "norm", "unif")
    #       - fargs: list of arguments to pass to the distribution
    # y0: initial values for each parameter
    # verbose: print progress
    # returns: list of data.tables with edges and nodes

    nodes <- network$nodes
    edges <- network$edges
    n_nodes <- nrow(nodes)
    n_edges <- nrow(nodes)
    initial_ms <- nodes$m
    initial_Ms <- nodes$M
    if (is.na(y0)) y0 <- initial_ms

    # get vector of sources and list of recipients
    recipients <- make_recipients_list(network)
    sources <- make_source_vector(network)

    # initialize output
    output <- matrix(NA_real_, n_nodes * iters, nrow = iters, ncol = n_nodes)
    colnames(output) <- nodes$ID

    current <- y0
    # output[1, ] <- y0


    for (iter in (1:iters) ) {
        if (verbose) {
            sprintf("iteration %d \n", iter) |> cat()
        }

        for (v in 1:n_nodes) {
            src <- sources[v]
            rec <- recipients[[v]]
            # prior   <- priors[[v]]
            y_s <- current[nodes$ID %in% src]
            y_r <- current[nodes$ID %in% rec]

            m <- max(c(initial_ms[v], y_s))
            M <- min(c(initial_Ms[v], y_r))

            draw_v <- one_step_arsp_gammagi(
                shape = shape, rate = rate,
                y_s = y_s, y_r = y_r, a = m, b = M
            )

            current[v] <- draw_v

            if (verbose) {
                sprintf(
                    "individual %d: %f (m = %f, M = %f)\n",
                    v, draw_v, m, M
                ) |> cat()
            }
        }
        output[iter, ] <- current

        if (verbose) cat("\n")
    }
    # remove first line, initialization
    return(output[-1, ])
}

plot_network_and_posterior <- function(net, post=NULL){

    tmp1 <- network::network(x=network$edges,directed=TRUE)

    p0 <- ggnet::ggnet2(tmp1, size = 8, label = TRUE, label.size = 5,
        arrow.size=5, arrow.gap=.05)
         
    p1 <- ggplot(network$nodes, aes(y=ID, color=ID)) +
        geom_point(aes(x=m)) +
        geom_point(aes(x=M)) +
        geom_linerange(aes(xmin=m, xmax=M)) +
        labs(color="individual", x="date of infection", y=NULL) +
        theme_minimal()

    p2 <- plot_posterior_doi(DT=posterior)

    if( is.null(post) ){
        ( p0 / p1 ) 
    }else{
        require(patchwork)
        (( p0 / p1 ) | p2 ) + plot_layout(guides="collect")
    }
}
