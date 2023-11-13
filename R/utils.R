verbose_print <- function(verbose=FALSE, ...){
    if( verbose ) print(...)
}

check_network <- function(network=NULL, edges, nodes, check_connected = TRUE ){

    if(! is.null(network)){
        edges <- network$edges
        nodes <- network$nodes
    }

    # check format
    stopifnot( is.data.table(edges) )
    stopifnot( names(edges) == c("SOURCE", "RECIPIENT") )
    stopifnot( is.data.table(nodes) )
    stopifnot( names(nodes) == c("ID", "m", "M") )

    # check times
    stopifnot( nodes[, all(m <= M)] )

    # check that individuals appear as recipients at most once
    if( edges[, uniqueN(RECIPIENT) != .N ]){
        stop("Individuals appear as recipients more than once")
    }

    # check that individuals are in nodes
    edges_ids <- edges[, unique(c(SOURCE, RECIPIENT))]
    if( ! all( edges_ids %in% nodes_dt$ID) & all( nodes_dt$ID %in% edges_ids) ){
        stop("IDs in edges and nodes do not match")
    }

    # check network is connected:
    if( check_connected ){
        stopifnot( "the transmission network is not connected" =  nrow(edges) == nrow(nodes) - 1)
    }

    return(TRUE)
}


propagate_minimum_infection_date <- function(network=NULL, edges, nodes, i0=NA_character_, verbose=FALSE){
    
    if( ! is.null(network) ){
        edges <- network$edges
        nodes <- network$nodes
    }

    # in the case of a connected network, the index case is the only non-recipient
    if( is.na(i0) ) i0 <- setdiff( nodes$ID, edges$RECIPIENT )
    stopifnot( "i0 is not in nodes"=i0 %in% nodes$ID )


    queue <- c(i0)
    output <- copy(network$nodes)

    while( length(queue) ){

        i <- queue[1]
        m_i <- output[ ID == i, m ]
        queue <- queue[-1]
        verbose_print(paste("Individual ", i), verbose=verbose)

        # get the recipients of i
        recipients <- edges[ SOURCE == i, RECIPIENT ]

        for ( rec in recipients ){
            output[ ID == rec, m :=  max(m_i, m)]
        }

        # add the recipients to the queue
        queue <- c(recipients, queue)

    }
    return(output)
}

propagate_maximum_infection_date <- function(network=NULL, edges, nodes, verbose=FALSE){
    
    if( ! is.null(network) ){
        edges <- network$edges
        nodes <- network$nodes
    }

    # find individuals who are not sources 
    queue <- setdiff( nodes$ID, edges$SOURCE )
    output <- copy(network$nodes)

    while( length(queue) ){

        i <- queue[1] 
        M_i <- output[ ID == i, M ]
        queue <- queue[-1]
        verbose_print(paste("Individual ", i), verbose=verbose)

        # get the source of i
        source <- edges[ RECIPIENT == i, SOURCE ]
        if( ! length(source) ){ return(output) }
        output[ ID == source, M := min(M_i, M)]

        # add the source to the queue
        queue <- c(source, queue)
    }
    return(output)
}

make_source_vector <- function(network){
    merge( network$nodes, network$edges, by.x="ID", by.y="RECIPIENT", all.x=TRUE)$SOURCE
}

make_recipients_list <- function(network){
    # make list of recipients for each individual in the network 
    nodes <- network$nodes
    edges <- network$edges
    sources <- merge( nodes, edges, by.x="ID", by.y="RECIPIENT", all.x=TRUE)$SOURCE
    tmp <- merge( nodes, edges, by.x="ID", by.y="SOURCE", all.x=TRUE)
    tmp <- tmp[, list(REC = list(REC = RECIPIENT)), by=ID]
    recipients <- setNames(tmp$REC, tmp$ID)
    return(recipients)
}

sample_truncated_distribution <- function(fname="unif", m, M, fargs){

    # fname <- prefix to distributoin from which we want to sample

    F <- function(q) {
      do.call(paste0("p", fname), c(list(q = q), fargs))
    }
    F_inv <- function(p) {
      do.call(paste0("q", fname), c(list(p = p), fargs))
    }

    # F <- function(q) punif(q=q, min=-2, max=5)
    # F_inv <- function(p) qunif(p=p, min=-2, max=5)
    # m <- 0; M <- 1

    u_sample <-  runif(1, min=F(m), max=F(M)) 
    return(F_inv(u_sample))
    
}

estimate_dates_network_gibbs <- function(network, iters, priors, y0=NA_real_, verbose = FALSE ){

    # network: list of data.tables with edges and nodes
    # iters: number of iterations
    # priors: list specifying priors for each individual. 
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
    if( is.na(y0)) y0 <- initial_ms

    # get vector of sources and list of recipients
    recipients <- make_recipients_list(network)
    sources <- make_source_vector(network)

    # initialize output
    output <- matrix(NA_real_, n_nodes * (iters + 1), nrow=iters + 1, ncol=n_nodes)
    colnames(output) <- nodes$ID
    output[1,  ] <- y0

    
    for( iter in (1:iters) + 1 ){

        if( verbose ){
            sprintf("iteration %d \n", iter)  |> cat()
        }
        

        for (v in 1:n_nodes){

            src     <- sources[v]
            rec     <- recipients[[v]]
            prior   <- priors[[v]]
            m <- max(c(initial_ms[v], output[ iter - 1, nodes$ID %in% src ]))
            M <- min(c(initial_Ms[v], output[ iter - 1, nodes$ID %in% rec ]))

            draw_v <- sample_truncated_distribution(
                fname = prior$fname,
                fargs = prior$fargs,
                m = m,
                M = M
            ) 
            output[iter, v] <- draw_v

            if( verbose ){
                sprintf(
                    "individual %d: %f (m = %f, M = %f)\n",
                    v, draw_v , m, M) |> cat()
            }

        }

        if( verbose ) cat('\n')
    }
    # remove first line, initialization
    return(output[-1, ])
}


make_normal_priors <- function(network, t_s, t_e){
    # make normal priors for each individual in the network
    #   centered at t_s + t_emd/2 with variance = length of interval / 4 .
    # network: list of data.tables with edges and nodes
    # returns: list of data.tables with edges and nodes

    nodes <- network$nodes
    n_nodes <- nrow(nodes)
    .make_fargs <- function(idx){
        list(
            fname = "norm",
            fargs = list(
                mean = (t_s[idx] + t_e[idx])/2,
                sd = (t_e[idx] - t_s[idx])/4
            )
        )
    }
    priors_specification <- lapply( 1:n_nodes, .make_fargs )
    return(priors_specification)

}

make_uniform_priors <- function(network){
    # make uniform priors for each individual in the network
    # network: list of data.tables with edges and nodes
    # returns: list of data.tables with edges and nodes

    nodes <- network$nodes
    n_nodes <- nrow(nodes)
    .make_fargs <- function(idx){
        list(
            fname = "unif",
            fargs = list(
                min = nodes$m[idx],
                max = nodes$M[idx]
            )
        )
    }
    priors_specification <- lapply( 1:n_nodes, .make_fargs )
    return(priors_specification)
}

plot_posterior_doi <- function(DT){

    posterior_dt <- as.data.table(DT)

    list_gghists <- lapply( 
        colnames(posterior_dt),
        function(c){
            force(c)
            geom_density( aes(x=!!sym(c), color=c, fill=c), alpha=.2) 
        })

    force(list_gghists)

    ggplot(posterior_dt, aes()) +
        list_gghists + 
        theme_minimal() +
        labs( x = "date of infection", y = "posterior", color="individual", fill="individual")

}

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
