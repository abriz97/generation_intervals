#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
#include <iostream>

// write a test function
// [[Rcpp::export]]
NumericVector test(NumericVector x) {
    return x;
}

// Function to perform Gibbs sampling with uniform priors
// [[Rcpp::export]]
NumericMatrix GibbsSamplerUniform(NumericVector initial_values,
        int n_iter,
        IntegerVector source,
        List recipients,
        NumericVector min_infection_date,
        NumericVector max_infection_date) {

    int n_individuals = initial_values.size();

    // Initialise chains
    NumericMatrix samples(n_iter, n_individuals);
    // we are not doing the proper conditioning
    NumericVector y(n_individuals);
    y = initial_values;

    // samples(0, _) = initial_values;

    
    //
    // Perform Gibbs sampling
    for (int iter = 0; iter < n_iter; iter++) {
        for (int i = 0; i < n_individuals; i++) {

            double a = min_infection_date[i];
            double b = max_infection_date[i]; 

            // Get the source of individual i
            int source_i = source[i];
            if ( source_i != -1 ){
                a = std::max(a, y(source_i));
            }

            // Get the recipients of individual i
            IntegerVector recipients_i = recipients[i];

            // check recipient size is larger than 0 before looping
            for (int j = 0; j < recipients_i.size(); j++) {
                if( recipients_i[j] != -1){
                    double recipient_value = y(recipients_i[j]);
                    b = std::min(b, recipient_value);
                }
            }

            // Sample from uniform distribution in [a, b]
            double sample = R::runif(a, b);
            y(i) = sample;
        }

        // append the new sample to the chain
        samples(iter, _) = y;
    }

    return samples;
}
