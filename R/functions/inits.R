# a function to generate initial values for the model

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = rep(
        1,
        constant_list$ndata
      ),
      x = matrix(
        1,
        ncol = constant_list$ncity,
        nrow = constant_list$nspecies
      ),
      a_community = rnorm(
        constant_list$nomega
      ),
      tau_community_omega = rgamma(
        constant_list$nomega,
        1,
        1
      ),
      a_species = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('nomega', 'nspecies')]
            )
          )
        ),
        nrow = constant_list$nomega,
        ncol = constant_list$nspecies
      ),
      b_community = rnorm(
        constant_list$npsi
      ),
      tau_community_psi = rgamma(
        constant_list$npsi,
        1,
        1
      ),
      tau_shape_psi = runif(
        constant_list$npsi,
        0.5,
        2
      ),
      tau_rate_psi = runif(
        constant_list$npsi,
        0.5,
        2
      ),
      b_species = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('npsi', 'npspecies')]
            )
          )
        ),
        nrow = constant_list$npsi,
        ncol = constant_list$nspecies
      ),
      tau_species_psi = matrix(
        rgamma(
          prod(
            unlist(
              constant_list[c('npsi', 'nspecies')]
            )
          ),
          1,
          1
        ),
        nrow = constant_list$npsi,
        ncol = constant_list$nspecies
      ),
      b_species_city = array(
        rnorm(
          prod(
            unlist(
              constant_list[c('npsi', 'nspecies', 'ncity')]
            )
          )
        ),
        dim = unlist(
          constant_list[c('npsi', 'nspecies', 'ncity')]
        )
      ),
      theta_community = rnorm(
        1
      ),
      tau_community_theta = rgamma(
        1,
        1,
        1
      ),
      tau_shape_theta = runif(
        1,
        0.5,
        2
      ),
      tau_rate_theta = runif(
        1,
        0.5,
        2
      ),
      theta_species = rnorm(
        constant_list$nspecies
      ),
      tau_species_theta = rgamma(
        constant_list$nspecies,
        1,
        1
      ),
      theta_psi = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('nspecies', 'ncity')]
            )
          )
        ),
        nrow = constant_list$nspecies,
        ncol = constant_list$ncity
      ),
      c_community = rnorm(
        constant_list$nrho
      ),
      tau_community_rho = rgamma(
        constant_list$nrho,
        1,
        1
      ),
      tau_shape_rho = runif(
        constant_list$nrho,
        0.5,
        2
      ),
      tau_rate_rho = runif(
        constant_list$nrho,
        0.5,
        2
      ),
      c_species = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('nrho', 'nspecies')]
            )
          )
        ),
        nrow = constant_list$nrho,
        ncol = constant_list$nspecies
      ),
      tau_species_rho = matrix(
        rgamma(
          prod(
            unlist(
              constant_list[c('nrho', 'nspecies')]
            )
          ),
          1,
          1
        ),
        nrow = constant_list$nrho,
        ncol = constant_list$nspecies
      ),
      c_species_city = array(
        rnorm(
          prod(
            unlist(
              constant_list[c('nrho', 'nspecies', 'ncity')]
            )
          )
        ),
        dim = unlist(
          constant_list[c('nrho', 'nspecies', 'ncity')]
        )
      ),
      city_shape_psi = runif(
        1,
        0.5,
        2
      ),
      city_shape_rho = runif(
        1,
        0.5,
        2
      ),
      city_rate_psi = runif(
        1,
        0.5,
        2
      ),
      city_rate_rho = runif(
        1,
        0.5,
        2
      ),
      city_tau_psi = rgamma(
        constant_list$ncity,
        1,
        1
      ),
      city_tau_rho = rgamma(
        constant_list$ncity,
        1,
        1
      ),
      ssc_psi = rnorm(
        constant_list$nseason_params
      ),
      ssc_rho = rnorm(
        constant_list$nseason_params
      )
    )
  }
  return(
    switch(
      chain,           
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
  )
}
