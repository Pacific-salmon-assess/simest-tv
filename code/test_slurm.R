#====== 
#Test rslurm functions
#catarina Wor 
#Nov 2022
#--------

packages <- c("rslurm")
install.packages(setdiff(packages, rownames(installed.packages()))) 

for(package.i in packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
      )
    )
}

test_func <- function(par_mu, par_sd) {
	par_sd_vec <- rep(par_sd,length(par_mu))
    samp <- rnorm(10^6, par_mu, par_sd_vec)
    c(s_mu = mean(samp), s_sd = sd(samp))
}

pars <- list(par_mu = 1:10,
             par_sd = 0.3)


sjob <- slurm_map(pars,
	test_func, 
	jobname = 'test_apply',
    nodes = 2, cpus_per_node = 2, submit = FALSE)


cleanup_files(sjob)