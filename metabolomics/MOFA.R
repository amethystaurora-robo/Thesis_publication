'''
This file runs Multi-omics Factor Analysis. 
It looks at the transcriptomic and metabolomic layers and determines latent factors
which account for shared variation at both omics levels.

'''

library(data.table)
library('MOFA2')
library(reticulate)
library(ggplot2)
py_config()

combined_neg <- read.csv('mofa_neg_omics.csv')
combined_pos <- read.csv('mofa_pos_omics.csv')


MOFAobject <- create_mofa(combined_neg)
plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
data_opts$scale_groups <- TRUE
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$save_interrupted <- TRUE
tail(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"model_new.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile,use_basilisk=TRUE)


model <- load_model('model_new.hdf5')
plot_variance_explained(model, x="view", y="factor")
head(model@cache$variance_explained$r2_per_factor[[1]])

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]

plot_factor(model, 
            factor = 1:3)