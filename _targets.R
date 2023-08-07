# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    # argparse
    "DSS",
    "magrittr",
    "dplyr",
    "stringr",
    "parallel",
    "tidyverse"), 
  format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # For distributed computing in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller with 2 workers which will run as local R processes:
  #
  #   controller = crew::crew_controller_local(workers = 2)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package. The following
  # example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     workers = 50,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.0".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "multicore")

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

tar_source("R/")


# Replace the target list below with your own:
list(
  tar_target(master_df, clean_master_samplesheet("../DataRaw/masterSamplesheet.csv")),
  tar_target(valid_ids, get_valid_ids(master_df, "../DataRaw/MCovRawSplit/")),
  tar_target(control_ids, get_control_samples(master_df, valid_ids)),
  tar_target(load_ids, get_load_samples(master_df, valid_ids))

)

# bootstrap_params <- tribble(
#   ~bootstrap_ix, ~seed, 
#   1,  "a",  1,    
#   1,  "b",  2,
#   1,  "c",  3  
# ) 

# mapped_pipeline <- tar_map(
#   values = bootstrap_params, 
#   names = "name",
#   tar_target(filtered, filter_M_Cov()),
#   tar_target(prin_comps,
#              f(arg1, arg2)),
#   tar_target(second_target,
#              g(first_target)
#   )
# )

# tar_pipeline(pipeline)