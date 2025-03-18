# # test summary table
# setClass("horseshoeocc_mcmc", slots = list())
#
# setGeneric("summary")
#
# setMethod("summary",
#           signature(object = "horseshoeocc_mcmc"),
#           jags_summary)
#
# sum_test <- summary(test$mcmc)
#
# setClass("horseshoeocc", slots = list())
#
# setGeneric("summary")
#
# setMethod("summary",
#           signature(object = "horseshoeocc"),
#           jags_summary)
#
# sum_test <- summary(test)
#
# # test trace plots
# setGeneric(
#   "plot",
#   function(x, which = c("beta"), ...) standardGeneric("plot")
# )
#
# setMethod("plot",
#           signature(x = "horseshoeocc"),
#           trace_plot)
#
# setMethod("plot",
#           signature(x = "horseshoeocc_mcmc"),
#           trace_plot)
#
# plot(test)
# plot(test, "meff")
# plot(test, "alpha")
# plot(test, c("alpha", "beta"))
# plot(test, 1)
#
# plot(test$mcmc)
#
# # test credible intervals
# setClass("horseshoeocc_sum")
#
# setGeneric(
#   "plot",
#   function(x, which = c("beta"), hdi = TRUE, equal = NULL, ...) standardGeneric("plot")
# )
#
# setMethod("plot",
#           signature(x = "horseshoeocc_sum"),
#           cred_plot)
#
# plot(sum_test, which = "beta", hdi = T)
# plot(sum_test, c("meff"))
# plot(sum_test, c("alpha", "beta"), equal = 0 )
# plot(sum_test, 1, equal = 0)
#
# # test summarize and plot derived parameters
# setClass("horseshoeocc_derived")
#
# setGeneric("summary")
#
# setMethod("summary",
#           signature(object = "horseshoeocc_derived"),
#           derived_summary)
#
# testderive_sum <- summary(test_derive)
#
# setClass("horseshoeocc_derivedsum")
# setGeneric(
#   "plot",
#   function(x, which = c("psi"), hdi = FALSE, equal = NULL, ...) standardGeneric("plot")
# )
#
# setMethod("plot",
#           signature(x = "horseshoeocc_derivedsum"),
#           cred_derive)
#
# plot(testderive_sum, "psi")
#
# plots_p <- plot(testderive_sum, "p")
#
