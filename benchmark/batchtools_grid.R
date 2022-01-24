source("_parameters.R")

# source functions

source("_functions.R")

## Define Algorithm

reg_dir =  "~/nobackup/bt_asd"
if (TESTMODE) {
  unlink(reg_dir, recursive = TRUE)
} else if (fs::dir_exists(reg_dir)) {
  reg = loadRegistry(reg_dir, writeable = TRUE)
}
reg = makeExperimentRegistry(file.dir = reg_dir, seed = 1, packages = c("asd.mod", "ParamHelpers", "smoof", "mlrMBO", "mlr3misc"))
if (reg$cluster.functions$name == "Interactive") {
  reg$cluster.functions = makeClusterFunctionsMulticore(3)
}

fun = addToEnvironment(fun, dictionaries = list(effects = EFFECTS))
addProblem("test1", data = fun)

addAlgorithm("eval", fun = function(job, data, instance, fun_x, n_cases, corr, nsim, ...) {
  dots = list(...)
  ps = getParamSet(data)
  des = do.call(data.frame, dots[getParamIds(ps)])
  fun_x = dfRowToList(df = des, par.set = ps, i = 1)
  dots = dots[setdiff(names(dots), names(des))]
  do.call(data, c(list(x = fun_x, n_cases = n_cases, corr = corr, nsim = nsim), dots))
})

addAlgorithm("mbo", fun = function(job, data, instance, corr, nsim, n_cases, mbo_iters, select_max_val, ...) {
  library(mlrMBO)
  fun = data
  ps = getParamSet(fun)
  param_d = map_int(ps$pars, function(x) max(length(x$values), 1L))
  mbo_des = generateDesign(n = 4L * sum(param_d) , par.set = ps)
  #mbo_des$y = vapply(ParamHelpers::dfRowsToList(df = mbo_des, par.set = ps), fun, numeric(1L))
  lrn.km = makeLearner("regr.km", covtype = "matern5_2", optim.method = "gen", control = list(trace = FALSE), scaling = FALSE, predict.type = "se", nugget.estim = TRUE)
  lrn.km = makeRemoveConstantFeaturesWrapper(lrn.km)
  lrn.km = makeDummyFeaturesWrapper(lrn.km)
  lrn.km = makeImputeWrapper(lrn.km, classes = list(numeric = imputeMax(2), factor = imputeConstant("<missing>")))

  impute_y = function(x, y, opt.path) {
    #mainly fails for stage_ratio close to 0 or 1 (no surpise)
    message(BBmisc::convertToShortString(x), "as x failed!", "\n")
    min(getOptPathY(opt.path), na.rm = TRUE)
  }

  mbo_ctrl = makeMBOControl(on.surrogate.error = "warn", final.method = "best.predicted", impute.y.fun = impute_y, final.evals = 20)
  mbo_ctrl = setMBOControlInfill(mbo_ctrl, crit = makeMBOInfillCritAEI(aei.use.nugget = TRUE))
  mbo_ctrl = setMBOControlTermination(mbo_ctrl, iters = mbo_iters)

  res = mbo(fun, control = mbo_ctrl, design = mbo_des, learner = lrn.km, more.args = list(corr = corr, nsim = nsim, n_cases = n_cases, ...))
  res$models = NULL
  res$final.opt.state = NULL
  res$control = NULL
  res
})

ades = generate_eval_mbo_design(n_cases = NCASES, nsim = NSIM, effect = names(EFFECTS))
addExperiments(algo.designs = ades, repls = REPLS)

# add 5times sim steps for one experiment
ades_nsim_high = generate_eval_mbo_design(n_cases = tail(NCASES,1), nsim = c(NSIM, 5 * NSIM), effect = "paper")
addExperiments(algo.designs = ades_nsim_high, repls = max(4 * REPLS, 100))

if (TESTMODE) {
  stop()
}

jdf = getJobTable()
jdf[algorithm == "eval", chunk := batchtools::chunk(job.id, chunk.size = 200)]
jdf[algorithm == "mbo", chunk := max(jdf$chunk, na.rm = TRUE) + (1:.N)]

submitJobs(jdf)
waitForJobs()

# submitJobs(jdf[algorithm == "mbo" & map(algo.pars, "effect") %in% names(EFFECTS),])

if (FALSE) {
  #find expired but reads
  ids_expired = jdf[findExpired(),]
  ids_expired[, is_ready := fs::file_exists(paste0(reg$file.dir, "/results/", job.id, ".rds"))]
  reg$status[ids_expired[is_ready == TRUE, .(job.id)], done := Sys.time()]
}

jdf_done = getJobTable(findDone())
res_eval = reduceResultsDataTable(jdf_done[algorithm == "eval", ], fun = function(res, ...) {
  extras = attr(res, "extras")
  extras = lapply(extras, function(x) {attributes(x) = NULL; x}) #clean extras
  c(y = res[[1]], extras)
})
res_eval = unwrap(res_eval, "result")
res_eval = ljoin(jdf_done[algorithm == "eval", ], res_eval)
saveRDS(res_eval, "../benchmark_results/batchtools_grid_res_eval.rds")

res_mbo = reduceResultsDataTable(jdf_done[algorithm == "mbo" & map(algo.pars, "effect") %in% names(EFFECTS),], fun = function(res) {
  opdf = as.data.table(res$opt.path)
  opdf[, c("error.model", "eol", "aei", "train.time", "propose.time", "se", "mean", "tau") := NULL]
  best_el = getOptPathEl(res$opt.path, index = res$best.ind)
  c(res$x, y = res$y, best_index = res$best.ind, best_el$extra[c("stage_1_arms", "stage_1_n","stage_2_arms", "stage_2_n")], opt.path = list(opdf))
})
res_mbo = ljoin(jdf_done[algorithm == "mbo", ], res_mbo)
saveRDS(res_mbo, "../benchmark_results/batchtools_grid_res_mbo.rds")

if (FALSE) {
  #debug calibration
  x = list(select = 6, stage_ratio = 2.622986e-05, epsilon = NA_real_, thresh = 0.00152911)
  n_cases = 1000
  effect = list(early=c(0,0.68,0.82,0.95,0.91),final=c(0,0.13,0.17,0.23,0.20))
  corr = 0.4
  nsim = 100
  
  opdf = as.data.table(res$opt.path)
  p = ggplot(data = opdf[select == "4",], mapping = aes(x = stage_ratio, y = epsilon, color = y))
  p + geom_point()
  
  p = ggplot(data = opdf[select == "6",], mapping = aes(x = stage_ratio, y = thresh, color = y))
  p + geom_point()
  
  p = ggplot(opdf, aes(x = dob, y = y, color = select))
  p + geom_point()
  
}
