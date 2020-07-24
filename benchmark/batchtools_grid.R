library(batchtools)
library(smoof)
library(asd.mod)
library(data.table)
library(mlr3misc)

TESTMODE = Sys.info()[["nodename"]] %in% c("JayBook2")

if (TESTMODE) {
  #try out params
  NSIM = 50
  MBOITERS = 10
  GRIDRES = 4
  REPLS = 1
  NCASES = c(1000,2000)
  EFFECTS = list(
    paper = list(early = c(0,0.68,0.82,0.95,0.91), final = c(0,0.13,0.17,0.23,0.20))
  )
} else {
  #real params
  NSIM = 1000
  MBOITERS = 100
  GRIDRES = 25
  REPLS = 10
  NCASES = c(500,1000,2000)
  EFFECTS = list(
    paper = list(early = c(0,0.68,0.82,0.95,0.91), final = c(0,0.13,0.17,0.23,0.20)),
    linear = {x = c(0,2,4,6,8)/10; list(early = x, final = x/4)},
    #threashold = {x = c(0,0,0,4,8)/10; list(early = x, final = x/4)},
    #saturation = {x = c(0,5,7,7.8,8)/10; list(early = x, final = x/4)},
    sigmoid = {x = c(0,1,2,7,8)/10; list(early = x, final = x/4)},
    #paper_mod = list(early = c(0,0.68,0.82,0.93,0.93), final = c(0,0.13,0.17,0.23,0.20)),
    paper2 = list(early = c(0,0.68,0.82,0.95,0.91), final = 2*c(0,0.13,0.17,0.23,0.20))# ,
    #paper_mod2 = list(early = c(0,0.68,0.82,0.93,0.93), final = 2*c(0,0.13,0.17,0.23,0.20))
  )
}

## Define Problem

fn = function(
  x, 
  n_cases = 1000, 
  effect = list(early=c(0,0.68,0.82,0.95,0.91),final=c(0,0.13,0.17,0.23,0.20)), 
  corr = 0.4, 
  nsim = 100) {

  if (is.character(effect)) {
    effect = dictionaries$effects[[effect]]
  }
  k_plus_1 = length(effect$early)
  p = x$stage_ratio

  calculate_n = function(n_select) { # calculate n for each arm in stage1 and stage2 given the number of arrms in stage2, and the ratio p
    const = n_cases/(k_plus_1 * p + (n_select+1) * (1-p))
    list(
      stage1 = round(p * const), 
      stage2 = round((1-p) * const)
    )
  }

  calculate_nselect = function(res) { #calculate the nselect from the results after the simulation ran
    if (!is.null(res$error)) {
      return(0)
    }
    as.numeric(sum(res$res.5$sim.res$count.total * seq_along(res$res.5$sim.res$count.total)) / res$res.1$simulations["n"])
  }

  # wrapper to safely run treatsel.sim (necessary for init design)
  treatsel_sim_safe = function(n, effect, nsim, corr, seed = sample(.Machine$integer.max, 1), select, level = 0.025, ptest = c(3, 4), epsilon = epsilon, thresh = thresh, ...) {
    capture.output({
      res = tryCatch({
        treatsel.sim(n = n, effect = effect, outcome = list(early = "N", final = "N"), nsim = nsim, corr = corr, seed = seed, select = select, level = level, ptest = ptest, epsilon = epsilon, thresh = thresh, ...)
        # on error return result that leads to power 0
      }, error = function(e) list(res.5 = list(sim.res = list(sim.reject = 0), sim.res.n = nsim), error = as.character(e)))
    })
    return(res)
  }
  if (x$select %in% c(0,1:3,5)) { # all go to second stage
    select2 = switch(as.character(x$select),
      "0" = k_plus_1-1,
      "1" = 1,
      "2" = 2,
      "3" = 3,
      "5" = 1) #1 random
  } else if (x$select %in% c(4,6)) { # select arms according to thresold of p value or epsilon # 4 and 6
    # first we have to calibrate / calculate the expected number of arms in stage 2
    calibrate_n = function(x) {
      calc_stage2 = function(stage1, p) { #calc n for stage 2 per branch
        ((1-p)/p) * stage1
      }
      # function to calculate how far off we are for the desired number of target treatments, given the patients per arm in stage 1 (x$stage1)
      fun_target_treatments = makeSingleObjectiveFunction(
        name = "target treatments", 
        has.simple.signature = FALSE, 
        fn = function(x, xx, target_treatments) {
          assert_int(target_treatments)
          n = list(stage1 = x$stage1, stage2 = 1) #stage2 is of no importance here
          res = treatsel_sim_safe(n = n, effect = effect, nsim = nsim, corr = corr, select = xx$select, epsilon = xx$epsilon, thresh = xx$thresh)
          mean_select = calculate_nselect(res)
          k_plus_1 = length(effect$early)
          # n = p * n + (1-p) * n mit p = xx$stage_ratio und p * n = x$stage1 => (1-p) * n = (1-p)/p * (p*n) 
          stage2 = ((1-xx$stage_ratio)/xx$stage_ratio) * x$stage1
          n_cases_expect = k_plus_1 * x$stage1 + (mean_select + 1) * calc_stage2(x$stage1, xx$stage_ratio)
          y = (n_cases_expect - target_treatments)^2
          attr(y, "extras") = list(mean_select = mean_select)
          return(y)
        },
        par.set = makeParamSet(makeIntegerParam("stage1", lower = ceiling(n_cases / length(effect$early) * 0.01), upper = ceiling(n_cases / length(effect$early)))),
        noisy = TRUE
      )
      
      ctrl = makeMBOControl(final.method = "best.predicted")
      ctrl = setMBOControlInfill(ctrl, crit = crit.aei)
      des = generateDesign(n = 8, getParamSet(fun_target_treatments))
      res = mbo(fun_target_treatments, more.args = list(xx = x, target_treatments = n_cases), control = ctrl, design = des)
      n = list(
        stage1 = res$x$stage1, 
        stage2 = round(calc_stage2(res$x$stage1, x$stage_ratio)), 
        mean_select = getOptPathEl(res$opt.path, getOptPathBestIndex(res$opt.path))$extra$mean_select)
      return(n)
    }
    calib_res = calibrate_n(x)
    select2 = calib_res$mean_select
  }
  n = calculate_n(select2)
  res = treatsel_sim_safe(n = n, effect = effect, nsim = nsim, corr = corr, select = x$select, epsilon = x$epsilon, thresh = x$thresh)
  y = as.numeric(res$res.5$sim.res$sim.reject / res$res.5$sim.res.n)
  attr(y, "extras") = list(
    stage_1_arms = k_plus_1, 
    stage_1_n = res$res.1$sample_sizes[["stage1"]] %??% NA_real_, 
    stage_2_arms = calculate_nselect(res) + 1, 
    stage_2_n = res$res.1$sample_sizes[["stage2"]] %??% NA_real_
  )
  return(y)
}

fun = makeSingleObjectiveFunction(
  "treatsel", 
  par.set = makeParamSet(
    makeDiscreteParam("select", values = c(0,1,2,3,4,6)),
    makeNumericParam("stage_ratio", lower = 0, upper = 1),
    makeNumericParam("epsilon", lower = 0, upper = 4, requires = quote(select == 4)),
    makeNumericParam("thresh", lower = 0, upper = 10, requires = quote(select == 6))
  ),
  fn = fn, 
  has.simple.signature = FALSE, 
  noisy = TRUE, 
  minimize = FALSE
)

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

addToEnvironment = function(f, ...) {
  dots = list(...)
  e = new.env()
  for (n in names(dots)) {
    e[[n]] = dots[[n]]
  }
  environment(f) = e
  return(f)
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

  mbo_ctrl = makeMBOControl(on.surrogate.error = "warn", final.method = "best.predicted", impute.y.fun = impute_y, final.evals = 10)
  mbo_ctrl = setMBOControlInfill(mbo_ctrl, crit = makeMBOInfillCritAEI(aei.use.nugget = TRUE))
  mbo_ctrl = setMBOControlTermination(mbo_ctrl, iters = mbo_iters)

  res = mbo(fun, control = mbo_ctrl, design = mbo_des, learner = lrn.km, more.args = list(corr = corr, nsim = nsim, n_cases = n_cases, ...))
  res$models = NULL
  res$final.opt.state = NULL
  res$control = NULL
  res
})

const_design = expand.grid(n_cases = NCASES, nsim = NSIM, corr = 0.4, effect = names(EFFECTS), stringsAsFactors = FALSE)
eval_design = generateGridDesign(getParamSet(fun), GRIDRES)
eval_design = subset(eval_design, stage_ratio != 1 & stage_ratio != 0)
eval_design = merge.data.frame(eval_design, const_design)
setDT(eval_design)

mbo_design = cbind(const_design, mbo_iters = MBOITERS)
ades = list(
  eval = eval_design,
  mbo = mbo_design
)

addExperiments(algo.designs = ades, repls = REPLS)

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

jdf = getJobTable(findDone())
res_eval = reduceResultsDataTable(jdf[algorithm == "eval", ], fun = function(res, ...) {
  extras = attr(res, "extras")
  extras = lapply(extras, function(x) {attributes(x) = NULL; x}) #clean extras
  c(y = res[[1]], extras)
})
res_eval = unwrap(res_eval, "result")
res_eval = ljoin(jdf[algorithm == "eval", ], res_eval)
saveRDS(res_eval, "batchtools_grid_res_eval.rds")

res_mbo = reduceResultsDataTable(jdf[algorithm == "mbo" & map(algo.pars, "effect") %in% names(EFFECTS),], fun = function(res) {
  opdf = as.data.table(res$opt.path)
  opdf[, c("error.model", "eol") := NULL]
  best_el = getOptPathEl(res$opt.path, index = res$best.ind)
  c(res$x, y = res$y, best_el$extra[c("stage_1_arms", "stage_1_n","stage_2_arms", "stage_2_n")], opt.path = list(opdf))
})
res_mbo = ljoin(jdf[algorithm == "mbo", ], res_mbo)
saveRDS(res_mbo, "batchtools_grid_res_mbo.rds")

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