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

addToEnvironment = function(f, ...) {
  dots = list(...)
  e = new.env()
  for (n in names(dots)) {
    e[[n]] = dots[[n]]
  }
  environment(f) = e
  return(f)
}

generate_eval_mbo_design = function(n_cases, nsim, effect, par_set = getParamSet(fun), grid_res = GRIDRES) {
  const_design = expand.grid(n_cases = n_cases, nsim = nsim, corr = 0.4, effect = effect, stringsAsFactors = FALSE)
  eval_design = generateGridDesign(par_set, grid_res)
  eval_design = subset(eval_design, stage_ratio != 1 & stage_ratio != 0)
  eval_design = merge.data.frame(eval_design, const_design)
  setDT(eval_design)
  
  mbo_design = cbind(const_design, mbo_iters = MBOITERS)
  list(
    eval = eval_design,
    mbo = mbo_design
  )
}
