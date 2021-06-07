set.seed(1) #actually nothing stochastic happens here

## ----setup, include=FALSE----------------------------------------------------- 
library(data.table)
library(ggplot2)
library(scales)
library(stringi)
library(knitr)
library(kableExtra)
library(dplyr)
library(mlr3misc)
library(grid)
library(gridExtra)
library(memoise)

root = rprojroot::find_root(rprojroot::is_git_root)
setwd(paste0(root,"/generator"))


FIG_WIDTH = 9
FIG_HEIGHT = 4.5
FNT_SMALL = 6
CACHE = cache_filesystem(path = "memoise")

plot_wrapper_uncached = function(name, fig.height = FIG_HEIGHT, fig.width = FIG_WIDTH, expr, crop = TRUE) {
  path = paste0("../generated/figures/", name, ".pdf")
  cairo_pdf(filename = path, width = fig.width, height = fig.height, onefile = TRUE)
  print(eval(expr))
  dev.off()
  #browser()
  if (crop) {
    knitr::plot_crop(path)
  }
}
plot_wrapper = memoise(plot_wrapper_uncached, cache = CACHE)

kable_to_text = function(text, name) {
  text = stri_replace_all_fixed(text, pattern = "\\label{tab:}", sprintf("\\label{tab:%s}", name))
  writeLines(text, con = paste0("../generated/tables/", name, ".tex"))
}

own_theme = function(...) {
  theme_light(...) + 
  theme(legend.margin=margin(0,0,0,0), strip.background = element_rect(fill = "white", color = "white"), strip.text = element_text(color = "grey10"))
}
theme_set(own_theme())


options(
  knitr.table.format = "latex",
  knitr.kable.NA = '', 
  knitr.kable.NaN = ''
)


## ----helpers------------------------------------------------------------------
# paper parameters
EFFECTS = list(
  paper = list(early = c(0,0.68,0.82,0.95,0.91), final = c(0,0.13,0.17,0.23,0.20)),
  linear = {x = c(0,2,4,6,8)/10; list(early = x, final = x/4)}, #yes
  # threshold = {x = c(0,0,0,4,8)/10; list(early = x, final = x/4)},
  # saturation = {x = c(0,5,7,7.8,8)/10; list(early = x, final = x/4)},
  sigmoid = {x = c(0,1,2,7,8)/10; list(early = x, final = x/4)}, #yes
  # paper_mod = list(early = c(0,0.68,0.82,0.93,0.93), final = c(0,0.13,0.17,0.23,0.20)), 
  paper2 = list(early = c(0,0.68,0.82,0.95,0.91), final = 2*c(0,0.13,0.17,0.23,0.20)), #yes
  # paper_mod2 = list(early = c(0,0.68,0.82,0.93,0.93), final = 2*c(0,0.13,0.17,0.23,0.20))
  paper_rev = list(early = c(0,0.13,0.17,0.23,0.20), final = c(0,0.68,0.82,0.95,0.91))
)
effect_names = setNames(names(EFFECTS), names(EFFECTS)) #to be populated with better names
main_effect_names = names(EFFECTS)[1:4]

select_labels = c("0" = "all", "1" = "1 best", "2" = "2 best", "3" = "3 best", "4" = "epsilon rule", "5" = "random", "6" = "threshold rule")
select_labels_colors = #http://paletton.com/#uid=74Z140kqdu7ghF3lowyu1ppwGk7
  c(
    "all" = "#A61E96",
    "1 best" = "#FF9B55",
    "2 best" = "#CA5B0D",
    "3 best" = "#A04200",
    "epsilon rule" = "#1B9786",
    "random" = "#000000",
    "threshold rule" = "#C0E62A"
  )
algorithm_labels = c("mbo" = "BO", "eval" = "Grid", "mbo grid" = "BO Grid", "grid" = "Grid", "grid7" = "Grid Small")
algorithm_labels_color = c("mbo" = "#F04A2B", "eval" = "#099438", "mbo grid" = "#690D86", "grid" = "#099438", "grid7" = "#44CA6D") #https://paletton.com/#uid=6060R0krPwadDMilcBuwcpxDti6

ncases_labels = c(500,1000,2000)
ncases_labels = paste0("n[cases] ", ncases_labels)
names(ncases_labels) = c(500,1000,2000)

names(effect_labels) = effect_labels = c("linear", "paper", "paper2", "sigmoid")



## ----table_effect_names-------------------------------------------------------
tmp = lapply(EFFECTS[main_effect_names], do.call, what = rbind)
tmp = lapply(tmp, function(x) cbind.data.frame(stage = rownames(x), x))
tmp = lapply(tmp, function(x) {colnames(x) = c("stage", seq_len(ncol(x)-1)-1); x}) # colnames start at 0
tmp = rbindlist(tmp, idcol = TRUE)
setnames(tmp, ".id", "Effect Set")
setnames(tmp, "stage", "Stage")
kable(tmp, booktabs = TRUE, caption = "Effect sizes used for simulation in the first (early) and second (final) stage", label = "table_effect_names") %>% 
  collapse_rows(columns = 1:2, latex_hline = "major") %>% 
  add_header_above(c(" " = 2, "Treatment" = 5)) %>% 
  kable_styling(position = "center") %>% 
  kable_to_text("table_effect_names")

## ----read_data----------------------------------------------------------------
read_results_uncached = function(select_labels, case = "default", effect_names = main_effect_names) {
  if (case == "default") {
    case_selector = quote(nsim == 1000 & repl <= 20)
  } else if (case == "nsim") {
    case_selector = quote(n_cases == 2000 & effect == "paper")
  } else {
    stop("not supported case")
  }
  
  read_res_eval_uncached = function() {
    readRDS("../benchmark_results/batchtools_grid_res_eval.rds")
  }
  read_res_eval = memoise(read_res_eval_uncached)
  
  read_res_mbo_uncached = function() {
    readRDS("../benchmark_results/batchtools_grid_res_mbo.rds")
  }
  read_res_mbo = memoise(read_res_mbo_uncached)
  
  res_eval = read_res_eval()
  res_eval = res_eval[,.(job.id, y, stage_1_arms, stage_1_n, stage_2_arms, stage_2_n, repl, problem, prob.pars, algorithm, algo.pars, time.running)]
  res_eval = res_eval[map_chr(algo.pars, "effect") %in% effect_names,]

  algo.par.names = c("algorithm", names(res_eval$algo.pars[[1]]))
  algo.par.names.meta = c("algorithm", "effect", "corr", "nsim", "n_cases")
  algo.par.names.optim = setdiff(algo.par.names, algo.par.names.meta)

  res_eval = unnest(res_eval, cols = c("prob.pars", "algo.pars"))
  res_eval = res_eval[eval(case_selector),]
  #res_eval = res_eval[repl %in% 1:10,] #to run faster #FIXME remove for real results

  res_mbo = read_res_mbo()
  res_mbo = res_mbo[,.(job.id, result, repl, problem, prob.pars, algorithm, algo.pars, time.running)]
  res_mbo = res_mbo[map_chr(algo.pars, "effect") %in% effect_names,]
  res_mbo = unnest(res_mbo, cols = c("result", "prob.pars", "algo.pars"))
  res_mbo = res_mbo[eval(case_selector),]
  
  res_eval[, select := factor(select, levels = names(select_labels), labels = select_labels)]
  res_mbo[, select := factor(select, levels = names(select_labels), labels = select_labels)]
  
  ## ----estimate res_grid--------------------------------------------------------
  # we calculate res_grid which emulates a grid search with an independent validation of the best found x value
  # 1. for each repl we look for the best y and take the x values (randomly 1 on ties)
  # 2. we obtain the y values of the given x-values of the other 19 repls
  # 3. we average the external y values
  # 4. we sum up the time that the grid needed on one replication
  calc_res_grid_internal = function(.SD) {
    optim_dt = copy(.SD)
    this_repl = optim_dt$repl[1]
    optim_dt = optim_dt[y == max(y),]
    if (nrow(optim_dt) > 1) optim_dt = optim_dt[sample.int(nrow(optim_dt), 1),] # in case of multiple best
    valid_dt = res_eval[repl != this_repl, ]
    res = merge(optim_dt[,-"y"], valid_dt[, c(algo.par.names, "y"), with = FALSE], by = algo.par.names)
    # write result into x values dt
    optim_dt$y_valid = mean(res$y)
    #strangely we have 3 to 4 missing values here, probably batchtools had a hick up
    times = .SD$time.running
    if (any(is.na(times)) && mean(is.na(times)) < 0.01) {
      times[is.na(times)] = mean(times, na.rm = TRUE)
    } else if (any(is.na(times))) {
      # investigate please!
      browser()
    }
    optim_dt$time.running = sum(times)
    optim_dt[, c(algo.par.names.meta, "repl") := NULL]
    return(optim_dt)
  }
  res_grid = res_eval[, calc_res_grid_internal(.SD), 
    by = c(algo.par.names.meta, "repl"), 
    .SDcols = c("y", "repl", algo.par.names, "time.running")]
  res_grid[, algorithm := "grid"]
  
  # calculate grid 7 result (subset of res_grid)
  grid7 = local({
    library(smoof)
    source("../benchmark/_functions.R")
    full_design = generateGridDesign(getParamSet(fun), 25) # instead of just saying 7, we create a set with the index using elem 2, 5, 9, 13, 17, 21 and 25. Not the first because it is often 0!
    num_cols = map_lgl(full_design, is.numeric)
    for (num_col in names(num_cols[num_cols])) {
      col = full_design[[num_col]]
      keep = sort(unique(col))[c(2, seq(1,25, by = 4)[-1])] #every 4th element but start with 2 instead of 1 because 1 often is unfeasible
      full_design = full_design[is.na(full_design[[num_col]]) | full_design[[num_col]] %in% keep, ]
    }
    full_design
  })
  setDT(grid7)
  grid7[, select := factor(select, levels = names(select_labels), labels = select_labels)]
  res_eval7 = merge(res_eval, grid7, by = colnames(grid7), all.y = FALSE) # NA actually matches NA
  res_grid7 = res_eval7[, calc_res_grid_internal(.SD), 
                      by = c(algo.par.names.meta, "repl"), 
                      .SDcols = c("y", "repl", algo.par.names, "time.running")]
  res_grid7[, algorithm := "grid7"]
  
  
  ## calculate y result on grid
  fct.vars = c("select", "effect", "n_cases", "nsim", "corr") # first match where fcts are equal
  num.vars = setdiff(algo.par.names, c(fct.vars, "algorithm")) # knn on numeric rest
  res_mbo[, y_grid := {
    mbo_eval = .SD
    matches = batchtools::rjoin(res_eval, mbo_eval[, fct.vars, with=FALSE]) # look for grid results with same factors
    num.matches = matches[, num.vars, with=FALSE] # just take numeric columns
    non.na.cols = map_lgl(num.matches, function(x) !all(is.na(x))) # just take not na + numeric columsn
    non.na.cols = names(non.na.cols)[non.na.cols]
    mbo_eval[, non.na.cols, with = FALSE]
    knn_res = FNN::get.knnx(data = num.matches[, non.na.cols, with = FALSE], query = mbo_eval[, non.na.cols, with = FALSE], k = 10) #find 10=max_repl next evals on grid
    if (length(unique(knn_res$nn.dist[1,])) != 1) {
      stop ("knn found different results. Should not happen!")
    }
    mean(matches[knn_res$nn.index[1,],]$y)
  }, by = list(rownames(res_mbo))]


  ## ----correct res_mbo----------------------------------------------------------
  # we remove the runtime of the validation (dob = 101) from the total "time.running"
  # determine max dob (probably 101) because it contains outside mbo evals
  max_dob = max(map_int(map(res_mbo$opt.path, "dob"), max))

  res_mbo[, time.running := time.running - as.difftime(sum(.SD$opt.path[[1]][prop.type == "final_eval", ]$exec.time), units = "secs") , by = list(rownames(res_mbo))]
  data.table::setnames(res_mbo, "y", "y_valid")
  res_mbo[, y := .SD$opt.path[[1]][.SD$best_index, ]$y, by = list(rownames(res_mbo))]
  

  list(res_mbo = res_mbo, res_grid = res_grid, res_grid7 = res_grid7, res_eval = res_eval, algo.par.names = algo.par.names, algo.par.names.meta = algo.par.names.meta, algo.par.names.optim = algo.par.names.optim, max_dob = max_dob)
}

read_results = memoise(read_results_uncached, cache = CACHE)
tmp = read_results(select_labels)
algo.par.names = tmp$algo.par.names
algo.par.names.meta = tmp$algo.par.names.meta
algo.par.names.optim = tmp$algo.par.names.optim
max_dob = tmp$max_dob

get_res_grid = function(case = "default", ...) {
  read_results(select_labels, case, ...)$res_grid
}

get_res_grid7 = function(case = "default", ...) {
  read_results(select_labels, case, ...)$res_grid7
}

get_res_mbo = function(case = "default", ...) {
  read_results(select_labels, case, ...)$res_mbo
}

get_res_eval = function(case = "default", ...) {
  read_results(select_labels, case, ...)$res_eval
}

get_res_mbogrid = function(case = "default", ...) {
  tmp_mbogrid = get_res_mbo(case = case, ...)
  tmp_mbogrid[, algorithm := "mbo grid"]
  tmp_mbogrid[, y_valid := NULL]
  setnames(tmp_mbogrid, "y_grid", "y_valid")
}

get_res_eval_average = function(case = "default", ...) {
  get_res_eval(case = case, ...)[, list(mean_y = mean(y)), by = algo.par.names]
}

## ----put_ext_eval_mbo_results, eval = FALSE-----------------------------------
## res_mbo[, y_tuning := y]
## map_dbl(res_mbo$opt.path, function(op) mean(op[dob == max(dob), y])) == res_mbo$y_tuning
## map_dtr(res_mbo$opt.path, function(op) {
##   y = op[dob == max(dob), y]
##   data.table(y = mean(y), ysd = sd(y))
## })
## res_mbo[, .(y,ysd) := map_dtr(opt.path, function(op) {
##   y = op[dob == max(dob), y]
##   c(mean(y), sd(y))
## })]
## res_mbo[, y - y_tuning]

## ----plot_allbest, fig.height = 1.5 * FIG_HEIGHT------------------------------
# find best performing solution for each select methods and across all epsilon and threshold choices
plot_wrapper(name = "plot_allbest", fig.height = 1.5 * FIG_HEIGHT, expr = {
  res_ave = get_res_eval_average()
  res_best = res_ave[, .SD[order(-mean_y)[1], ], by = c("select", algo.par.names.meta)]
  res_best = res_best[, !c("mean_y", "stage_ratio")]
  # merge so the plot only contains the best curves (applies for epsilon and threshold rule)
  df = merge(get_res_eval(), res_best, all.x = FALSE, all.y = TRUE, by = colnames(res_best))
  dfmean = df[, lapply(.SD, mean),by = c(algo.par.names), .SDcols = c("y", "stage_1_arms", "stage_1_n")] #stage_1_n can vary a little
  g = ggplot(df, aes(x = (stage_1_arms * stage_1_n), y = y, color = select, group = paste(select, epsilon, thresh)))
  g = g + geom_point(aes(size = algorithm, shape = algorithm), data = get_res_mbo(), shape = 16)
  g = g + geom_line(data = dfmean)
  g = g + geom_point(aes(size = algorithm, shape = algorithm), alpha = 0.15, shape = 16)
  g = g + scale_color_manual(values = select_labels_colors)
  g = g + scale_size_manual(labels = algorithm_labels[c("mbo" ,"eval")], values = c("mbo" = 2, "eval" = 0.5))
  g = g + facet_grid(effect~n_cases, scales = "free", labeller = label_bquote(cols = {n[total]==.(n_cases)}))
  g = g + labs(x = expression(k[1] %.% n[stage1]))
  g = g + theme(legend.position = "bottom")
  g
})

## ----plot_appendix_allbest_paper_rev------------------------------
# find best performing solution for each select methods and across all epsilon and threshold choices
plot_wrapper(name = "plot_appendix_allbest_paper_rev", fig.height = FIG_HEIGHT, expr = {
  res_ave = get_res_eval_average(effect_names = c("paper", "paper_rev"))
  res_best = res_ave[, .SD[order(-mean_y)[1], ], by = c("select", algo.par.names.meta)]
  res_best = res_best[, !c("mean_y", "stage_ratio")]
  # merge so the plot only contains the best curves (applies for epsilon and threshold rule)
  df = merge(get_res_eval(effect_names = c("paper", "paper_rev")), res_best, all.x = FALSE, all.y = TRUE, by = colnames(res_best))
  dfmean = df[, lapply(.SD, mean),by = c(algo.par.names), .SDcols = c("y", "stage_1_arms", "stage_1_n")] #stage_1_n can vary a little
  g = ggplot(df, aes(x = (stage_1_arms * stage_1_n), y = y, color = select, group = paste(select, epsilon, thresh)))
  g = g + geom_point(aes(size = algorithm, shape = algorithm), data = get_res_mbo(effect_names = c("paper", "paper_rev")), shape = 16)
  g = g + geom_line(data = dfmean)
  g = g + geom_point(aes(size = algorithm, shape = algorithm), alpha = 0.15, shape = 16)
  g = g + scale_color_manual(values = select_labels_colors)
  g = g + scale_size_manual(labels = algorithm_labels[c("mbo" ,"eval")], values = c("mbo" = 2, "eval" = 0.5))
  g = g + facet_grid(effect~n_cases, scales = "free", labeller = label_bquote(cols = {n[total]==.(n_cases)}))
  g = g + labs(x = expression(k[1] %.% n[stage1]))
  g = g + theme(legend.position = "bottom")
  g
})


## ----plot_best_x, fig.height=1.5*FIG_HEIGHT-----------------------------------
# find best from grid 
#g = ggplot(data = res_mbo, aes(x = (stage_1_arms * stage_1_n), y = y, color = select))
plot_wrapper(name = "plot_best_x", fig.height = 1.5 * FIG_HEIGHT, expr = {
  tmp = get_res_grid()
  tmp = rbind(tmp, get_res_mbo()[, colnames(tmp), with = FALSE])
  g = ggplot(data = tmp, aes(x = stage_ratio, y = y_valid, color = select, shape = algorithm))
  g = g + geom_point(size = 3)
  #g = g + geom_text(data = res_mbo[nsim == 1000 & repl <= 10 & select == "epsilon rule", ], aes(label = round(epsilon,2)), hjust = 0, vjust = 1, show.legend = FALSE)
  #g = g + geom_text(data = res_mbo[nsim == 1000 & repl <= 10 & select == "threshold rule", ], aes(label = round(thresh,2)), hjust = 0, vjust = 1, show.legend = FALSE)
  g = g + facet_wrap(effect~n_cases, scales = "free", ncol = 3) # labels not used
  g = g + scale_x_continuous(expand = expansion(mult = 0.2))
  g = g + scale_y_continuous(expand = expansion(mult = 0.2))
  g = g + scale_shape_manual(labels = algorithm_labels, values = c("mbo" = 16, "grid" = 21))
  g = g + scale_color_manual(labels = select_labels, values = select_labels_colors)
  g = g + labs(x = expression(r), y = expression(y[valid]))
  g = g + theme(legend.position = "bottom", strip.background = element_blank(), strip.text.x = element_blank())
  #grid headlines
  col_heads = lapply(unique(tmp$n_cases), function(x) bquote(n[total]==.(x))) %>% do.call("expression",.) %>% lapply(textGrob, gp = gpar(fontsize = 10, color = "grey10"))
  row_heads = levels(factor(tmp$effect)) %>% lapply(textGrob, gp = gpar(fontsize = 10, col = "grey10"), rot=90*3)
  layout_mat = matrix(c(
    1, 2, 3, NA,
    8, 8, 8,  4,
    8, 8, 8,  5,
    8, 8, 8,  6,
    8, 8, 8,  7,
    8, 8, 8, NA
  ), byrow = TRUE, ncol = 4)
  grid.arrange(grobs=c(col_heads, row_heads, list(g)),layout_matrix=layout_mat, widths=c(10,10,10,1), heights=c(1,5,5,5,5,3))
})

## ----plot_appendix_best_x_paper_rev -----------------------------------
# find best from grid 
#g = ggplot(data = res_mbo, aes(x = (stage_1_arms * stage_1_n), y = y, color = select))
plot_wrapper(name = "plot_appendix_best_x_paper_rev", fig.height = FIG_HEIGHT, expr = {
  tmp = get_res_grid(effect_names = c("paper", "paper_rev"))
  tmp = rbind(tmp, get_res_mbo(effect_names = c("paper", "paper_rev"))[, colnames(tmp), with = FALSE])
  g = ggplot(data = tmp, aes(x = stage_ratio, y = y_valid, color = select, shape = algorithm))
  g = g + geom_point(size = 3)
  g = g + facet_wrap(effect~n_cases, scales = "free", ncol = 3) # labels not used
  g = g + scale_x_continuous(expand = expansion(mult = 0.2))
  g = g + scale_y_continuous(expand = expansion(mult = 0.2))
  g = g + scale_shape_manual(labels = algorithm_labels, values = c("mbo" = 16, "grid" = 21))
  g = g + scale_color_manual(labels = select_labels, values = select_labels_colors)
  g = g + labs(x = expression(r), y = expression(y[valid]))
  g = g + theme(legend.position = "bottom", strip.background = element_blank(), strip.text.x = element_blank())
  #grid headlines
  col_heads = lapply(unique(tmp$n_cases), function(x) bquote(n[total]==.(x))) %>% do.call("expression",.) %>% lapply(textGrob, gp = gpar(fontsize = 10, color = "grey10"))
  row_heads = levels(factor(tmp$effect)) %>% lapply(textGrob, gp = gpar(fontsize = 10, col = "grey10"), rot=90*3)
  layout_mat = matrix(c(
    1, 2, 3, NA,
    6, 6, 6,  4,
    6, 6, 6,  5,
    6, 6, 6, NA
  ), byrow = TRUE, ncol = 4)
  grid.arrange(grobs=c(col_heads, row_heads, list(g)),layout_matrix=layout_mat, widths=c(10,10,10,1), heights=c(1,5,5,3))
})


## ---- plot_opt_path -----------------------------

# calculate y perf of mbo runs
plot_wrapper(name = "plot_opt_path", fig.height = 1.6 * FIG_HEIGHT, expr = {
  df = get_res_mbo(effect_names = names(EFFECTS))
  common_names = intersect(colnames(df), colnames(df$opt.path[[1]]))
  data.table::setnames(df, common_names, paste0("opt.", common_names))
  df = tidyr::unnest(df, "opt.path")
  setDT(df)
  df = df[prop.type != "final_eval", cummax_y := cummax(y), by = c(algo.par.names.meta, "repl")]
  
  # calculate theoretical best y from grid
  res_ave = get_res_eval_average(effect_names = names(EFFECTS))
  df_best = res_ave[,.SD[order(-mean_y)[1]], by = c("effect", "n_cases", "algorithm")][,.(effect, n_cases, algorithm, best_y = mean_y)]
  df_best[algorithm == "eval", algorithm := "grid"]
  df = merge(df_best[, -"algorithm"], df, by = c("effect", "n_cases"))
  
  dfmean = df[, list(cummax_y_mean = mean(cummax_y), y_mean = mean(y), y_sd = sd(y)), by = c(algo.par.names.meta, "dob", "best_y")]
  
  g = ggplot(df[dob > 0 & prop.type != "final_eval", ], aes(x = dob, y = cummax_y-best_y, color = algorithm))
  g = g + geom_line(alpha = 0.3, aes(group = paste0(repl)))
  g = g + geom_line(data = dfmean, aes(y = cummax_y_mean-best_y))
  g = g + facet_wrap(effect~n_cases, scales = "free_y", labeller = label_both, ncol = 3)
  g = g + geom_hline(yintercept = 0 , color = colorspace::darken(algorithm_labels_color[["grid"]], amount = 0.6))
  g = g + geom_label(data = df_best, aes(label = formatC(best_y, 4)), y = 0, x = 90, vjust = 1.5, show.legend = FALSE)
  g = g + scale_color_manual(labels = algorithm_labels, values = algorithm_labels_color)
  g = g + theme(legend.position = "bottom", strip.background = element_blank(), strip.text.x = element_blank())
  g = g + labs(x = "iteration", y = expression(y[iter]*"*" - y[grid]*"*"), color = NULL)
  g = g + guides(colour = guide_legend(override.aes = list(size=2)))
  #grid headlines
  col_heads = lapply(unique(df$n_cases), function(x) bquote(n[total]==.(x))) %>% do.call("expression",.) %>% lapply(textGrob, gp = gpar(fontsize = 10, color = "grey10"))
  row_heads = levels(factor(df$effect)) %>% lapply(textGrob, gp = gpar(fontsize = 10, col = "grey10"), rot=90*3)
  layout_mat = matrix(c(
    1, 2, 3, NA,
    9, 9, 9,  4,
    9, 9, 9,  5,
    9, 9, 9,  6,
    9, 9, 9,  7,
    9, 9, 9,  8,
    9, 9, 9, NA
  ), byrow = TRUE, ncol = 4)
  grid.arrange(grobs=c(col_heads, row_heads, list(g)),layout_matrix=layout_mat, widths=c(10,10,10,1), heights=c(1,5,5,5,5,5,2.5))
})

create_boxplot_df = function(case = "default", ...) {
  tmp = get_res_grid(case = case, ...)
  rbind(
    tmp, 
    get_res_grid7(case = case, ...)[, colnames(tmp), with = FALSE],
    get_res_mbo(case = case, ...)[, colnames(tmp), with = FALSE], 
    get_res_mbogrid(case = case, ...)[, colnames(tmp), with = FALSE]
  )
}
## ----plot_boxplot_valid_y-----------------------------------------------------
plot_wrapper(name = "plot_boxplot_valid_y", fig.height = FIG_HEIGHT * 0.9, expr = {
  tmp = create_boxplot_df()
  g = ggplot(tmp, aes(x = algorithm, y = y_valid, color = algorithm, fill = algorithm))
  # mylabels = function(labels) {
  #   do.call(map, args = c(list(paste), labels))
  #   map(paste, labels[[1]], labels[[2]])
  #   paste(value, variable)
  # }
  
  g = g + facet_wrap(effect~n_cases, scales = "free_y", ncol = 6, labeller = label_bquote({atop(.(effect), n[total]==.(n_cases))}))
  darker_colors = colorspace::darken(algorithm_labels_color, amount = 0.6)
  names(darker_colors) = names(algorithm_labels_color)
  g = g + scale_fill_manual(labels = algorithm_labels, values = algorithm_labels_color) + scale_color_manual(labels = algorithm_labels, values = darker_colors)
  g = g + scale_x_discrete(labels = algorithm_labels)
  g = g + geom_boxplot() + theme(legend.position = "right")
  g = g + labs(x = NULL, y = expression(y[valid]), color = NULL, fill = NULL)
  g = g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  g
})

## ----plot_boxplot_valid_y_5000-----------------------------------------------------
plot_wrapper(name = "plot_boxplot_valid_y_5000", fig.height = 1.6 * FIG_HEIGHT * 0.35, fig.width = 0.35 * FIG_WIDTH, expr = {
  tmp = create_boxplot_df(case = "nsim")
  g = ggplot(tmp, aes(x = as.factor(nsim), y = y_valid, color = algorithm, fill = algorithm))
  #g = g + facet_grid(.~nsim, scales = "free", labeller = label_both)
  darker_colors = colorspace::darken(algorithm_labels_color, amount = 0.6)
  names(darker_colors) = names(algorithm_labels_color)
  g = g + scale_fill_manual(labels = algorithm_labels, values = algorithm_labels_color) + scale_color_manual(labels = algorithm_labels, values = darker_colors)
  g = g + theme(legend.position = "bottom")
  g = g + geom_boxplot()
  g = g + labs(x = expression(n[sim]), y = expression(y[valid]), color = NULL, fill = NULL)
  g = g + guides(fill=guide_legend(nrow=2,byrow=FALSE))
  g
})

## ----plot_appendix_boxplot_paper_rev-----------------------------------------------------
plot_wrapper(name = "plot_appendix_boxplot_paper_rev", fig.height = FIG_HEIGHT * 0.7, fig.width = FIG_WIDTH, expr = {
  tmp = create_boxplot_df(case = "default", effect_names = c("paper", "paper_rev"))
  g = ggplot(tmp, aes(x = algorithm, y = y_valid, color = algorithm, fill = algorithm))
  g = g + facet_wrap(effect~n_cases, scales = "free_y", ncol = 6, labeller = label_bquote({atop(.(effect), n[total]==.(n_cases))}))
  darker_colors = colorspace::darken(algorithm_labels_color, amount = 0.6)
  names(darker_colors) = names(algorithm_labels_color)
  g = g + scale_fill_manual(labels = algorithm_labels, values = algorithm_labels_color) + scale_color_manual(labels = algorithm_labels, values = darker_colors)
  g = g + scale_x_discrete(labels = algorithm_labels)
  g = g + geom_boxplot()
  g = g + labs(x = NULL, y = expression(y[valid]), color = NULL, fill = NULL)
  g = g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  g
})

## ----plot_opt_path_5000-----------------------------

# calculate y perf of mbo runs
plot_wrapper(name = "plot_opt_path_5000", fig.height = 1.6 * FIG_HEIGHT * 0.35, fig.width = 0.64 * FIG_WIDTH, expr = {
  df = get_res_mbo(case = "nsim")
  common_names = intersect(colnames(df), colnames(df$opt.path[[1]]))
  data.table::setnames(df, common_names, paste0("opt.", common_names))
  df = tidyr::unnest(df, "opt.path")
  setDT(df)
  df = df[prop.type != "final_eval", cummax_y := cummax(y), by = c(algo.par.names.meta, "repl")]
  
  # calculate theoretical best y from grid
  res_ave = get_res_eval_average(case = "nsim")
  df_best = res_ave[,.SD[order(-mean_y)[1]], by = algo.par.names.meta][, c(algo.par.names.meta, "mean_y"), with = FALSE]
  data.table::setnames(df_best, old = "mean_y", new = "best_y")
  df_best[algorithm == "eval", algorithm := "grid"]
  df = merge(df_best[, -"algorithm"], df, by = setdiff(algo.par.names.meta, "algorithm"))
  
  dfmean = df[, list(cummax_y_mean = mean(cummax_y), y_mean = mean(y), y_sd = sd(y)), by = c(algo.par.names.meta, "dob", "best_y")]
  
  g = ggplot(df[dob > 0 & prop.type != "final_eval", ], aes(x = dob, y = cummax_y-best_y, color = algorithm))
  g = g + geom_line(alpha = 0.1, aes(group = paste0(repl)), size = 0.5)
  g = g + geom_line(data = dfmean, aes(y = cummax_y_mean-best_y), color = colorspace::darken(algorithm_labels_color[["mbo"]], amount = 0.1))
  g = g + facet_grid(.~nsim, scales = "free")
  g = g + geom_hline(yintercept = 0 , color = colorspace::darken(algorithm_labels_color[["grid"]], amount = 0.1))
  g = g + geom_label(data = df_best, aes(label = formatC(best_y, 4)), y = 0, x = 90, vjust = 1.5, show.legend = FALSE)
  g = g + scale_color_manual(labels = algorithm_labels, values = algorithm_labels_color)
  g = g + theme(legend.position = "bottom")
  g = g + labs(x = "iteration", y = expression(y[iter]*"*" - y[grid]*"*"), color = NULL)
  g = g + guides(colour = guide_legend(override.aes = list(size=1, alpha = 1)))
  g
})


## ----table_best---------------------------------------------------------------
#best of grid
#FIXME: calculate time of complete grid, devide through 10
if (FALSE) { # this table does not make much sense
  res_ave = get_res_eval_average()
  df = res_ave[,.SD[order(-mean_y)[1:3]], by = c("effect", "n_cases")][,.(effect, n_cases, select, stage_ratio, epsilon, mean_y)]
  #mbo average
  res_ave_mbo = get_res_mbo()[, list(mean_y = mean(y)), by = algo.par.names.meta]
  df = rbind(df, res_ave_mbo[, .(effect, n_cases, mean_y, select = algorithm)], fill = TRUE)
  setkey(df, effect, n_cases, mean_y)
  knitr::kable(df, booktabs = TRUE, caption = "Best configurations per ncases and effects", longtable = FALSE, digits = 3, label = "table_best") %>% # linesep = c(rep("",4), "\\addlinespace") 
    kable_styling(position = "center", font_size = FNT_SMALL) %>% 
    collapse_rows(1:2, latex_hline = "custom", custom_latex_hline = 1:2) %>% 
    kable_to_text("table_best")  
}



## ----table_time---------------------------------------------------------------
#table(res_eval$n_cases, res_eval$effect)
tmp = get_res_grid()
tmp = rbind(tmp, get_res_mbo()[, colnames(tmp), with = FALSE], get_res_grid7()[, colnames(tmp), with = FALSE])
tmp = tmp[, list(time.running = mean(as.numeric(time.running, unit = "hours"), na.rm = TRUE)), by = c("algorithm", "effect", "n_cases")]
setkeyv(tmp, c("algorithm", "effect", "n_cases"))
tmp[, time.running := round(time.running, 1)]
tmp2 = tidyr::pivot_wider(tmp, id_cols = "algorithm", names_from = c("effect", "n_cases"), values_from = "time.running")
tmp2$algorithm = dplyr::recode(tmp2$algorithm, !!!algorithm_labels)
tmp2 = cbind(tmp2, data.frame("evals" = c(1350,126,116)))
kable(tmp2, 
  booktabs = TRUE, 
  col.names = c("", as.character(tmp[algorithm == "mbo"]$n_cases), "\\emph{evals}"),
  caption = "Average runtime in hours, for evaluating one grid and one optimization run of MBO",
  escape = FALSE,
  label = "table_time"
) %>% 
  add_header_above(c(" " = 1, table(tmp$effect)/3, " " = 1)) %>% 
  kable_styling(position = "center", latex_options = "scale_down") %>% 
  kable_to_text("table_time")

# res_grid[effect == "sigmoid" & n_cases == 1000, ]

# 
# g = ggplot(tmp, aes(x = algorithm, y = as.numeric(time.running, unit = "hours"), color = algorithm, fill = algorithm))
# g = g + facet_wrap(effect~n_cases, scales = "free", labeller = label_both, ncol = 3)
# darker_colors = colorspace::darken(algorithm_labels_color, amount = 0.6)
# names(darker_colors) = names(algorithm_labels_color)
# g = g + scale_fill_manual(values = algorithm_labels_color) + scale_color_manual(values = darker_colors)
# #g = g + scale_y_log10()
# g + geom_boxplot()


## ----debug, eval = FALSE, include=FALSE---------------------------------------
## res_eval[effect == "sigmoid" & n_cases == 1000,
##   {
##     optim_dt = copy(.SD)
##     this_repl = optim_dt$repl[1]
##     optim_dt = optim_dt[sample(which(y == max(y)), 1),]
##     valid_dt = res_eval[repl != this_repl, ]
##     res = merge(optim_dt[,-"y"], valid_dt[, c(algo.par.names, "y"), with = FALSE], by = algo.par.names)
##     # write result into x values dt
##     optim_dt$y = mean(res$y)
##     #strangely we have 3 to 4 missing values here, probably batchtools had a hick up
##     times = .SD$time.running
##     if (any(is.na(times)) && mean(is.na(times)) < 0.004) {
##       times[is.na(times)] = mean(times, na.rm = TRUE)
##     }
##     optim_dt$time.running = sum(times)
##     optim_dt[, c(algo.par.names.meta, "repl") := NULL]
##     optim_dt
##   },
##   by = c(algo.par.names.meta, "repl"),
##   .SDcols = c("y", "repl", algo.par.names, "time.running")]

# --- optimal epsilon values ------

plot_wrapper(name = "plot_appendix_opt_epsilon", fig.height = 1.5 * FIG_HEIGHT, fig.width = 0.5 * FIG_WIDTH, crop = FALSE, expr = {
  tmp = get_res_eval_average(effect_names = names(EFFECTS))
  tmp = tmp[select == "epsilon rule",]
  tmp[, mean_y := BBmisc::normalize(x = mean_y, method = "range"), by = .(effect, n_cases)]
  best = tmp[, .SD[mean_y == max(mean_y)], by = .(effect, n_cases)]
  ggplot(tmp, aes(x = stage_ratio, y = epsilon)) + geom_raster(aes(fill = mean_y)) + facet_grid(effect~n_cases) + geom_point(data = best, color = "white", size = 0.7) + scale_fill_continuous(guide = FALSE) + scale_fill_gradientn(colours = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"), guide = FALSE) + labs(title = "Optimal ϵ values for epsilon rule", y = "ϵ") + theme(axis.text = element_text(size = 7)) + scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE))
})

# --- optimal tau values ------

plot_wrapper(name = "plot_appendix_opt_thresh", fig.height = 1.5 * FIG_HEIGHT, fig.width = 0.5 * FIG_WIDTH, crop = FALSE, expr = {
  tmp = get_res_eval_average(effect_names = names(EFFECTS))
  tmp = tmp[select == "threshold rule",]
  tmp[, mean_y := BBmisc::normalize(x = mean_y, method = "range"), by = .(effect, n_cases)]
  best = tmp[, .SD[mean_y == max(mean_y)], by = .(effect, n_cases)]
  ggplot(tmp, aes(x = stage_ratio, y = thresh)) + geom_raster(aes(fill = mean_y)) + facet_grid(effect~n_cases) + geom_point(data = best, color = "white", size = 0.7) + scale_fill_continuous(guide = FALSE) + scale_fill_gradientn(colours = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"), guide = FALSE) + labs(title = "Optimal τ values for threshold rule", x = "r", y = expression(tau)) + theme(axis.text = element_text(size = 7)) + scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE))
})
