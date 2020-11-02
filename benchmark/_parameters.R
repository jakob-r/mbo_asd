library(batchtools)
library(smoof)
library(asd.mod)
library(data.table)
library(mlr3misc)

TESTMODE = Sys.info()[["nodename"]] %in% c("JayBook2", "HERKULES5ARCH")

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
  REPLS = 20
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
