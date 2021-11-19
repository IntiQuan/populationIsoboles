# -------------------------------------------------------------------------#
# Equations ----
# -------------------------------------------------------------------------#
modelObjects <-
  list(
    odes = c(Adx1 = "-kax1*Adx1+Fabs1x1*factorUnits*INPUT1",
             Acx1 = "kax1*Adx1-Q1x1/Vcx1*Acx1+Q1x1/Vp1x1*Ap1x1-Q2x1/Vcx1*Acx1+
              Q2x1/Vp2x1*Ap2x1-CLx1*(Acx1/Vcx1)/((1+(Acx2/Vcx2)/Ki21))-
              (VMAXx1*(Acx1/Vcx1))/(KMx1*((1+(Acx2/Vcx2)/Ki21))+(Acx1/Vcx1))",
             Ap1x1 = "Q1x1/Vcx1*Acx1-Q1x1/Vp1x1*Ap1x1",
             Ap2x1 = "Q2x1/Vcx1*Acx1-Q2x1/Vp2x1*Ap2x1",
             Adx2 = "-kax2*Adx2+Fabs1x2*factorUnits*INPUT2",
             Acx2 = "kax2*Adx2-Q1x2/Vcx2*Acx2+Q1x2/Vp1x2*Ap1x2-Q2x2/Vcx2*Acx2+Q2x2/Vp2x2*Ap2x2-
              CLx2*(Acx2/Vcx2)/((1+(Acx1/Vcx1)/Ki12))-
              (VMAXx2*(Acx2/Vcx2))/(KMx2*((1+(Acx1/Vcx1)/Ki12))+(Acx2/Vcx2))",
             Ap1x2 = "Q1x2/Vcx2*Acx2-Q1x2/Vp1x2*Ap1x2",
             Ap2x2 = "Q2x2/Vcx2*Acx2-Q2x2/Vp2x2*Ap2x2",
             PL = "GR-((EMAXx1*((((Acx1/Vcx1)+yps)^hillx1/(EC50x1^hillx1+((Acx1/Vcx1)+yps)^hillx1))))+
              (EMAXx2*((((Acx2/Vcx2)+yps)^hillx2/(EC50x2^hillx2+((Acx2/Vcx2)+yps)^hillx2))))-
              Gamma*((((Acx1/Vcx1)+yps)^hillx1/(EC50x1^hillx1+((Acx1/Vcx1)+yps)^hillx1)))*
              ((((Acx2/Vcx2)+yps)^hillx2/(EC50x2^hillx2+((Acx2/Vcx2)+yps)^hillx2)))*
              ((1/(1+(EMAXx1/EMAXx2)^20))*EMAXx1+(1-(1/(1+(EMAXx1/EMAXx2)^20)))*EMAXx2))",
             INPUT1 = "0", INPUT2 = "0"),

    observables = c(OUTPUT1 = "PL"),

    trafo = c(Adx1 = "0", Acx1 = "0", Ap1x1 = "0",
              Ap2x1 = "0", Adx2 = "0", Acx2 = "0", Ap1x2 = "0",
              Ap2x2 = "0", PL = "PLbase",
              INPUT1 = "0", INPUT2 = "0",
              kax1 = "kax1",
              Fabs1x1 = "1", factorUnits = "1",
              Q1x1 = "Q1x1",
              Vcx1 = "Vcx1",
              Vp1x1 = "Vp1x1",
              Q2x1 = "0", Vp2x1 = "1",
              CLx1 = "CLx1",
              Vcx2 = "Vcx2",
              Ki21 = "1", VMAXx1 = "0", KMx1 = "1",
              kax2 = "kax2", Fabs1x2 = "1",
              Q1x2 = "Q1x2",
              Vp1x2 = "Vp1x2",
              Q2x2 = "0", Vp2x2 = "1", CLx2 = "CLx2",
              Ki12 = "1", VMAXx2 = "0", KMx2 = "1",
              GR = "GR", EMAXx1 = "EMAXx1", yps = "1e-6",
              hillx1 = "hillx1", EC50x1 = "EC50x1",
              EMAXx2 = "EMAXx2", hillx2 = "hillx2",
              EC50x2 = "EC50x2", Gamma = "Gamma",
              ton_INPUT1_1 = "(0)", toff_INPUT1_1 = "(0.0001)",
              ton_INPUT2_1 = "(0)", toff_INPUT2_1 = "(0.0001)",
              xon_INPUT1_1 = "(AMTx1/0.0001)", xon_INPUT2_1 = "(AMTx2/0.0001)",
              AMTx1 = "AMTx1",
              AMTx2 = "AMTx2",
              PLbase = "PLbase"
    ),

    events = structure(list(
      var = c("INPUT1", "INPUT1", "INPUT2", "INPUT2"),
      time = c("ton_INPUT1_1","toff_INPUT1_1", "ton_INPUT2_1", "toff_INPUT2_1"),
      value = c("xon_INPUT1_1","0", "xon_INPUT2_1", "0"),
      root = c(NA, NA, NA, NA),
      method = c("replace", "replace", "replace", "replace"
      )), class = "data.frame", row.names = c(NA, -4L))
  )


# -------------------------------------------------------------------------#
# Compile the model ----
# -------------------------------------------------------------------------#
# For more information on the model syntax and combining multiple functions to a prediction function, read ?dMod::prdfn
.tempdir = tempdir()
.currentdir <- getwd()
setwd(.tempdir)
x <- Xs(odemodel(modelObjects$odes, events = modelObjects$events, estimate = c("kax1")))
g <- Y(modelObjects$observables, x, attach.input = FALSE)
p <- P(modelObjects$trafo)
setwd(.currentdir)
# .. "return" the model prediction function(times, pars) -----
prd <- (g*x*p)
