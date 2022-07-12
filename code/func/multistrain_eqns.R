# Model 1: SEIR ----------------------------------------------------------------
multi_eqn_mod1 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR)*(1 - ((fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet3*fS*fI3) -(fbet4*fS*fI4) -(fbet5*fS*fI5) -(fbet6*fS*fI6) -(fbet7*fS*fI7) -(fbet8*fS*fI8) -(fbet9*fS*fI9) -(fbet10*fS*fI10) -
      nat_mort*fS
    dfE1 = (fbet1*fS*fI1) -sig*fE1 -nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) -sig*fE2 -nat_mort*fE2
    dfE3 = (fbet3*fS*fI3) -sig*fE3 -nat_mort*fE3
    dfE4 = (fbet4*fS*fI4) -sig*fE4 -nat_mort*fE4
    dfE5 = (fbet5*fS*fI5) -sig*fE5 -nat_mort*fE5
    dfE6 = (fbet6*fS*fI6) -sig*fE6 -nat_mort*fE6
    dfE7 = (fbet7*fS*fI7) -sig*fE7 -nat_mort*fE7
    dfE8 = (fbet8*fS*fI8) -sig*fE8 -nat_mort*fE8
    dfE9 = (fbet9*fS*fI9) -sig*fE9 -nat_mort*fE9
    dfE10 = (fbet10*fS*fI10) -sig*fE10 -nat_mort*fE10
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2
    dfI3 = sig*fE3 -gamm*(1-p_3)*fI3 -mort*p_3*fI3 -nat_mort*fI3
    dfI4 = sig*fE4 -gamm*(1-p_4)*fI4 -mort*p_4*fI4 -nat_mort*fI4
    dfI5 = sig*fE5 -gamm*(1-p_5)*fI5 -mort*p_5*fI5 -nat_mort*fI5
    dfI6 = sig*fE6 -gamm*(1-p_6)*fI6 -mort*p_6*fI6 -nat_mort*fI6
    dfI7 = sig*fE7 -gamm*(1-p_7)*fI7 -mort*p_7*fI7 -nat_mort*fI7
    dfI8 = sig*fE8 -gamm*(1-p_8)*fI8 -mort*p_8*fI8 -nat_mort*fI8
    dfI9 = sig*fE9 -gamm*(1-p_9)*fI9 -mort*p_9*fI9 -nat_mort*fI9
    dfI10 = sig*fE10 -gamm*(1-p_10)*fI10 -mort*p_10*fI10 -nat_mort*fI10
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*(1-p_3)*fI3 +gamm*(1-p_4)*fI4 +gamm*(1-p_5)*fI5 +gamm*(1-p_6)*fI6 +gamm*(1-p_7)*fI7 +gamm*(1-p_8)*fI8 +gamm*(1-p_9)*fI9 +gamm*(1-p_10)*fI10-
      nat_mort*fR
    return(list(c(dfS, 
                  dfE1, dfE2, dfE3, dfE4, dfE5, dfE6, dfE7, dfE8, dfE9, dfE10, 
                  dfI1, dfI2, dfI3, dfI4, dfI5, dfI6, dfI7, dfI8, dfI9, dfI10,
                  dfR)))})}

# Model 2: SEIR with vaccination -----------------------------------------------
eqn_mod2 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR + 
               fV + fVE1 + fVE2 + fVE3 + fVE4 + fVE5 + fVE6 + fVE7 + fVE8 + fVE9 + fVE10 + fVI1 + fVI2 + fVI3 + fVI4 + fVI5 + fVI6 + fVI7 + fVI8 + fVI9 + fVI10)*(1 - ((fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR + 
                                                                                                                                                                          fV + fVE1 + fVE2 + fVE3 + fVE4 + fVE5 + fVE6 + fVE7 + fVE8 + fVE9 + fVE10 + fVI1 + fVI2 + fVI3 + fVI4 + fVI5 + fVI6 + fVI7 + fVI8 + fVI9 + fVI10) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet3*fS*fI3) -(fbet4*fS*fI4) -(fbet5*fS*fI5) -(fbet6*fS*fI6) -(fbet7*fS*fI7) -(fbet8*fS*fI8) -(fbet9*fS*fI9) -(fbet10*fS*fI10) -
      (fbet1*fS*fVI1) -(fbet2*fS*fVI2) -(fbet3*fS*fVI3) -(fbet4*fS*fVI4) -(fbet5*fS*fVI5) -(fbet6*fS*fVI6) -(fbet7*fS*fVI7) -(fbet8*fS*fVI8) -(fbet9*fS*fVI9) -(fbet10*fS*fVI10)-
      v*fS +v_hat*fV -nat_mort*fS
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fVI1) -sig*fE1 -nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fVI2) -sig*fE2 -nat_mort*fE2
    dfE3 = (fbet3*fS*fI3) +(fbet2*fS*fVI3) -sig*fE3 -nat_mort*fE3
    dfE4 = (fbet4*fS*fI4) +(fbet2*fS*fVI4) -sig*fE4 -nat_mort*fE4
    dfE5 = (fbet5*fS*fI5) +(fbet2*fS*fVI5) -sig*fE5 -nat_mort*fE5
    dfE6 = (fbet6*fS*fI6) +(fbet2*fS*fVI6) -sig*fE6 -nat_mort*fE6
    dfE7 = (fbet7*fS*fI7) +(fbet2*fS*fVI7) -sig*fE7 -nat_mort*fE7
    dfE8 = (fbet8*fS*fI8) +(fbet2*fS*fVI8) -sig*fE8 -nat_mort*fE8
    dfE9 = (fbet9*fS*fI9) +(fbet2*fS*fVI9) -sig*fE9 -nat_mort*fE9
    dfE10 = (fbet10*fS*fI10) +(fbet2*fS*fVI10) -sig*fE10 -nat_mort*fE10
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2
    dfI3 = sig*fE3 -gamm*(1-p_3)*fI3 -mort*p_3*fI3 -nat_mort*fI3
    dfI4 = sig*fE4 -gamm*(1-p_4)*fI4 -mort*p_4*fI4 -nat_mort*fI4
    dfI5 = sig*fE5 -gamm*(1-p_5)*fI5 -mort*p_5*fI5 -nat_mort*fI5
    dfI6 = sig*fE6 -gamm*(1-p_6)*fI6 -mort*p_6*fI6 -nat_mort*fI6
    dfI7 = sig*fE7 -gamm*(1-p_7)*fI7 -mort*p_7*fI7 -nat_mort*fI7
    dfI8 = sig*fE8 -gamm*(1-p_8)*fI8 -mort*p_8*fI8 -nat_mort*fI8
    dfI9 = sig*fE9 -gamm*(1-p_9)*fI9 -mort*p_9*fI9 -nat_mort*fI9
    dfI10 = sig*fE10 -gamm*(1-p_10)*fI10 -mort*p_10*fI10 -nat_mort*fI10
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*(1-p_3)*fI3 +gamm*(1-p_4)*fI4 +gamm*(1-p_5)*fI5 +gamm*(1-p_6)*fI6 +gamm*(1-p_7)*fI7 +gamm*(1-p_8)*fI8 +gamm*(1-p_9)*fI9 +gamm*(1-p_10)*fI10 +
      gamm*fVI1 +gamm*fVI2 +gamm*fVI3 +gamm*fVI4 +gamm*fVI5 +gamm*fVI6 +gamm*fVI7 +gamm*fVI8 +gamm*fVI9 +gamm*fVI10 -nat_mort*fR
    dfV = -(fbet1*fV*fI1) -(fbet2*fV*fI2) -(fbet1*fV*fVI1) -(fbet2*fV*fVI2) +v*fS -v_hat*fV -nat_mort*fV
    dfVE1 = (fbet1*fV*fI1) +(fbet1*fV*fVI1) -sig*fVE1 -nat_mort*fVE1
    dfVE2 = (fbet2*fV*fI2) +(fbet2*fV*fVI2) -sig*fVE2 -nat_mort*fVE2
    dfVE3 = (fbet3*fV*fI3) +(fbet3*fV*fVI3) -sig*fVE3 -nat_mort*fVE3
    dfVE4 = (fbet4*fV*fI4) +(fbet4*fV*fVI4) -sig*fVE4 -nat_mort*fVE4
    dfVE5 = (fbet5*fV*fI5) +(fbet5*fV*fVI5) -sig*fVE5 -nat_mort*fVE5
    dfVE6 = (fbet6*fV*fI6) +(fbet6*fV*fVI6) -sig*fVE6 -nat_mort*fVE6
    dfVE7 = (fbet7*fV*fI7) +(fbet7*fV*fVI7) -sig*fVE7 -nat_mort*fVE7
    dfVE8 = (fbet8*fV*fI8) +(fbet8*fV*fVI8) -sig*fVE8 -nat_mort*fVE8
    dfVE9 = (fbet9*fV*fI9) +(fbet9*fV*fVI9) -sig*fVE9 -nat_mort*fVE9
    dfVE10 = (fbet10*fV*fI10) +(fbet10*fV*fVI10) -sig*fVE10 -nat_mort*fVE10
    dfVI1 = sig*fVE1 -gamm*fVI1 -nat_mort*fVI1
    dfVI2 = sig*fVE2 -gamm*fVI2 -nat_mort*fVI2
    dfVI3 = sig*fVE3 -gamm*fVI3 -nat_mort*fVI3
    dfVI4 = sig*fVE4 -gamm*fVI4 -nat_mort*fVI4
    dfVI5 = sig*fVE5 -gamm*fVI5 -nat_mort*fVI5
    dfVI6 = sig*fVE6 -gamm*fVI6 -nat_mort*fVI6
    dfVI7 = sig*fVE7 -gamm*fVI7 -nat_mort*fVI7
    dfVI8 = sig*fVE8 -gamm*fVI8 -nat_mort*fVI8
    dfVI9 = sig*fVE9 -gamm*fVI9 -nat_mort*fVI9
    dfVI10 = sig*fVE10 -gamm*fVI10 -nat_mort*fVI10
    return(list(c(dfS, 
                  dfE1, dfE2, dfE3, dfE4, dfE5, dfE6, dfE7, dfE8, dfE9, dfE10,
                  dfI1, dfI2, dfI3, dfI4, dfI5, dfI6, dfI7, dfI8, dfI9, dfI10,
                  dfR, dfV, 
                  dfVE1, dfVE2, dfVE3, dfVE4, dfVE5, dfVE6, dfVE7, dfVE8, dfVE9, dfVE10,
                  dfVI1, dfVI2, dfVI3, dfVI4, dfVI5, dfVI6, dfVI7, dfVI8, dfVI9, dfVI10)))})}





# Model 3: SEIR with vaccination and migration ---------------------------------