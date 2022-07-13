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
multi_eqn_mod2 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR + 
               fV + fVE1 + fVE2 + fVE3 + fVE4 + fVE5 + fVE6 + fVE7 + fVE8 + fVE9 + fVE10 + fVI1 + fVI2 + fVI3 + fVI4 + fVI5 + fVI6 + fVI7 + fVI8 + fVI9 + fVI10)*(1 - ((fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR + 
                                                                                                                                                                          fV + fVE1 + fVE2 + fVE3 + fVE4 + fVE5 + fVE6 + fVE7 + fVE8 + fVE9 + fVE10 + fVI1 + fVI2 + fVI3 + fVI4 + fVI5 + fVI6 + fVI7 + fVI8 + fVI9 + fVI10) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet3*fS*fI3) -(fbet4*fS*fI4) -(fbet5*fS*fI5) -(fbet6*fS*fI6) -(fbet7*fS*fI7) -(fbet8*fS*fI8) -(fbet9*fS*fI9) -(fbet10*fS*fI10) -
      (fbet1*fS*fVI1) -(fbet2*fS*fVI2) -(fbet3*fS*fVI3) -(fbet4*fS*fVI4) -(fbet5*fS*fVI5) -(fbet6*fS*fVI6) -(fbet7*fS*fVI7) -(fbet8*fS*fVI8) -(fbet9*fS*fVI9) -(fbet10*fS*fVI10) -
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
multi_eqn_mod3 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR + 
               fV + fVE1 + fVE2 + fVE3 + fVE4 + fVE5 + fVE6 + fVE7 + fVE8 + fVE9 + fVE10 + fVI1 + fVI2 + fVI3 + fVI4 + fVI5 + fVI6 + fVI7 + fVI8 + fVI9 + fVI10)*(1 - ((fS + fE1 + fE2 + fE3 + fE4 + fE5 + fE6 + fE7 + fE8 + fE9 + fE10 + fI1 + fI2 + fI3 + fI4 + fI5 + fI6 + fI7 + fI8 + fI9 + fI10 + fR + 
                                                                                                                                                                          fV + fVE1 + fVE2 + fVE3 + fVE4 + fVE5 + fVE6 + fVE7 + fVE8 + fVE9 + fVE10 + fVI1 + fVI2 + fVI3 + fVI4 + fVI5 + fVI6 + fVI7 + fVI8 + fVI9 + fVI10) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet3*fS*fI3) -(fbet4*fS*fI4) -(fbet5*fS*fI5) -(fbet6*fS*fI6) -(fbet7*fS*fI7) -(fbet8*fS*fI8) -(fbet9*fS*fI9) -(fbet10*fS*fI10) -
      (fbet1*fS*fVI1) -(fbet2*fS*fVI2) -(fbet3*fS*fVI3) -(fbet4*fS*fVI4) -(fbet5*fS*fVI5) -(fbet6*fS*fVI6) -(fbet7*fS*fVI7) -(fbet8*fS*fVI8) -(fbet9*fS*fVI9) -(fbet10*fS*fVI10) -
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
    dmS = b*(mS + mE1 + mE2 + mE3 + mE4 + mE5 + mE6 + mE7 + mE8 + mE9 + mE10 + mI1 + mI2 + mI3 + mI4 + mI5 + mI6 + mI7 + mI8 + mI9 + mI10 + mR + 
               mV + mVE1 + mVE2 + mVE3 + mVE4 + mVE5 + mVE6 + mVE7 + mVE8 + mVE9 + mVE10 + mVI1 + mVI2 + mVI3 + mVI4 + mVI5 + mVI6 + mVI7 + mVI8 + mVI9 + mVI10)*(1 - ((mS + mE1 + mE2 + mE3 + mE4 + mE5 + mE6 + mE7 + mE8 + mE9 + mE10 + mI1 + mI2 + mI3 + mI4 + mI5 + mI6 + mI7 + mI8 + mI9 + mI10 + mR + 
                                                                                                                                                                          mV + mVE1 + mVE2 + mVE3 + mVE4 + mVE5 + mVE6 + mVE7 + mVE8 + mVE9 + mVE10 + mVI1 + mVI2 + mVI3 + mVI4 + mVI5 + mVI6 + mVI7 + mVI8 + mVI9 + mVI10) / ((b / (b - nat_mort)) * (pop_size)))) -
      (mbet1*mS*mI1) -(mbet2*mS*mI2) -(mbet3*mS*mI3) -(mbet4*mS*mI4) -(mbet5*mS*mI5) -(mbet6*mS*mI6) -(mbet7*mS*mI7) -(mbet8*mS*mI8) -(mbet9*mS*mI9) -(mbet10*mS*mI10) -
      (mbet1*mS*mVI1) -(mbet2*mS*mVI2) -(mbet3*mS*mVI3) -(mbet4*mS*mVI4) -(mbet5*mS*mVI5) -(mbet6*mS*mVI6) -(mbet7*mS*mVI7) -(mbet8*mS*mVI8) -(mbet9*mS*mVI9) -(mbet10*mS*mVI10) +
      v_hat*mV -nat_mort*mS
    dmE1 = (mbet1*mS*mI1) +(mbet1*mS*mVI1) -sig*mE1 -nat_mort*mE1
    dmE2 = (mbet2*mS*mI2) +(mbet2*mS*mVI2) -sig*mE2 -nat_mort*mE2
    dmE3 = (mbet3*mS*mI3) +(mbet2*mS*mVI3) -sig*mE3 -nat_mort*mE3
    dmE4 = (mbet4*mS*mI4) +(mbet2*mS*mVI4) -sig*mE4 -nat_mort*mE4
    dmE5 = (mbet5*mS*mI5) +(mbet2*mS*mVI5) -sig*mE5 -nat_mort*mE5
    dmE6 = (mbet6*mS*mI6) +(mbet2*mS*mVI6) -sig*mE6 -nat_mort*mE6
    dmE7 = (mbet7*mS*mI7) +(mbet2*mS*mVI7) -sig*mE7 -nat_mort*mE7
    dmE8 = (mbet8*mS*mI8) +(mbet2*mS*mVI8) -sig*mE8 -nat_mort*mE8
    dmE9 = (mbet9*mS*mI9) +(mbet2*mS*mVI9) -sig*mE9 -nat_mort*mE9
    dmE10 = (mbet10*mS*mI10) +(mbet2*mS*mVI10) -sig*mE10 -nat_mort*mE10
    dmI1 = sig*mE1 -gamm*(1-p_1)*mI1 -mort*p_1*mI1 -nat_mort*mI1
    dmI2 = sig*mE2 -gamm*(1-p_2)*mI2 -mort*p_2*mI2 -nat_mort*mI2
    dmI3 = sig*mE3 -gamm*(1-p_3)*mI3 -mort*p_3*mI3 -nat_mort*mI3
    dmI4 = sig*mE4 -gamm*(1-p_4)*mI4 -mort*p_4*mI4 -nat_mort*mI4
    dmI5 = sig*mE5 -gamm*(1-p_5)*mI5 -mort*p_5*mI5 -nat_mort*mI5
    dmI6 = sig*mE6 -gamm*(1-p_6)*mI6 -mort*p_6*mI6 -nat_mort*mI6
    dmI7 = sig*mE7 -gamm*(1-p_7)*mI7 -mort*p_7*mI7 -nat_mort*mI7
    dmI8 = sig*mE8 -gamm*(1-p_8)*mI8 -mort*p_8*mI8 -nat_mort*mI8
    dmI9 = sig*mE9 -gamm*(1-p_9)*mI9 -mort*p_9*mI9 -nat_mort*mI9
    dmI10 = sig*mE10 -gamm*(1-p_10)*mI10 -mort*p_10*mI10 -nat_mort*mI10
    dmR = gamm*(1-p_1)*mI1 +gamm*(1-p_2)*mI2 +gamm*(1-p_3)*mI3 +gamm*(1-p_4)*mI4 +gamm*(1-p_5)*mI5 +gamm*(1-p_6)*mI6 +gamm*(1-p_7)*mI7 +gamm*(1-p_8)*mI8 +gamm*(1-p_9)*mI9 +gamm*(1-p_10)*mI10 +
      gamm*mVI1 +gamm*mVI2 +gamm*mVI3 +gamm*mVI4 +gamm*mVI5 +gamm*mVI6 +gamm*mVI7 +gamm*mVI8 +gamm*mVI9 +gamm*mVI10 -nat_mort*mR
    dmV = -(mbet1*mV*mI1) -(mbet2*mV*mI2) -(mbet1*mV*mVI1) -(mbet2*mV*mVI2) -v_hat*mV -nat_mort*mV
    dmVE1 = (mbet1*mV*mI1) +(mbet1*mV*mVI1) -sig*mVE1 -nat_mort*mVE1
    dmVE2 = (mbet2*mV*mI2) +(mbet2*mV*mVI2) -sig*mVE2 -nat_mort*mVE2
    dmVE3 = (mbet3*mV*mI3) +(mbet3*mV*mVI3) -sig*mVE3 -nat_mort*mVE3
    dmVE4 = (mbet4*mV*mI4) +(mbet4*mV*mVI4) -sig*mVE4 -nat_mort*mVE4
    dmVE5 = (mbet5*mV*mI5) +(mbet5*mV*mVI5) -sig*mVE5 -nat_mort*mVE5
    dmVE6 = (mbet6*mV*mI6) +(mbet6*mV*mVI6) -sig*mVE6 -nat_mort*mVE6
    dmVE7 = (mbet7*mV*mI7) +(mbet7*mV*mVI7) -sig*mVE7 -nat_mort*mVE7
    dmVE8 = (mbet8*mV*mI8) +(mbet8*mV*mVI8) -sig*mVE8 -nat_mort*mVE8
    dmVE9 = (mbet9*mV*mI9) +(mbet9*mV*mVI9) -sig*mVE9 -nat_mort*mVE9
    dmVE10 = (mbet10*mV*mI10) +(mbet10*mV*mVI10) -sig*mVE10 -nat_mort*mVE10
    dmVI1 = sig*mVE1 -gamm*mVI1 -nat_mort*mVI1
    dmVI2 = sig*mVE2 -gamm*mVI2 -nat_mort*mVI2
    dmVI3 = sig*mVE3 -gamm*mVI3 -nat_mort*mVI3
    dmVI4 = sig*mVE4 -gamm*mVI4 -nat_mort*mVI4
    dmVI5 = sig*mVE5 -gamm*mVI5 -nat_mort*mVI5
    dmVI6 = sig*mVE6 -gamm*mVI6 -nat_mort*mVI6
    dmVI7 = sig*mVE7 -gamm*mVI7 -nat_mort*mVI7
    dmVI8 = sig*mVE8 -gamm*mVI8 -nat_mort*mVI8
    dmVI9 = sig*mVE9 -gamm*mVI9 -nat_mort*mVI9
    dmVI10 = sig*mVE10 -gamm*mVI10 -nat_mort*mVI10
    return(list(c(dfS,
                  dfE1, dfE2, dfE3, dfE4, dfE5, dfE6, dfE7, dfE8, dfE9, dfE10,
                  dfI1, dfI2, dfI3, dfI4, dfI5, dfI6, dfI7, dfI8, dfI9, dfI10,
                  dfR, dfV, 
                  dfVE1, dfVE2, dfVE3, dfVE4, dfVE5, dfVE6, dfVE7, dfVE8, dfVE9, dfVE10,
                  dfVI1, dfVI2, dfVI3, dfVI4, dfVI5, dfVI6, dfVI7, dfVI8, dfVI9, dfVI10,
                  dmS,
                  dmE1, dmE2, dmE3, dmE4, dmE5, dmE6, dmE7, dmE8, dmE9, dmE10,
                  dmI1, dmI2, dmI3, dmI4, dmI5, dmI6, dmI7, dmI8, dmI9, dmI10,
                  dmR, dmV, 
                  dmVE1, dmVE2, dmVE3, dmVE4, dmVE5, dmVE6, dmVE7, dmVE8, dmVE9, dmVE10,
                  dmVI1, dmVI2, dmVI3, dmVI4, dmVI5, dmVI6, dmVI7, dmVI8, dmVI9, dmVI10)))})}
