# Model 1: SEIR ----------------------------------------------------------------
eqn_mod1 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -nat_mort*fS
    dfE1 = (fbet1*fS*fI1) -sig*fE1 -nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) -sig*fE2 -nat_mort*fE2
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 -nat_mort*fR
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR)))})}

# Model 1 modification: SEIR from Ian's paper ----------------------------------
eqn_mod1_ian <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -nat_mort*fS
    dfE1 = (fbet1*fS*fI1) -sig*fE1 -nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) -sig*fE2 -nat_mort*fE2
    dfI1 = sig*fE1 -gamm*fI1 -p_1*100*fI1 -nat_mort*fI1
    dfI2 = sig*fE2 -gamm*fI2 -p_2*100*fI2 -nat_mort*fI2
    dfR = gamm*fI1 +gamm*fI2 -nat_mort*fR
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR)))})}

# Model 2: SEIR with vaccination -----------------------------------------------
eqn_mod2 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet1*fS*fVI1) -(fbet2*fS*fVI2) -v*fS +v_hat*fV -nat_mort*fS
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fVI1) -sig*fE1 -nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fVI2) -sig*fE2 -nat_mort*fE2
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*fVI1 +gamm*fVI2 -nat_mort*fR
    dfV = -(fbet1*fV*fI1) -(fbet2*fV*fI2) -(fbet1*fV*fVI1) -(fbet2*fV*fVI2) +v*fS -v_hat*fV -nat_mort*fV
    dfVE1 = (fbet1*fV*fI1) +(fbet1*fV*fVI1) -sig*fVE1 -nat_mort*fVE1
    dfVE2 = (fbet2*fV*fI2) +(fbet2*fV*fVI2) -sig*fVE2 -nat_mort*fVE2
    dfVI1 = sig*fVE1 -gamm*fVI1 -nat_mort*fVI1
    dfVI2 = sig*fVE2 -gamm*fVI2 -nat_mort*fVI2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR, dfV, dfVE1, dfVE2, dfVI1, dfVI2)))})}

# Model 3: SEIR with vaccination and migration, no slaughter ----------------------
eqn_mod3 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet1*fS*fVI1) -(fbet2*fS*fVI2) -v*fS +v_hat*fV -nat_mort*fS +m_mf*mS -m_fm*fS
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fVI1) -sig*fE1 -nat_mort*fE1 +m_mf*mE1 -m_fm*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fVI2) -sig*fE2 -nat_mort*fE2 +m_mf*mE2 -m_fm*fE2
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1 +m_mf*mI1 -m_fm*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2 +m_mf*mI2 -m_fm*fI2
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*fVI1 +gamm*fVI2 -nat_mort*fR +m_mf*mR -m_fm*fR
    dfV = -(fbet1*fV*fI1) -(fbet2*fV*fI2) -(fbet1*fV*fVI1) -(fbet2*fV*fVI2) +v*fS -v_hat*fV -nat_mort*fV +m_mf*mV -m_fm_vax*fV
    dfVE1 = (fbet1*fV*fI1) +(fbet1*fV*fVI1) -sig*fVE1 -nat_mort*fVE1 +m_mf*mVE1 -m_fm_vax*fVE1
    dfVE2 = (fbet2*fV*fI2) +(fbet2*fV*fVI2) -sig*fVE2 -nat_mort*fVE2 +m_mf*mVE2 -m_fm_vax*fVE2
    dfVI1 = sig*fVE1 -gamm*fVI1 -nat_mort*fVI1 +m_mf*mVI1 -m_fm_vax*fVI1
    dfVI2 = sig*fVE2 -gamm*fVI2 -nat_mort*fVI2 +m_mf*mVI2 -m_fm_vax*fVI2
    dmS = b*(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)*(1 - ((mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (mbet1*mS*mI1) -(mbet2*mS*mI2) -(mbet1*mS*mVI1) -(mbet2*mS*mVI2) +v_hat*mV -nat_mort*mS -m_mf*mS +m_fm*fS
    dmE1 = (mbet1*mS*mI1) +(mbet1*mS*mVI1) -sig*mE1 -nat_mort*mE1 -m_mf*mE1 +m_fm*fE1
    dmE2 = (mbet2*mS*mI2) +(mbet2*mS*mVI2) -sig*mE2 -nat_mort*mE2 -m_mf*mE2 +m_fm*fE2
    dmI1 = sig*mE1 -gamm*(1-p_1)*mI1 -mort*p_1*mI1 -nat_mort*mI1 -m_mf*mI1 +m_fm*fI1
    dmI2 = sig*mE2 -gamm*(1-p_2)*mI2 -mort*p_2*mI2 -nat_mort*mI2 -m_mf*mI2 +m_fm*fI2
    dmR = gamm*(1-p_1)*mI1 +gamm*(1-p_2)*mI2 +gamm*mVI1 +gamm*mVI2 -nat_mort*mR -m_mf*mR +m_fm*fR
    dmV = -(mbet1*mV*mI1) -(mbet2*mV*mI2) -(mbet1*mV*mVI1) -(mbet2*mV*mVI2) -v_hat*mV -nat_mort*mV -m_mf*mV +m_fm_vax*fV
    dmVE1 = (mbet1*mV*mI1) +(mbet1*mV*mVI1) -sig*mVE1 -nat_mort*mVE1 -m_mf*mVE1 +m_fm_vax*fVE1
    dmVE2 = (mbet2*mV*mI2) +(mbet2*mV*mVI2) -sig*mVE2 -nat_mort*mVE2 -m_mf*mVE2 +m_fm_vax*fVE2
    dmVI1 = sig*mVE1 -gamm*mVI1 -nat_mort*mVI1 -m_mf*mVI1 +m_fm_vax*fVI1
    dmVI2 = sig*mVE2 -gamm*mVI2 -nat_mort*mVI2 -m_mf*mVI2 +m_fm_vax*fVI2
    dM_FM = -M_FM +m_fm*fS +m_fm*fE1 +m_fm*fE2 +m_fm*fI1 +m_fm*fI2 +m_fm*fR +m_fm_vax*fV +m_fm_vax*fVE1 +m_fm_vax*fVE2 +m_fm_vax*fVI1 +m_fm_vax*fVI2
    dM_MF = -M_MF +m_mf*mS +m_mf*mE1 +m_mf*mE2 +m_mf*mI1 +m_mf*mI2 +m_mf*mR +m_mf*mV +m_mf*mVE1 +m_mf*mVE2 +m_mf*mVI1 +m_mf*mVI2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR, dfV, dfVE1, dfVE2, dfVI1, dfVI2,
                  dmS, dmE1, dmE2, dmI1, dmI2, dmR, dmV, dmVE1, dmVE2, dmVI1, dmVI2,
                  dM_FM, dM_MF)))})}

# Model 4: SEIR with vaccination and migration and slaughter -------------------
eqn_mod4 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet1*fS*fVI1) -(fbet2*fS*fVI2) -v*fS +v_hat*fV -nat_mort*fS +m_mf*mS*(1-p_s) -m_fm*fS
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fVI1) -sig*fE1 -nat_mort*fE1 +m_mf*mE1*(1-p_s) -m_fm*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fVI2) -sig*fE2 -nat_mort*fE2 +m_mf*mE2*(1-p_s) -m_fm*fE2
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1 +m_mf*mI1*(1-p_s) -m_fm*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2 +m_mf*mI2*(1-p_s) -m_fm*fI2
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*fVI1 +gamm*fVI2 -nat_mort*fR +m_mf*mR*(1-p_s) -m_fm*fR
    dfV = -(fbet1*fV*fI1) -(fbet2*fV*fI2) -(fbet1*fV*fVI1) -(fbet2*fV*fVI2) +v*fS -v_hat*fV -nat_mort*fV +m_mf*mV*(1-p_s) -m_fm_vax*fV
    dfVE1 = (fbet1*fV*fI1) +(fbet1*fV*fVI1) -sig*fVE1 -nat_mort*fVE1 +m_mf*mVE1*(1-p_s) -m_fm_vax*fVE1
    dfVE2 = (fbet2*fV*fI2) +(fbet2*fV*fVI2) -sig*fVE2 -nat_mort*fVE2 +m_mf*mVE2*(1-p_s) -m_fm_vax*fVE2
    dfVI1 = sig*fVE1 -gamm*fVI1 -nat_mort*fVI1 +m_mf*mVI1*(1-p_s) -m_fm_vax*fVI1
    dfVI2 = sig*fVE2 -gamm*fVI2 -nat_mort*fVI2 +m_mf*mVI2*(1-p_s) -m_fm_vax*fVI2
    dmS = b*(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)*(1 - ((mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (mbet1*mS*mI1) -(mbet2*mS*mI2) -(mbet1*mS*mVI1) -(mbet2*mS*mVI2) +v_hat*mV -nat_mort*mS -m_mf*mS +m_fm*fS
    dmE1 = (mbet1*mS*mI1) +(mbet1*mS*mVI1) -sig*mE1 -nat_mort*mE1 -m_mf*mE1 +m_fm*fE1
    dmE2 = (mbet2*mS*mI2) +(mbet2*mS*mVI2) -sig*mE2 -nat_mort*mE2 -m_mf*mE2 +m_fm*fE2
    dmI1 = sig*mE1 -gamm*(1-p_1)*mI1 -mort*p_1*mI1 -nat_mort*mI1 -m_mf*mI1 +m_fm*fI1
    dmI2 = sig*mE2 -gamm*(1-p_2)*mI2 -mort*p_2*mI2 -nat_mort*mI2 -m_mf*mI2 +m_fm*fI2
    dmR = gamm*(1-p_1)*mI1 +gamm*(1-p_2)*mI2 +gamm*mVI1 +gamm*mVI2 -nat_mort*mR -m_mf*mR +m_fm*fR
    dmV = -(mbet1*mV*mI1) -(mbet2*mV*mI2) -(mbet1*mV*mVI1) -(mbet2*mV*mVI2) -v_hat*mV -nat_mort*mV -m_mf*mV +m_fm_vax*fV
    dmVE1 = (mbet1*mV*mI1) +(mbet1*mV*mVI1) -sig*mVE1 -nat_mort*mVE1 -m_mf*mVE1 +m_fm_vax*fVE1
    dmVE2 = (mbet2*mV*mI2) +(mbet2*mV*mVI2) -sig*mVE2 -nat_mort*mVE2 -m_mf*mVE2 +m_fm_vax*fVE2
    dmVI1 = sig*mVE1 -gamm*mVI1 -nat_mort*mVI1 -m_mf*mVI1 +m_fm_vax*fVI1
    dmVI2 = sig*mVE2 -gamm*mVI2 -nat_mort*mVI2 -m_mf*mVI2 +m_fm_vax*fVI2
    dM_FM = -M_FM +(1-p_s)*(m_fm*fS +m_fm*fE1 +m_fm*fE2 +m_fm*fI1 +m_fm*fI2 +m_fm*fR +m_fm_vax*fV +m_fm_vax*fVE1 +m_fm_vax*fVE2 +m_fm_vax*fVI1 +m_fm_vax*fVI2)
    dM_MF = -M_MF +m_mf*mS +m_mf*mE1 +m_mf*mE2 +m_mf*mI1 +m_mf*mI2 +m_mf*mR +m_mf*mV +m_mf*mVE1 +m_mf*mVE2 +m_mf*mVI1 +m_mf*mVI2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR, dfV, dfVE1, dfVE2, dfVI1, dfVI2,
                  dmS, dmE1, dmE2, dmI1, dmI2, dmR, dmV, dmVE1, dmVE2, dmVI1, dmVI2,
                  dM_FM, dM_MF)))})}

# Model 5: SEIR with vaccination and migration, no migration infectious --------
eqn_mod5 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR + fV + fVE1 + fVE2 + fVI1 + fVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -(fbet1*fS*fVI1) -(fbet2*fS*fVI2) -v*fS +v_hat*fV -nat_mort*fS +m_mf*mS*(1-p_s) -m_fm*fS
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fVI1) -sig*fE1 -nat_mort*fE1 +m_mf*mE1*(1-p_s) -m_fm*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fVI2) -sig*fE2 -nat_mort*fE2 +m_mf*mE2*(1-p_s) -m_fm*fE2
    dfI1 = sig*fE1 -gamm*(1-p_1)*fI1 -mort*p_1*fI1 -nat_mort*fI1
    dfI2 = sig*fE2 -gamm*(1-p_2)*fI2 -mort*p_2*fI2 -nat_mort*fI2
    dfR = gamm*(1-p_1)*fI1 +gamm*(1-p_2)*fI2 +gamm*fVI1 +gamm*fVI2 -nat_mort*fR +m_mf*mR*(1-p_s) -m_fm*fR
    dfV = -(fbet1*fV*fI1) -(fbet2*fV*fI2) -(fbet1*fV*fVI1) -(fbet2*fV*fVI2) +v*fS -v_hat*fV -nat_mort*fV +m_mf*mV*(1-p_s) -m_fm_vax*fV
    dfVE1 = (fbet1*fV*fI1) +(fbet1*fV*fVI1) -sig*fVE1 -nat_mort*fVE1 +m_mf*mVE1*(1-p_s) -m_fm_vax*fVE1
    dfVE2 = (fbet2*fV*fI2) +(fbet2*fV*fVI2) -sig*fVE2 -nat_mort*fVE2 +m_mf*mVE2*(1-p_s) -m_fm_vax*fVE2
    dfVI1 = sig*fVE1 -gamm*fVI1 -nat_mort*fVI1 +m_mf*mVI1*(1-p_s) -m_fm_vax*fVI1
    dfVI2 = sig*fVE2 -gamm*fVI2 -nat_mort*fVI2 +m_mf*mVI2*(1-p_s) -m_fm_vax*fVI2
    dmS = b*(mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2)*(1 - ((mS + mE1 + mE2 + mI1 + mI2 + mR + mV + mVE1 + mVE2 + mVI1 + mVI2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (mbet1*mS*mI1) -(mbet2*mS*mI2) -(mbet1*mS*mVI1) -(mbet2*mS*mVI2) +v_hat*mV -nat_mort*mS -m_mf*mS +m_fm*fS
    dmE1 = (mbet1*mS*mI1) +(mbet1*mS*mVI1) -sig*mE1 -nat_mort*mE1 -m_mf*mE1 +m_fm*fE1
    dmE2 = (mbet2*mS*mI2) +(mbet2*mS*mVI2) -sig*mE2 -nat_mort*mE2 -m_mf*mE2 +m_fm*fE2
    dmI1 = sig*mE1 -gamm*(1-p_1)*mI1 -mort*p_1*mI1 -nat_mort*mI1
    dmI2 = sig*mE2 -gamm*(1-p_2)*mI2 -mort*p_2*mI2 -nat_mort*mI2
    dmR = gamm*(1-p_1)*mI1 +gamm*(1-p_2)*mI2 +gamm*mVI1 +gamm*mVI2 -nat_mort*mR -m_mf*mR +m_fm*fR
    dmV = -(mbet1*mV*mI1) -(mbet2*mV*mI2) -(mbet1*mV*mVI1) -(mbet2*mV*mVI2) -v_hat*mV -nat_mort*mV -m_mf*mV +m_fm_vax*fV
    dmVE1 = (mbet1*mV*mI1) +(mbet1*mV*mVI1) -sig*mVE1 -nat_mort*mVE1 -m_mf*mVE1 +m_fm_vax*fVE1
    dmVE2 = (mbet2*mV*mI2) +(mbet2*mV*mVI2) -sig*mVE2 -nat_mort*mVE2 -m_mf*mVE2 +m_fm_vax*fVE2
    dmVI1 = sig*mVE1 -gamm*mVI1 -nat_mort*mVI1 -m_mf*mVI1 +m_fm_vax*fVI1
    dmVI2 = sig*mVE2 -gamm*mVI2 -nat_mort*mVI2 -m_mf*mVI2 +m_fm_vax*fVI2
    dM_FM = -M_FM +(1-p_s)*(m_fm*fS +m_fm*fE1 +m_fm*fE2 +m_fm*fI1 +m_fm*fI2 +m_fm*fR +m_fm_vax*fV +m_fm_vax*fVE1 +m_fm_vax*fVE2 +m_fm_vax*fVI1 +m_fm_vax*fVI2)
    dM_MF = -M_MF +m_mf*mS +m_mf*mE1 +m_mf*mE2 +m_mf*mI1 +m_mf*mI2 +m_mf*mR +m_mf*mV +m_mf*mVE1 +m_mf*mVE2 +m_mf*mVI1 +m_mf*mVI2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR, dfV, dfVE1, dfVE2, dfVI1, dfVI2,
                  dmS, dmE1, dmE2, dmI1, dmI2, dmR, dmV, dmVE1, dmVE2, dmVI1, dmVI2,
                  dM_FM, dM_MF)))})}

