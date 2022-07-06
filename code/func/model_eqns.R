# Model 1: SEIR model ----------------------------------------------------------
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

# Model 2: SEIR model with vaccination -----------------------------------------
# Model 3: SEIR model with vaccination and space -------------------------------
eqn_mod3 <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    # Farm patch
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR1 + fR2 + fV + fV_E1 + fV_I1 + fV_E2 + fV_I2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR1 + fR2 + fV + fV_E1 + fV_I1 + fV_E2 + fV_I2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -
      (fbet1*fS*fV_I1) -(fbet2*fS*fV_I2) -
      v*fS +v_hat*fV -
      m_fm*fS +m_mf*mS -
      nat_mort*fS +
      theta*fR1 + theta*fR2
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fV_I1) -
      sig*fE1 -
      m_fm*fE1 +m_mf*mE1 -
      nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fV_I2) -
      sig*fE2 -
      m_fm*fE2 +m_mf*mE2 -
      nat_mort*fE2
    dfI1 = sig*fE1 -
      gamm*(1-p_1)*fI1 -mort*p_1*fI1 -
      m_fm*fI1 +m_mf*mI1 -
      nat_mort*fI1
    dfI2 = sig*fE2 -
      gamm*(1-p_2)*fI2 -mort*p_2*fI2 -
      m_fm*fI2 + m_mf*mI2 -
      nat_mort*fI2
    dfR1 = gamm*(1-p_1)*fI1 -
      theta*fR1 -
      m_fm*fR1 +m_mf*mR1 -
      nat_mort*fR1
    dfR2 = gamm*(1-p_2)*fI2 -
      theta*fR2 -
      m_fm*fR2 +m_mf*mR2 -
      nat_mort*fR2
    dfV = v*(fS) -v_hat*fV -
      (fbet1*fV*fI1) -(fbet2*fV*fI2) -
      (fbet1*fV*fV_I1) -(fbet2*fV*fV_I2) +
      gamm*fV_I1 +gamm*fV_I2 -
      m_fm_vax*fV + m_mf*mV -
      nat_mort*fV
    dfV_E1 = (fbet1*fV*fI1) +(fbet1*fV*fV_I1) -
      sig*fV_E1 -
      m_fm_vax*fV_E1 +m_mf*mV_E1 -
      nat_mort*fV_E1
    dfV_I1 = sig*fV_E1 -
      gamm*fV_I1 -
      m_fm_vax*fV_I1 +m_mf*mV_I1 -
      nat_mort*fV_I1
    dfV_E2 = (fbet2*fV*fI2) +(fbet2*fV*fV_I2) -
      sig*fV_E2 -
      m_fm_vax*fV_E2 +m_mf*mV_E2 -
      nat_mort*fV_E2
    dfV_I2 = sig*fV_E2 -
      gamm*fV_I2 -
      m_fm_vax*fV_I2 +m_mf*mV_I2 -
      nat_mort*fV_I2
    # Market patch
    dmS = b*(mS + mE1 + mE2 + mI1 + mI2 + mR1 + mR2 + mV + mV_E1 + mV_I1 + mV_E2 + mV_I2)*(1 - ((mS + mE1 + mE2 + mI1 + mI2 + mR1 + mR2 + mV + mV_E1 + mV_I1 + mV_E2 + mV_I2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (mbet1*mS*mI1) -(mbet2*mS*mI2) -
      (mbet1*mS*mV_I1) -(mbet2*mS*mV_I2) -
      m_mf*mS +m_fm*fS -
      nat_mort*mS +
      theta*mR1 + theta*mR2
    dmE1 = (mbet1*mS*mI1) +(mbet1*mS*mV_I1) -
      sig*mE1 -
      m_mf*mE1 +m_fm*fE1 -
      nat_mort*mE1
    dmE2 = (mbet2*mS*mI2) +(mbet2*mS*mV_I2)  -
      sig*mE2 -
      m_mf*mE2 +m_fm*fE2 -
      nat_mort*mE2
    dmI1 = sig*mE1 -
      gamm*(1-p_1)*mI1 -mort*p_1*mI1 -
      m_mf*mI1 +m_fm*fI1 -
      nat_mort*mI1
    dmI2 = sig*mE2 -
      gamm*(1-p_2)*mI2 -mort*p_2*mI2 -
      m_mf*mI2 +m_fm*fI2 -
      nat_mort*mI2
    dmR1 = gamm*(1-p_1)*mI1 -
      theta*mR1 -
      m_mf*mR1 +m_fm*fR1 -
      nat_mort*mR1
    dmR2 = gamm*(1-p_2)*mI2 -
      theta*mR2 -
      m_mf*mR2 +m_fm*fR2 -
      nat_mort*mR2
    dmV = -v_hat*mV -
      (mbet1*mV*mI1) -(mbet2*mV*mI2) -
      (mbet1*mV*mV_I1) -(mbet2*mV*mV_I2) +
      gamm*mV_I1 + gamm*mV_I2 +
      m_fm*fV -m_mf*mV -
      nat_mort*mV
    dmV_E1 = (mbet1*mV*mI1) +(mbet1*mV*mV_I1) -
      sig*mV_E1 +
      m_fm*fV_E1 -m_mf*mV_E1 -
      nat_mort*mV_E1
    dmV_I1 = sig*mV_E1 -
      gamm*mV_I1 +
      m_fm*fV_I1 -m_mf*mV_I1 -
      nat_mort*mV_I1
    dmV_E2 = (mbet2*mV*mI2) +(mbet2*mV*mV_I2) -
      sig*mV_E2 +
      m_fm*fV_E2 -m_mf*mV_E2 -
      nat_mort*mV_E2
    dmV_I2 = sig*mV_E2 -
      gamm*mV_I2 +
      m_fm*fV_I2 -m_mf*mV_I2 -
      nat_mort*mV_I2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR1, dfR2,
                  dfV, dfV_E1, dfV_I1, dfV_E2, dfV_I2, 
                  dmS, dmE1, dmE2, dmI1, dmI2, dmR1, dmR2, 
                  dmV, dmV_E1, dmV_I1, dmV_E2, dmV_I2)))})}


