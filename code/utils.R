
library(ggpackets)
library(tidyverse)

get_ggColorHue <- function(n) {
  # get_ggColorHue.R
  # 
  # Function will spit out color codes for a given n, just like ggplot2. 
  # Helpful when you want to manually set colors in plots.
  #
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} # function get_ggColorHue


geom_scatterbox <- ggpacket() +
  geom_jitter() +
  geom_boxplot(fill = NA, colour = "#000000", outlier.shape = NA)


cohens_d <- function(x, y, DIM=1, SIGN=TRUE, na.rm=TRUE) {
  #
  # Function will compute cohen's d effect size.
  # Generalized to work on either matrices, data frames, vectors
  #
  # INPUT
  #	x <- matrix or numeric vector
  #	y <- matrix or numeric vector
  #	DIM <- specify the dimension which samples are along
  #
  # Example usage
  #
  # x <- cbind(rnorm(100,0,1), rnorm(100,0,1), rnorm(100,0,1))
  # y <- cbind(rnorm(100,1,1), rnorm(100,2,1), rnorm(100,3,1))
  # d <- cohens_d(x, y, 1)
  #
  # written by mvlombardo - 28.08.2015
  #
  
  library(matrixStats)
  
  # if x and y are vectors, coerce them into matrices
  if (class(x)=="numeric" | class(x)=="integer") {
    x <- as.matrix(x)
  } # if
  
  if (class(y)=="numeric" | class(y)=="integer") {
    y <- as.matrix(y)
  }# if
  
  if (na.rm==TRUE){
    missingValDecision = TRUE
  } else {
    missingValDecision = FALSE
  }
  
  # n-1 for x and y
  lx <- dim(x)[DIM]-1
  ly <- dim(y)[DIM]-1
  
  # if samples are along the rows
  if (DIM==1){
    if (SIGN){
      # mean difference (numerator)
      md <- colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision)
    } else{
      # mean difference (numerator)
      md <- abs(colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision))
    }# if (SIGN)
    # pooled variance (denominator), but before any sqrt is done
    csd <- (lx * rowVars(t(x),na.rm = missingValDecision)) + (ly * rowVars(t(y), na.rm = missingValDecision))
    
    # else if samples are along the columns
  } else if (DIM==2){
    if (SIGN){
      # mean difference (numerator)
      md <- rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision)
    } else{
      # mean difference (numerator)
      md <- abs(rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision))
    }# if (SIGN)
    # pooled variance (denominator), but before any sqrt is done
    csd <- lx * rowVars(x, na.rm = missingValDecision) + ly * rowVars(y, na.rm = missingValDecision)
  }# end if
  
  # divide pooled variance by sum of n-1 for x and y and then square root it
  csd <- sqrt(csd/(lx + ly))
  # compute cohen's d
  cd  <- md/csd
}# end cohens_d <- function(x, y, DIM)



prep_dataset <- function(data_fname){

  require(here)
  require(tidyverse)
  require(foreign)
  
  # read in data
  data = read.spss(file = data_fname, to.data.frame = TRUE)
  
  # clean up
  na_rows = rowSums(is.na(data))==dim(data)[2]
  data = data %>% filter(!na_rows)
  data$Responder_Status_VABS_ABC[is.na(data$Responder_Status_VABS_ABC)] = NA
  data$Responder_Status_VABS_ABC[data$Responder_Status_VABS_ABC==999] = NA
  
  # define variables to label which dataset and site the data comes from
  data$subid = factor(as.character(data$MIRA_ID))
  
  # define variables to label which dataset and site the data comes from
  data$dataset = factor(data$Intervention_Study)
  data$site = factor(data$Intervention_Site)
  
  # define a variable to label which treatment type child was given
  data$tx = factor(data$Child_Tx_Type, levels = c("ASAP","LEAP","ESDM","NDBI","NDBI + VF","TEACCH","EIBI","PRT","STAR/IMPACT","EIBI+PRT","DTT","STAR","ABA"))
  data$tx_broad = factor(data$Child_Tx_Type_Broad_Categories, levels = c("Other","EIBI","NDBI","ESDM"))
  # custom Tx variable where ESDM and Other are the only levels
  data$tx_mvl = as.character(data$tx)
  data$tx_mvl[!is.element(data$tx,"ESDM") & !is.na(data$tx)] = "Other"
  data$tx_mvl = factor(data$tx_mvl, levels = c("Other","ESDM"))
  
  # define a variable that conveys whether the child was a Responder or Non-Responder 
  # based on whether they moved from 18 mo age equivalence pre-treatment to >18 mo
  # age equivalence post-treatment.
  data$responder18 = NA
  data$responder18[data$SM_Moved_from_below_18_mo_to_more_advanced_phase=="Yes moved from 18 to higher phase"] = "Responder"
  data$responder18[data$SM_Moved_from_below_18_mo_to_more_advanced_phase=="Did not move from 18 to higher phase"] = "Non-Responder"
  data$responder18 = factor(data$responder18, levels = c("Non-Responder","Responder"))

  data$responder15 = NA
  data$responder15[data$SM_Moved_from_below_15_mo_to_more_advanced_phase=="Yes moved from below 15 to higher phase"] = "Responder"
  data$responder15[data$SM_Moved_from_below_15_mo_to_more_advanced_phase=="Did not move from below 15 to higher phase"] = "Non-Responder"
  data$responder15 = factor(data$responder15, levels = c("Non-Responder","Responder"))
  
  data$responder12 = NA
  data$responder12[data$SM_Moved_from_below_12_mo_to_more_advanced_phase=="Yes moved from below 12 to higher phase"] = "Responder"
  data$responder12[data$SM_Moved_from_below_12_mo_to_more_advanced_phase=="Did not move from below 12 to higher phase"] = "Non-Responder"
  data$responder12 = factor(data$responder12, levels = c("Non-Responder","Responder"))

  data$responder24 = NA
  data$responder24[data$SM_Moved_from_below_24_mo_to_more_advanced_phase=="Yes moved from 24 to higher phase"] = "Responder"
  data$responder24[data$SM_Moved_from_below_24_mo_to_more_advanced_phase=="Did not move from 24 to higher phase"] = "Non-Responder"
  data$responder24 = factor(data$responder24, levels = c("Non-Responder","Responder"))
  
  # define a variable that shows whether the child was minimally verbal (<18 months) pre-treatment
  data$pre_mv18 = NA
  data$pre_mv18[data$Pre_SM_First_Words_Minimally_Verbal_Below_18=="Yes below 18 months"] = "Minimally Verbal"
  data$pre_mv18[data$Pre_SM_First_Words_Minimally_Verbal_Below_18=="Above 18 months"] = "Verbal"
  data$pre_mv18 = factor(data$pre_mv18, levels = c("Minimally Verbal","Verbal"))
  
  # define a variable that shows whether the child was minimally verbal (<18 months) post-treatment
  data$post_mv18 = NA
  data$post_mv18[data$Post_SM_First_Words_Minimally_Verbal_Below_18=="Yes below 18 months"] = "Minimally Verbal"
  data$post_mv18[data$Post_SM_First_Words_Minimally_Verbal_Below_18=="Above 18 months"] = "Verbal"
  data$post_mv18 = factor(data$post_mv18, levels = c("Minimally Verbal","Verbal"))
  
  # define a variable that conveys treatment intensity
  data$intensity = data$Child_Tx_HrsPerWeek
  
  # define a variable that conveys the child's sex
  data$sex = factor(data$Child_Sex)
  data$sex = factor(data$sex, levels = c("Female","Male"))
  
  # define a variable that conveys the child's race 
  data$race = factor(data$Child_Race)
  
  # define a variable that conveys age at treatment start in months
  data$age_at_start = data$Pre_Age_Months
  data$age_at_start_grp = factor(data$Pre_Age_Group)
  
  # define a variable that conveys age at treatment end in months
  data$age_at_end = data$Post_Age_Months
  
  # define a variable that conveys treatment duration
  # as well as variables that show whether a child's treatment was >24 months or 
  # <6 months
  data$tx_duration = data$Duration_from_Pre_to_Post_assessment
  data$tx_duration_beyond24 = data$tx_duration>24
  data$tx_duration_lessthan6 = data$tx_duration<6
  
  # define variables that convey what the child's pre-treatment developmental quotient was
  # and which measure the DQ was derived from
  data$pretx_dq = data$Pre_IQ_Overall
  data$pretx_verbal_dq = data$Pre_IQ_Verbal
  data$pretx_nonverbal_dq = data$Pre_IQ_NonVerbal
  data$pretx_dq_measure = factor(data$Pre_IQ_Overall_Measure_Used)
  
  # define variables that convey what the child's post-treatment developmental quotient was
  # and which measure the DQ was derived from
  data$posttx_dq = data$Post_IQ_Overall
  data$posttx_dq_measure = factor(data$Post_IQ_Overall_Measure_Used)
  
  # define variables that show whether the child was intellectually disabled pre- and post-treatment
  data$pre_id = factor(data$Pre_Child_IntelDis_Range)
  data$post_id = factor(data$Post_Child_IntelDis_Range)
  
  # define measure used to define pre- and post-treatment expressive language age-equivalent
  # to establish verbal status
  data$pre_verbalstatus_measure = factor(data$Pre_Standardized_Measure_used_for_Verbal_Status)
  data$post_verbalstatus_measure = factor(data$Post_Standardized_Measure_used_for_Verbal_Status)
  
  # define variable for pre-treatment imitation skills
  data$pre_imitation = data$Pre_Imitation_Percent
  
  # define variables for VABS
  # pre-treatment
  data$pre_vabs_explang_ae = data$Pre_VABS_ExpLang_AE
  data$pre_vabs_reclang_ae = data$Pre_VABS_RecLang_AE
  data$pre_vabs_personal_ae = data$Pre_VABS_Personal_AE
  data$pre_vabs_domestic_ae = data$Pre_VABS_Domestic_AE
  data$pre_vabs_community_ae = data$Pre_VABS_Community_AE
  data$pre_vabs_interpersonal_ae = data$Pre_VABS_IntRel_AE
  data$pre_vabs_playleisure_ae = data$Pre_VABS_PlayandLeisure_AE
  data$pre_vabs_coping_ae = data$Pre_VABS_Coping_AE
  data$pre_vabs_grossmotor_ae = data$Pre_VABS_GrossMotor_AE
  data$pre_vabs_finemotor_ae = data$Pre_VABS_FineMotor_AE
  data$pre_vabs_abc = data$Pre_VABS_ABC
  data$pre_vabs_comm_ss = data$Pre_VABS_Communication_SS
  data$pre_vabs_dls_ss = data$Pre_VABS_DailyLiving_SS
  data$pre_vabs_soc_ss = data$Pre_VABS_Socialization_SS
  data$pre_vabs_motor_ss = data$Pre_VABS_Motor_SS
  # post-treatment
  data$post_vabs_explang_ae = data$Post_VABS_ExpLang_AE
  data$post_vabs_reclang_ae = data$Post_VABS_RecLang_AE
  data$post_vabs_personal_ae = data$Post_VABS_Personal_AE
  data$post_vabs_domestic_ae = data$Post_VABS_Domestic_AE
  data$post_vabs_community_ae = data$Post_VABS_Community_AE
  data$post_vabs_interpersonal_ae = data$Post_VABS_IntRel_AE
  data$post_vabs_playleisure_ae = data$Post_VABS_PlayandLeisure_AE
  data$post_vabs_coping_ae = data$Post_VABS_Coping_AE
  data$post_vabs_grossmotor_ae = data$Post_VABS_GrossMotor_AE
  data$post_vabs_finemotor_ae = data$Post_VABS_FineMotor_AE
  data$post_vabs_abc = data$Post_VABS_ABC
  data$post_vabs_comm_ss = data$Post_VABS_Comm_SS
  data$post_vabs_dls_ss = data$Post_VABS_DailyLiving_SS
  data$post_vabs_soc_ss = data$Post_VABS_Socialization_SS
  data$post_vabs_motor_ss = data$Post_VABS_Motor_SS
  
  # define variables for ADOS
  data$pre_ados_css = data$Pre_ADOS_CSS_Total
  data$pre_ados_css_sa = data$Pre_ADOS_CSS_SA
  data$pre_ados_css_rrb = data$Pre_ADOS_CSS_RRB
  
  # define variables for MSEL
  # pre-treatment
  data$pre_msel_vr_t = data$Pre_MSEL_VisRec_TS
  data$pre_msel_fm_t = data$Pre_MSEL_FineMot_TS
  data$pre_msel_el_t = data$Pre_MSEL_ExpLang_TS
  data$pre_msel_rl_t = data$Pre_MSEL_RecLang_TS
  data$pre_msel_vr_ae = data$Pre_MSEL_VisRec_AE
  data$pre_msel_fm_ae = data$Pre_MSEL_FineMot_AE
  data$pre_msel_el_ae = data$Pre_MSEL_ExpLang_AE
  data$pre_msel_rl_ae = data$Pre_MSEL_RecLang_AE
  data$pre_msel_elc = data$Pre_MSEL_ELC_SS
  # post-treatment
  data$post_msel_vr_t = data$Post_MSEL_VisRec_TS
  data$post_msel_fm_t = data$Post_MSEL_FineMot_TS
  data$post_msel_el_t = data$Post_MSEL_ExpLang_TS
  data$post_msel_rl_t = data$Post_MSEL_RecLang_TS
  data$post_msel_vr_ae = data$Post_MSEL_VisRec_AE
  data$post_msel_fm_ae = data$Post_MSEL_FineMot_AE
  data$post_msel_el_ae = data$Post_MSEL_ExpLang_AE
  data$post_msel_rl_ae = data$Post_MSEL_RecLang_AE
  data$post_msel_elc = data$Post_MSEL_ELC
  
  
  vars2use = c("subid","dataset","site","sex","race","tx","tx_mvl","tx_broad","intensity",
               "responder18","responder15","responder12","responder24","pre_mv18","post_mv18",
               "tx_duration","tx_duration_beyond24","tx_duration_lessthan6",
               "age_at_start","age_at_start_grp","age_at_end",
               "pretx_dq","pretx_verbal_dq","pretx_nonverbal_dq","pretx_dq_measure",
               "posttx_dq","posttx_dq_measure",
               "pre_id","post_id",
               "pre_verbalstatus_measure","post_verbalstatus_measure",
               "pre_imitation",
               "pre_vabs_explang_ae","pre_vabs_reclang_ae",
               "pre_vabs_personal_ae","pre_vabs_domestic_ae",
               "pre_vabs_community_ae","pre_vabs_interpersonal_ae",
               "pre_vabs_playleisure_ae","pre_vabs_coping_ae",
               "pre_vabs_grossmotor_ae","pre_vabs_finemotor_ae",
               "pre_vabs_abc",
               "pre_vabs_comm_ss","pre_vabs_dls_ss","pre_vabs_soc_ss","pre_vabs_motor_ss",
               "post_vabs_explang_ae","post_vabs_reclang_ae",
               "post_vabs_personal_ae","post_vabs_domestic_ae",
               "post_vabs_community_ae","post_vabs_interpersonal_ae",
               "post_vabs_playleisure_ae","post_vabs_coping_ae",
               "post_vabs_grossmotor_ae","post_vabs_finemotor_ae",
               "post_vabs_abc",
               "post_vabs_comm_ss","post_vabs_dls_ss","post_vabs_soc_ss","post_vabs_motor_ss",
               "pre_ados_css","pre_ados_css_sa","pre_ados_css_rrb",
               "pre_msel_vr_t","pre_msel_fm_t","pre_msel_el_t","pre_msel_rl_t",
               "pre_msel_vr_ae","pre_msel_fm_ae","pre_msel_el_ae","pre_msel_rl_ae","pre_msel_elc",
               "post_msel_vr_t","post_msel_fm_t","post_msel_el_t","post_msel_rl_t",
               "post_msel_vr_ae","post_msel_fm_ae","post_msel_el_ae","post_msel_rl_ae","post_msel_elc")
  
  # make a mask to remove anyone that has >24 mo or <6 mo of treatment duration
  mask = !data$tx_duration_beyond24 | data$tx_duration_lessthan6 
  
  # apply mask and the select key columns of interest
  tidy_data = data %>% filter(mask) %>% select(vars2use)
  
  # make rownames the subid
  rownames(tidy_data) = tidy_data$subid
  
  return(tidy_data)
  
} # function prep_dataset


scale_data <- function(data, vars2scale){
  data[,vars2scale] = scale(data[,vars2scale])
  return(data)
} # function scale_data


logistic_regression <- function(data, 
                                formula, 
                                vars2scale, responder_var, 
                                pcavars,
                                do_pca = FALSE, 
                                pcs2use = 5,
                                recon_pcs = NULL,
                                print_output = TRUE){
  
  require(lme4)

  # grab subset of data for just responders or non-responders
  # mask = !is.na(data$responder)
  # data2use = data %>% filter(mask)

  mask2 = rowSums(is.na(data[,c(responder_var, vars2scale)]))==0
  data2use = data %>% filter(mask2)
  rownames(data2use) = data2use$subid
  
  # principal component analysis on the continuous predictors
  if (do_pca){
    
    # pcs2use = 5 # how many PCs to use as predictors
    pc_names = character(length=pcs2use)
    for (ipc in 1:pcs2use){
      pc_names[ipc] = sprintf("PC%d",ipc)
    }
    
    # pca_res = prcomp(data2use[,pcavars])
    pca_res = prcomp(data2use[,pcavars], scale=TRUE)
    pca_scores = data.frame(pca_res$x)
    pca_scores$subid = rownames(pca_scores)
    
    data2use = merge(data2use, pca_scores, by = "subid")
  } # if (do_pca)
  
  # scale data
  scaled_data2use = scale_data(data2use, vars2scale)
  
  # fit model
  mod2use = glmer(formula = formula, 
                  data = scaled_data2use, 
                  family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
  
  if(print_output){
    print(mod2use, corr=FALSE)
  }
  
  # compute standard errors 
  se = sqrt(diag(vcov(mod2use)))
  
  # compute 95% CIs
  tab = cbind(Est = fixef(mod2use), 
              LL = fixef(mod2use) - 1.96 * se, 
              UL = fixef(mod2use) + 1.96 *se)
  
  # convert to odds ratios
  res = data.frame(cbind(tab,exp(tab)))
  colnames(res) = c("beta","ci.low95","ci.high85","OR","ci.OR.low95","ci.OR.high95")
  a = summary(mod2use)
  
  # grab p-values and FDR
  res$pval = as.numeric(a$coefficients[,"Pr(>|z|)"])
  res$fdr = p.adjust(res$pval, method = "fdr")
  if(print_output){
    print(res)
  }
  
  # export results and model
  result = list()
  result[["data"]] = data2use
  result[["scaled_data"]] = scaled_data2use
  result[["model"]] = mod2use
  result[["table"]] = res
  if (do_pca){
    result[["pca_res"]] = pca_res
  }
  
  if (!is.null(recon_pcs)){
    mu = colMeans(scaled_data2use[,pcavars])
    tmp_pca_res = prcomp(scaled_data2use[,pcavars])
    Xhat = tmp_pca_res$x[,recon_pcs] %*% t(tmp_pca_res$rotation[,recon_pcs])
    Xhat = scale(Xhat, center = -mu, scale = FALSE)
    recon_data = data.frame(Xhat)
    recon_data$subid = factor(data2use$subid)
    rownames(recon_data) = recon_data$subid
    recon_data = merge(recon_data, data2use[,c("subid",responder_var)], by="subid")
    recon_data[,responder_var] = factor(recon_data[,responder_var])
    result[["pca_recon_data"]] = recon_data
  }

  return(result)
} # function logistic_regression
