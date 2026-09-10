#' Simulate a time based on a provided survival function.
#'
#' @param time vector of times corresponding to the survival function
#' @param surv vector of survival probabilities at each time.
#' @return A data frame with "time" and "status" elements, where status=1 indicates the event occurred.
#' @export
gettime<-function(time,surv){

  unif1=runif(n=1)
  cdf=1-surv

  eventtime=approx(x=cdf, y=time, xout=unif1, rule=2)$y

  eventstatus=ifelse(eventtime==max(time),0,1)

  return(data.frame(time = eventtime, status = eventstatus))
}

###########################################################################################
#' Simulate a time of cancer death for a specific cancer
#' Simulate a time of cancer death for a specific cancer in a specific stage based on provided survival distributions.
#' Uses the inverse CDF method.
#' @param the_stage stage at clinical diagnosis
#' @param the_cancer_site cancer site (one of "Anus","Breast","Bladder","Colorectal","Esophagus","Headandneck","Gastric","Liver" ,"Lung","Pancreas", "Prostate", "Renal", "Ovary", "Uterine")
#' @param the_sex Male or Female
#' @param the_model_type Weibull or  Loglogistic
#' @param cancer_survival_dist: a data frame with columns: surv, time, site, stage, sex, model_type
#' @return A data frame with "time" and "status" elements, where status=1 indicates the event occurred.
#' @export
# OUTPUTS: A numeric value representing the simulated time of death due to a specific cancer
#################################################################################################
sim_cancer_death <- function(the_stage, the_cancer_site, the_sex,the_model_type, cancer_survival_dist){

  # Filter the survival distribution based on the type and stage
  survival_dist_indiv = filter(cancer_survival_dist,cancer_site == paste(the_cancer_site), stage==paste(the_stage),
                               sex==paste(the_sex),model_type==paste(the_model_type))


  # Get the survival time based on the distribution
  death_info = gettime(time = survival_dist_indiv$time, surv = survival_dist_indiv$surv)
  death_time = death_info$time

  return(death_time)
}

#' Simulate a time of cancer death for a specific cancer in a specific stage based on provided parametric survival distributions.
#'
#' @param the_stage stage at clinical diagnosis
#' @param the_cancer_site cancer site (one of "Anus","Breast","Bladder","Colorectal","Esophagus","Headandneck","Gastric","Liver" ,"Lung","Pancreas", "Prostate", "Renal", "Ovary", "Uterine")
#' @param the_sex Male or Female
#' @param the_model_type Weibull or  Loglogistic
#' @param param_table a data frame with columns: intercept (based on survreg fit), scale (based on survreg fit) site, stage, sex, model_type
#' @param ID optional ID to set seed
#' @param novel_tx_received NEW (Step 5): 0/1/NA. If 1, death time is drawn from treatment_surv_table instead of the standard parametric distribution.
#'   Defaults to 0 so all existing callers are unaffected.
#' @param treatment_surv_table NEW (Step 5): pre-computed empirical curves (pred_treated / treatment_parametric_surv_fits), used only when
#'   novel_tx_received == 1. Defaults to NULL.
#' @return A data frame with "time" and "status" elements, where status=1 indicates the event occurred.
#' @export
#' @import survobj
#sim_cancer_death_param <- function(the_stage, the_cancer_site, the_sex,ID=NA,the_model_type, param_table){
sim_cancer_death_param <- function(the_stage, the_cancer_site, the_sex, the_race = NA, ID=NA, the_model_type,
                                   param_table, novel_tx_received = 0, treatment_surv_table = NULL){

 if(!is.na(ID)){set.seed(ID)}

  # ==============================================================================
  # NEW novel-treatment branch. If this person was assigned novel
  # treatment, draw a death time from treatment_surv_table)
  # ==============================================================================
  if(!is.na(novel_tx_received) && novel_tx_received == 1){

  # NEW checkpoint - only stops here when someone actually gets novel treatment
  #  browser()

    if(is.null(treatment_surv_table)){
      stop("novel_tx_received == 1 but no treatment_surv_table was provided.")
    }

    curve <- treatment_surv_table %>%
      filter(stage == paste(the_stage),
             cancer_site == paste(the_cancer_site),
             sex == paste(the_sex),
             race == paste(the_race))

    if(nrow(curve) == 0){
      stop("No matching novel-treatment survival curve for: ",
           the_stage, " / ", the_cancer_site, " / ", the_sex, " / ", the_race)
    }

    death_info <- gettime(time = curve$time, surv = curve$surv)
    return(death_info$time)
  }
  # ==============================================================================
  # END NEW branch. Everything below is the same.
  # ==============================================================================
 # browser()

  # Filter the survival distribution based on the type and stage
  survival_dist_indiv = filter(param_table,cancer_site == paste(the_cancer_site), stage==paste(the_stage),
                               sex==paste(the_sex),model_type==paste(the_model_type))


 # browser()

  if(length(survival_dist_indiv$model_type == "Loglogistic")==0){browser()}

  if(survival_dist_indiv$model_type=="Loglogistic"){
    the_survobj=s_loglogistic(intercept = survival_dist_indiv$intercept, scale =survival_dist_indiv$scale)
  }

  if(survival_dist_indiv$model_type=="Weibull"){
    the_survobj=s_weibull(intercept = survival_dist_indiv$intercept, scale =survival_dist_indiv$scale)
  }

  death_time=rsurv(the_survobj,n=1)

 # browser()

  return(death_time)
}


#' @export
sim_cancer_deaths_screen_no_screen<-function(clinical_diagnosis_time,clinical_diagnosis_stage,cancer_site,sex,ID,
                                             screen_diagnosis_stage,surv_param_table, optimistic_surv_param_table=NULL){

  if(!is.na(ID)){
    set.seed(ID)
  }

  #Set surv param table for screen-detected individuals to optimistic if required
  if(!is.null(optimistic_surv_param_table)){
    screen_surv_param_table=optimistic_surv_param_table
  }else{
    screen_surv_param_table=surv_param_table
  }

  #cancer death time in absence of screening
  cancer_death_time_no_screen=clinical_diagnosis_time+sim_cancer_death_param(the_stage=clinical_diagnosis_stage,
                                                                             the_cancer_site=cancer_site,
                                                                             the_sex=sex,
                                                                             the_model_type="Loglogistic",
                                                                             param_table=surv_param_table,ID=ID)



  #for non-optimistic scenario, generate cancer death under screening as before
  if(is.null(optimistic_surv_param_table)){

    cancer_death_time_screen=ifelse((screen_diagnosis_stage!=clinical_diagnosis_stage&!is.na(screen_diagnosis_stage))&!
                                      is.na(clinical_diagnosis_stage),
                                    clinical_diagnosis_time+
                                      sim_cancer_death_param(the_stage="Early",
                                                             the_cancer_site=cancer_site,
                                                             the_sex=sex,
                                                             the_model_type="Loglogistic",
                                                             param_table=screen_surv_param_table,
                                                             ID=ID),
                                    cancer_death_time_no_screen)
  }else{#optimistic scenario--use optimistic surv param tables

    cancer_death_time_screen=ifelse(!is.na(screen_diagnosis_stage)&!
                                      is.na(clinical_diagnosis_stage),
                                    clinical_diagnosis_time+
                                      sim_cancer_death_param(the_stage=screen_diagnosis_stage,
                                                             the_cancer_site=cancer_site,
                                                             the_sex=sex,
                                                             the_model_type="Loglogistic",
                                                             param_table=screen_surv_param_table,
                                                             ID=ID),
                                    cancer_death_time_no_screen)


  }


  return(data.frame(cancer_death_time_screen=cancer_death_time_screen, cancer_death_time_no_screen=cancer_death_time_no_screen,ID=ID,cancer_site=cancer_site))
}

