###########################################################
#' Simulate a cohort of individuals with and without Multicancer Early Detection (MCED) screening
#'
#' @description
#' This function simulates cancer outcomes in a population with a designated starting age using a "parallel universe" approach.
#' That is, cancer outcomes in the population are simulated with and without MCED screening.  The natural history (i.e., the times of cancer onset and clinical diagnosis) for
#' each individual is the same in both screening and no-screening scenarios.
#'
#' The user specifies cancer sites in the MCED screening test.
#' The user also provides sensitivity of the tests for early and late-stage disease, where early refers
#' to AJCC 7 stages I-II and late, III-IV, except for pancreas cancer, where early is stage I and late, II-IV.
#'
#' Natural history models are based on built-in fitted models that are calibrated to SEER 2015-2021 data by age and sex and can be specified based on
#' user-provided inputs about the overall mean sojourn time (OMST) and the late mean sojourn time (LMST) for each cancer site.
#' The function tracks only the first cancer diagnosis based on pre-clinical onset.
#'
#' Other-cause mortality is based on all cause mortality tables from the Human Mortality Database that have been adjusted to remove the mortality due to
#' cancers included in the MCED tests.  Cancer-specific mortality in the screen arm assumes a stage-shift benefit of screening. That is, individuals
#' who are diagnosed in early stage under screening but who would have been diagnosed in late stage clinically, are assume to remain in early stage post lead-time.
#' For both screened and unscreened individuals, cancer mortality is projected from the point of clinical diagnosis (i.e., post lead-time) to prevent
#' lead-time bias.
#'
#' @param cancer_sites Vector of cancer sites (allowable values include
#'   "Anus",  "Bladder", "Esophagus", "Gastric", "Headandneck",
#'  "Liver", "Lung", "Lymphoma", "Ovary", "Pancreas", "Renal", "Uterine")
#'
#' @param LMST_vec Numeric vector of late mean sojourn times (years) for each cancer site.
#' @param OMST_vec Numeric vector of overall mean sojourn times (years) for each cancer site.
#' @param test_performance_dataframe Data frame with test sensitivity/specificity info.
#' @param MCED_specificity Overall specificity for the MCED test (not specific to cancer site)
#' @param starting_age Numeric starting age for simulation
#' @param ending_age Numeric ending age for simulation
#' @param num_screens Number of screening rounds.
#' @param screen_interval Interval between screening rounds.
#' @param num_males Number of male individuals to simulate.
#' @param num_females Number of female individuals to simulate.
#' @param all_rates_male List of transition matrices for males.
#' @param all_rates_female List of transition matrices for females.
#' @param all_meta_data_female Metadata for female cancer sites.
#' @param all_meta_data_male Metadata for male cancer sites.
#' @param cdc_data CDC mortality data.
#' @param hmd_data Human Mortality Database data.
#' @param MCED_cdc CDC data for MCED.
#' @param surv_param_table Data frame of cancer-specific survival parameters.
#' @param CRC_data CRC data.
#' @param simulation_seed Set seed for each simulation
#' @export
#'
#' @return A data frame with combined simulated results for all individuals.
#' The function returns each individual's first cancer site,  age and stage of clinical diagnosis in
#' absence of screening, age and stage at screen diagnosis, the time of other-cause mortality, and the time of cancer-specific mortality in absence
#' and presence of screening.   Cancer diagnosis and death times are presented both with and without competing other-cause mortality.
#'
#' @examples
#'library(MCEDsimCarolyn)
#'
#'# Load the other-cause mortality tables
#'data("cdc_hmd_data")
#'# Load the prefitted natural history models
#'data("combined_fits")
#'# Load the prefitted cause-specific survival models
#'data("parametric_surv_fits")
#'#load("/home/groups/CEDAR/MCED_sim/parametric_surv_fits.rda")
#'
#'theseed      <- 1
#'scenario_no  <- 3
#'
#'cancer_sites_vec <- c(
#'  "Anus",  "Bladder", "Esophagus", "Gastric", "Headandneck",
#'  "Liver", "Lung", "Lymphoma", "Ovary", "Pancreas", "Renal", "Uterine")
#'
#'OMST_vec <- rep(2, 12)
#'LMST_vec <- rep(0.5, 12)
#'
#'early_sens <- c(
#'  0.5,  0.18, 0.48, 0.33, 0.72, 0.81,
#'  0.40, 0.61, 0.60, 0.61, 0.07, 0.18)
#'
#'late_sens <- c(
#'  1.00, 0.83, 0.97, 0.94, 0.93, 1.00,
#'  0.93, 0.94, 0.90, 0.94, 0.45, 0.81)
#'
#'test_performance_dataframe <- data.frame(early_sens  = early_sens,
#'                                         late_sens   = late_sens,
#'                                         cancer_site = cancer_sites_vec)
#'
#'set.seed(123)
#'results <- sim_MCED_parallel_universe_before_CRC(cancer_sites          = cancer_sites_vec,
#'                                                 LMST_vec              = LMST_vec,
#'                                                 OMST_vec              = OMST_vec,
#'                                                 test_performance_dataframe = test_performance_dataframe,
#'                                                 starting_age          = 45,
#'                                                 ending_age            = 500,
#'                                                 num_screens           = 30,
#'                                                 screen_interval       = 1,
#'                                                 num_males             = 50,
#'                                                 num_females           = 50,
#'                                                 all_rates_male        = all_rates_male,
#'                                                 all_rates_female      = all_rates_female,
#'                                                 all_meta_data_female  = all_meta_data_female,
#'                                                 all_meta_data_male    = all_meta_data_male,
#'                                                 cdc_data              = all_cause_cdc,
#'                                                 hmd_data              = hmd_data,
#'                                                 MCED_cdc              = MCED_cdc,
#'                                                 surv_param_table      = param_table,
#'                                                 MCED_specificity      = 0.995,
#'                                                 simulation_seed       = theseed)
sim_MCED_parallel_universe_before_CRC <- function(cancer_sites,
                                                  LMST_vec,
                                                  OMST_vec,
                                                  test_performance_dataframe,
                                                  MCED_specificity,
                                                  starting_age,
                                                  ending_age,
                                                  num_screens,
                                                  screen_interval,
                                                  num_males,
                                                  num_females,
                                                  all_rates_male,
                                                  all_rates_female,
                                                  all_meta_data_female,
                                                  all_meta_data_male,
                                                  cdc_data,
                                                  hmd_data,
                                                  MCED_cdc,
                                                  surv_param_table,
                                                  optimistic_surv_param_table=NULL,
                                                  simulation_seed){



 total_individuals=num_males+num_females


 start_male=(simulation_seed-1)*(total_individuals)+1

 end_male=start_male+num_males-1
   # Create a vector of IDs
  IDs_male <- start_male:end_male

  # Female IDs: continue sequentially after males
  IDs_female <-(end_male+1):(num_females+end_male)


  # ---- Extract Sex-Specific Rate Matrices ----
  # Extract rate matrices matrices based on OMST and LMST specs (Male)
  rates_list_male = get_filtered_rates(the_omsts = OMST_vec, the_lmsts = LMST_vec,
                                       all_meta_data = all_meta_data_male,
                                       all_rates = all_rates_male, the_cancer_sites = cancer_sites)

  sites_male = rates_list_male$cancer_sites
  rates_list_male = rates_list_male$rates_list


  # Extract rate matrices matrices based on OMST and LMST specs (Female)
  rates_list_female = get_filtered_rates(the_omsts = OMST_vec, the_lmsts = LMST_vec,
                                         all_meta_data = all_meta_data_female,
                                         all_rates = all_rates_female, the_cancer_sites = cancer_sites)
  sites_female = rates_list_female$cancer_sites
  rates_list_female = rates_list_female$rates_list

  # ---- Extract Test Performance Parameters ----
  # Extract sensitivities and specificity based on selected cancer sites
  test_performance_male = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_male))
  test_performance_female = test_performance_dataframe %>% filter(cancer_site %in% as.vector(sites_female))

  #Get the other-cause death tables for men and women
  other_cause_death_male=make_othercause_death_table(cdc_data=cdc_data,
                                                     MCED_cdc=MCED_cdc,
                                                     hmd_data=hmd_data,
                                                     the_starting_age = starting_age,
                                                     the_sex="Male",
                                                     selected_cancers=sites_male,
                                                     the_year=2018)

  other_cause_death_female=make_othercause_death_table(cdc_data=cdc_data,
                                                       MCED_cdc=MCED_cdc,
                                                       hmd_data=hmd_data,
                                                       the_starting_age = starting_age,
                                                       the_sex="Female",
                                                       selected_cancers=sites_female,
                                                       the_year=2018)

  # ---- Simulate Individual Outcomes ----
  # Use mapply to apply the sim_individual_MCED function to each ID (males)
  results_list_male <- mapply(sim_individual_MCED,
                              ID = IDs_male,
                              MoreArgs = list(rates_list=rates_list_male,
                                              cancer_sites=sites_male,
                                              test_performance=test_performance_male,
                                              other_cause_death_dist=other_cause_death_male,
                                              starting_age=starting_age,
                                              num_screens=num_screens,
                                              screen_interval=screen_interval,
                                              end_time=ending_age,
                                              surv_param_table=surv_param_table,
                                              optimistic_surv_param_table=optimistic_surv_param_table,

                                              sex="Male",MCED_specificity=MCED_specificity),
                              SIMPLIFY = FALSE)

  # Use mapply to apply the sim_individual_MCED function to each ID (females)
  results_list_female <- mapply(sim_individual_MCED,
                                ID = IDs_female,
                                MoreArgs = list(rates_list=rates_list_female,
                                                cancer_sites=sites_female,
                                                test_performance=test_performance_female,
                                                other_cause_death_dist=other_cause_death_female,
                                                starting_age=starting_age,
                                                num_screens=num_screens,
                                                screen_interval=screen_interval,
                                                end_time=ending_age,
                                                surv_param_table=surv_param_table,
                                                optimistic_surv_param_table=optimistic_surv_param_table,
                                                sex="Female",
                                                MCED_specificity=MCED_specificity),
                                SIMPLIFY = FALSE)


  #Get the first cancer and additional cancers for all individuals (female)
  first_site_female=lapply(results_list_female,"[[","first_result")
  additional_sites_female=lapply(results_list_female,"[[","stored_result")

  #Get the first cancer and additional cancers for all individuals (male)
  first_site_male=lapply(results_list_male,"[[","first_result")
  additional_sites_male=lapply(results_list_male,"[[","stored_result")


  # Combine all individual results (first cancers)
  combined_first_results_males <- do.call(rbind, first_site_male)%>%mutate(sex="Male")
  combined_first_results_females <- do.call(rbind, first_site_female)%>%mutate(sex="Female")
  combined_first_results=bind_rows(combined_first_results_males,combined_first_results_females)%>%
    mutate(start_age=starting_age,end_time=ending_age)

  # Combine all individual results (additional cancers)
  combined_additional_results_males <- do.call(rbind, additional_sites_male)%>%mutate(sex="Male")
  combined_additional_results_females <- do.call(rbind, additional_sites_female)%>%mutate(sex="Female")
  combined_additional_results=bind_rows(combined_additional_results_males,combined_additional_results_females)%>%
    mutate(start_age=starting_age,end_time=ending_age)

  return(list(
    combined_additional_results=combined_additional_results,
    combined_first_results=combined_first_results
  ))
}

