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

###########################################################
#' Combine MCED Results with Colorectal Cancer Data
#'
#' @description
#' Integrates colorectal cancer (CRC) outcomes with MCED simulation results. The function:
#' (1) filters additional cancers to those clinically diagnosable before other-cause death,
#' (2) simulates cancer-specific survival for these additional cancers,
#' (3) reassigns additional cancers to individuals without a primary cancer diagnosis using
#'     5-year other-cause-death age strata, and
#' (4) computes competing-risk outcomes (diagnosis ages/events and death ages/events) for
#'     both screening and no-screening scenarios, including an overdiagnosis indicator.
#'
#' @param merged_CRC_MCED_results  Data frame that combine the CRC with MCEDsim data for a particular scenario.
#' @param combined_additional_results Data frame of additional (non-first) cancer diagnoses
#' @param starting_age Numeric starting age for cohort
#' @param ending_age Numeric ending age for simulation
#' @param surv_param_table Data frame of cancer-specific survival parameters
#' @param optimistic_surv_param_table Data frame of cancer-specific optimistic survival parameters
#'
#' @details
#' The reassignment algorithm attempts to match each additional cancer to an individual
#'  without a primary cancer diagnosis in the same 5-year other-cause-death age stratum.
#'  If a match is found, the additional cancer is reassigned to that individual; unmatched
#'  additional cancers increment no_match_counter.
#' @export
#'
#' @return A list with two elements: results(A data frame for all individuals after reassignment and
#'                                   competing-risk processing)
#'                                   no_match_counter (total number of additional cancers that could
#'                                   not be reassigned to a no-primary-cancer individual within the same age stratum)
#'
#' @examples
#' out <- combine_MCED_CRC(merged_CRC_MCED_results = merged_CRC_MCED_results,
#'                        combined_additional_results = combined_additional_results,
#'                         starting_age = 45,
#'                        ending_age   = 500,
#'                        CRC_data = CRC_data,
#'                        surv_param_table = param_table
#' )
combine_MCED_CRC<- function(merged_CRC_MCED_results,
                            starting_age,
                            ending_age,
                            surv_param_table,
                            optimistic_surv_param_table=NULL,
                            pairs_table = NULL,
                            screen_interval=1,
                            num_screens=30,
                            fp_rate = 0.005) {


  # Pull merged MCED+CRC cohort
  combined_first_results=merged_CRC_MCED_results$joined_data

# ################ verison 6 ######
# Add eligibility column if it doesn't exist
 #   if (!"clinical_diagnosis_time_eligibility" %in% names(combined_first_results)) {
#     combined_first_results <- combined_first_results %>% mutate(clinical_diagnosis_time_eligibility = NA_real_)
#   }

   combined_first_results <- combined_first_results %>%
     mutate(elig_time = ifelse(!is.na(clinical_diagnosis_time_eligibility), clinical_diagnosis_time_eligibility, clinical_diagnosis_time),

  # if elig_time is still NA, set to Inf to keep in no_primary_cancer
  elig_time = ifelse(is.na(elig_time), Inf, elig_time))

   primary_cancer <- combined_first_results %>% filter(elig_time <= other_cause_death_time)

    no_primary_cancer <- combined_first_results %>% filter(elig_time > other_cause_death_time) %>%
      mutate(age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, by = 5)))

  # Build the "cancer pool" we want to assign to unaffected people:
  #  additional cancers (MCED additional cancers)
  #   excess cancers (MCED primaries displaced when CRC became primary)
  combined_additional_results <- bind_rows(merged_CRC_MCED_results$combined_additional_results, merged_CRC_MCED_results$excess_cancers)

  # Keep only cancers that could be clinically diagnosed before OC death
  # If there are any additional/excess cancers eligible for reassignment, proceed
  combined_additional_results <- combined_additional_results %>% filter(clinical_diagnosis_time <= other_cause_death_time)

  at_least_one_additional_result=dim(combined_additional_results)[1]>0
  no_match_counter=0

  # ensure pairs always exists for returning
  pairs <- pairs_table


  # If there are cancers to reassign, run matching & reassignment
  if(at_least_one_additional_result){



    #   set.seed(12345)
    # Simulate cancer-specific deaths for additional cancers
    addtl_cancer_deaths=mapply(
      FUN="sim_cancer_deaths_screen_no_screen",
      clinical_diagnosis_time  = combined_additional_results$clinical_diagnosis_time,
      clinical_diagnosis_stage = combined_additional_results$clinical_diagnosis_stage,
      cancer_site = combined_additional_results$cancer_site,
      sex = combined_additional_results$sex,
      ID = combined_additional_results$ID,
      screen_diagnosis_stage = combined_additional_results$screen_diagnosis_stage,
      MoreArgs = list(surv_param_table=surv_param_table,
                      optimistic_surv_param_table=optimistic_surv_param_table),
      SIMPLIFY = F)

    # Join cancer-specific deaths with cancer diagnoses for additional cancers
    combined_additional_results = data.frame(do.call(rbind,addtl_cancer_deaths))%>%
      inner_join(combined_additional_results, by=c("ID","cancer_site"))


    # Keep only additional/excess cancers that are clinically diagnosable before other-cause death (eligible for reassignment)
    combined_additional_results <- combined_additional_results %>%
      mutate(age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0,150,by=5))) %>%
      filter(clinical_diagnosis_time <= other_cause_death_time)

    # counts before matching (for validation)
    N_cancers_to_assign <- nrow(combined_additional_results)
    N_primary_before    <- nrow(primary_cancer)
    N_unaffected_before <- nrow(no_primary_cancer)


# cancers_to_reassign: all additional/excess cancers that need to be assigned to unaffected people
# cancers_to_reassign <- combined_additional_results %>% mutate(treat = 1L, row_id = paste0("T_", row_number()))

# using cancer_site for row_id can reduce the number of p_id for whom the LSD additivity are not holding.
   cancers_to_reassign <- combined_additional_results %>% mutate(treat = 1L, row_id = paste0("T_", ID, "_", cancer_site))


   unaffected_people <- no_primary_cancer %>% mutate(treat = 0L, row_id = paste0("C_", ID))

    # ================================
    # Reassignment using MatchIt
    # ================================


    # if pairs_table is NULL, create pairs (do matching)
    # else, use existing pairs_table (skip matching)
    if (is.null(pairs_table)) {

      # no pairs provided- do matching


      # Combine both groups for matching algorithm
      combined_for_matching <- bind_rows(cancers_to_reassign, unaffected_people)


        # Run matching: 1-to-1, without replacement
      m.out <- matchit( treat ~ other_cause_death_time,
                        data    = combined_for_matching,
                        method  = "nearest",
                        ratio   = 1,
                        replace = FALSE,
                        exact   = ~ sex
      )

      # Extract matched data
      matched_data <- match.data(m.out)
      matched_cancers  <- matched_data %>% filter(treat == 1)
      matched_unaffected  <- matched_data %>% filter(treat == 0)


      # pairs table with original IDs included
      pairs <- matched_cancers %>% rename(cancer_id = row_id) %>% select(subclass, cancer_id, cancer_original_ID = ID) %>%
        inner_join(matched_unaffected %>% rename(unaffected_id = row_id) %>% select(subclass, unaffected_id, unaffected_original_ID = ID),
                   by = "subclass")

    }
    if (nrow(pairs) > 0) {


      reassigned_cancers <- cancers_to_reassign %>% inner_join(pairs, by = c("row_id" = "cancer_id")) %>%
        left_join(unaffected_people %>% select(row_id, ID_new = ID, sex_new = sex, death_time_new = other_cause_death_time, p_id_new = p_id),
          by = c("unaffected_id" = "row_id")) %>%
        mutate(ID = ID_new,
               sex = sex_new,
               other_cause_death_time = death_time_new,
               p_id = p_id_new,
          age_OC_death_cat = cut(death_time_new, breaks = seq(0,150,by=5)))%>%
                select(-treat, -row_id, -unaffected_id, -ID_new, -p_id_new, -sex_new, -death_time_new)

    #  browser()

      # Add reassigned cancers to the primary cancer dataset
      primary_cancer <- bind_rows(primary_cancer, reassigned_cancers)


      # Remove matched individuals from the pool of unaffected people
      unaffected_who_received_cancers <- unique(reassigned_cancers$ID)
      no_primary_cancer <- no_primary_cancer %>% filter(!(ID %in% unaffected_who_received_cancers))

      # Counters
      N_matched <- nrow(reassigned_cancers)
      no_match_counter <- nrow(cancers_to_reassign) - N_matched
      N_unmatched <- no_match_counter

      # Sanity check: after reassignment
      N_primary_after    <- nrow(primary_cancer)
      N_unaffected_after <- nrow(no_primary_cancer)
    }

  }


  # Final data with all cancers combined (first cancers and reassigned cancers)
  combined_results=bind_rows(primary_cancer,no_primary_cancer)

  end_time <- ending_age

  #Process data with other cause death as a censoring event

  #Ascertain age at at screen and clinical diagnosis in presence of other cause death
  #Ascertain age at death under screening scenarios and no screening scenarios in presence of other cause death
  #Ascertain age at diagnosis under screening scenario (can be either screen or clinical)
  #Ascertain mode of diagnosis under screening scenario
  #Ascertain stage at diagnsosis under screening scenario
  #Ascertain if individual was overdiagnosed (screen detected but died due to other causes prior to clinical diagnosis)
  combined_results=combined_results %>% mutate(clin_dx_age = pmin(other_cause_death_time,clinical_diagnosis_time,end_time,na.rm = T),
                                               clin_dx_event = case_when(
                                                 clin_dx_age == other_cause_death_time ~ "other_cause_death",
                                                 clin_dx_age == end_time ~ "censor",
                                                 clin_dx_age == clinical_diagnosis_time ~ "clin_cancer_diagnosis",
                                                 .default = NA
                                               ),
                                               clin_dx_event_stage = case_when(clin_dx_event == "clin_cancer_diagnosis" & clinical_diagnosis_stage == "Early"~1,
                                                                               clin_dx_event == "clin_cancer_diagnosis" & clinical_diagnosis_stage == "Late"~2,
                                                                               .default = 3),
                                               screen_dx_age = pmin(other_cause_death_time,screen_diagnosis_time,end_time,na.rm = T),
                                               screen_dx_event = case_when(
                                                 screen_dx_age == other_cause_death_time ~ "other_cause_death",
                                                 screen_dx_age == end_time ~ "censor",
                                                 screen_dx_age == screen_diagnosis_time ~ "screen_cancer_diagnosis",
                                                 .default = NA
                                               ),
                                               screen_dx_event_stage = case_when(screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Early"~1,
                                                                                 screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Late"~2,
                                                                                 .default = 3),
                                               death_age_no_screen=pmin(other_cause_death_time,cancer_death_time_no_screen,end_time,na.rm = T),
                                               death_age_screen=pmin(other_cause_death_time,cancer_death_time_screen,end_time,na.rm = T),
                                               death_event_no_screen=case_when(
                                                 death_age_no_screen == other_cause_death_time ~ "other_cause_death",
                                                 death_age_no_screen == end_time ~ "censor",
                                                 death_age_no_screen == cancer_death_time_no_screen ~ "cancer_death",
                                                 .default = NA
                                               ),
                                               death_event_screen=case_when(
                                                 death_age_screen == other_cause_death_time ~ "other_cause_death",
                                                 death_age_screen == end_time ~ "censor",
                                                 death_age_screen == cancer_death_time_screen ~ "cancer_death",
                                                 .default = NA
                                               ),
                                               diagnosis_age_screen_scenario=pmin(clin_dx_age,screen_dx_age,na.rm=T),
                                               diagnosis_event_screen_scenario=ifelse(screen_dx_age<=clin_dx_age, screen_dx_event,
                                                                                      clin_dx_event),
                                               diagnosis_event_stage_screen_scenario=case_when(screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Early"~1,
                                                                                               screen_dx_event == "screen_cancer_diagnosis" & screen_diagnosis_stage == "Late" ~2,
                                                                                               (screen_dx_event!="screen_cancer_diagnosis" & clin_dx_event=="clin_cancer_diagnosis") & clinical_diagnosis_stage=="Early"~1,
                                                                                               (screen_dx_event !="screen_cancer_diagnosis" & clin_dx_event=="clin_cancer_diagnosis") & clinical_diagnosis_stage=="Late"~2,
                                                                                               .default = 3),
                                               life_years_diff=death_age_screen-death_age_no_screen,
                                               overdiagnosis=ifelse(screen_dx_event=="screen_cancer_diagnosis"&clin_dx_event=="other_cause_death",1,0)

    # Late-stage reduction: cancer would have been Late clinically but caught Early by screening
#    late_stage_reduction = ifelse(clin_dx_event_stage == 2 & diagnosis_event_stage_screen_scenario == 1, 1, 0)

  )

# test..
#  combined_results <- combined_results %>% select(ID, onset_time) %>% sample_n(50)
#  test_data <- test_data %>% mutate(n_false_positives = rbinom(n(), size = n_screens_before_onset, prob = fp_rate))

# Calculate the number of screens each person received before cancer onset, capped at the maximum number of screens (num_screens).
#   Uses floor() to count only fully completed screening intervals between starting_age and onset_time.
 combined_results <- combined_results %>%
   mutate(n_screens_before_onset = pmin(floor((onset_time - starting_age) / screen_interval), num_screens))

#
# n_screens_before_onset = pmin(floor((onset_time - starting_age) / screen_interval),
#                               floor(num_screens/screen_interval))


# Helper function to draw the number of false positives per person using binomial sampling.
#   Each persons draw is seeded by their ID to ensure reproducibility across runs.

 draw_false_positives <- function(ID, n_screens, fp_rate) {
   if (is.na(n_screens) || is.na(ID)) return(NA)
   set.seed(as.integer(ID))
   rbinom(1, size = n_screens, prob = fp_rate)
 }

# Apply the false positive draw row-wise across all individuals.
# Each person's number of false positives is drawn from Binomial(n_screens, fp_rate), where n_screens is their number of screens
#  before onset and fp_rate is the per-screen false positive rate.
 combined_results <- combined_results %>%
   mutate(n_false_positives = mapply(draw_false_positives,
                                     ID = ID,
                                     n_screens = n_screens_before_onset,
                                     fp_rate = fp_rate))


  return(list(combined_results = combined_results, pairs = pairs,no_match_counter = no_match_counter))
}

#########################################
#' Combine MCED and CRC data for specific age category and sex

#' @description
#' Integrates CRC screening outcomes with MCED results within a specified 5-year
#'  other-cause-death age category, processing males and females separately and
#'  then combining the outputs.
#'
#' @param the_CRC_data Processed CRC data with age categories
#' @param combined_additional_results MCED additional cancer outcomes
#' @param combined_first_results MCED first cancer diagnoses
#' @param the_age_cat Character string: age category (e.g., "(45,50]")
#'
#' @export
#' @return A list with with results_all: Combined MCED + CRC outcomes for age category
#'                  no_match_all: Total number of unmatched reassignments summed
#'                                 over males and females
#'
#'
#' @examples
#' out_age_bin <- combine_by_age_cat_sex(
#'                            the_CRC_data = crc_processed,
#'             combined_additional_results = combined_addtl,
#'                  combined_first_results = combined_first,
#'                             the_age_cat = "(60,65]")
#'
combine_by_age_cat_sex<-function(the_CRC_data,
                                 combined_additional_results,
                                 combined_first_results,
                                 the_age_cat){

  # Filter CRC data by sex and age category
  the_CRC_data_female=filter(the_CRC_data,sex=="Female"&age_OC_death_cat==the_age_cat)
  the_CRC_data_male=filter(the_CRC_data,sex=="Male"&age_OC_death_cat==the_age_cat)

  # Add age categories to MCED datasets
  combined_additional_results<-combined_additional_results%>%
    mutate(age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0,150,by=5)))
  combined_first_results<-combined_first_results%>%
    mutate(age_OC_death_cat=cut(other_cause_death_time,breaks=seq(0,150,by=5)))

  # Filter MCED data by sex and age category.
  combined_additional_results_male=combined_additional_results%>%filter(sex=="Male"&age_OC_death_cat==the_age_cat)
  combined_first_results_male=combined_first_results%>%filter(sex=="Male"&age_OC_death_cat==the_age_cat)

  combined_additional_results_female=combined_additional_results%>%filter(sex=="Female"&age_OC_death_cat==the_age_cat)
  combined_first_results_female=combined_first_results%>%filter(sex=="Female"&age_OC_death_cat==the_age_cat)

  # Combine MCED and CRC results for males
  results_CRC_male <- combine_MCED_CRC(combined_first_results=combined_first_results_male,
                                       combined_additional_results=combined_additional_results_male,
                                       ending_age=500,
                                       starting_age=45,
                                       CRC_data=the_CRC_data_male,
                                       surv_param_table=param_table)

  # Combine MCED and CRC results for females
  results_CRC_female <- combine_MCED_CRC(combined_first_results=combined_first_results_female,
                                         combined_additional_results=combined_additional_results_female,
                                         ending_age=500,
                                         starting_age=45,
                                         CRC_data=the_CRC_data_female,
                                         surv_param_table=param_table)

  # Combine results across sexes
  results_all=bind_rows(results_CRC_female$results,results_CRC_male$results)
  # Sum no-match counts across sexes
  no_match_all=results_CRC_female$no_match_counter+results_CRC_male$no_match_counter

  return(list(results_all=results_all,no_match_all=no_match_all))
}












