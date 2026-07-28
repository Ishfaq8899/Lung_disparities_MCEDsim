#' Build revised CRC-SPIN person-level data for a screening scenario
#'
#' Reads CRC-SPIN person-level simulation output for a specified screening scenario
#' and returns a standardized dataset with harmonized variable names used in
#' downstream MCED/CRC merging workflows.
#'
#' Supported scenarios include:
#' \itemize{
#'   \item \code{"control"} (no screening)
#'   \item \code{"colonoscopy"}
#'   \item \code{"annual_fit"}
#'   \item \code{"mced_annual"}, \code{"mced_biennial"}, \code{"mced_triennial"} (MCED including CRC;
#'   returns screened trajectory and joins control CRC trajectory by \code{p_id})
#'   \item \code{"CRC_single_cancer_annual"}, \code{"CRC_single_cancer_biennial"},
#'   \code{"CRC_single_cancer_triennial"} (CRC-only MCED; implemented via CRC-SPIN scenarioIDs 4/5/6)
#'   \item \code{"MCED_annual_no_CRC"}, \code{"MCED_biennial_no_CRC"},
#'   \code{"MCED_triennial_no_CRC"} (MCED used for non-CRC cancers only; CRC remains control trajectory,
#'   joined by \code{p_id})
#' }
#'
#' Output columns depend on scenario type:
#' \itemize{
#'   \item Clinical scenarios return \code{clinical_diagnosis_time},
#'   \code{cancer_death_time_no_screen}, \code{clinical_diagnosis_stage}.
#'   \item Screen scenarios return \code{screen_diagnosis_time},
#'   \code{cancer_death_time_screen}, \code{screen_diagnosis_stage}, and may also include
#'   control variables joined by \code{p_id}.
#' }
#'
#' @param crc_spin_scenario Character scalar specifying the CRC-SPIN scenario name. Must be one of
#' \code{"control"}, \code{"colonoscopy"}, \code{"annual_fit"},
#' \code{"mced_annual"}, \code{"mced_biennial"}, \code{"mced_triennial"},
#' \code{"CRC_single_cancer_annual"}, \code{"CRC_single_cancer_biennial"},
#' \code{"CRC_single_cancer_triennial"},
#' \code{"MCED_annual_no_CRC"}, \code{"MCED_biennial_no_CRC"}, \code{"MCED_triennial_no_CRC"}.
#'
#' @param CRC_SPIN_dir Path to directory containing CRC-SPIN person-level \code{.rds}
#' files named like \code{scr_parmID_0_blockID_0_scenarioID_<ID>_person.rds}.
#'
#'
#' @return A data.frame containing harmonized CRC-SPIN person-level data. At minimum includes
#' \code{sex}, \code{p_id}, \code{start_age}, \code{other_cause_death_time}, \code{cancer_site}, and
#' \code{age_OC_death_cat}, plus scenario-dependent diagnosis/death/stage fields.
#'
#' @details
#' Sex is standardized to \code{"Male"}/\code{"Female"} using the CRC-SPIN \code{fem} indicator.
#' \code{cancer_site} is set to \code{"Colorectal"}. Stage is collapsed to \code{"Early"} for stages
#' 1--2 and \code{"Late"} otherwise. Other-cause death age is binned into 5-year intervals using
#' \code{cut(other_cause_death_time, breaks = seq(0, 150, 5))}.
#'
#' @examples
#' \dontrun{
#' crc_control <- build_crc_data_revised("control", CRC_SPIN_dir = "path/to/crc_spin")
#' crc_fit     <- build_crc_data_revised("annual_fit", CRC_SPIN_dir = "path/to/crc_spin")
#'
#' crc_list <- list(control = crc_control, annual_fit = crc_fit)
#' }
#'
#' @export
build_crc_data_revised <- function(crc_spin_scenario, CRC_SPIN_dir) {

    if (!crc_spin_scenario %in% c("control", "colonoscopy", "annual_fit",
                                  "mced_annual", "mced_biennial", "mced_triennial",
                                  "CRC_single_cancer_annual", "CRC_single_cancer_biennial", "CRC_single_cancer_triennial",
                                  "MCED_annual_no_CRC", "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")) {
      stop("Invalid crc_spin_scenario")
    }

    if (crc_spin_scenario == "control") {
      control_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))

      # CONTROL (no screening)
      crc_data_revised <- control_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc,
               clinical_diagnosis_time = age_clindet,
               cancer_death_time_no_screen = age_death_crc_noscreen) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               clinical_diagnosis_stage = case_when(
                 stage_clindet %in% c(1, 2) ~ "Early",
                 stage_clindet %in% c(3, 4) ~ "Late",
                 TRUE ~ NA_character_),
               cancer_site = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
        select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
               cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site,
               age_OC_death_cat, age_death, n_preclin)

      # =========================
      # Colonoscopy (scenarioID 2)
      # =========================
    } else if (crc_spin_scenario == "colonoscopy") {
      colonoscopy_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_2_person.rds"))


      crc_data_revised <- colonoscopy_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               clinical_diagnosis_time     = pmin(age_screendet, age_clindet, na.rm=TRUE),
               cancer_death_time_no_screen = pmin(age_death_crc_screendet, age_death_crc_noscreen, na.rm=TRUE),
               clinical_diagnosis_stage    = ifelse(stage_det %in% c(1,2), "Early", "Late"),
               cancer_site      = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0,150,5))) %>%
        select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
               cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site,
               age_OC_death_cat, age_death, n_preclin)

      #==================================================
      # Scenario 3: Annual FIT screening (crc_spin_scenarioID = 3 )
      #   - CRC-SPIN FIT strategy
      #   - Screen diagnosis time = age_screendet
      #   - Cancer death time = age_death_crc_screendet
      #==================================================
    } else if (crc_spin_scenario == "annual_fit") {
      annual_fit_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_3_person.rds"))

      crc_data_revised <- annual_fit_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               clinical_diagnosis_time     = pmin(age_screendet, age_clindet, na.rm=TRUE),
               cancer_death_time_no_screen = pmin(age_death_crc_screendet, age_death_crc_noscreen, na.rm=TRUE),
               clinical_diagnosis_stage    = ifelse(stage_det %in% c(1,2), "Early", "Late"),
               cancer_site      = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0,150,5))) %>%
        select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
               cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site, age_OC_death_cat, age_death, n_preclin)

      #==================================================
      # Scenarios 4–6: MCED including CRC
      #   mced_annual / mced_biennial / mced_triennial
      #==================================================
    } else if (crc_spin_scenario %in% c("mced_annual", "mced_biennial", "mced_triennial")) {

      screen_id <- case_when(crc_spin_scenario == "mced_annual" ~ 4L,
                             crc_spin_scenario == "mced_biennial" ~ 5L,
                             TRUE ~ 6L)

      control_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))
      screen_data  <- readRDS(file.path(CRC_SPIN_dir, paste0("scr_parmID_0_blockID_0_scenarioID_", screen_id, "_person.rds")))

      control_data_revised <- control_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc,
               clinical_diagnosis_time = age_clindet,
               cancer_death_time_no_screen = age_death_crc_noscreen) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               clinical_diagnosis_stage = case_when(
                 stage_clindet %in% c(1, 2) ~ "Early",
                 stage_clindet %in% c(3, 4) ~ "Late",
                 TRUE ~ NA_character_),
               cancer_site = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
        select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time,
               cancer_death_time_no_screen, clinical_diagnosis_stage, cancer_site,
               age_OC_death_cat, age_death, n_preclin)

      crc_data_revised <- screen_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc,
               screen_diagnosis_time = age_screendet,
               cancer_death_time_screen = age_death_crc_screendet) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               screen_diagnosis_stage = case_when(
                 stage_screendet %in% c(1, 2) ~ "Early",
                 stage_screendet %in% c(3, 4) ~ "Late",
                 TRUE ~ NA_character_),
               cancer_site = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
        select(sex, p_id, start_age, other_cause_death_time, screen_diagnosis_time,
               cancer_death_time_screen, screen_diagnosis_stage, cancer_site,
               age_OC_death_cat, age_death, n_preclin) %>%
        left_join(select(control_data_revised, p_id, clinical_diagnosis_time,
                         clinical_diagnosis_stage, cancer_death_time_no_screen), by = "p_id")

      #==================================================
      # Scenarios 7–9: CRC_single_cancer (annual/biennial/triennial)
      #==================================================
    } else if (crc_spin_scenario %in% c("CRC_single_cancer_annual", "CRC_single_cancer_biennial", "CRC_single_cancer_triennial")) {

      screen_id <- case_when(crc_spin_scenario == "CRC_single_cancer_annual" ~ 4L,
                             crc_spin_scenario == "CRC_single_cancer_biennial" ~ 5L,
                             TRUE ~ 6L)

      control_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))
      CRC_single_cancer_data <- readRDS(file.path(CRC_SPIN_dir, paste0("scr_parmID_0_blockID_0_scenarioID_", screen_id, "_person.rds")))

      # version 6
      crc_data_revised <- CRC_single_cancer_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               clinical_diagnosis_time = pmin(age_clindet, age_screendet, na.rm = TRUE),
               clinical_diagnosis_stage = case_when(
                 stage_det %in% c(1, 2) ~ "Early",
                 stage_det %in% c(3, 4) ~ "Late",
                 TRUE ~ NA_character_),
               cancer_death_time_no_screen = age_death_crc_noscreen,
               cancer_site = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
        select(sex, p_id, start_age, other_cause_death_time,
               clinical_diagnosis_time,
               cancer_death_time_no_screen, clinical_diagnosis_stage,
               cancer_site, age_OC_death_cat, age_death, n_preclin)



      #==================================================
      # Scenarios 10–12: MCED_no_CRC (annual/biennial/triennial)
      #==================================================
    } else if (crc_spin_scenario %in% c("MCED_annual_no_CRC", "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")) {

      screen_id <- case_when(crc_spin_scenario == "MCED_annual_no_CRC" ~ 4L,
                             crc_spin_scenario == "MCED_biennial_no_CRC" ~ 5L,
                             TRUE ~ 6L)

      control_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))
      MCED_no_CRC_data <- readRDS(file.path(CRC_SPIN_dir, paste0("scr_parmID_0_blockID_0_scenarioID_", screen_id, "_person.rds")))

      control_data_revised <- control_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc,
               clinical_diagnosis_time = age_clindet,
               cancer_death_time_no_screen = age_death_crc_noscreen) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
               clinical_diagnosis_stage = case_when(
                 stage_clindet %in% c(1, 2) ~ "Early",
                 stage_clindet %in% c(3, 4) ~ "Late",
                 TRUE ~ NA_character_),
               cancer_site = "Colorectal",
               screen_diagnosis_stage=clinical_diagnosis_stage,
               screen_diagnosis_time=clinical_diagnosis_time,
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
        select(sex, p_id, start_age, other_cause_death_time, clinical_diagnosis_time, screen_diagnosis_time,
               cancer_death_time_no_screen, clinical_diagnosis_stage, screen_diagnosis_stage, cancer_site,
               age_OC_death_cat, age_death, n_preclin)

      crc_data_revised <- MCED_no_CRC_data %>%
        rename(start_age = age_at_index,
               other_cause_death_time = age_death_oc,
            #   screen_diagnosis_time = age_clindet,
               cancer_death_time_screen = age_death_crc_noscreen) %>%
        mutate(sex = ifelse(fem == FALSE, "Male", "Female"),
       #        screen_diagnosis_stage = case_when(
      #           stage_clindet %in% c(1, 2) ~ "Early",
      #           stage_clindet %in% c(3, 4) ~ "Late",
       #          TRUE ~ NA_character_),
               cancer_site = "Colorectal",
               age_OC_death_cat = cut(other_cause_death_time, breaks = seq(0, 150, 5))) %>%
        select(sex, p_id, start_age, other_cause_death_time,
               cancer_death_time_screen, cancer_site,age_OC_death_cat, age_death, n_preclin) %>%
        left_join(select(control_data_revised, p_id, clinical_diagnosis_time, screen_diagnosis_time, screen_diagnosis_stage,
                         clinical_diagnosis_stage, cancer_death_time_no_screen), by = "p_id")
    }

    final_control_data <- readRDS(file.path(CRC_SPIN_dir, "scr_parmID_0_blockID_0_scenarioID_1_person.rds"))

    final_control_data <- final_control_data %>% mutate(clinical_diagnosis_time_eligibility = age_clindet) %>%
      select(p_id, clinical_diagnosis_time_eligibility)

    crc_data_revised <- crc_data_revised %>% left_join(final_control_data, by ="p_id")


    return(crc_data_revised)
  }

# ==============================================================
#' Merge CRC-SPIN and MCED simulation outputs for a specified scenario
#'
#' Loads MCED simulation output (annual/biennial/triennial) and merges it with a
#' preprocessed CRC-SPIN scenario dataset (stored in \code{crc_list} inside an \code{.RData} file).
#' Individuals are aligned by sorting by sex, assigning a synthetic \code{new_ID} within each dataset,
#' and joining MCED and CRC records on \code{new_ID}. Sex-mismatched pairs are removed after the join.
#'
#'
#' The function supports two main classes of scenarios:
#' \itemize{
#'   \item \strong{Clinical scenarios} (\code{control}, \code{colonoscopy}, \code{annual_fit},
#'     \code{CRC_single_cancer_*}): CRC diagnosis/death/stage are treated as "clinical"
#'     variables in CRC-SPIN and merged into MCED outputs for non-CRC cancers.
#'   \item \strong{Screen scenarios} (\code{mced_*}, \code{MCED_*_no_CRC}): CRC diagnosis/death/stage
#'     are treated as "screen" variables and merged into MCED outputs in a manner consistent with
#'     MCED screening assumptions.
#' }
#'
#' After merging, the function overwrites MCED \code{other_cause_death_time} with the CRC-SPIN value,
#' determines whether CRC becomes the primary cancer site, moves displaced MCED primaries to an
#' \code{excess_cancers} dataset, and updates \code{combined_additional_results} to match the final cohort.
#'
#' After merging, the function:
#' \itemize{
#'   \item uses CRC-SPIN \code{other_cause_death_time} to overwrite MCED other-cause death time
#'   \item determines whether CRC becomes the primary cancer site (if CRC diagnosis precedes
#'     other-cause death and no earlier primary cancer exists per MCED timing)
#'   \item flags individuals with both a primary CRC assignment and an MCED primary cancer
#'     prior to other-cause death, exporting the MCED cancer as an "excess/additional" cancer
#'   \item updates \code{combined_additional_results} to use CRC-SPIN other-cause death time and
#'     filters additional cancers to match the final joined cohort
#' }
#'
#' @param scenario_name Must be one of:
#' \code{"control"}, \code{"colonoscopy"}, \code{"annual_fit"},
#' \code{"CRC_single_cancer_annual"}, \code{"CRC_single_cancer_biennial"},
#' \code{"CRC_single_cancer_triennial"},
#' \code{"mced_annual"}, \code{"mced_biennial"}, \code{"mced_triennial"},
#' \code{"MCED_annual_no_CRC"}, \code{"MCED_biennial_no_CRC"}, \code{"MCED_triennial_no_CRC"}.
#'
#' @param mced_data_dir Path to directory containing MCED combined result files named like
#' \code{mced_scenario_<id>_combined_results.Rdata}, where \code{id} is:
#' \itemize{
#'   \item 1 for annual scenarios
#'   \item 2 for biennial scenarios
#'   \item 3 for triennial scenarios
#' }
#'
#' @param crc_data_path Path to an \code{.RData} file that contains \code{crc_list}, a named list
#' of CRC-SPIN processed scenario data frames keyed by \code{scenario_name}.
#'
#' @return A named \code{list} with components:
#' \itemize{
#'   \item \code{joined_data}: merged person-level dataset combining MCED and CRC-SPIN outputs
#'   \item \code{excess_cancers}: MCED primary cancers for individuals whose primary assignment
#'     becomes CRC (treated as additional cancers)
#'   \item \code{combined_additional_results}: MCED additional cancer dataset filtered/updated
#'     to align with \code{joined_data} and updated other-cause death time
#'   \item \code{scenario}: the input \code{scenario_name}
#' }
#'
#' @details
#' This function expects the MCED \code{.Rdata} file to load objects named
#' \code{combined_first_results} and \code{combined_additional_results}.
#' It also expects the CRC \code{.RData} to load an object named \code{crc_list}.
#'
#' The merge uses \code{sex} ordering and assigns \code{new_ID = seq_len(n)} separately
#' to MCED and CRC datasets after sorting by \code{sex}. Individuals with sex mismatches
#' between MCED and CRC after joining are removed for internal consistency.
#'
#' @examples
#' \dontrun{
#' res <- process_mced_crc_data(
#'   scenario_name = "mced_annual",
#'   mced_data_dir  = "path/to/MCED_data",
#'   crc_data_path  = "path/to/crc_spin_initial_process_data.RData"
#' )
#'
#' head(res$joined_data)
#' }
#'
#'
#' @export
process_mced_crc_data <- function(scenario_name, mced_data_dir, crc_data_path,
                                  window_type="OMST_1_LMST_0.5") {

  # Map window type to base file index
  window_map <- c(OMST_1_LMST_0.5 = 0L,
                  OMST_2_LMST_1   = 3L,
                  OMST_2_LMST_0.5 = 6L,
                  Litrature_LMST_1   = 9L,
                  Litrature_LMST_0.5 = 12L)

  window_base <- unname(window_map[window_type])

  if (is.na(window_base)) {
    stop("Invalid window_type: ", window_type)
  }

  # Pick which MCED results file to load based on scenario name
  # Each window_type corresponds to a block of 3 MCED files:
  # annual    -> window_base + 1
  # biennial  -> window_base + 2
  # triennial -> window_base + 3

  if (scenario_name %in% c("mced_biennial", "MCED_biennial_no_CRC", "CRC_single_cancer_biennial")) {
    mced_id <- window_base + 2L
  } else if (scenario_name %in% c("mced_triennial", "MCED_triennial_no_CRC", "CRC_single_cancer_triennial")) {
    mced_id <- window_base + 3L
  } else {
    mced_id <-  window_base + 1L
  }

  # Load MCED annual/biennial/triennial results
  load(file.path(mced_data_dir, paste0("mced_scenario_", mced_id, "_combined_results.Rdata")))
  if (!exists("combined_first_results")) stop("combined_first_results not found in MCED file")

  # Load CRC list (all 12 scenarios)
  load(crc_data_path)

  # Define which scenarios use clinical vs screen variables
  clinical_scenarios <- c("control", "colonoscopy", "annual_fit","CRC_single_cancer_annual", "CRC_single_cancer_biennial", "CRC_single_cancer_triennial")
  screen_scenarios   <- c("mced_annual", "mced_biennial", "mced_triennial", "MCED_annual_no_CRC", "MCED_biennial_no_CRC", "MCED_triennial_no_CRC")

  # Determine scenario type
  if (scenario_name %in% clinical_scenarios) {
    scenario_type <- "clinical"
  } else if (scenario_name %in% screen_scenarios) {
    scenario_type <- "screen"
  } else {
    stop("Invalid scenario_name: ", scenario_name)
  }

  #============================================
  # Merge processed CRC-SPIN data with MCED data for scenarios without MCED screening for non-CRC cancers
  # Scenarios include control, colonoscopy, annual_fit, CRC_single_cancer
  # The processed CRC-SPIN specifies CRC cancer diagnoses with "clinical" labeled variables, because these are used
  # when converting to MCEDsim variables in absence of MCED screening for non-CRC cancers.
  # Natural history assumption: fast-fast-optimistic
  #At the end of the script, we have two files: excess cancers that includes non-CRC cancers (storing time of clinical diagnosis)
  #     joined_data that consists of the merged CRC-SPIN and MCEDsim files
  #============================================
  if (scenario_type == "clinical") {

    # Get CRC data from crc_list
    CRC_data <- crc_list[[scenario_name]]


######## version 6 ####################
    CRC_data <- CRC_data %>% mutate(fem = ifelse(sex == "Female", TRUE, FALSE)) %>%
      select(c("fem","sex","p_id","start_age","other_cause_death_time",
                      "clinical_diagnosis_time",
                      "clinical_diagnosis_time_eligibility",
                      "cancer_death_time_no_screen","clinical_diagnosis_stage",
                      "cancer_site","n_preclin","age_death"))

    CRC_data <- CRC_data %>% arrange(p_id, sex) %>% mutate(new_ID = seq(1, nrow(CRC_data)))
    combined_first_results <- combined_first_results %>% arrange(ID, sex) %>% mutate(new_ID = seq(1, nrow(combined_first_results)))

    joined_data <- combined_first_results %>%
      left_join(select(CRC_data, any_of(c("fem","sex","p_id","new_ID",
                                          "other_cause_death_time",
                                          "clinical_diagnosis_time",
                                          "clinical_diagnosis_time_eligibility",
                                          "cancer_death_time_no_screen",
                                          "clinical_diagnosis_stage",
                                          "cancer_site","n_preclin","age_death"))),
                by = c("new_ID" = "new_ID"))


    #If sex allocations are not aligned, remove individuals with mismatched sex in joined data--may also have to remove these people from combined additional cancers
    #Also set other-cause death time to be CRC-SPIN OC death time
    joined_data <- joined_data %>% filter(sex.x == sex.y) %>% mutate(other_cause_death_time = other_cause_death_time.y, sex = sex.x)


    ###### New
    # If individual has CRC pre-cursor lesions that would become clinical (n_preclin>0) and no primary cancer before OC death, assign CRC to be the cancer site of primary cancer.
    joined_data <- joined_data %>% mutate(cancer_site = ifelse((!is.na(n_preclin) & n_preclin > 0), cancer_site.y, as.character(cancer_site.x)))

    # If individual has both pre-cursor lesions for CRC and one or more MCED cancers before age_death: move the MCED cancer to an additional cancer file
    joined_data <- joined_data %>% mutate(MCED_and_CRC = ifelse(cancer_site == "Colorectal" & onset_time < age_death, TRUE, FALSE))


    # Extract individuals with both MCED and CRC cancers to a separate dataset
    # These MCED cancers will be treated as "additional cancers" rather than primary
    excess_cancers <- subset(joined_data, MCED_and_CRC == TRUE) %>%
      mutate(clinical_diagnosis_time     = clinical_diagnosis_time.x,
             cancer_death_time_no_screen = cancer_death_time_no_screen.x,
             clinical_diagnosis_stage    = clinical_diagnosis_stage.x,
             cancer_site = cancer_site.x) %>%
      select(-ends_with(".x"), -ends_with(".y")) %>%
      select(c(names(combined_first_results)))



###### version 4-5 ######
# Set variables in joined_data depending on if they have CRC or MCED as their primary cancer
# If CRC is the primary cancer, use CRC values (.y), otherwise use MCED values (.x)
   joined_data <- joined_data %>%
     mutate(clinical_diagnosis_time     = ifelse(cancer_site == "Colorectal", clinical_diagnosis_time.y,     clinical_diagnosis_time.x),
            cancer_death_time_no_screen = ifelse(cancer_site == "Colorectal", cancer_death_time_no_screen.y, cancer_death_time_no_screen.x),
            clinical_diagnosis_stage    = ifelse(cancer_site == "Colorectal", clinical_diagnosis_stage.y,    clinical_diagnosis_stage.x)) %>%
     select(-ends_with(".x"), -ends_with(".y"))


    #============================================
    # Merge processed CRC-SPIN data with MCED data for scenarios with MCED screening for non-CRC cancers
    # Scenarios include mced, MCED_no_CRC
    # The processed CRC-SPIN specifies CRC cancer diagnoses with "screen" labeled variables, because these are used
    # when converting to MCEDsim variables in presence of MCED screening for non-CRC cancers.
    # Natural history assumption: fast-fast-optimistic
    #At the end of the script, we have two files: excess cancers that includes non-CRC cancers (storing time of screened diagnosis)
    #                             joined_data that consists of the merged CRC-SPIN and MCEDsim files
    #============================================
  } else if (scenario_type == "screen") {

    # Get CRC data from crc_list
    CRC_data             <- crc_list[[scenario_name]]

    control_process_data <- crc_list$control
    CRC_data <- CRC_data %>% left_join(control_process_data %>% select(p_id), by="p_id")

    # Load processed CRC data
    CRC_data <- CRC_data %>% mutate(fem = ifelse(sex == "Female", TRUE, FALSE)) %>%
      select(c("fem","sex","p_id","start_age","other_cause_death_time","screen_diagnosis_time",
               "cancer_death_time_screen","screen_diagnosis_stage","cancer_site",
               "clinical_diagnosis_time",
               "clinical_diagnosis_time_eligibility",                                # added
               "clinical_diagnosis_stage","n_preclin","age_death"))

    # sort by sex and set new ids in both datasets
    CRC_data               <- CRC_data %>% arrange(p_id, sex) %>% mutate(new_ID = seq(1, nrow(CRC_data)))
    combined_first_results <- combined_first_results %>% arrange(ID, sex) %>% mutate(new_ID = seq(1, nrow(combined_first_results)))

    # Join datasets by new id
    # Note: variable names ending in .x are from MCEDsim .y are from CRC-spin
    joined_data <- combined_first_results %>%
      left_join(CRC_data %>% select(fem, sex, p_id, new_ID, other_cause_death_time,
                                    screen_diagnosis_time, cancer_death_time_screen, screen_diagnosis_stage,
                                    cancer_site, clinical_diagnosis_time, clinical_diagnosis_stage,
                                    clinical_diagnosis_time_eligibility,                   # added
                                    n_preclin, age_death), by="new_ID")

    #If sex allocations are not aligned, remove individuals with mismatched sex in joined data--may also have to remove these people from combined additional cancers
    #Also set other-cause death time to be CRC-SPIN OC death time
    joined_data <- joined_data %>% filter(sex.x == sex.y) %>% mutate(other_cause_death_time = other_cause_death_time.y, sex = sex.x)

    # If individual has CRC pre-cursor lesions that would become clinical (n_preclin>0) and no primary cancer before OC death, assign CRC to be the cancer site of primary cancer.
    joined_data <- joined_data %>% mutate(cancer_site = ifelse((!is.na(n_preclin) & n_preclin > 0), cancer_site.y, as.character(cancer_site.x)))

    # If individual has both pre-cursor lesions for CRC and one or more MCED cancers before age_death: move the MCED cancer to an additional cancer file
    joined_data <- joined_data %>% mutate(MCED_and_CRC = ifelse(cancer_site == "Colorectal" & onset_time < age_death, TRUE, FALSE))


    ########
    # Extract individuals with both MCED and CRC cancers to a separate dataset
    # These MCED cancers will be treated as "additional cancers" rather than primary
    excess_cancers <- subset(joined_data, MCED_and_CRC == TRUE) %>%
      mutate(screen_diagnosis_time    = screen_diagnosis_time.x,
             cancer_death_time_screen = cancer_death_time_screen.x,
             screen_diagnosis_stage   = screen_diagnosis_stage.x,
             clinical_diagnosis_time  = clinical_diagnosis_time.x,
             clinical_diagnosis_stage = clinical_diagnosis_stage.x,
             cancer_site = cancer_site.x) %>%
      select(-ends_with(".x"), -ends_with(".y")) %>%
      select(c(names(combined_first_results)))



    # Set variables in joined_data depending on if they have CRC or MCED as their primary cancer
    # If CRC is the primary cancer, use CRC values (.y), otherwise use MCED values (.x)
    joined_data <- joined_data %>%
      mutate(screen_diagnosis_time    = ifelse(cancer_site == "Colorectal", screen_diagnosis_time.y,    screen_diagnosis_time.x),
             cancer_death_time_screen = ifelse(cancer_site == "Colorectal", cancer_death_time_screen.y, cancer_death_time_screen.x),
             screen_diagnosis_stage   = ifelse(cancer_site == "Colorectal", screen_diagnosis_stage.y,   screen_diagnosis_stage.x),
             clinical_diagnosis_time  = ifelse(cancer_site == "Colorectal", clinical_diagnosis_time.y,  clinical_diagnosis_time.x),
             clinical_diagnosis_stage = ifelse(cancer_site == "Colorectal", clinical_diagnosis_stage.y, clinical_diagnosis_stage.x)) %>%
      select(-ends_with(".x"), -ends_with(".y"))
  }

  # Join combined_additional_results with joined_data to update other_cause_death_time with CRC-SPIN data.
  #   Filter to remove sex-mismatched individuals (previously excluded from joined_data) to maintain consistency
  #     between primary and additional cancer datasets.
  combined_additional_results <- combined_additional_results %>%
    filter(ID %in% joined_data$ID) %>%
    select(-other_cause_death_time) %>%
    left_join(joined_data %>% select(ID, other_cause_death_time), by="ID")

  return(list(joined_data                 = joined_data,
              excess_cancers              = excess_cancers,
              scenario                    = scenario_name,
              combined_additional_results = combined_additional_results))
}



