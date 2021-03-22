study_metadata =
  read_excel(
    path    = metadata_file,
    sheet   = "main",
    skip    = 0,
    trim_ws = TRUE,
    na      = "n/a",
    .name_repair = janitor::make_clean_names
  ) %>%
  select(
    -ord,
    -pos,
    -i7_index,
    -index,
    -i5_index,
    -index2,
    -correction_needed,
    -study_group,
    -rin,
    -mess,
    -abc_or_mess_or_control
  ) %>%
  rename(
    sample_name = nova_seq_sample_id,
    ethnicity = race_code,
  ) %>%
  filter(
    !is.na(sample_name),
    project %in% (project_groups_to_include %||% unique(.data$project_group)),
    project %nin% project_groups_to_exclude
  ) %>%
  mutate(
    sample_name =
      janitor::make_clean_names(
        string = sample_name,
        case = "all_caps"
        ),
    across(
      .cols =
        c(
          sex,
          ethnicity,
          run_id,
          disease_class
          ),
      .fns = as_factor
      ),
    age = as.numeric(age)
    ) %>%
  distinct()

non_project_controls =
  read_excel(
    path = main_sample_list,
    sheet = main_sample_sheet,
    .name_repair = janitor::make_clean_names
    ) %>%
  mutate(
    disease_class = tolower(disease_class),
    project_group = "control"
    ) %>%
  filter(disease_class == "control") %>%
  #Select the portions of the metadata that are useful:
  select(
    sample_name = nova_seq_sample_id,
    disease_class,
    project_group,
    sex,
    ethnicity = race_code,
    visit_ref,
    subject_ref,
    sample_alias,
    project,
    initial_concentration_ng_ul,
    final_concentration_ng_ul,
    study,
    age,
    run_id
    ) %>%
  mutate(
    across(
      .cols =
        c(
          sex,
          ethnicity,
          run_id
        ),
      .fns = as_factor
      ),
    age = as.numeric(age)
    ) %>%
  distinct()

md =
  bind_rows(
    study_metadata,
    non_project_controls
    )

tx_sample_names =
  dir(
    path = seq_file_directory,
    pattern = "quant.sf.gz",
    recursive = TRUE,
    full.name = TRUE
  ) %>%
  grep(
    pattern = "Undetermined|NONE",
    invert = TRUE,
    value = TRUE) %>%
  str_split(pattern = "/") %>%
  map_chr(~pluck(.x, length(.x)-1)) %>%
  str_remove(pattern = '(_[L|S][[:digit:]]+)+') %>%
  janitor::make_clean_names(case = "all_caps")

tx_files =
  dir(
    path = seq_file_directory,
    pattern = "quant.sf.gz",
    recursive = TRUE,
    full.name = TRUE
  ) %>%
  grep(
    pattern = "Undetermined|NONE",
    invert = TRUE,
    value = TRUE
  ) %>%
  set_names(nm = tx_sample_names) %>%
  purrr::discard(
    .p =
      is.na(
        match(
          names(.),
          md[["sample_name"]]
        )
      )
  )

final_md = filter(
  .data = md,
  sample_name %in% names(tx_files),
  sample_name %nin% samples_to_manually_remove,
  str_detect(
    string = sample_name,
    pattern = "_2$",
    negate = TRUE
    )
  ) %>%
  filter(run_id != "S4_004_2") %>%
  mutate(
    disease_class = as_factor(disease_class) %>% fct_relevel({{control_group}}),
    across(
      .cols = where(is.factor),
      .fns = fct_drop
    )
  ) %>%
  column_to_rownames('sample_name')

#Inspect the metadata:
md_cat_data = inspect_cat(final_md)

md_num_data = inspect_num(final_md)

# samples = target({
#   message("A count file was not found for the following selected samples:")
#   print(md$sample_name[(md$sample_name %nin% names(tx_files))])
pruned_samples = tx_files[rownames(final_md)]
# ),

counts =
  tximport(
    pruned_samples,
    type = "salmon",
    txIn = TRUE,
    txOut = FALSE,
    tx2gene = annot,
    importer = fread
  )

dds_import =
  DESeqDataSetFromTximport(
    txi = counts,
    colData = final_md,
    design = study_design
  )
