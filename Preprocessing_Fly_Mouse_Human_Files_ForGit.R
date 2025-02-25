library(tidyverse)

## Cutting down the mouse data set size (saved initially on local drive)
uniprot_mouse <- read_csv("Galaxy9-UCSC Mouse_unipAliSwissprot_2.csv")

# Step 1: Remove unwanted columns first
message("Step 1: Removing columns")
step1_data <- uniprot_mouse %>%
  dplyr::select(-thickStart, -thickEnd, -oChromStart, -oCDS,
                -protShortNames, -geneSynonyms, -score, -seqType,
                -misMatch, -repMatch, -nCount, -blockCount,
                -strand, -oStrand, -protAltFullNames, -protAltShortNames, -isoIds, -transList, -ensGene,
                -ensProt, -ensTrans, -refSeq, -match, -uniprotName, -reserved, -blockSizes, -functionText)
glimpse(step1_data)

# Step 2: Calculate metrics
message("Step 2: Calculating metrics")
step2_data <- step1_data %>%
  mutate(
    aa_length = str_length(oSequence),
    gene_size = chromEnd - chromStart,
    intergenic_distance = chromStart - lag(chromEnd),
    relative_pos = chromStart / chromSize
  )
glimpse(step2_data)

# Step 3: Remove duplicates
message("Step 3: Removing duplicates")
step3_data <- step2_data %>%
  distinct(oSequence, .keep_all = TRUE)
glimpse(step3_data)

# Step 4: Apply filters
message("Step 4: Applying filters")
clean_mouse <- step3_data %>%
  filter(
    status == "Manually reviewed (Swiss-Prot)",
    aa_length >= 50,
    aa_length <= 10000
  )
glimpse(clean_mouse)

# Write the final cleaned file and loaded this to GitHub
write_csv(clean_mouse, "Galaxy9_UCSC_Mouse_unipAliSwissprot_clean.csv")

#####

## Cutting down the human data set - stored on local drive
uniprot_human <- read_csv("Galaxy7_UCSC_Human_ unipAliSwissprot_genome_original.csv")

# Step 1: Remove unwanted columns first
message("Step 1: Removing columns")
step1_data <- uniprot_human %>%
  dplyr::select(-thickStart, -thickEnd, -oChromStart, -oCDS,
                -protShortNames, -geneSynonyms, -score, -seqType,
                -misMatch, -repMatch, -nCount, -blockCount,
                -strand, -oStrand, -protAltFullNames, -protAltShortNames, -isoIds, -transList, -ensGene,
                -ensProt, -ensTrans, -refSeq, -match, -uniprotName, -reserved, -blockSizes, -functionText)
glimpse(step1_data)

# Step 2: Calculate metrics
message("Step 2: Calculating metrics")
step2_data <- step1_data %>%
  mutate(
    aa_length = str_length(oSequence),
    gene_size = chromEnd - chromStart,
    intergenic_distance = chromStart - lag(chromEnd),
    relative_pos = chromStart / chromSize
  )
glimpse(step2_data)

# Step 3: Remove duplicates
message("Step 3: Removing duplicates")
step3_data <- step2_data %>%
  distinct(oSequence, .keep_all = TRUE)
glimpse(step3_data)

# Step 4: Apply filters
message("Step 4: Applying filters")
clean_human <- step3_data %>%
  filter(
    status == "Manually reviewed (Swiss-Prot)",
    aa_length >= 50,
    aa_length <= 10000
  )
glimpse(clean_human)

# Write the final cleaned file uploaded to GitHub
write_csv(clean_human, "Galaxy7-UCSC_Human_ unipAliSwissprot_genome_clean.csv")


# 3. Process unique proteins for fly

## Cutting down the mouse data set size (saved initially on local drive)
uniprot_fly <- read_csv("Galaxy6_UCSC_D_melanogaster_unipAliSwissprot_genome.csv")

# Step 1: Remove unwanted columns first
message("Step 1: Removing columns")
step1_data <- uniprot_fly %>%
  dplyr::select(-thickStart, -thickEnd, -oChromStart, -oCDS,
                -protShortNames, -geneSynonyms, -score, -seqType,
                -misMatch, -repMatch, -nCount, -blockCount,
                -strand, -oStrand, -protAltFullNames, -protAltShortNames, -isoIds, -transList, -ensGene,
                -ensProt, -ensTrans, -refSeq, -match, -uniprotName, -reserved, -blockSizes, -functionText)
glimpse(step1_data)

# Step 2: Calculate metrics
message("Step 2: Calculating metrics")
step2_data <- step1_data %>%
  mutate(
    aa_length = str_length(oSequence),
    gene_size = chromEnd - chromStart,
    intergenic_distance = chromStart - lag(chromEnd),
    relative_pos = chromStart / chromSize
  )
glimpse(step2_data)

# Step 3: Remove duplicates
message("Step 3: Removing duplicates")
step3_data <- step2_data %>%
  distinct(oSequence, .keep_all = TRUE)
glimpse(step3_data)

# Step 4: Apply filters
message("Step 4: Applying filters")
clean_fly <- step3_data %>%
  filter(
    status == "Manually reviewed (Swiss-Prot)",
    aa_length >= 50,
    aa_length <= 10000
  )
glimpse(clean_fly)

# Write the final cleaned file uploaded to GitHub
write_csv(clean_fly, "Galaxy6_UCSC_D_melanogaster_unipAliSwissprot_genome.csv")


