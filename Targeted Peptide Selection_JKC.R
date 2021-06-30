install.packages("OrgMassSpecR")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(kableExtra)
library(OrgMassSpecR)




#---------INTRO---------------------------------------------------------------------------

# Before starting you'll want to download the study organism's proteome from
# UniProt as an Excel file. Then convert this excel file to .csv format and
# ensure that your column headers are appropriate for what information you
# want to include. Save the proteome to your working directory in the .csv
# format.

# Here you will read in your proteome file to a data frame:
Proteome <- read.csv('C:/Users/jkcan/Desktop/Julia/R and RStudio Work/Targeted Peptide Selectrion - Bertrand Lab/Peptide Selection Code/Frag_Proteome.csv', stringsAsFactors=FALSE)

# Look over your data frame to check that data matches up with UniProt results
# (e.g. Frag has 18,075 proteins...this should be equivalent to nrow)
names(Proteome)
str(Proteome)
dim(Proteome)
nrow(Proteome)
summary(Proteome)

#---------END INTRO-------------------------------------------------------------------------


#---------PROTEIN OF INTEREST------------------------------------------------------------

# Now that the proteome is read in, you can choose a protein of interest to inspect for
# Trypsin digestion.

# Example: Rhodopsin A0A1E7G037

# Filter the data to select only for Rhodopsin by its accession code:

Protein_of_interest <- filter(Proteome, Accession == "A0A1E7G037")

sequence <- select(Protein_of_interest, Sequence)


# Check out the info for the filtered protein of interest and its sequence:

head(Protein_of_interest)


#-------END PROTEIN OF INTEREST-----------------------------------------------------------



#-------PERFORM TRYPSIN DIGESTION-------------------------------------------------------

# Cleave the full protein sequence after K and R to mimic Trypsin digestion:

# To do this we'll use the Digest function in the OrgMassSpecR package.

# Description: Cleave an amino acid sequence (a protein or peptide) according to
# enzyme specific rules and calculate the precursor ion m/z values.
# Supporting info found at: https://cran.r-project.org/web/packages/OrgMassSpecR/OrgMassSpecR.pdf

# Put the output of the tryptic digest into a data frame:
Cleaved_protein <- Digest(as.character(sequence), enzyme = "trypsin.strict", missed = 0, IAA = TRUE,
       N15 = FALSE, custom = list())


#-------END TRYPSIN DIGESTION---------------------------------------------------------



#------APPLY THE GUIDELINES FOR CHOOSING PEPTIDES (HOOFNAGLE ET AL. 2016)------

# Create lists of character/numeric values that meet all of the guidelines
# for choosing peptides:

peptides <- c(Cleaved_protein$peptide)

length <- c(nchar(Cleaved_protein$peptide))

missed_cleavage <- c(grepl(pattern = "K|R", x = Cleaved_protein$peptide))

three_or_more_Cys <- c(str_count(Cleaved_protein$peptide, fixed("C")) < 3)

Cys_or_Gln_Nterminus <- c(grepl(pattern = "C|Q", x = Cleaved_protein$peptide[1]))

Avoid_Met <- c(str_count(Cleaved_protein$peptide, fixed("M")) == 0)

Avoid_Dibasic_aa <- c(!grepl(pattern = "KK|RR|KR", x = Cleaved_protein$peptide))

Avoid_Hys <- c(str_count(Cleaved_protein$peptide, fixed("H")) == 0)


# Create a data frame that gives the result of the rule check for each peptide:
Peptides <- data.frame(peptides, length, missed_cleavage, three_or_more_Cys,
                       Cys_or_Gln_Nterminus, Avoid_Met, Avoid_Dibasic_aa,
                       Avoid_Hys, stringsAsFactors=FALSE)


# Check the rules and append a new column to the Peptides data frame
# that says YES if the peptide is a candidate.

Peptides$Candidate <- with(Peptides, ifelse(peptides %in% c(Cleaved_protein$peptide) & length >= 8 & length <= 21 &
                                      missed_cleavage == TRUE & three_or_more_Cys == TRUE & Cys_or_Gln_Nterminus == TRUE &
                                        Avoid_Met == TRUE & Avoid_Dibasic_aa == TRUE & Avoid_Hys == TRUE,
                                  "YES", "No"))

#------END APPLY RULES----------------------------------------------------------



#-----SELECT CANDIDATE PEPTIDES AND EXPORT TO FINAL TABLE-----------------------

# Combine the peptide info with the protein of interest data:

combined_data <- merge(Peptides, Protein_of_interest, all.x = TRUE, all.y = TRUE)

# Now use 'select' to select only the columns with peptides, whether or not the peptide is a candidate, and associated protein information:

Selection <- select(combined_data, peptides, Candidate, Accession, Protein_name)

# Now use 'filter' to filter through the selection for candidate peptides ONLY:

Candidate_Peptides <- filter(Selection, Candidate == "YES")


#-----END SELECT CANDIDATES AND EXPORT------------------------------------------------------




#-------RESULTS TABLE FORMATTING------------------------------------------------------

# Present resulting candidate peptides in a finalized table:
Candidate_Peptides %>%
  kbl(caption = "Candidate Peptide Information") %>%
  kable_material(c("striped", "hover")) %>%
  column_spec(3, width = "10em")




# Other table format options:
# Table format option 1
Candidate_Peptides %>%
  kbl(caption = "Candidate Peptide Information") %>%
  kable_classic(full_width = F, html_font = "Cambria")
save_kable("My_Candidate_Peptides.png")

# Table format option 2
Candidate_Peptides %>%
  kbl(caption = "Candidate Peptide Information") %>%
  kable_material(c("striped", "hover")) %>%
  column_spec(3, width = "10em")

# Table format option 3
kbl(Candidate_Peptides) %>%
  kable_paper(full_width = F) %>%
  column_spec(3, width = "30em", background = "yellow")



#-----END RESULTS TABLE FORMATTING-----------------------------------------------------










