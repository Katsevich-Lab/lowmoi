% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mann_whitney_perm_2.R
\name{mann_whitney_perm}
\alias{mann_whitney_perm}
\title{Mann-whitney test (with permutations)}
\usage{
mann_whitney_perm(
  response_odm,
  grna_odm,
  response_grna_group_pairs,
  B = 1000,
  progress = TRUE,
  full_output = FALSE
)
}
\arguments{
\item{response_odm}{an expression ODM of responses (typically genes)}

\item{grna_odm}{an ODM of either grna expressions (counts) or grna
assignments (logical)}

\item{response_grna_group_pairs}{a data frame with columns \code{response_id} and
\code{grna_group} giving the response ID / grna group pairs to analyze.}

\item{B}{number of permutation replicates}

\item{progress}{print progress messages?}

\item{full_output}{return the full output (including the resampled statistics; TRUE) or a reduced output (FALSE)?}
}
\value{
a data frame with columns \code{response_id}, \code{grna_group}, and \code{p_value}.
}
\description{
Runs a Mann-Whitney test using permutations.
}
\details{
The \code{grna_group} column of \code{response_grna_group_pairs} contains the names of grna groups.
These grna groups are assumed to be present in the \code{target} column of the feature
covariate matrix of \code{grna_odm}. This column should contain entries "non-targeting"
indicating the non-targeting grnas. If \code{target_type} is present as a column of the
feature covariate matrix of \code{grna_odm}, \code{target_type} is ignored.
}
\examples{
\dontrun{
response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
response_grna_group_pairs <-
 expand.grid(grna_group = c("CUL3", "CMTM6"),
             response_id = sample(ondisc::get_feature_ids(response_odm), 5))
}
}
