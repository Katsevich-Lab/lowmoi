% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liscovitch_method.R
\name{liscovitch_method}
\alias{liscovitch_method}
\title{Liscovitch method}
\usage{
liscovitch_method(response_odm, grna_odm, response_grna_group_pairs)
}
\arguments{
\item{response_odm}{an expression ODM of responses (typically genes)}

\item{grna_odm}{an ODM of either grna expressions (counts) or grna
assignments (logical)}

\item{response_grna_group_pairs}{a data frame with columns \code{response_id} and
\code{grna_group} giving the response ID / grna group pairs to analyze.}
}
\value{
a data frame with columns \code{response_id}, \code{grna_group}, and \code{p_value}.
}
\description{
Implements the differential expression method of Liscovitch
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
response_odm <- load_dataset_modality("liscovitch/experiment_small/chromatin")
grna_odm <- load_dataset_modality("liscovitch/experiment_small/grna_assignment")
response_grna_group_pairs <- expand.grid(response_id = response_odm |> ondisc::get_feature_ids(),
grna_group = c("ARID1A", "ATRX"))
res <- liscovitch_method(response_odm, grna_odm, response_grna_group_pairs)
}
}
