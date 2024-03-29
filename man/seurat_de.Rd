% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_de.R
\name{seurat_de}
\alias{seurat_de}
\title{Seurat DE}
\usage{
seurat_de(response_odm, grna_odm, response_grna_group_pairs)
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
Implements the Seurat differential expression method. The (default) wilcox test is used,
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
# a schraivogel example
response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
grna_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
response_grna_group_pairs <-
 expand.grid(grna_group = c("OXR1-TSS", "LRRCC1-TSS"),
             response_id = sample(ondisc::get_feature_ids(response_odm), 50))

# a frangieh example
response_odm <- load_dataset_modality("frangieh/control/gene")
grna_odm <- load_dataset_modality("frangieh/control/grna_assignment")
response_grna_group_pairs <- data.frame(grna_group = "A2M",
response_id = sample(ondisc::get_feature_ids(response_odm), 1))

# a papalexi example
response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
response_grna_group_pairs <- data.frame(grna_group = "CUL3",
response_id = sample(ondisc::get_feature_ids(response_odm), 4))
}
}
