#' Summarise modality communities in a catmodgraph
#'
#' Summarises the structure of communities detected in a \code{catmodgraph}
#' object. The output reports which modalities belong to each community,
#' which variables contribute to each community, and the internal edge
#' cohesion. This is descriptive output for interpreting the community
#' structure of a modality co-association graph.
#'
#' @param x A \code{catmodgraph} object that has already been processed by
#'   \code{\link{cluster_modalities}()}.
#'
#' @return A list of class \code{"catmodcommunity"} with components:
#'   \describe{
#'     \item{\code{community_summary}}{Data frame with one row per modality
#'       community, reporting size, variables represented, and mean internal
#'       edge weight.}
#'     \item{\code{community_members}}{Data frame listing all modality
#'       nodes and their assigned community.}
#'     \item{\code{variable_composition}}{Table of variable counts by
#'       community, useful for seeing which variables span which
#'       communities.}
#'   }
#'
#' @details
#' Communities are interpreted as \strong{groups of modalities (factor
#' levels) that co-associate across different variables}, not as
#' respondent segments or latent classes. For respondent-level
#' segmentation, use the \pkg{poLCA} or \pkg{FactoMineR} packages.
#'
#' @examples
#' df <- expand_table(Titanic)
#' mg <- build_modality_graph(df)
#' mg <- cluster_modalities(mg)
#' comm <- summarise_modality_communities(mg)
#' comm
#'
#' @importFrom igraph as_data_frame induced_subgraph E ecount
#' @export
summarise_modality_communities <- function(x) {
  
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  }
  
  if (is.null(x$membership)) {
    stop("Run `cluster_modalities()` before summarising modality communities.",
         call. = FALSE)
  }
  
  verts <- igraph::as_data_frame(x$graph, what = "vertices")
  verts$community <- x$membership[verts$name]
  
  community_ids <- sort(unique(verts$community))
  
  community_summary <- do.call(
    rbind,
    lapply(community_ids, function(k) {
      sub_nodes <- verts$name[verts$community == k]
      sub_vars  <- verts$variable[verts$community == k]
      
      subg <- igraph::induced_subgraph(x$graph, vids = sub_nodes)
      ecount_sub <- igraph::ecount(subg)
      
      mean_weight <- if (ecount_sub > 0L) {
        mean(igraph::E(subg)$weight, na.rm = TRUE)
      } else {
        NA_real_
      }
      
      data.frame(
        community            = k,
        n_modalities         = length(sub_nodes),
        n_variables          = length(unique(sub_vars)),
        variables            = paste(sort(unique(sub_vars)), collapse = ", "),
        mean_internal_weight = mean_weight,
        stringsAsFactors     = FALSE
      )
    })
  )
  
  community_members <- verts[, c("name", "variable", "modality", "community"),
                             drop = FALSE]
  
  variable_composition <- with(
    community_members,
    table(community = community, variable = variable)
  )
  
  out <- list(
    community_summary    = community_summary,
    community_members    = community_members,
    variable_composition = variable_composition
  )
  
  class(out) <- "catmodcommunity"
  out
}

#' @export
print.catmodcommunity <- function(x, ...) {
  cat("catmodcommunity object\n")
  cat("  Modality communities:", nrow(x$community_summary), "\n")
  cat("  Community sizes     :",
      paste(x$community_summary$n_modalities, collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.catmodcommunity <- function(object, ...) {
  object$community_summary
}