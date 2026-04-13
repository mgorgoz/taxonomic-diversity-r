#' Anatolian Forest Trees: Multi-Site Species Data
#'
#' A data frame containing 20 tree species from Anatolian forests,
#' distributed across three sample plots with varying community
#' compositions. Species abundances follow the Westhoff & van der Maarel
#' (1973) scale (1--9). Taxonomic classification includes seven ranks
#' from species to kingdom.
#'
#' The three sites represent different forest types:
#' \describe{
#'   \item{Karisik_Orman}{Mixed forest -- both conifers and broadleaves (12 species)}
#'   \item{Yaprakli_Orman}{Broadleaf-dominated forest (13 species)}
#'   \item{Konifer_Orman}{Conifer-dominated forest (8 species)}
#' }
#'
#' @format A data frame with 33 rows and 9 columns:
#' \describe{
#'   \item{Site}{Sample plot name (character)}
#'   \item{Species}{Binomial species name with underscore separator (character)}
#'   \item{Genus}{Genus (character)}
#'   \item{Family}{Family (character)}
#'   \item{Order}{Order (character)}
#'   \item{Class}{Class (character)}
#'   \item{Phylum}{Phylum / Division (character)}
#'   \item{Kingdom}{Kingdom (character)}
#'   \item{Abundance}{Westhoff abundance value, integer 1--9 (numeric)}
#' }
#'
#' @details
#' This dataset can be used directly with \code{\link{batch_analysis}}
#' for multi-site analysis:
#'
#' \preformatted{batch_analysis(anatolian_trees)}
#'
#' To extract a single community for use with \code{\link{ozkan_pto}}
#' or \code{\link{compare_indices}}:
#'
#' \preformatted{
#' site1 <- anatolian_trees[anatolian_trees$Site == "Karisik_Orman", ]
#' community <- setNames(site1$Abundance, site1$Species)
#' tax_tree  <- site1[, c("Species", "Genus", "Family", "Order",
#'                         "Class", "Phylum", "Kingdom")]
#' ozkan_pto(community, tax_tree)
#' }
#'
#' @references
#' Westhoff, V. & van der Maarel, E. (1973). The Braun-Blanquet
#' approach. In: R.H. Whittaker (ed.), Ordination and classification
#' of communities. Handbook of Vegetation Science 5, 617--726.
#'
#' @seealso \code{\link{batch_analysis}} for multi-site analysis,
#'   \code{\link{gazi_comm}} and \code{\link{gazi_gytk}} for a
#'   single-community example.
#'
#' @examples
#' data(anatolian_trees)
#' head(anatolian_trees)
#'
#' # Multi-site analysis
#' batch_analysis(anatolian_trees)
#'
#' # Single site extraction
#' site1 <- anatolian_trees[anatolian_trees$Site == "Karisik_Orman", ]
#' comm  <- setNames(site1$Abundance, site1$Species)
#' tax   <- site1[, c("Species", "Genus", "Family", "Order",
#'                     "Class", "Phylum", "Kingdom")]
#' ozkan_pto(comm, tax)
#'
"anatolian_trees"


#' Example Community Vector: 8 Anatolian Tree Species
#'
#' A named numeric vector of species abundances for a single forest
#' community with 8 Anatolian tree species. Abundance values follow
#' the Westhoff & van der Maarel (1973) scale (1--9). This vector
#' mirrors the hypothetical example in Ozkan (2018).
#'
#' @format A named numeric vector with 8 elements.
#'   Names are species binomials (underscore-separated); values are
#'   integer abundances (1--4).
#'
#' @details
#' The species include 3 genera from Pinaceae, 2 from Fagaceae,
#' 1 each from Cupressaceae and Betulaceae, spanning 2 orders
#' (Pinales, Fagales).
#'
#' Pair with \code{\link{gazi_gytk}} for analysis:
#'
#' \preformatted{
#' ozkan_pto(gazi_comm, gazi_gytk)
#' compare_indices(gazi_comm, gazi_gytk)
#' }
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336--346.
#'
#' Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based
#' taxonomic diversity measures with the other diversity measures and
#' introduction to the new proposed (reinforced) estimators. FORESTIST,
#' 72(2). DOI: 10.5152/forestist.2021.21025
#'
#' @seealso \code{\link{gazi_gytk}} for the matching taxonomy,
#'   \code{\link{anatolian_trees}} for a multi-site dataset.
#'
#' @examples
#' data(gazi_comm)
#' data(gazi_gytk)
#'
#' # Ozkan pTO
#' ozkan_pto(gazi_comm, gazi_gytk)
#'
#' # All indices at once
#' compare_indices(gazi_comm, gazi_gytk)
#'
"gazi_comm"


#' Example Taxonomy: 8 Anatolian Tree Species
#'
#' A data frame containing the taxonomic hierarchy for the 8 species
#' in \code{\link{gazi_comm}}, with 3 taxonomic ranks (Genus, Family,
#' Order). This compact taxonomy table is designed for quick
#' demonstrations and unit testing.
#'
#' @format A data frame with 8 rows and 4 columns:
#' \describe{
#'   \item{Species}{Binomial species name (character)}
#'   \item{Genus}{Genus (character)}
#'   \item{Family}{Family (character)}
#'   \item{Order}{Order (character)}
#' }
#'
#' @details
#' The taxonomy represents:
#' \itemize{
#'   \item 8 genera: Pinus, Cedrus, Quercus, Fagus, Juniperus, Carpinus
#'   \item 4 families: Pinaceae (3 spp), Fagaceae (3 spp), Cupressaceae (1), Betulaceae (1)
#'   \item 2 orders: Pinales (4 spp), Fagales (4 spp)
#' }
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336--346.
#'
#' Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based
#' taxonomic diversity measures with the other diversity measures and
#' introduction to the new proposed (reinforced) estimators. FORESTIST,
#' 72(2). DOI: 10.5152/forestist.2021.21025
#'
#' @seealso \code{\link{gazi_comm}} for the matching community vector,
#'   \code{\link{build_tax_tree}} for building custom taxonomies.
#'
#' @examples
#' data(gazi_gytk)
#' gazi_gytk
#'
#' # Compute taxonomic distance matrix
#' tax_distance_matrix(gazi_gytk)
#'
"gazi_gytk"
