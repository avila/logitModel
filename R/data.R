#' Donner Party Dataset
#'
#' "During the summer of 1846, the 87 members of the Donner Party travelled an
#' untested route between Fort Bridger, Wyoming, and the Humboldt River, Nevada.
#' As a result of delays caused by this choice, they became stranded by heavy snows
#' on the east flank of the Sierra Nevada. By the time the last member of the party
#' was rescued, 40 had died" [Donald K. Grayson, 1990].
#'
#' The main families in the Donner party were: Donner, Graves, Breen and Reed.
#' The families of Murphy, Foster and Pike are grouped as 'MurFosPik', those of
#' Fosdick and Wolfinger are coded as 'FosdWolf', and all others as 'Other'.
#' This dataset was retrieved from the vcdExtra package.
#'
#' @format A data frame with 90 obs. of 5 variables:
#' \describe{
#' \item{family}{Factor w/ 10 levels. Family name.}
#' \item{age}{int. Age.}
#' \item{sex}{Factor w/ 2 levels "Female","Male".}
#' \item{survived}{int. 1 = Suvived, 0 = Decesead.}
#' \item{death}{POSIXct, format: "1846-12-29". }
#' }
#' @usage data(DonnerData)
#' @source \url{https://en.wikipedia.org/wiki/Donner_Party}
#' @source \url{https://cran.r-project.org/web/packages/vcdExtra/index.html}
#' @references Donner Party Deaths: A Demographic Assessment, Donald K. Grayson,
#' Journal of Anthropological Research, Vol. 46, No. 3 (Autumn, 1990), pp. 223-242
#' @examples
#' require(logitModel)
#' data(DonnerData)
#' pairs(DonnerData)
"DonnerData"
