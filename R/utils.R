# Utility Functions
logist <- function(eta) return(exp(eta) / (1 + exp(eta)))

nextPlot <- function() {
  prompt <- "Enter anything to continue or [q] to quit: "
  UInput <- readline(prompt=prompt)
  if (UInput %in% c("q", "Q")) return(invisible("no"))
  return("yes")
}
