# Utility Functions

logist <- function(eta) return(exp(eta) / (1 + exp(eta)))

nextPlot <- function() {
  # this functions asks the user if wants to continue to
  # cycle through the diagnostics plots or quit the precedure
  prompt <- "Enter anything to continue or [q] to quit: "
  UInput <- readline(prompt=prompt)
  if (UInput %in% c("q", "Q")) return(invisible("no"))
  return("yes")
}
