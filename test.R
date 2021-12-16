get_yesterday <- function() {
  x <- (as.character(Sys.Date() -1))
  return (x)
}