layers2list <- function(r)
    lapply(layer.apply(r, function(i) c(i)), function(i) i[[1]])
