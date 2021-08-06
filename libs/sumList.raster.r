sumList.raster <- function(rs) {
    r = rs[[1]]
    for (i in rs[-1]) r = r + i
    r
}
