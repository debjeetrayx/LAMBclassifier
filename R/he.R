# ---- Optional H&E layer for Xenium results ----
# Explorer's matrix.csv maps full-resolution H&E pixels to Xenium pixels; cells
# are placed on the H&E with its inverse.

.LAMB_HE_TARGET <- 6000   # long side, in pixels, of the H&E read (whole scan)
.LAMB_HE_DRAW <- 2500     # long side, in pixels, of the section as drawn

.lamb_he_layer <- function(he, xy_um, pixel_size, say) {
    if (is.null(he$alignment)) {
        stop("he_image needs he_alignment: the matrix.csv that Xenium Explorer exports when the H&E is aligned.",
             call. = FALSE)
    }
    if (!file.exists(he$image)) stop("H&E image not found: ", he$image, call. = FALSE)
    if (!file.exists(he$alignment)) stop("Alignment file not found: ", he$alignment, call. = FALSE)
    M <- as.matrix(suppressWarnings(utils::read.csv(he$alignment, header = FALSE)))   # Explorer omits the final newline
    if (!identical(dim(M), c(3L, 3L)) || !is.numeric(M)) {
        stop("he_alignment must be Xenium Explorer's 3 x 3 matrix.csv.", call. = FALSE)
    }
    im <- .lamb_he_read(he$image, he$downsample)
    X <- cbind(xy_um$x / pixel_size, xy_um$y / pixel_size, 1)
    H <- X %*% t(solve(M))                  # Xenium pixels -> full-resolution H&E pixels
    xy <- cbind(H[, 1] / H[, 3], H[, 2] / H[, 3]) / im$downsample
    inside <- xy[, 1] >= 0 & xy[, 2] >= 0 & xy[, 1] <= ncol(im$raster) & xy[, 2] <= nrow(im$raster)
    say("H&E: ", he$image, " read at 1/", round(im$downsample, 2), " of full resolution (",
        ncol(im$raster), " x ", nrow(im$raster), " px); alignment ", he$alignment, "; ",
        round(100 * mean(inside, na.rm = TRUE), 1), "% of cells fall on the image.")
    if (mean(inside, na.rm = TRUE) < 0.9) {
        say("Warning: fewer than 90% of cells land on the H&E. Check that the alignment file belongs to ",
            "this section and that he_downsample is right.")
    }
    # crop to the cells' extent (a scan can hold more than one section)
    ok <- inside & !is.na(inside)
    if (!any(ok)) stop("No cell lands on the H&E image: the alignment file probably belongs to another image or section.", call. = FALSE)
    pad <- 0.03 * max(diff(range(xy[ok, 1])), diff(range(xy[ok, 2])))
    x0 <- max(1, floor(min(xy[ok, 1]) - pad)); x1 <- min(ncol(im$raster), ceiling(max(xy[ok, 1]) + pad))
    y0 <- max(1, floor(min(xy[ok, 2]) - pad)); y1 <- min(nrow(im$raster), ceiling(max(xy[ok, 2]) + pad))
    ras <- im$raster[y0:y1, x0:x1]
    xy <- cbind(xy[, 1] - x0 + 1, xy[, 2] - y0 + 1)
    step <- max(1L, ceiling(max(dim(ras)) / .LAMB_HE_DRAW))   # keep the PDF small
    if (step > 1) {
        ras <- ras[seq(1, nrow(ras), step), seq(1, ncol(ras), step)]
        xy <- xy / step
    }
    list(raster = ras, xy = xy)
}

.lamb_he_read <- function(path, downsample) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("png", "jpg", "jpeg")) {
        pkg <- if (ext == "png") "png" else "jpeg"
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the '", pkg, "' package to read this image.", call. = FALSE)
        if (is.null(downsample)) {
            stop("For a PNG / JPEG H&E, say how many times smaller it is than the aligned full-resolution ",
                 "image, for example he_downsample = 16.", call. = FALSE)
        }
        a <- if (ext == "png") png::readPNG(path) else jpeg::readJPEG(path)
        return(list(raster = grDevices::as.raster(a[, , seq_len(min(3, dim(a)[3])), drop = FALSE]), downsample = downsample))
    }
    if (!requireNamespace("terra", quietly = TRUE)) {
        stop("Install the 'terra' package to read a TIFF H&E, or export it as PNG.", call. = FALSE)
    }
    # noflip: H&E files carry no map coordinates; without it terra returns them upside down
    r <- tryCatch(suppressWarnings(terra::rast(path, noflip = TRUE)), error = function(e) NULL)
    if (is.null(r)) {
        stop("This H&E file cannot be read (often JPEG2000-compressed OME-TIFF). Export a smaller copy as PNG, ",
             "for example in QuPath (File > Export images, downsample 16), and run again with ",
             "he_image = \"that.png\", he_downsample = 16.", call. = FALSE)
    }
    full_w <- terra::ncol(r)
    f <- max(1, max(terra::ncol(r), terra::nrow(r)) / .LAMB_HE_TARGET)
    s <- terra::spatSample(r[[seq_len(min(3, terra::nlyr(r)))]], size = ceiling(terra::ncell(r) / f^2),
                           method = "regular", as.raster = TRUE)   # uses the file's overviews when present
    a <- terra::as.array(s)
    if (max(a, na.rm = TRUE) > 1) a <- a / 255
    a[is.na(a)] <- 1
    if (dim(a)[3] == 1) a <- array(rep(a, 3), c(dim(a)[1:2], 3))
    list(raster = grDevices::as.raster(a), downsample = full_w / terra::ncol(s))
}

.lamb_draw_he <- function(img, rx, ry) {
    graphics::plot(NA, xlim = c(0, ncol(img$raster)), ylim = c(nrow(img$raster), 0), asp = 1, axes = FALSE,
                   xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    graphics::rasterImage(img$raster, 0, nrow(img$raster), ncol(img$raster), 0, interpolate = TRUE)
}
