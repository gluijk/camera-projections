# Linear to Cylindrical, Equirectangular, Spherical and Stereographic camera projections
# www.overfitting.net
# https://www.overfitting.net/

library(Rcpp)
library(tiff)


# Add grid to input image
add_grid <- function(img,
                     n_gridx = 12,
                     linewidth = 4,
                     colour = c(1, 1, 0)) {
    
    # For 3:2 format images, n_gridx = 12 -> 12x8 squares grid
    # For 4:3 format images, n_gridx = 16 -> 16x12 squares grid
    
    stopifnot(length(dim(img)) == 3, dim(img)[3] == 3)
    
    M <- dim(img)[1]
    N <- dim(img)[2]
    step <- N / n_gridx
    
    # Determinar si la cuadrícula debe tener una línea cruzando el centro exacto.
    # Si hay un número par de celdas, hay un número impar de líneas -> Línea central.
    has_center_line <- (n_gridx %% 2 == 0)
    
    # Función pura: genera posiciones basadas en el paso y el modo de alineación
    make_lines <- function(L, step, center_line = TRUE) {
        centre <- (L + 1) / 2
        
        # Calculamos suficientes pasos para cubrir la mitad de la imagen
        max_steps <- ceiling((L / 2) / step)
        
        if (center_line) {
            # Multiplicadores enteros: 0, ±1, ±2... (cae exactamente en el centro)
            multipliers <- seq(-max_steps, max_steps, by = 1)
        } else {
            # Multiplicadores fraccionarios: ±0.5, ±1.5... (el centro queda vacío)
            multipliers <- seq(-max_steps, max_steps, by = 1) + 0.5
        }
        
        # Proyectar posiciones y filtrar las que caen fuera de los límites de la imagen
        pos <- centre + multipliers * step
        pos <- pos[pos >= 1 & pos <= L]
        
        return(sort(unique(round(pos))))
    }
    
    # Ambas dimensiones usan ahora la misma lógica de paridad geométrica
    cols <- make_lines(N, step, center_line = has_center_line)
    rows <- make_lines(M, step, center_line = has_center_line)
    
    # Ensanchar líneas según linewidth
    expand <- function(v, max_val) {
        expanded <- unlist(lapply(v, function(x) {
            x - floor((linewidth - 1) / 2) + 0:(linewidth - 1)
        }))
        unique(sort(expanded[expanded >= 1 & expanded <= max_val]))
    }
    
    cols <- expand(cols, N)
    rows <- expand(rows, M)
    
    # Aplicar color (asignación vectorizada)
    img[, cols, 1] <- colour[1]
    img[, cols, 2] <- colour[2]
    img[, cols, 3] <- colour[3]
    
    img[rows, , 1] <- colour[1]
    img[rows, , 2] <- colour[2]
    img[rows, , 3] <- colour[3]
    
    return(img)
}


# Quick visualization wrapper using base R
plot_array <- function(arr, title) {
    # Convert array layout into an R native raster representation
    # Coordinates might need a transpose depending on how you view matrix directions
    grid::grid.newpage()
    grid::grid.raster(aperm(arr, c(1, 2, 3)))
}


#######################################
# Compile all reproyection functions

# Conversions from Linear to other projections
sourceCpp("linear2cylindrical.cpp")
sourceCpp("linear2equirectangular.cpp")
sourceCpp("linear2spherical.cpp")
sourceCpp("linear2stereographic.cpp")

# Back to Linear projections (useful to check if perfect circles become ellipses)
sourceCpp("spherical2linear.cpp")
sourceCpp("stereographic2linear.cpp")



#######################################
# Examples



# Example 1: neighbourhood park with smartphone UWA (13mm FF eq.)
img=readTIFF("linear_13mm.tif")
img=add_grid(img, n_gridx = 16, colour = c(1, 0, 0), linewidth = 3)
writeTIFF(img, "linear_13mm_grid.tif")

cylindrical_13mm <- linear_to_cylindrical_rcpp(img, fl_FF_mm = 13)
# plot_array(cylindrical_13mm)
writeTIFF(cylindrical_13mm, "cylindrical_13mm.tif")



# Example 2: landscape with Laowa FF 12mm
img=readTIFF("laowa12mm.tif")
img=add_grid(img, n_gridx = 22)
writeTIFF(img, "laowa12mm_grid.tif")
    
cylindrical_12mm <- linear_to_cylindrical_rcpp(img, fl_FF_mm = 12)
equirectangular_12mm <- linear_to_equirectangular_rcpp(img, fl_FF_mm = 12)
spherical_12mm <- linear_to_spherical_perspective_rcpp(img, fl_FF_mm = 12)
stereographic_12mm <- linear_to_stereographic_rcpp(img, fl_FF_mm = 13.3)

writeTIFF(cylindrical_12mm, "cylindrical_12mm.tif")
writeTIFF(equirectangular_12mm, "equirectangular_12mm.tif")
writeTIFF(spherical_12mm, "spherical_12mm.tif")
writeTIFF(stereographic_12mm, "stereographic_12mm_laowa.tif")



# Example 3: synthetic circles for spherical 2 linear transform and back to spherical
img=readTIFF("circles.tif")
img=add_grid(img, n_gridx = 12, colour = c(1,0,0), linewidth = 1)
writeTIFF(img, "circles_grid.tif")

linear_12mm <- spherical_perspective_to_linear_rcpp(img, fl_FF_mm = 12)
spherical_12mm <- linear_to_spherical_perspective_rcpp(linear_12mm, fl_FF_mm = 12)

writeTIFF(linear_12mm, "linear_12mm.tif")
writeTIFF(spherical_12mm, "spherical_12mm_back.tif")



# Example 4: real balls
img=readTIFF("linear_12mm_real2.tif")
img=add_grid(img, n_gridx = 12, colour = c(0,1,1))
writeTIFF(img, "linear_12mm_real2_grid.tif")

cylindrical_12mm <- linear_to_cylindrical_rcpp(img, fl_FF_mm = 12)
writeTIFF(cylindrical_12mm, "cylindrical_12mm_real2.tif")

spherical_12mm <- linear_to_spherical_perspective_rcpp(img, fl_FF_mm = 16)  # 16 ~corrects spheres
writeTIFF(spherical_12mm, "spherical_12mm_real2.tif")
linear_12mm <- spherical_perspective_to_linear_rcpp(spherical_12mm, fl_FF_mm = 16)  # 16 ~corrects spheres
writeTIFF(linear_12mm, "linear_12mm_real2_back.tif")  # OK, image back is perfect

stereographic_12mm <- linear_to_stereographic_rcpp(img, fl_FF_mm = 13.3)  # 13.3 corrects spheres better
writeTIFF(stereographic_12mm, "stereographic_12mm_real2.tif")
linear_12mm <- stereographic_to_linear_rcpp(stereographic_12mm, fl_FF_mm = 13.3)  # 13.3 corrects spheres better
writeTIFF(linear_12mm, "linear_12mm_real2_back.tif")  # OK, image back is perfect



# Example 5: synthetic circles for stereographic 2 linear transform and back to stereographic
img=readTIFF("circles.tif")
img=add_grid(img, n_gridx = 12, colour = c(1,0,0), linewidth = 1)
writeTIFF(img, "circles_grid.tif")

linear_12mm <- stereographic_to_linear_rcpp(img, fl_FF_mm = 12)
stereographic_12mm <- linear_to_stereographic_rcpp(linear_12mm, fl_FF_mm = 12)

writeTIFF(linear_12mm, "linear_12mm.tif")
writeTIFF(stereographic_12mm, "stereographic_12mm_back.tif")
