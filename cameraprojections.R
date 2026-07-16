# Linear (Pinhole camera) to Cylindrical, Equirectangular, Mercator, Panini, Orthographic, Spherical
# and Stereographic camera reprojections
# www.overfitting.net
# https://www.overfitting.net/

library(Rcpp)
library(tiff)


# Auxiliar functions (the projection transformation equations are implemented in cpp files)

# Generate equispaced circles
circles_matrix <- function(NUMX, NUMY, DIMX, R = 1, value.circle = 0.7, value.background = 0.3)
{
    stopifnot(NUMX >= 1, NUMY >= 1, DIMX >= NUMX, R >= 0, R <= 1)
    
    # Resolución vertical
    DIMY <- round(DIMX * NUMY / NUMX)
    # Distancia entre centros
    step <- DIMX / NUMX
    # Radio en píxeles
    radius <- R * step / 2
    # Imagen
    img <- matrix(value.background, DIMY, DIMX)
    # Coordenadas de los píxeles
    x <- seq_len(DIMX)
    y <- seq_len(DIMY)
    # Centros de los círculos
    cx <- (seq_len(NUMX) - 0.5) * step
    cy <- (seq_len(NUMY) - 0.5) * step
    
    for (yc in cy)
        for (xc in cx)
        {
            mask <- outer(
                y, x,
                function(yy, xx)
                    (xx - xc)^2 + (yy - yc)^2 <= radius^2
            )
            img[mask] <- value.circle
        }
    
    # Convertir a RGB
    out <- array(0, dim = c(DIMY, DIMX, 3))
    out[,,1] <- img
    out[,,2] <- img
    out[,,3] <- img
    out
}


# Add grid to input image
add_grid_old <- function(img, n_gridx = 12, linewidth = 4, colour = c(1, 1, 0)) {
    
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



bresenham_line <- function(x0,y0,x1,y1){
    x0<-as.integer(x0);y0<-as.integer(y0);x1<-as.integer(x1);y1<-as.integer(y1)
    dx<-abs(x1-x0); sx<-if(x0<x1)1L else -1L
    dy<--abs(y1-y0); sy<-if(y0<y1)1L else -1L
    err<-dx+dy
    xs<-integer(); ys<-integer()
    repeat{
        xs<-c(xs,x0); ys<-c(ys,y0)
        if(x0==x1 && y0==y1) break
        e2<-2*err
        if(e2>=dy){ err<-err+dy; x0<-x0+sx}
        if(e2<=dx){ err<-err+dx; y0<-y0+sy}
    }
    cbind(x=xs,y=ys)
}

draw_bresenham <- function(img, x0, y0, x1, y1, linewidth = 1, colour = c(1,1,0)) {
    M <- dim(img)[1]
    N <- dim(img)[2]
    pts <- bresenham_line(x0, y0, x1, y1)
    
    radius <- (linewidth - 1) / 2
    rmax <- ceiling(radius)
    
    for(i in seq_len(nrow(pts))) {
        x <- pts[i,1]
        y <- pts[i,2]
        
        for(dx in -rmax:rmax) {
            for(dy in -rmax:rmax) {
                # Disco
                if(dx*dx + dy*dy <= radius*radius) {
                    xx <- x + dx
                    yy <- y + dy
                    
                    if(xx >= 1 && xx <= N &&
                       yy >= 1 && yy <= M) {
                        img[yy, xx, 1] <- colour[1]
                        img[yy, xx, 2] <- colour[2]
                        img[yy, xx, 3] <- colour[3]
                    }
                }
            }
        }
    }
    
    img
}

# Add grid to input image
add_grid <- function(img, n_gridx = 12, linewidth = 4, colour = c(1, 1, 0), diagonal = TRUE) {
    # For 3:2 format images, n_gridx = 12 -> 12x8 squares grid
    # For 4:3 format images, n_gridx = 16 -> 16x12 squares grid
    
    stopifnot(length(dim(img)) == 3, dim(img)[3] == 3)
    M <- dim(img)[1]
    N <- dim(img)[2]
    step <- N / n_gridx
    
    # Si hay un número par de celdas, una línea pasa exactamente por el centro
    has_center_line <- (n_gridx %% 2 == 0)
    
    make_lines <- function(L, step, center_line = TRUE) {
        
        centre <- (L + 1) / 2
        max_steps <- ceiling((L / 2) / step)
        
        if (center_line) {
            multipliers <- seq(-max_steps, max_steps, by = 1)
        } else {
            multipliers <- seq(-max_steps, max_steps, by = 1) + 0.5
        }
        
        pos <- centre + multipliers * step
        pos <- pos[pos >= 1 & pos <= L]
        sort(unique(round(pos)))
    }
    
    expand <- function(v, max_val) {
        expanded <- unlist(lapply(v, function(x) {
            x - floor((linewidth - 1) / 2) + 0:(linewidth - 1)
        }))
        unique(sort(expanded[expanded >= 1 & expanded <= max_val]))
    }
    
    cols <- make_lines(N, step, center_line = has_center_line)
    rows <- make_lines(M, step, center_line = has_center_line)

    cols <- expand(cols, N)
    rows <- expand(rows, M)
    
    # Líneas verticales
    img[, cols, 1] <- colour[1]
    img[, cols, 2] <- colour[2]
    img[, cols, 3] <- colour[3]
    
    # Líneas horizontales
    img[rows, , 1] <- colour[1]
    img[rows, , 2] <- colour[2]
    img[rows, , 3] <- colour[3]
    
    # Diagonales
    if (diagonal) {
        img <- draw_bresenham(img, x0 = 1, y0 = 1, x1 = N, y1 = M, linewidth = linewidth, colour = colour)
        img <- draw_bresenham(img, x0 = N, y0 = 1, x1 = 1, y1 = M, linewidth = linewidth, colour = colour)
    }
    
    img
}


#######################################
# Compile all reprojection functions

# Conversions from Linear (Pinhole camera) to other projections...
sourceCpp("linear2cylindrical.cpp")  # preserves vertical lines
sourceCpp("linear2equirectangular.cpp")  # preserves vertical lines
sourceCpp("linear2mercator.cpp")  # preserves vertical lines, conformal vs 3D sphere, restores spheres as circles
sourceCpp("linear2panini.cpp")  # preserves vertical and radial lines
sourceCpp("linear2orthographic.cpp")  # preserves radial lines
sourceCpp("linear2spherical.cpp")  # preserves radial lines
sourceCpp("linear2stereographic.cpp")  # preserves radial lines, conformal vs 3D sphere, restores spheres as circles


# Back to Linear projections (useful to check if perfect circles become ellipses)
sourceCpp("mercator2linear.cpp")
sourceCpp("spherical2linear.cpp")
sourceCpp("stereographic2linear.cpp")



#######################################
# Exercises (4)


# Example 1: landscape with Laowa 12mm FF (Dpreview)
img=readTIFF("laowa12mm.tif")
img=add_grid(img, n_gridx = 12)
writeTIFF(img, "laowa12mm_grid.tif")
    
cylindrical_12mm <- linear_to_cylindrical_rcpp(img, fl_FF_mm = 12)
equirectangular_12mm <- linear_to_equirectangular_rcpp(img, fl_FF_mm = 12)
mercator_12mm <- linear_to_mercator_rcpp(img, fl_FF_mm = 12)
panini_12mm <- linear_to_panini_rcpp(img, fl_FF_mm = 12)
orthographic_12mm <- linear_to_orthographic_rcpp(img, fl_FF_mm = 12)
spherical_12mm <- linear_to_spherical_rcpp(img, fl_FF_mm = 12)
stereographic_12mm <- linear_to_stereographic_rcpp(img, fl_FF_mm = 12)

writeTIFF(cylindrical_12mm, "laowa12mm_cylindrical.tif")
writeTIFF(equirectangular_12mm, "laowa12mm_equirectangular.tif")
writeTIFF(mercator_12mm, "laowa12mm_mercator.tif")
writeTIFF(panini_12mm, "laowa12mm_panini.tif")
writeTIFF(orthographic_12mm, "laowa12mm_orthographic.tif")
writeTIFF(spherical_12mm, "laowa12mm_spherical.tif")
writeTIFF(stereographic_12mm, "laowa12mm_stereographic.tif")



# Example 2: synthetic circles for Mercator/Spherical/Stereographic to Linear transformation and back
img=circles_matrix(12, 8, 1920*4, R = 0.5)
img=add_grid(img, n_gridx = 12, colour = c(0,1,1), linewidth = 3)
writeTIFF(img, "circles12mm_grid.tif")

linear_12mm <- mercator_to_linear_rcpp(img, fl_FF_mm = 12)
mercator_12mm <- linear_to_mercator_rcpp(linear_12mm, fl_FF_mm = 12)
writeTIFF(linear_12mm, "circles12mm_mercator2linear.tif")
writeTIFF(mercator_12mm, "circles12mm_mercator_back.tif")

linear_12mm <- spherical_to_linear_rcpp(img, fl_FF_mm = 12)
spherical_12mm <- linear_to_spherical_rcpp(linear_12mm, fl_FF_mm = 12)
writeTIFF(linear_12mm, "circles12mm_spherical2linear.tif")
writeTIFF(spherical_12mm, "circles12mm_spherical_back.tif")

linear_12mm <- stereographic_to_linear_rcpp(img, fl_FF_mm = 12)
stereographic_12mm <- linear_to_stereographic_rcpp(linear_12mm, fl_FF_mm = 12)
writeTIFF(linear_12mm, "circles12mm_stereographic2linear.tif")
writeTIFF(stereographic_12mm, "circles12mm_stereographic_back.tif")



# Example 3: real spheres
img=readTIFF("spheres12mm.tif")
img=add_grid(img, n_gridx = 12, colour = c(0,1,1))
writeTIFF(img, "spheres12mm_grid.tif")

# Spherical with 16mm is almost equivalent to stereographic with 13.3mm
# but only Mercator (locally) and stereographic restore spheres as circles

cylindrical_12mm <- linear_to_cylindrical_rcpp(img, fl_FF_mm = 12)
equirectangular_12mm <- linear_to_equirectangular_rcpp(img, fl_FF_mm = 12)
mercator_12mm <- linear_to_mercator_rcpp(img, fl_FF_mm = 13.3)  # 13.3 restores spheres as circles better
panini_12mm <- linear_to_panini_rcpp(img, fl_FF_mm = 12)
orthographic_12mm <- linear_to_orthographic_rcpp(img, fl_FF_mm = 12)
spherical_12mm <- linear_to_spherical_rcpp(img, fl_FF_mm = 12)  # 16mm would restore spheres as circles
stereographic_12mm <- linear_to_stereographic_rcpp(img, fl_FF_mm = 13.3)  # 13.3 restores spheres as circles better

writeTIFF(cylindrical_12mm, "spheres12mm_cylindrical.tif")
writeTIFF(equirectangular_12mm, "spheres12mm_equirectangular.tif")
writeTIFF(mercator_12mm, "spheres12mm_mercator.tif")
writeTIFF(panini_12mm, "spheres12mm_panini.tif")
writeTIFF(orthographic_12mm, "spheres12mm_orthographic.tif")
writeTIFF(spherical_12mm, "spheres12mm_spherical.tif")
writeTIFF(stereographic_12mm, "spheres12mm_stereographic.tif")



# Example 4: Cervino 18mm (picture by Javier Camacho Gimeno)
img=readTIFF("cervino18mm.tif")
img=add_grid(img, n_gridx = 12, colour = c(1,1,0), linewidth = 2)
writeTIFF(img, "cervino18mm_grid.tif")

cylindrical_18mm <- linear_to_cylindrical_rcpp(img, fl_FF_mm = 18)
writeTIFF(cylindrical_18mm, "cervino18mm_cylindrical.tif")
