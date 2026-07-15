#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linear_to_stereographic_rcpp(NumericVector img, double fl_FF_mm) {
	// Mercator/Stereographic son proyecciones conformes = tienen la propiedad de conservación de ángulos locales respecto a la esfera 3D
	// En una proyección conforme, cualquier forma infinitesimal (como una pequeña esfera o un círculo en la superficie terrestre)
	// se mantiene como una forma similar en el mapa
    
    IntegerVector dims = img.attr("dim");
    if (dims.size() != 3) {
        stop("Input image must be a 3D array (M x N x 3).");
    }
    
    int H_in = dims[0]; 
    int W_in = dims[1]; 
    int channels = dims[2];
    
    double diag_mm = std::sqrt(36.0 * 36.0 + 24.0 * 24.0); 
    double diag_pixel = std::sqrt((double)H_in * H_in + (double)W_in * W_in);
    double f_pixel = fl_FF_mm * (diag_pixel / diag_mm);
    
    int H_out = H_in;
    int W_out = W_in;
    
    double cy_in = H_in / 2.0;
    double cx_in = W_in / 2.0;
    double cy_out = H_out / 2.0;
    double cx_out = W_out / 2.0;
    
    NumericVector out_img(Dimension(H_out, W_out, channels));
    
	// Precalculamos las distancias Y al cuadrado para acelerar la hipotenusa en el bucle interior
    std::vector<double> y_out_vals(H_out);
    std::vector<double> y_out_sq(H_out);
    for (int r_out = 0; r_out < H_out; ++r_out) {
        double y = r_out - cy_out;
        y_out_vals[r_out] = y;
        y_out_sq[r_out] = y * y;
    }
    
    double* p_out = out_img.begin();
    const double* p_in = img.begin();
    
    int in_channel_stride = H_in * W_in;
    int out_channel_stride = H_out * W_out;
    
	// Bucle de columnas
    for (int c_out = 0; c_out < W_out; ++c_out) {
        
		// --- Invariantes de la columna ---
        double x_out = c_out - cx_out;
        double x_out_sq = x_out * x_out;
        int out_col_idx = H_out * c_out;
        
		// Bucle de filas
        for (int r_out = 0; r_out < H_out; ++r_out) {
            int out_pixel_idx = r_out + out_col_idx;
            double y_out = y_out_vals[r_out];
            
            // Distancia radial desde el centro en la imagen Estereográfica (salida)
            double R_out = std::sqrt(x_out_sq + y_out_sq[r_out]);
            
            // Ángulo incidente desde el eje óptico inverso usando el modelo Estereográfico:
                // R = 2 * f * tan(theta / 2) -> theta = 2 * atan(R / (2 * f))
            double theta = 2.0 * std::atan(R_out / (2.0 * f_pixel));
            
            // Límite de campo lineal (no podemos buscar datos más allá de 90 grados)
            if (theta >= 1.5707963267948966) {
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_pixel_idx + ch * out_channel_stride] = 0.0;
                }
                continue;
            }
            
            // En una imagen de cámara oscura rectilínea, R_in responde a f * tan(theta)
            double R_in = f_pixel * std::tan(theta);
            
            // Escala radial
            double scale = (R_out > 1e-6) ? (R_in / R_out) : 1.0;
            
            // Deformación radial acoplada
            double c_in = cx_in + x_out * scale;
            double r_in = cy_in + y_out * scale;
            
            int c0 = std::floor(c_in);
            int c1 = c0 + 1;
            int r0 = std::floor(r_in);
            int r1 = r0 + 1;
            
            if (c0 < 0 || c1 >= W_in || r0 < 0 || r1 >= H_in) {
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_pixel_idx + ch * out_channel_stride] = 0.0;
                }
                continue;
            }
            
            double delta_c = c_in - c0;
            double one_minus_delta_c = 1.0 - delta_c;
            double delta_r = r_in - r0;
            double one_minus_delta_r = 1.0 - delta_r;
            
            int in_c0_stride = c0 * H_in;
            int in_c1_stride = c1 * H_in;
            
            for (int ch = 0; ch < channels; ++ch) {
                int ch_offset = ch * in_channel_stride;
                
                double p00 = p_in[r0 + in_c0_stride + ch_offset];
                double p10 = p_in[r0 + in_c1_stride + ch_offset];
                double p01 = p_in[r1 + in_c0_stride + ch_offset];
                double p11 = p_in[r1 + in_c1_stride + ch_offset];
                
                double interpolated_val = one_minus_delta_r * (one_minus_delta_c * p00 + delta_c * p10) +
                    delta_r           * (one_minus_delta_c * p01 + delta_c * p11);
                
                p_out[out_pixel_idx + ch * out_channel_stride] = interpolated_val;
            }
        }
    }
    
    return out_img;
}