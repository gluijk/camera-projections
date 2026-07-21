#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linear_to_panini_rcpp(NumericVector img, double fl_FF_mm, double d = 1.0, double zoom = 1.0) {
	// Parámetro d: distancia al centro de proyección (normalmente entre 0.0 y 1.0)
	// 		d = 0 : proyección lineal/rectilínea
	// 		d = 1 : proyección estereográfica cilíndrica o Panini estándar
	// 		d >> 1: proyección cilíndrica ortográfica
	// Esta proyección cilíndrica/panorámica mantiene rectas las líneas verticales y las diagonales que pasan por el centro
    IntegerVector dims = img.attr("dim");
    if (dims.size() != 3) {
        stop("Input image must be a 3D array (M x N x 3).");
    }
    
    int H_in = dims[0]; 
    int W_in = dims[1]; 
    int channels = dims[2];
    
    double diag_mm = std::sqrt(36.0 * 36.0 + 24.0 * 24.0); 
    double diag_pixel = std::sqrt((double)H_in * H_in + (double)W_in * W_in);
    
    // Diferenciamos la focal de la imagen original plana (in) de la focal proyectada con zoom (out)
    // Utilizamos f_pixel_out para establecer las coordenadas geométricas de salida
    double f_pixel_in = fl_FF_mm * (diag_pixel / diag_mm);
    double f_pixel_out = f_pixel_in * zoom;
    
    int H_out = H_in;
    int W_out = W_in;
    
    double cy_in = H_in / 2.0;
    double cx_in = W_in / 2.0;
    double cy_out = H_out / 2.0;
    double cx_out = W_out / 2.0;
    
    NumericVector out_img(Dimension(H_out, W_out, channels));
    
    std::vector<double> Y_vals(H_out);
    for (int r_out = 0; r_out < H_out; ++r_out) {
        Y_vals[r_out] = (r_out - cy_out) / f_pixel_out;
    }
    
    double* p_out = out_img.begin();
    const double* p_in = img.begin();
    
    int in_channel_stride = H_in * W_in;
    int out_channel_stride = H_out * W_out;
    
    // Constantes de la proyección Panini pre-calculadas
    double d_plus_1 = d + 1.0;
    double d_plus_1_sq = d_plus_1 * d_plus_1;
    double d_sq_minus_1 = (d * d) - 1.0;
    
    for (int c_out = 0; c_out < W_out; ++c_out) {
        
        // Coordenada horizontal normalizada (u)
        double u = (c_out - cx_out) / f_pixel_out;
        double u_sq = u * u;
        
        // Discriminante de la ecuación cuadrática
        double inner = d_plus_1_sq - u_sq * d_sq_minus_1;
        
        if (inner < 0.0) {
            for (int r_out = 0; r_out < H_out; ++r_out) {
                int out_idx = r_out + H_out * c_out;
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_idx + ch * out_channel_stride] = 0.0;
                }
            }
            continue;
        }
        
        // Calculamos C = cos(phi)
        double C = (-u_sq * d + d_plus_1 * std::sqrt(inner)) / (u_sq + d_plus_1_sq);
        
        // Si C <= 0, el ángulo es mayor a 90 grados (el rayo apunta hacia atrás)
        if (C <= 0.0) {
            for (int r_out = 0; r_out < H_out; ++r_out) {
                int out_idx = r_out + H_out * c_out;
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_idx + ch * out_channel_stride] = 0.0;
                }
            }
            continue; 
        }
        
        // Calculamos S = sin(phi)
        double S = u * (d + C) / d_plus_1;
        
        // Retorno al plano rectilíneo: X = tan(phi) = sin(phi)/cos(phi)
        double X_in = S / C;
        double c_in = f_pixel_in * X_in + cx_in;
        
        int c0 = std::floor(c_in);
        int c1 = c0 + 1;
        double delta_c = c_in - c0;
        double one_minus_delta_c = 1.0 - delta_c;
        
        if (c0 < 0 || c1 >= W_in) {
            for (int r_out = 0; r_out < H_out; ++r_out) {
                int out_idx = r_out + H_out * c_out;
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_idx + ch * out_channel_stride] = 0.0;
                }
            }
            continue;
        }
        
        // Factor Y exacto para esta columna: contrarresta el estiramiento rectilíneo
        double Y_col_factor = f_pixel_in * (d + C) / (d_plus_1 * C); 
        
        int in_c0_stride = c0 * H_in;
        int in_c1_stride = c1 * H_in;
        int out_col_idx = H_out * c_out;
        
        for (int r_out = 0; r_out < H_out; ++r_out) {
            
            double r_in = Y_col_factor * Y_vals[r_out] + cy_in;
            
            int r0 = std::floor(r_in);
            int r1 = r0 + 1;
            
            int out_pixel_idx = r_out + out_col_idx;
            
            if (r0 < 0 || r1 >= H_in) {
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_pixel_idx + ch * out_channel_stride] = 0.0;
                }
                continue;
            }
            
            double delta_r = r_in - r0;
            double one_minus_delta_r = 1.0 - delta_r;
            
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
