#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linear_to_cylindrical_rcpp(NumericVector img, double fl_FF_mm, double zoom = 1.0) {
    IntegerVector dims = img.attr("dim");
    if (dims.size() != 3) {
        stop("Input image must be a 3D array (M x N x 3).");
    }
    
    int H_in = dims[0]; 
    int W_in = dims[1]; 
    int channels = dims[2];
    
    double diag_mm = std::sqrt(36.0 * 36.0 + 24.0 * 24.0); 
    double diag_pixel = std::sqrt((double)H_in * H_in + (double)W_in * W_in);
    
    // Diferenciamos la focal de la imagen de entrada y la de salida aplicando el zoom
    double f_pixel_in = fl_FF_mm * (diag_pixel / diag_mm);
    double f_pixel_out = f_pixel_in * zoom;
    
    int H_out = H_in;
    int W_out = W_in;
    
    double cy_in = H_in / 2.0;
    double cx_in = W_in / 2.0;
    double cy_out = H_out / 2.0;
    double cx_out = W_out / 2.0;
    
    NumericVector out_img(Dimension(H_out, W_out, channels));
    
    // Pre-calculate row-dependent Y values to avoid recalculating in the outer loop
    // Utilizamos f_pixel_out para establecer las coordenadas geométricas de salida
    std::vector<double> Y_vals(H_out);
    for (int r_out = 0; r_out < H_out; ++r_out) {
        Y_vals[r_out] = (r_out - cy_out) / f_pixel_out;
    }
    
    // Raw pointers for fast, direct contiguous memory access (bypasses Rcpp operator overhead)
    double* p_out = out_img.begin();
    const double* p_in = img.begin();
    
    // Pre-calculate channel memory block strides
    int in_channel_stride = H_in * W_in;
    int out_channel_stride = H_out * W_out;
    
    // Loop columns first (Outer)
    for (int c_out = 0; c_out < W_out; ++c_out) {
        
        // --- Column Invariants Hoisted Here ---
        // Utilizamos f_pixel_out para determinar los ángulos visuales en la proyección cilíndrica
        double theta = (c_out - cx_out) / f_pixel_out;
        double X = std::sin(theta);
        double Z = std::cos(theta);
        
        // If Z <= 0, the ray points backwards. Entire column is out of bounds.
        if (Z <= 0) {
            for (int r_out = 0; r_out < H_out; ++r_out) {
                int out_idx = r_out + H_out * c_out;
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_idx + ch * out_channel_stride] = 0.0;
                }
            }
            continue; 
        }
        
        double inv_Z = 1.0 / Z;
        
        // Trazamos el rayo de vuelta hacia la imagen original plana usando f_pixel_in
        double c_in = f_pixel_in * (X * inv_Z) + cx_in;
        
        int c0 = std::floor(c_in);
        int c1 = c0 + 1;
        double delta_c = c_in - c0;
        double one_minus_delta_c = 1.0 - delta_c;
        
        // Check if column is completely out of bounds
        if (c0 < 0 || c1 >= W_in) {
            for (int r_out = 0; r_out < H_out; ++r_out) {
                int out_idx = r_out + H_out * c_out;
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_idx + ch * out_channel_stride] = 0.0;
                }
            }
            continue;
        }
        
        // Pre-calculate column memory strides for the input array
        int in_c0_stride = c0 * H_in;
        int in_c1_stride = c1 * H_in;
        int out_col_idx = H_out * c_out;
        
        // Loop rows (Inner)
        for (int r_out = 0; r_out < H_out; ++r_out) {
            
            // Compute vertical projection using the pre-calculated Y data
            // Trazamos el rayo vertical en la imagen de entrada con f_pixel_in
            double r_in = f_pixel_in * (Y_vals[r_out] * inv_Z) + cy_in;
            
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
            
            // Loop channels (Innermost - only 3 iterations, easy for compiler to unroll)
            for (int ch = 0; ch < channels; ++ch) {
                int ch_offset = ch * in_channel_stride;
                
                // Direct pointer offsets replace the lambda function
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
