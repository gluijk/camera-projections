#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linear_to_orthographic_rcpp(NumericVector img, double fl_FF_mm, double zoom = 1.0) {
    IntegerVector dims = img.attr("dim");
    if (dims.size() != 3) {
        stop("Input image must be a 3D array (M x N x 3).");
    }
    
    int H_in = dims[0]; 
    int W_in = dims[1]; 
    int channels = dims[2];
    
    double diag_mm = std::sqrt(36.0 * 36.0 + 24.0 * 24.0); 
    double diag_pixel = std::sqrt((double)H_in * H_in + (double)W_in * W_in);
    
    // The base focal length for sampling the input rectilinear image
    double f_pixel = fl_FF_mm * (diag_pixel / diag_mm);
    
    // The scaled focal length for the output orthographic mapping
    double f_pixel_zoom = f_pixel * zoom;
    
    int H_out = H_in;
    int W_out = W_in;
    
    double cy_in = H_in / 2.0;
    double cx_in = W_in / 2.0;
    double cy_out = H_out / 2.0;
    double cx_out = W_out / 2.0;
    
    NumericVector out_img(Dimension(H_out, W_out, channels));
    
    // Pre-calculate row-dependent Y values scaled by the zoomed focal length
    std::vector<double> Y_vals(H_out);
    for (int r_out = 0; r_out < H_out; ++r_out) {
        Y_vals[r_out] = (r_out - cy_out) / f_pixel_zoom;
    }
    
    // Raw pointers for fast, direct contiguous memory access
    double* p_out = out_img.begin();
    const double* p_in = img.begin();
    
    // Pre-calculate channel memory block strides
    int in_channel_stride = H_in * W_in;
    int out_channel_stride = H_out * W_out;
    
    // Loop columns first (Outer)
    for (int c_out = 0; c_out < W_out; ++c_out) {
        
        // Compute normalized X on the orthographic plane
        double X = (c_out - cx_out) / f_pixel_zoom;
        double X_sq = X * X;
        
        // If X^2 >= 1.0, the ray misses the unit sphere entirely. 
        // The entire column is out of bounds.
        if (X_sq >= 1.0) {
            for (int r_out = 0; r_out < H_out; ++r_out) {
                int out_idx = r_out + H_out * c_out;
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_idx + ch * out_channel_stride] = 0.0;
                }
            }
            continue; 
        }
        
        int out_col_idx = H_out * c_out;
        
        // Loop rows (Inner)
        for (int r_out = 0; r_out < H_out; ++r_out) {
            
            int out_pixel_idx = r_out + out_col_idx;
            
            double Y = Y_vals[r_out];
            double R_sq = X_sq + Y * Y;
            
            // Check if ray is outside the orthographic hemisphere.
            // Using 0.999999 to protect against divide-by-zero when Z approaches 0.
            if (R_sq >= 0.999999) {
                for (int ch = 0; ch < channels; ++ch) {
                    p_out[out_pixel_idx + ch * out_channel_stride] = 0.0;
                }
                continue;
            }
            
            // Reconstruct Z on the unit sphere
            double Z = std::sqrt(1.0 - R_sq);
            double inv_Z = 1.0 / Z;
            
            // Compute intersection on the flat pinhole image plane using original f_pixel
            double c_in = f_pixel * (X * inv_Z) + cx_in;
            double r_in = f_pixel * (Y * inv_Z) + cy_in;
            
            int c0 = std::floor(c_in);
            int c1 = c0 + 1;
            int r0 = std::floor(r_in);
            int r1 = r0 + 1;
            
            // Check if the mapped pixel is completely out of the input image bounds
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
            
            // Pre-calculate column memory strides for the input array
            int in_c0_stride = c0 * H_in;
            int in_c1_stride = c1 * H_in;
            
            // Loop channels (Innermost - easily unrolled)
            for (int ch = 0; ch < channels; ++ch) {
                int ch_offset = ch * in_channel_stride;
                
                // Direct pointer offsets
                double p00 = p_in[r0 + in_c0_stride + ch_offset];
                double p10 = p_in[r0 + in_c1_stride + ch_offset];
                double p01 = p_in[r1 + in_c0_stride + ch_offset];
                double p11 = p_in[r1 + in_c1_stride + ch_offset];
                
                // Bilinear interpolation
                double interpolated_val = one_minus_delta_r * (one_minus_delta_c * p00 + delta_c * p10) +
                                          delta_r           * (one_minus_delta_c * p01 + delta_c * p11);
                
                p_out[out_pixel_idx + ch * out_channel_stride] = interpolated_val;
            }
        }
    }
    
    return out_img;
}
