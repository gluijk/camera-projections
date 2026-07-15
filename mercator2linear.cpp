#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mercator_to_linear_rcpp(NumericVector img, double fl_FF_mm) {

    IntegerVector dims = img.attr("dim");
    if (dims.size() != 3)
        stop("Input image must be a 3D array (M x N x 3).");

    int H_in = dims[0];
    int W_in = dims[1];
    int channels = dims[2];

    double diag_mm = std::sqrt(36.0 * 36.0 + 24.0 * 24.0);
    double diag_pixel = std::sqrt((double)H_in * H_in + (double)W_in * W_in);
    double f_pixel = fl_FF_mm * (diag_pixel / diag_mm);

    int H_out = H_in;
    int W_out = W_in;

    double cy_in  = H_in  / 2.0;
    double cx_in  = W_in  / 2.0;
    double cy_out = H_out / 2.0;
    double cx_out = W_out / 2.0;

    NumericVector out_img(Dimension(H_out, W_out, channels));

    double *p_out = out_img.begin();
    const double *p_in = img.begin();

    int in_channel_stride  = H_in * W_in;
    int out_channel_stride = H_out * W_out;

    // Precalcular coordenadas Y
    std::vector<double> y_out(H_out);
    std::vector<double> y_out_sq(H_out);

    for(int r=0;r<H_out;r++) {
        double y=r-cy_out;
        y_out[r]=y;
        y_out_sq[r]=y*y;
    }

    for(int c_out=0;c_out<W_out;c_out++) {

        double x = c_out-cx_out;
        double x2 = x*x;

        double theta = std::atan(x/f_pixel);

        double c_in = f_pixel*theta + cx_in;

        int c0 = std::floor(c_in);
        int c1 = c0+1;

        if(c0<0 || c1>=W_in) {

            for(int r=0;r<H_out;r++) {

                int out_idx=r+H_out*c_out;

                for(int ch=0;ch<channels;ch++)
                    p_out[out_idx+ch*out_channel_stride]=0.0;
            }

            continue;
        }

        double dc = c_in-c0;
        double omdc = 1.0-dc;

        int in_c0_stride=c0*H_in;
        int in_c1_stride=c1*H_in;
        int out_col_idx=H_out*c_out;

        for(int r_out=0;r_out<H_out;r_out++) {

            int out_pixel_idx=r_out+out_col_idx;

            double y=y_out[r_out];

            double denom = std::sqrt(x2+y_out_sq[r_out]+f_pixel*f_pixel);

            double sin_phi = y/denom;

            // Protección numérica
            if(sin_phi>=1.0) sin_phi=0.999999999;
            if(sin_phi<=-1.0) sin_phi=-0.999999999;

            double v = std::atanh(sin_phi);

            double r_in = f_pixel*v + cy_in;

            int r0=std::floor(r_in);
            int r1=r0+1;

            if(r0<0 || r1>=H_in) {

                for(int ch=0;ch<channels;ch++)
                    p_out[out_pixel_idx+ch*out_channel_stride]=0.0;

                continue;
            }

            double dr=r_in-r0;
            double omdr=1.0-dr;

            for(int ch=0;ch<channels;ch++) {

                int ch_offset=ch*in_channel_stride;

                double p00=p_in[r0+in_c0_stride+ch_offset];
                double p10=p_in[r0+in_c1_stride+ch_offset];
                double p01=p_in[r1+in_c0_stride+ch_offset];
                double p11=p_in[r1+in_c1_stride+ch_offset];

                double value=
                    omdr*(omdc*p00+dc*p10)+
                    dr  *(omdc*p01+dc*p11);

                p_out[out_pixel_idx+ch*out_channel_stride]=value;
            }
        }
    }

    return out_img;
}
