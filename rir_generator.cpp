/*
Program     : StIM Room Impulse Response Generator
 
Description : Computes the response of an acoustic source to a microphone
              residing in an adjacent room using the Structuralimage method [1,2].
 
              [1] J.B. Allen and D.A. Berkley,
              Image method for efficiently simulating small-room acoustics,
              Journal Acoustic Society of America, 65(4), April 1979, p 943.
 
              [2] E.Shalev I.Cohen and D.Levov,
              ndoors audio classification with structure image methodfor simulating multi-room acoustics'
              The Journal of the Acoustical Society of America, 150(4):3059–3073, 2021.
*/

#define _USE_MATH_DEFINES
#include "mex.h"
#include "math.h"

#ifndef M_PI 
    #define M_PI 3.14159265358979323846 
#endif

#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

int box_ray(const double L[], double xv[], const double xf[]) {
 // Same as Step 4

 
    double ndx = xv[0] - xf[0];
    double ndy = xv[1] - xf[1];
    double ndz = xv[2] - xf[2];


    double sxy = ndx * L[1];
    double sxz = ndx * L[2];
    double syx = ndy * L[0];
    double syz = ndy * L[2];
    double szx = ndz * L[0];
    double szy = ndz * L[1];

    double cxy = xf[0] * xv[1] - xf[1] * xv[0];
    double cxz = xf[0] * xv[2] - xf[2] * xv[0];
    double cyz = xf[1] * xv[2] - xf[2] * xv[1];


    double axy = abs(ndx * ndy);
    double axz = abs(ndx * ndz);
    double ayz = abs(ndy * ndz);

    int face_num = 0;
    double face_tau = abs(ndz * axy);
    double tau;

    if (xv[0] < 0 && xf[0] > 0) {
        tau = -xv[0] * ayz;
        if (tau < face_tau && cxy >= 0 && cxz >= 0 && cxy <= -sxy && cxz <= -sxz) { 
         face_tau = tau;
         face_num = 1;
        }
    } else if (xf[0] < L[0] && xv[0] > L[0]) {
        tau = (xv[0] - L[0]) * ayz;
        if (tau < face_tau && cxy <= syx && cxz <= szx && cxy >= syx - sxy && cxz >= szx - sxz) {
         face_tau = tau;
         face_num = 2;
        }
    }
    if (xv[1] < 0 && xf[1] > 0) {
        tau = -xv[1] * axz;
        if (tau < face_tau && cxy <= 0 && cyz >= 0 && cxy >= syx && cyz <= -syz) {
         face_tau = tau;
         face_num = 3;
    }
    } else if (xv[1] > L[1] && xf[1] < L[1]) {
        tau = (xv[1] - L[1]) * axz;
        if (tau < face_tau && cxy >= -sxy && cyz <= szy && cxy <= syx - sxy && cyz >= szy - syz) {
         
         face_tau = tau;
         face_num = 4;
        }
    }
    if (xv[2] < 0 && xf[2] > 0) {
        tau = -xv[2] * axy;
        if (tau < face_tau && cxz <= 0 && cyz <= 0 && cxz >= szx && cyz >= szy) {
         face_tau = tau;
         face_num = 5;
        }
    } else if (xv[2] > L[2] && xf[2] < L[2]) {
        tau = (xv[2] - L[2]) * axy;
        if (tau < face_tau && cxz >= -sxz && cyz >= -syz && cxz <= szx - sxz && cyz <= szy - syz) {
         face_tau = tau;
         face_num = 6;
        }
    }

    return face_num;
}


double sinc(double x) {
    if (x == 0) return 1.0;
    return sin(x) / x;
}

double sim_microphone(double x, double y, double z, double* angle, char mtype) {
 // Same as Step 4
 if (mtype == 'b' || mtype == 'c' || mtype == 's' || mtype == 'h') {
     double gain, vartheta, varphi, rho;

     switch (mtype) {
     case 'b': rho = 0; break;
     case 'h': rho = 0.25; break;
     case 'c': rho = 0.5; break;
     case 's': rho = 0.75; break;
     }
  
  
  vartheta = acos(z / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
  varphi = atan2(y, x);

  gain = sin(M_PI / 2 - angle[1]) * sin(vartheta) * cos(angle[0] - varphi) +
               cos(M_PI / 2 - angle[1]) * cos(vartheta);
  gain = rho + (1 - rho) * gain;

  return gain;
 }
 return 1.0;
}


void ReciverIM(double c, double fs, const double* rr, int nMicrophones, const double* ss,
               const double* LL, const double* beta_input, int nSamples, char* microphone_type,
               double* angle, double* imp, int dim) {
    double* beta = new double[6];
    double reverberation_time = 0;

    if (mxGetN(beta_input) == 1) {
        double V = LL[0] * LL[1] * LL[2];
        double S = 2 * (LL[0] * LL[2] + LL[1] * LL[2] + LL[0] * LL[1]);
        reverberation_time = beta_input[0];
        if (reverberation_time != 0) {
            double alfa = 24 * V * log(10.0) / (c * S * reverberation_time);
            if (alfa > 1)
                mexErrMsgTxt("Error: Invalid reverberation time for room parameters.");
            for (int i = 0; i < 6; i++)
                beta[i] = sqrt(1 - alfa);
        } else {
            for (int i = 0; i < 6; i++)
                beta[i] = 0;
        }
    } else {
        for (int i = 0; i < 6; i++)
            beta[i] = beta_input[i];
    }
     if (dim == 2) {
        beta[4] = 0;
        beta[5] = 0;
    }


    double* r = new double[3];
    double* s = new double[3];
    double* L = new double[3];
    double* xp = new double[3];
    double Rm[3], Rp_plus_Rm[3], refl[3];
    double dist, fdist, gain;
    int n1, n2, n3, mx, my, mz, q, j, k, n, face_num;

    const double cTs = c / fs;
    const double Fc = 1.0;
    const int Tw = 2 * ROUND(0.004 * fs);
    double* LPI = new double[Tw];

    s[0] = ss[0] / cTs; s[1] = ss[1] / cTs; s[2] = ss[2] / cTs;
    L[0] = LL[0] / cTs; L[1] = LL[1] / cTs; L[2] = LL[2] / cTs;

    for (int idxMicrophone = 0; idxMicrophone < nMicrophones; idxMicrophone++) {
        r[0] = rr[idxMicrophone] / cTs;
        r[1] = rr[idxMicrophone + nMicrophones] / cTs;
        r[2] = rr[idxMicrophone + 2 * nMicrophones] / cTs;

        n1 = (int)ceil(nSamples / (2 * L[0]));
        n2 = (int)ceil(nSamples / (2 * L[1]));
        n3 = dim == 2 ? 0 : (int)ceil(nSamples / (2 * L[2]));

        for (mx = -n1; mx <= n1; mx++) {
            Rm[0] = 2 * mx * L[0];
            for (my = -n2; my <= n2; my++) {
                Rm[1] = 2 * my * L[1];
                for (mz = dim == 2 ? 0 : -n3; mz <= (dim == 2 ? 0 : n3); mz++) {
                    Rm[2] = 2 * mz * L[2];
                    for (q = 0; q <= 1; q++) {
                        Rp_plus_Rm[0] = (1 - 2 * q) * r[0] - s[0] + Rm[0];
                        xp[0] = 2 * mx * LL[0] + (1 - 2 * q) * rr[idxMicrophone];
                        refl[0] = pow(beta[0], abs(mx - q)) * pow(beta[1], abs(mx));
                        for (j = 0; j <= 1; j++) {
                            Rp_plus_Rm[1] = (1 - 2 * j) * r[1] - s[1] + Rm[1];
                            xp[1] = 2 * my * LL[1] + (1 - 2 * j) * rr[idxMicrophone + nMicrophones];
                            refl[1] = pow(beta[2], abs(my - j)) * pow(beta[3], abs(my));
                            for (k = dim == 2 ? 0 : 0; k <= (dim == 2 ? 0 : 1); k++) {
                                Rp_plus_Rm[2] = (1 - 2 * k) * r[2] - s[2] + Rm[2];
                                xp[2] = 2 * mz * LL[2] + (1 - 2 * k) * rr[idxMicrophone + 2 * nMicrophones];
                                refl[2] = pow(beta[4], abs(mz - k)) * pow(beta[5], abs(mz));

                                dist = sqrt(pow(Rp_plus_Rm[0], 2) + pow(Rp_plus_Rm[1], 2) + pow(Rp_plus_Rm[2], 2));
                                face_num = box_ray(LL, xp, ss);

                                if (face_num == 0) continue;

                                fdist = floor(dist);
                                if (fdist < nSamples) {
                                    gain = sim_microphone(Rp_plus_Rm[0], Rp_plus_Rm[1], Rp_plus_Rm[2], angle, microphone_type[0])
                                         * refl[0] * refl[1] * refl[2] / (4 * M_PI * dist * cTs);
                                    gain = gain * (1 - beta[face_num - 1]) / beta[face_num - 1];

                                    for (n = 0; n < Tw; n++)
                                        LPI[n] = 0.5 * (1 - cos(2 * M_PI * ((n + 1 - (dist - fdist)) / Tw))) * Fc * sinc(M_PI * Fc * (n + 1 - (dist - fdist) - (Tw / 2)));

                                    int startPosition = (int)fdist - (Tw / 2) + 1;
                                    for (n = 0; n < Tw; n++) {
                                        if (startPosition + n >= 0 && startPosition + n < nSamples)
                                            imp[idxMicrophone + nMicrophones * (startPosition + n)] += gain * LPI[n];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] LPI;
    delete[] r;
    delete[] s;
    delete[] L;
    delete[] xp;
    delete[] beta;

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 0) {
        mexPrintf("--------------------------------------------------------------------\n"
                  "| Room Impulse Response Generator                                  |\n"
                  "| Computes the response of an acoustic source to one or more       |\n"
                  "| microphones in a reverberant room using the image method [1,2].  |\n"
                  "| [1] J.B. Allen and D.A. Berkley, J. Acoust. Soc. Am., 1979.      |\n"
                  "| [2] E. Shalev, I. Cohen, D. Levov, J. Acoust. Soc. Am., 2021.    |\n"
                  "--------------------------------------------------------------------\n\n"
                  "function h = rir_generator(c, fs, r, s, Lr, Ls, beta_r, beta_s, nsample, mtype, orientation, dim);\n\n"
                  "Input parameters:\n"
                  " c        : sound velocity in m/s.\n"
                  " fs       : sampling frequency in Hz.\n"
                  " r        : M x 3 array of receiver coordinates (x,y,z) in m.\n"
                  " s        : 1 x 3 vector of source coordinates (x,y,z) in m.\n"
                  " Lr       : 1 x 3 vector of receiver room dimensions (x,y,z) in m.\n"
                  " Ls       : 1 x 3 vector of source room dimensions (x,y,z) in m.\n"
                  " beta_r   : 1 x 6 vector of reflection coefficients or single T60 value (s).\n"
                  " beta_s   : 1 x 6 vector of reflection coefficients or single T60 value (s).\n"
                  " nsample  : number of samples to calculate.\n\n"
                  " mtype       : [omnidirectional, subcardioid, cardioid, hypercardioid, bidirectional], default is omnidirectional.\n"
                  " orientation : microphone direction [azimuth, elevation] in radians, default is [0 0].\n\n"
                  "Output parameters:\n"
                  " h        : M x nsample matrix of impulse responses.\n\n");
        return;
    }
    if (nrhs < 9) mexErrMsgTxt("At least 9 inputs required.");
    if (nrhs > 12) mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");

    if (!(mxGetN(prhs[0]) == 1) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Invalid input c arguments!");
    if (!(mxGetN(prhs[1]) == 1) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Invalid input fs arguments!");
    if (!(mxGetN(prhs[2]) == 3) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
        mexErrMsgTxt("Invalid input r arguments!");
    if (!(mxGetN(prhs[3]) == 3) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
        mexErrMsgTxt("Invalid input s arguments!");
    if (!(mxGetN(prhs[4]) == 3) || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
        mexErrMsgTxt("Invalid input Lr arguments!");
    if (!(mxGetN(prhs[5]) == 3) || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]))
        mexErrMsgTxt("Invalid input Ls arguments!");
    if (!(mxGetN(prhs[6]) == 1 && mxGetM(prhs[6]) == 1 || mxGetN(prhs[6]) == 6) || !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]))
        mexErrMsgTxt("Invalid input beta_r arguments!");
    if (!(mxGetN(prhs[7]) == 1 && mxGetM(prhs[7]) == 1 || mxGetN(prhs[7]) == 6) || !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]))
        mexErrMsgTxt("Invalid input beta_s arguments!");
    if (!(mxGetN(prhs[8]) == 1) || !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]))
        mexErrMsgTxt("Invalid input nsample arguments!");

    double c = mxGetScalar(prhs[0]);
    double fs = mxGetScalar(prhs[1]);
    const double* rr = mxGetPr(prhs[2]);
    int nMicrophones = (int)mxGetM(prhs[2]);
    const double* ss = mxGetPr(prhs[3]);
    const double* LL_r = mxGetPr(prhs[4]);
    const double* LL_s = mxGetPr(prhs[5]);
    const double* beta_r = mxGetPr(prhs[6]);
    const double* beta_s = mxGetPr(prhs[7]);
    int nSamples = (int)mxGetScalar(prhs[8]);

    char* microphone_type;
    if (nrhs > 9 && !mxIsEmpty(prhs[9])) {
     microphone_type = new char[mxGetN(prhs[9]) + 1];
        mxGetString(prhs[9], microphone_type, mxGetN(prhs[9]) + 1);
    } else {
        microphone_type = new char[1];
        microphone_type[0] = 'o';
    }

    double angle[2];
    if (nrhs > 10 && !mxIsEmpty(prhs[10])) {
        const double* orientation = mxGetPr(prhs[10]);
        if (mxGetN(prhs[10]) == 1) {
            angle[0] = orientation[0];
            angle[1] = 0;
        } else {
            angle[0] = orientation[0];
            angle[1] = orientation[1];
        }
    } else {
        angle[0] = 0;
        angle[1] = 0;
    }
     int dim = 3;
    if (nrhs > 11 && !mxIsEmpty(prhs[11])) {
        dim = (int)mxGetScalar(prhs[11]);
        if (dim != 2 && dim != 3)
            mexErrMsgTxt("Invalid input dim: must be 2 or 3.");
    }

    plhs[0] = mxCreateDoubleMatrix(nMicrophones, nSamples, mxREAL);
    double* imp = mxGetPr(plhs[0]);

    double* r = new double[3];
    double* s = new double[3];
    double* L = new double[3];
    double* xp = new double[3];
    double Rm[3], Rp_plus_Rm[3], refl[3];
    double dist, fdist, gain;
    int n1, n2, n3, mx, my, mz, q, j, k, n, face_num;

    const double cTs = c / fs;
    const double Fc = 1.0;
    const int Tw = 2 * ROUND(0.004 * fs);
    double* LPI = new double[Tw];

    s[0] = ss[0] / cTs; s[1] = ss[1] / cTs; s[2] = ss[2] / cTs;
    L[0] = LL_s[0] / cTs; L[1] = LL_s[1] / cTs; L[2] = LL_s[2] / cTs;

    for (int idxMicrophone = 0; idxMicrophone < nMicrophones; idxMicrophone++) {
        r[0] = rr[idxMicrophone] / cTs;
        r[1] = rr[idxMicrophone + nMicrophones] / cTs;
        r[2] = rr[idxMicrophone + 2 * nMicrophones] / cTs;

        n1 = (int)ceil(nSamples / (2 * L[0]));
        n2 = (int)ceil(nSamples / (2 * L[1]));
        n3 = dim == 2 ? 0 : (int)ceil(nSamples / (2 * L[2]));

        for (mx = -n1; mx <= n1; mx++) {
            Rm[0] = 2 * mx * L[0];
            for (my = -n2; my <= n2; my++) {
                Rm[1] = 2 * my * L[1];
                for (mz = dim == 2 ? 0 : -n3; mz <= (dim == 2 ? 0 : n3); mz++) {
                    Rm[2] = 2 * mz * L[2];
                    for (q = 0; q <= 1; q++) {
                        Rp_plus_Rm[0] = (1 - 2 * q) * s[0] - r[0] + Rm[0];
                        xp[0] = 2 * mx * LL_s[0] + (1 - 2 * q) * ss[0];
                        refl[0] = pow(beta_s[0], abs(mx - q)) * pow(beta_s[1], abs(mx));
                        for (j = 0; j <= 1; j++) {
                            Rp_plus_Rm[1] = (1 - 2 * j) * s[1] - r[1] + Rm[1];
                            xp[1] = 2 * my * LL_s[1] + (1 - 2 * j) * ss[1];
                            refl[1] = pow(beta_s[2], abs(my - j)) * pow(beta_s[3], abs(my));
                            for (k = dim == 2 ? 0 : 0; k <= (dim == 2 ? 0 : 1); k++) {
                                Rp_plus_Rm[2] = (1 - 2 * k) * s[2] - r[2] + Rm[2];
                                xp[2] = 2 * mz * LL_s[2] + (1 - 2 * k) * ss[2];
                                refl[2] = pow(beta_s[4], abs(mz - k)) * pow(beta_s[5], abs(mz));

                                dist = sqrt(pow(Rp_plus_Rm[0], 2) + pow(Rp_plus_Rm[1], 2) + pow(Rp_plus_Rm[2], 2));
                                face_num = box_ray(LL_s, xp, rr);

                                if (face_num == 0) continue;

                                ReciverIM(c, fs, rr, nMicrophones, ss, LL_r, beta_r, nSamples, microphone_type, angle, imp);

                                fdist = floor(dist);
                                if (fdist < nSamples) {
                                    gain = sim_microphone(Rp_plus_Rm[0], Rp_plus_Rm[1], Rp_plus_Rm[2], angle, microphone_type[0])
                                         * refl[0] * refl[1] * refl[2] / (4 * M_PI * dist * cTs);

                                    for (n = 0; n < Tw; n++)
                                        LPI[n] = 0.5 * (1 - cos(2 * M_PI * ((n + 1 - (dist - fdist)) / Tw))) * Fc * sinc(M_PI * Fc * (n + 1 - (dist - fdist) - (Tw / 2)));

                                    int startPosition = (int)fdist - (Tw / 2) + 1;
                                    for (n = 0; n < Tw; n++) {
                                        if (startPosition + n >= 0 && startPosition + n < nSamples)
                                            imp[idxMicrophone + nMicrophones * (startPosition + n)] += gain * LPI[n];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] LPI;
    delete[] r;
    delete[] s;
    delete[] L;
    delete[] xp;
    delete[] microphone_type;
}
