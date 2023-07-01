// PVOCLIB_H

#ifndef PVOCLIB_H
#define PVOCLIB_H

#include <cmath>

const double TWOPI = M_PI * 2.;

void make_Window (float* out, int N, float a0, float a1, float a2) {
    // .5, .5, 0     --> hanning
    // .54, .46, 0   --> hamming
    // .42, .5, 0.08 --> blackman

	for (int i = 0; i < N; ++i) {
        out[i] = a0 - a1 * cos ((TWOPI * (float) i) / (N - 1)) + a2 *
        	cos ((2 * TWOPI * (float) i) / (N - 1)); // hamming, hann or blackman
    }
}

void fft (float *fftBuffer, long fftFrameSize, long sign) {
    float wr, wi, arg, *p1, *p2, temp;
    float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
    long i, bitm, j, le, le2, k;

    for (i = 2; i < 2*fftFrameSize-2; i += 2) {
        for (bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {
            if (i & bitm) j++;
            j <<= 1;
        }
        if (i < j) {
            p1 = fftBuffer+i;
            p2 = fftBuffer+j;
            temp = *p1;
            *(p1++) = *p2;
            *(p2++) = temp;
            temp = *p1;
            *p1 = *p2;
            *p2 = temp;
        }
    }
    for (k = 0, le = 2; k < (long)(log(fftFrameSize)/log(2.)); k++) {
        le <<= 1;
        le2 = le>>1;
        ur = 1.0;
        ui = 0.0;
        arg = M_PI / (le2>>1);
        wr = cos(arg);
        wi = sign*sin(arg);
        for (j = 0; j < le2; j += 2) {
            p1r = fftBuffer+j;
            p1i = p1r+1;
            p2r = p1r+le2;
            p2i = p2r+1;
            for (i = j; i < 2*fftFrameSize; i += le) {
                tr = *p2r * ur - *p2i * ui;
                ti = *p2r * ui + *p2i * ur;
                *p2r = *p1r - tr;
                *p2i = *p1i - ti;
                *p1r += tr;
                *p1i += ti;
                p1r += le;
                p1i += le;
                p2r += le;
                p2i += le;
            }
            tr = ur*wr - ui*wi;
            ui = ur*wi + ui*wr;
            ur = tr;
        }
    }
}

// rect to amp-freq
void convert (const float* cbuffer, float* amp, float* freq, float* oldPhi,
              int N, int hop, float R) {
    float osamp = N / hop;
    float freqPerBin = R / (float) N;
    float expct = TWOPI * (float) hop / (float) N;
    for (int i = 0; i < N; ++i) {
        amp[i] = 2. * sqrt (cbuffer[2 * i] * cbuffer[2 * i]
                            + cbuffer[2 * i + 1] * cbuffer[2 * i + 1]);

        float phase = atan2 (cbuffer[2 * i + 1], cbuffer[2 * i]);
        float tmp = phase - oldPhi[i];
        oldPhi[i] = phase;
        tmp -= (float) i * expct;
        int qpd = (int) (tmp / M_PI);
        if (qpd >= 0) qpd += qpd & 1;
        else qpd -= qpd & 1;
        tmp -= M_PI * (float) qpd;
        tmp = osamp * tmp / (2. * M_PI);
        tmp = (float) i * freqPerBin + tmp * freqPerBin;
        freq[i] = tmp;
    }
}

// amp-freq to rect
void unconvert (const float* amp, const float* freq, float* oldPhi,
                float* cbuffer, int N, int hop, float R) {
    float osamp = N / hop;
    float freqPerBin = R / (float) N;
    float expct = TWOPI * (float) hop / (float) N;

    for (int i = 0; i < N; ++i) {
        float tmp = freq[i];

        tmp -= (float)i * freqPerBin;
        tmp /= freqPerBin;
        tmp = TWOPI * tmp / osamp;
        tmp += (float) i * expct;
        oldPhi[i] += tmp;

        float phase = (oldPhi[i]); // - oldPhi[i - 1] - oldPhi[i + 1]);

        cbuffer[2 * i] = amp[i] * cos (phase);
        cbuffer[2 * i + 1] = amp[i] * sin (phase);
    }
}

#endif // PVOCLIB_H
