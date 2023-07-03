// nanocod.cpp

#include "WavFile.h"
#include "pvoclib.h"

#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;

int main (int argc, char* argv[]) {
    cout << "nanocod, ver. 0.1" << endl << endl;

    try {
        if (argc != 9) {
            throw runtime_error ("syntax: nancod input.wav output.wav nfft hop tstretch pshift denoise tau");
            return 0;
        }

        // read audio input file
        WavInFile input (argv[1]);
        int samps = input.getNumSamples ();
        int nchans = input.getNumChannels ();
        int nbits = input.getNumBits ();
        float sr = input.getSampleRate ();

        int nfft = atol (argv[3]);
        float ohop = atol (argv[4]);
        float oolap = nfft / ohop;

        float tstretch = atof (argv[5]);
        float hop = ohop / tstretch;
        float olap = nfft / hop;

        float pshift = atof (argv[6]);
        float threshold = atof (argv[7]);
        float tau = atof (argv[8]);
        float feedback = pow (10.0, (-3.0 * (float) ohop / (tau * sr)));

        // error check
        if ((((~nfft + 1) & nfft) != nfft) || nfft < 2) { // check pow 2
           throw runtime_error ("fft size must be a power of two");
        }
        if (hop <= 0) {
            throw runtime_error ("invalid hop size");
        }
        if (tstretch <= 0) {
            throw runtime_error ("invalid time strecth");
        }
        if (pshift <= 0) {
            throw runtime_error ("invalid pitch shift");
        }        
        if (feedback < 0 || feedback >= 1) {
            throw runtime_error ("invalid feedaback");
        }                

        cout << "sr      : " << sr << endl;
        cout << "samples : " << samps << endl;
        cout << "channels: " << nchans << endl << endl;
        
        cout << "stretch : " << tstretch << endl;
        cout << "shift   : " << pshift << endl;
        cout << "denoise : " << threshold << endl;
        cout << "tau     : " << tau << " (feedback: " << feedback << ")" << endl;
        cout << "nfft    : " << nfft << endl;
        cout << "in olap : " << olap << endl;
        cout << "out olap: " << oolap << endl << endl;

        float* window = new float[nfft];
        make_Window (window, nfft, .5, .5, 0);

        float* buffer = new float[samps * nchans + nfft];
        input.read (buffer, samps * nchans);

        int osamps = (int) ((float)  samps * nchans * tstretch) + nfft + (int) (tau * sr / nchans);
        float* obuffer = new float[osamps];
        memset (obuffer, 0, sizeof (float) * osamps);

        float* wksp = new float[2 * nfft];
        float* amp = new float[nfft];
        float* freq = new float[nfft];
        float* oamp = new float[nfft];
        float* ofreq = new float[nfft];
        float* phi = new float[nfft];
        float* ophi = new float[nfft];

        float* delay_amp = new float[nfft];
        float* delay_freq = new float[nfft];

        float imax = 0;
        for (int i = 0; i < samps * nchans; ++i) {
            if (buffer[i] > imax) imax = buffer[i];
        }

        // for each channel
        //      for each window of input file:
        //          - compute fft
        //          - process fft data
        //          - update phases
        //          - compute inverse fft
        cout << "processing..."; cout.flush ();
        for (int j = 0; j < nchans; ++j) {
            float pointer = 0;
            float opointer = 0;
            memset (phi, 0, sizeof (float) * nfft);
            memset (ophi, 0, sizeof (float) * nfft);
            memset (delay_amp, 0, sizeof (float) * nfft);
            memset (delay_freq, 0, sizeof (float) * nfft);

            while (pointer < osamps) {
                memset (wksp, 0, sizeof (float) * 2 * nfft);
                if (pointer < samps) {
                    for (int i = (int) pointer; i < (int) pointer + nfft; ++i) {
                        int rpos = (int) (nchans * i + j);
                        wksp[2 * (i - (int) pointer)] = buffer[rpos] * window[i - (int) pointer];
                        wksp[2 * (i - (int) pointer) + 1] = 0.;
                    }
                }

                fft (wksp, nfft, -1);
                convert (wksp, amp, freq, phi, nfft, hop, sr);
                
                // denoise
                for (int i = 0; i < nfft; ++i) {
                    if (amp[i] < threshold) amp[i] = 0;
                }

                memset (oamp, 0, sizeof (float) * nfft);
                memset (ofreq, 0, sizeof (float) * nfft);

                // pitch shift
                int hspect = nfft * pshift;
                for (int i = 0; i < nfft; ++i) {
                    int idx = (int) ((float) i / pshift);
                    if (idx < hspect) {
                        oamp[i] += amp[idx];
                        ofreq[i] = freq[idx] * pshift;
                    }
                }

                for (int i = 0; i < nfft; ++i) {
                    oamp[i] += delay_amp[i] * feedback;
                    ofreq[i] += delay_freq[i] * feedback;
                    delay_amp[i] = oamp[i];
                    delay_freq[i] = ofreq[i];
                }

                unconvert (oamp, ofreq, ophi, wksp, nfft, ohop, sr);
                fft (wksp, nfft, 1);

                for (int i = (int) opointer; i < (int) (opointer + nfft); ++i) {
                    int wpos = nchans * i + j;
                    if (wpos > osamps) break;
                    obuffer[wpos] += (wksp[2 * (i - (int) opointer)] * window[i - (int) (opointer)]);
                }

                pointer += hop;
                opointer += ohop;
            }
        }
        cout << "done" << endl;

        float omax = 0;
         for (int i = 0; i < samps * nchans; ++i) {
            if (obuffer[i] > omax) omax = obuffer[i];
        }

        // normalization
        for (int i = 0; i < osamps; ++i) {
            obuffer[i] *= (imax / omax);
        }

        // save audio output file
        WavOutFile output (argv[2], sr, nbits, nchans);
        output.write (obuffer, osamps);

        delete [] buffer;
        delete [] obuffer;
        delete [] wksp;
        delete [] window;
        delete [] amp;
        delete [] freq;
        delete [] oamp;
        delete [] ofreq;
        delete [] phi;
        delete [] ophi;

    } catch (exception& e) {
        cout << "error: " << e.what () << endl;
    } catch (...) {
        cout << "fatal error: unknown problem" << endl;
    }
    return 0;
}