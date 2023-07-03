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
        if (argc != 6) {
            throw runtime_error ("syntax: nancod input.wav output.wav tstretch pshift denoise");
            return 0;
        }

        // read audio input file
        WavInFile input (argv[1]);
        int samps = input.getNumSamples ();
        int nchans = input.getNumChannels ();
        int nbits = input.getNumBits ();
        float sr = input.getSampleRate ();

        int nfft = 4096;
        float ohop = 256;
        float oolap = nfft / ohop;

        float tstretch = atof (argv[3]);
        float hop = ohop / tstretch;
        float olap = nfft / hop;

        float pshift = atof (argv[4]);
        float threshold = atof (argv[5]);

        cout << "sr      : " << sr << endl;
        cout << "samples : " << samps << endl;
        cout << "channels: " << nchans << endl << endl;
        
        cout << "stretch : " << tstretch << endl;
        cout << "shift   : " << pshift << endl;
        cout << "denoise : " << threshold << endl;
        cout << "nfft    : " << nfft << endl;
        cout << "in olap : " << olap << endl;
        cout << "out olap: " << oolap << endl << endl;

        float* window = new float[nfft];
        make_Window (window, nfft, .5, .5, 0);

        float* buffer = new float[samps * nchans + nfft];
        input.read (buffer, samps * nchans);

        float* obuffer = new float[(int) ((float) samps * nchans * tstretch) + nfft];
        int v = (int) ((float) samps * nchans * tstretch) + nfft;
        cout << "v = "  << v << endl;
        memset (obuffer, 0, sizeof (float) * v);

        float* wksp = new float[2 * nfft];
        float* amp = new float[nfft];
        float* freq = new float[nfft];
        float* oamp = new float[nfft];
        float* ofreq = new float[nfft];
        float* phi = new float[nfft];
        float* ophi = new float[nfft];

        float norm = 0;
        for (int i = 0; i < samps * nchans; ++i) {
            if (buffer[i] > norm) norm = buffer[i];
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

            while (pointer < samps) {
                for (int i = (int) pointer; i < (int) pointer + nfft; ++i) {
                    int rpos = (int) (nchans * i + j);
                    wksp[2 * (i - (int) pointer)] = buffer[rpos] * window[i - (int) pointer];
                    wksp[2 * (i - (int) pointer) + 1] = 0.;
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

                unconvert (oamp, ofreq, ophi, wksp, nfft, ohop, sr);
                fft (wksp, nfft, 1);

                for (int i = (int) opointer; i < (int) (opointer + nfft); ++i) {
                    int wpos = nchans * i + j;
                    if (wpos > v) break;
                    obuffer[wpos] += (wksp[2 * (i - (int) opointer)] * window[i - (int) (opointer)]);
                }

                pointer += hop;
                opointer += ohop;
            }
        }
        cout << "done" << endl;

        float onorm = 0;
         for (int i = 0; i < samps * nchans; ++i) {
            if (obuffer[i] > onorm) onorm = obuffer[i];
        }

        // normalization
        for (int i = 0; i < (int) (samps * nchans * tstretch); ++i) {
            obuffer[i] *= (norm / onorm);
        }

        // save audio output file
        WavOutFile output (argv[2], sr, nbits, nchans);
        output.write (obuffer, (int) (samps * nchans * tstretch));
        delete [] buffer;
        delete [] obuffer;
        delete [] wksp;
        delete [] window;
        
    } catch (exception& e) {
        cout << "error: " << e.what () << endl;
    } catch (...) {
        cout << "fatal error: unknown problem" << endl;
    }
    return 0;
}