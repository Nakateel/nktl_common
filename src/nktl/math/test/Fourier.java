package nktl.math.test;

import nktl.math.Cnum;

/**
 * Created by Zheka Grushevskiy, NAKATEEL, 26.09.2017.
 */
public class Fourier {

    static double[][] test(Cnum[] signal, double sig_shift, double spec_shift){
        int N = signal.length;
        double[] signal_im = new double[N];
        for (int i = 0; i < N; i++)
            signal_im[i] = signal[i].im;

        Cnum[] spectrum = forward(signal, sig_shift, Cnum.array(N), spec_shift);
        double[] spectrum_mg = new double[N];
        for (int i = 0; i < N; i++){
            spectrum_mg[i] = spectrum[i].abs();
        }

        Cnum[] signal2 = backward(spectrum, spec_shift, Cnum.array(N), sig_shift);
        double[] signal2_im = new double[N];
        for (int i = 0; i < N; i++){
            signal2_im[i] = signal2[i].im;
        }

        return new double[][]{
                signal_im, spectrum_mg, signal2_im
        };
    }

    static double[][] test(int N, int frequencies, double lambda, double sig_shift, double spec_shift){
        Cnum[][] signals = new Cnum[frequencies][];
        Cnum[] signal = new Cnum[N];

        double allMult = Math.PI/lambda;

        for (int i = 0; i < frequencies; i++) {
            double mult = (1 + 2*i);
            signals[i] = new Cnum[N];
            for (int j = 0; j < N; j++) {
                signals[i][j] = new Cnum(allMult*j*mult)
                        .multIn(1/mult);
                if (signal[j] == null) signal[j] = new Cnum();
                signal[j].plusIn(signals[i][j]);
            }
        }
        return test(signal, sig_shift, spec_shift);
    }

    public static Cnum[] forward(Cnum[] src, double src_shift, Cnum[] dst, double dst_shift){
        int N = src.length;

        double constPow = -2*Math.PI/N, powK;

        for (int k = 0; k < N; k++) {
            powK = constPow*(k+src_shift);
            for (int n = 0; n < N; n++){
                dst[k].plusIn(src[n].mult(new Cnum(powK*(n+dst_shift))));
            }
            dst[k].multIn(1./N);
        }
        return dst;
    }

    public static Cnum[] backward(Cnum[] src, double src_shift, Cnum[] dst, double dst_shift){
        int N = src.length;

        double constPow = 2*Math.PI/N, powN;
        for (int n = 0; n < N; n++) {
            powN = constPow * (n + src_shift);
            for (int k = 0; k < N; k++) {
                dst[n].plusIn(src[k].mult(new Cnum(powN * (k + dst_shift))));
            }
        }
        return dst;
    }



}