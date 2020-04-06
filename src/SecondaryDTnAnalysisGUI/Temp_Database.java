package SecondaryDTnAnalysisGUI;

public class Temp_Database {

    public static Data[] data = {
            new Data(1.4642e+13, 8.5325e+11, 1.6100e+11, 5.4600e+09, 4.1498e+00, 1.6544e-01, 1.1282e+03, 7.7910e+01, 3.9887e+00)
    };

    public static class Data {

        double Y1n, Y1n_unc;
        double Y2n, Y2n_unc;
        double Te , Te_unc;

        double Ro;
        double t;
        double rho;

        public Data(double y1n, double y1n_unc, double y2n, double y2n_unc, double te, double te_unc, double ro, double t, double rho) {
            Y1n = y1n;
            Y1n_unc = y1n_unc;
            Y2n = y2n;
            Y2n_unc = y2n_unc;
            Te = te;
            Te_unc = te_unc;
            Ro = ro;
            this.t = t;
            this.rho = rho;
        }
    }
}
