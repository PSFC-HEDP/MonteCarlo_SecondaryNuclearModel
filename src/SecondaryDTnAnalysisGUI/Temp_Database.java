package SecondaryDTnAnalysisGUI;

public class Temp_Database {

    public static Data[] data = {
            new Data(3.1020e+11, 1.4774e+10, 2.0600e+09, 1.4600e+08, 1.5159e+00, 1.4120e-01, 1.1280e+03, 7.8208e+01, 4.6091e+00)
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
