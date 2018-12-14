package SecondaryDTnAnalysisGUI;

public class Temp_Database {

    public static Data[] data = {
            new Data(6.7575e+12, 3.1916e+11, 6.6400e+10, 2.2500e+09, 2.7700e+00-0.5, 1.6166e-01, 9.0968e+02, 6.5102e+01, 3.9889e+00)
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
