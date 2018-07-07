package SecondaryDTnAnalysisGUI;

public class Temp_Database {

    public static Data[] data = {
            new Data(1.5400E+11, 6.0200E+09, 6.9800E+08, 4.2700E+07, 1.85, 0.14, 1124.39808561, 189.91857633, 6.69),
            new Data(8.2355E+10, 3.3216E+09, 5.0100E+08, 3.7500E+07, 1.6529, 0.16686, 1125.22805634, 190.15545359, 6.67),
            new Data(2.9800E+11, 1.1800E+10, 1.9100E+09, 8.3500E+07, 2.1, 0.14, 1118.40064934, 189.43024987, 6.68),
            new Data(3.3449E+11, 1.3064E+10, 2.2100E+09, 9.4900E+07, 2.0675, 0.14, 1122.41248993, 190.57060902, 6.71),

            new Data(6.6500E+11, 2.5700E+10, 5.0700E+09, 1.8700E+08, 2.69, 0.14, 1012.4, 175.25, 6.67),
            new Data(4.6800E+11, 1.8200E+10, 3.1800E+09, 1.2400E+08, 2.48, 0.14, 1006.98852635, 174.65809795, 6.68),


            new Data(1.0104E+12, 3.9120E+10, 2.3100E+09, 9.8700E+07, 3.741516, 0.143574, 1073.965, 74.3069999999999, 3.980587),
            new Data(1.7411E+12, 7.8170E+10, 6.6700E+09, 2.7100E+08, 3.736592, 0.1692915, 1073.892, 74.2734, 3.98573),


            new Data(7.9383E+11, 3.0693E+10, 2.2700E+09, 9.5900E+07, 3.1702, 0.1412, 1068.3, 68.309, 3.99),

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
