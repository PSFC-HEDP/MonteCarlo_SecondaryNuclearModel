/**
 * Created by lahmann on 2016-10-04.
 */
public class SahaIonizationModel {

    private final static int[] degeneracyWeights = {
            2, 1,    4, 10, 16, 19, 16, 10, 4, 1,    4, 10, 16, 19, 16, 10, 4, 1
    };

    private final static int[] groundStates = {
            1, 1,    2, 2, 2, 2, 2, 2, 2, 2, 2,    3, 3, 3, 3, 3, 3, 3, 3
    };

    public static double getIonization(int Z, double Ti, double Te, double n){

        double[] epsilons = getEpsilons(Z);
        double[] ratios_ToPreviousState = new double[Z];            // Ratios of n(i+1)/n(i)
        double[] ratios_ToNeutralState  = new double[Z];            // Ratios of n(i+1)/n0
        double[] ratios_ToTotalDensity  = new double[Z+1];          // Ratios of n(i) / n


        double   lambda   = 2.1881e-9 / Math.sqrt(Te);      // DeBroglie Wavelength in cm
        Ti *= 1000;                                         // keV -> eV


        // As an initial guess, let's assume Zbar = Z
        double Zbar = Z;
        double relError = Double.MAX_VALUE;


        // Loop until we converge on Zbar
        while(relError > 1e-3) {

            System.out.println(Zbar);
            double ne = Zbar * n;                             // Free electron density;
            double inverseNeutralFraction = 1.0;              // Ratio of n/n0 (First value of sum n0/n0 is guaranteed to be 1)

            for (int i = 0; i < Z; i++) {
                double g0 = degeneracyWeights[i];
                double e0 = epsilons[i];

                double g1 = degeneracyWeights[i + 1];
                double e1 = epsilons[i + 1];

                double dE = (e1 - e0);

                // Saha equation for n(i+1) / n(i)
                ratios_ToPreviousState[i] = (g1 / g0) * Math.exp(-dE / Ti);
                ratios_ToPreviousState[i] *= 2 / Math.pow(lambda, 3);
                ratios_ToPreviousState[i] /= ne;


                // n(i+1) / n0 is simply the product of all the ratios up to this point
                ratios_ToNeutralState[i] = 1.0;
                for (int j = i; j >= 0; j--) {
                    ratios_ToNeutralState[i] *= ratios_ToPreviousState[j];
                }


                // The inverse of the neutral fraction is simply the sum (n0 + n1 + ...) / n0
                inverseNeutralFraction += ratios_ToNeutralState[i];
            }

            ratios_ToTotalDensity[0] = 1 / inverseNeutralFraction;
            double newZbar = 0;
            for (int i = 1; i < Z + 1; i++) {
                ratios_ToTotalDensity[i] = ratios_ToNeutralState[i - 1] / inverseNeutralFraction;
                newZbar += i * ratios_ToTotalDensity[i];
            }

            relError = Math.abs(newZbar - Zbar) / newZbar;
            Zbar = newZbar;
        }

        return Zbar;
    }

    private static double[] getEpsilons(int Z){
        double[] epsilons = new double[Z+1];
        epsilons[0] = 0;

        for (int i = 1; i < Z+1; i++){
            epsilons[i] = 13.6*Math.pow(i, 2);
            epsilons[i] /= groundStates[Z-i];
            for (int j = i-1; j >= 0; j--){
                epsilons[i] += epsilons[j];
            }
        }

        return epsilons;
    }

}
