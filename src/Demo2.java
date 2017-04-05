import MonteCarloParticleTracer.*;
import cStopPow.DoubleVector;
import cStopPow.StopPow;
import cStopPow.StopPow_LP;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by lahmann on 2017-02-27.
 */
public class Demo2 {

    private static double N170212_003_FillPressure_Torr = 7600.0;
    private static double N170212_003_FillTemperature_K = 293.15;
    private static double N170212_003_InnerDiameter_cm  = 0.3;

    private static double N170212_003_meanEnergy = 14.194;
    private static double N170212_003_inferredTion = 17.9;

    public static void main(String... args) throws Exception {


        //File cStopPowFile = new File("cStopPow/libcStopPow.so");        // Linux
        File cStopPowFile = new File("src/cStopPow/cStopPow.DLL");       // Windows
        System.load(cStopPowFile.getAbsolutePath());


        double[] P0s = new double[] {147};
        double[] burnTs = new double[] {6.8};

        for (double P0 : P0s) {
            for (double burnT : burnTs) {
                System.out.println("Starting at " + System.currentTimeMillis());
                testModel(1e-4 * P0, burnT, 1.0);
                System.out.println("Done at " + System.currentTimeMillis());
            }
        }

    }

    public static void testModel(double P0, double burnT, double TeFraction) throws Exception {

        final int NUM_PARTICLES = (int) 5e5;
        final double TOTAL_MASS = getN170212_003TotalMass();
        final double[] energyNodes = Utils.linspace(10.0, 16.0, 301);

        // Relation from Maria (03-02-2017) using basic reactivity scaling based on the measured yield
        //double P0 = 6.83088 * burnT * burnT - 28.7257 * burnT + 25.7562;        // in um
        //P0 *= 1e-4;     // um -> cm

        FileWriter w = new FileWriter("data.csv", true);
        w.write(String.format("%.4e,%.4e,", 1e4*P0, burnT));

        /**
         * Make the helium-3 plasma
         */
        PlasmaLayer plasma = PlasmaLayer.uniformPlasma(P0, TOTAL_MASS, burnT, TeFraction);
        plasma.addHelium3Species(1.0);



        /**
         * Make the D3He proton source
         */
        double temperatureSigmaMeV = 1e-3*Math.sqrt(burnT*5856.0);

        Distribution uniformDistribution = plasma.getSpatialDDnBurnDistribution();
        Distribution energyDistribution  = Distribution.normDistribution(14.7, temperatureSigmaMeV);

        ParticleDistribution particleDistribution = new ParticleDistribution(1, 1,
                uniformDistribution, energyDistribution);



        /**
         * Make the model
         */
        Model model = new Model("FiniteSource", particleDistribution, null, plasma);
        model.runSimulation(NUM_PARTICLES, null, energyNodes);



        /**
         * Fit a gaussian to the resulting tally and infer a Tion
         */
        WeightedObservedPoints points = new WeightedObservedPoints();
        for (int i = 0; i < model.energyNodes.length; i++){
            points.add(model.energyNodes[i], model.energyTallies[i]);
        }

        double[] parameters = GaussianCurveFitter.create().fit(points.toList());
        double meanEnergy_MeV = parameters[1];
        double sigma_keV = 1e3*parameters[2];

        double inferredTion = sigma_keV*sigma_keV / 5856.0;

        w.write(String.format("%.4e,%.4e,", meanEnergy_MeV, inferredTion));


        /**
         * Range that spectrum through CH to make meanEnergy consistent
         */
        Spectrum spectrum = new Spectrum(model.energyNodes, model.energyTallies);
        StopPow_LP stopPowLp = getCH_StopPow(1.0, 1.06);
        double thickness;
        try {
            thickness = stopPowLp.Thickness(meanEnergy_MeV, N170212_003_meanEnergy);
        }catch (IllegalArgumentException e){
            w.write("\n");
            w.close();
            return;
        }

        spectrum.range(stopPowLp, thickness);

        /**
         * Fit a gaussian to this new spectrum
         */
        points = new WeightedObservedPoints();
        for (int i = 0; i < spectrum.energies.size(); i++){
            points.add(spectrum.energies.get(i), spectrum.frequencies.get(i));
        }
        parameters = GaussianCurveFitter.create().fit(points.toList());
        meanEnergy_MeV = parameters[1];
        sigma_keV = 1e3*parameters[2];

        inferredTion = sigma_keV*sigma_keV / 5856.0;

        w.write(String.format("%.4e,%.4e,%.4e\n", thickness, meanEnergy_MeV, inferredTion));
        w.close();

        System.out.println(spectrum);
    }

    private static double getN170212_003TotalMass(){

        final double k  = 1.381e-23;        // Boltzmann's constant (J / K)
        final double mp = 1.673e-24;        // Proton mass (g)

        double P = N170212_003_FillPressure_Torr;
        P /= 7500.6;        // Torr -> J / cm^3

        double T = N170212_003_FillTemperature_K;

        double numberDensity = P / (k * T);                 // per cm^2
        double massDensity = 3.0 * mp * numberDensity;      // 3.0 because Helium-3 [g / cm^3]

        double volume = (4.0/3.0) * Math.PI * Math.pow(N170212_003_InnerDiameter_cm/2.0, 3); // cm^3

        double totalMass =  massDensity * volume;

        return totalMass;
    }

    private static PlasmaLayer getN170212_003Plasma(double P0, double burnT, double TeFraction){

        try {
            ArrayList<Double> massDensities = new ArrayList<>();
            ArrayList<Double> temperatures = new ArrayList<>();

            Scanner s = new Scanner(new File(".\\data\\LilacProfiles5.2atm70kJ_complete.txt"));

            while (s.hasNext()){
                s.next();
                massDensities.add(s.nextDouble());
                temperatures.add(1e-3*s.nextDouble());
            }

            double[] Ti = new double[temperatures.size()];
            double[] Te = new double[temperatures.size()];
            double[] rho = new double[massDensities.size()];
            for (int i = 0; i < massDensities.size(); i++){
                Ti[i] = temperatures.get(i);
                Te[i] = Ti[i];
                rho[i] = massDensities.get(i);
            }

            PlasmaLayer plasma = new PlasmaLayer(Ti, Te, rho);
            plasma.setP0(P0);
            plasma.setD3HepBurnAveragedIonTemperature(burnT);
            plasma.setElectronTemperatureFraction(TeFraction);
            plasma.setTotalMass(getN170212_003TotalMass());

            return plasma;

        }catch (Exception e){
            e.printStackTrace();
            System.exit(-1);
            return null;
        }
    }

    public static StopPow_LP getCH_StopPow(double Te, double rho){

        double mt = 1.0;
        double Zt = 1.0;
        double n  = rho / 1.67e-24 / (6 + 0.5);

        DoubleVector mf = new DoubleVector();
        mf.add(1.0);
        mf.add(12.0);
        mf.add(1.0 / 1836.0);


        DoubleVector Zf = new DoubleVector();
        Zf.add(1.0);
        Zf.add(6.0);
        Zf.add(-1.0);


        DoubleVector Tf = new DoubleVector();
        Tf.add(Te);
        Tf.add(Te);
        Tf.add(Te);


        DoubleVector nf = new DoubleVector();
        nf.add(0.5*n);
        nf.add(0.5*n);
        nf.add(0.5*n + 6.0*n);      // Assuming fully ionized


        StopPow_LP stopPowLp = new StopPow_LP(mt, Zt, mf, Zf, Tf, nf);
        stopPowLp.set_mode(StopPow.getMODE_RHOR());
        return stopPowLp;
    }

}
