import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

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

    public static void main(String... args) throws Exception {


        //File cStopPowFile = new File("cStopPow/libcStopPow.so");        // Linux
        File cStopPowFile = new File("src/cStopPow/cStopPow.DLL");       // Windows
        System.load(cStopPowFile.getAbsolutePath());

        testModel(1e-4*300, 8.5, 1.0);



    }

    public static void testModel(double P0, double burnT, double TeFraction) {

        final int NUM_PARTICLES = (int) 5e4;
        final double TOTAL_MASS = getN170212_003TotalMass();
        final Vector3D LOCATION = new Vector3D(52.0, 0.0, 0.0);

        // Relation from Maria (03-02-2017) using basic reactivity scaling based on the measured yield
        //double P0 = 6.83088 * burnT * burnT - 28.7257 * burnT + 25.7562;        // in um
        //P0 *= 1e-4;     // um -> cm


        /**
         * Make the helium-3 plasma
         */
        //Plasma plasma = getN170212_003Plasma(P0, burnT, TeFraction);
        Plasma plasma = Plasma.uniformPlasma(P0, TOTAL_MASS, burnT, TeFraction);
        plasma.addHelium3Species(1.0);

        //System.out.println(plasma);


        /**
         * Make the D3He proton source
         */
        double temperatureSigmaMeV = 1e-3*Math.sqrt(burnT*5856.0);

        Distribution uniformDistribution = plasma.getSpatialDDnBurnDistribution();
        Distribution energyDistribution  = Distribution.normDistribution(14.7, temperatureSigmaMeV);

        ParticleDistribution particleDistribution = new ParticleDistribution(1, 1,
                uniformDistribution, energyDistribution);

        //System.out.println(particleDistribution);


        /**
         * Make the model
         */
        FiniteSourceModel model = new FiniteSourceModel(particleDistribution, plasma, LOCATION);
        Particle[] particles = model.rangeParticles(NUM_PARTICLES);


        /**
         * Print results
         */
        try {
            FileWriter w = new FileWriter(".\\temp.dat");
            for (Particle particle : particles) {
                w.write(particle.getE() + "\n");
            }
            w.close();
        }
        catch (Exception e){

        }

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

    private static Plasma getN170212_003Plasma(double P0, double burnT, double TeFraction){

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

            Plasma plasma = new Plasma(Ti, Te, rho);
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

}
