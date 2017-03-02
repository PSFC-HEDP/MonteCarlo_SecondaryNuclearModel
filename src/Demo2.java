import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.File;

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

        testModel(1e-4*300, 8.5, 0.05);



    }

    public static void testModel(double P0, double burnT, double TeFraction) {

        final int NUM_PARTICLES = (int) 2e4;
        final double TOTAL_MASS = getN170212_003TotalMass();
        final Vector3D LOCATION = new Vector3D(52.0, 0.0, 0.0);


        /**
         * Make the helium-3 plasma
         */
        Plasma uniformPlasma = Plasma.uniformPlasma(P0, TOTAL_MASS, burnT, TeFraction);
        uniformPlasma.addHelium3Species(1.0);

        System.out.println(1000* P0 * uniformPlasma.getMassDensity(new Vector3D(0, 0 ,0)));


        /**
         * Make the D3He proton source
         */
        double temperatureSigmaMeV = 1e-3*Math.sqrt(burnT*5856.0);

        Distribution hotSpotDistribution = Distribution.deltaFunction(1e-32);
        Distribution uniformDistribution = uniformPlasma.getSpatialDDnBurnDistribution();
        Distribution energyDistribution  = Distribution.normDistribution(14.7, temperatureSigmaMeV);

        ParticleDistribution particleDistribution = new ParticleDistribution(1, 1,
                uniformDistribution, energyDistribution);


        /**
         * Make the model
         */
        FiniteSourceModel model = new FiniteSourceModel(particleDistribution, uniformPlasma, LOCATION);
        Particle[] particles = model.rangeParticles(NUM_PARTICLES);


        /**
         * Print results
         */
        for (Particle particle : particles){
            System.out.println(particle.getE());
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

}
