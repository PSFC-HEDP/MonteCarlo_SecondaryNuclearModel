import java.io.File;

public class Demo {

    public static void main(String... args) throws Exception {


	    //File cStopPowFile = new File("cStopPow/libcStopPow.so");        // Linux
        File cStopPowFile = new File("src/cStopPow/cStopPow.DLL");       // Windows
	    System.load(cStopPowFile.getAbsolutePath());


        //generateYieldRatioCurve(1.0, 1.0, false);
        testConvergence(1.0, 1.0, false);

    }

    public static void generateYieldRatioCurve(double Tion, double massDensity, boolean hotSpot) {

        final int NUM_PARTICLES = (int) 1e3;

        /**
         * Make base plasma object
         */
        Plasma uniformPlasma = Plasma.uniformPlasma(1, 1, Tion);
        uniformPlasma.setCenterMassDensity(massDensity);
        uniformPlasma.addDeuteriumSpecies(1.0);

        double[] plasmaRhoRs = Utils.logspace(-4, -1, 100);      // g per cm^2
        for (double rhoR : plasmaRhoRs) {
            // Set the radius
            double plasmaRadius = rhoR / massDensity;
            uniformPlasma.setP0(plasmaRadius);

            // Run the model
            double[] yieldRatios;
            if (hotSpot){
                yieldRatios = runHotSpotModel(uniformPlasma, NUM_PARTICLES);
            }else{
                yieldRatios = runUniformModel(uniformPlasma, NUM_PARTICLES);
            }

            // Print the results
            System.out.printf("%.4e %.4e %.4e\n", rhoR, yieldRatios[0], yieldRatios[1]);
        }

    }

    public static void testConvergence(double Tion, double massDensity, boolean hotSpot) {

        final double rhoR = 0.01;
        //final double protonReference = 2.3501e-04;      // Determined at 1 keV, 1 g/cc, 1 mg/cm^2 (Hot Spot)
        //final double neutronReference = 1.2743e-04;     // Determined at 1 keV, 1 g/cc, 1 mg/cm^2 (Hot Spot)
        final double protonReference = 4.5779e-04;      // Determined at 1 keV, 1 g/cc, 10 mg/cm^2 (Uniform)
        final double neutronReference = 2.0326e-03;     // Determined at 1 keV, 1 g/cc, 10 mg/cm^2 (Uniform)

        /**
         * Make base plasma object
         */
        Plasma uniformPlasma = Plasma.uniformPlasma(1, 1, Tion);
        uniformPlasma.setCenterMassDensity(massDensity);
        uniformPlasma.addDeuteriumSpecies(1.0);

        double plasmaRadius = rhoR / massDensity;
        uniformPlasma.setP0(plasmaRadius);

        double[] particles = Utils.logspace(2, 5, 1000);
        for (double numParticles : particles) {

            // Run the model
            double[] yieldRatios;
            if (hotSpot){
                yieldRatios = runHotSpotModel(uniformPlasma, (int) numParticles);
            }else{
                yieldRatios = runUniformModel(uniformPlasma, (int) numParticles);
            }

            // Print the results
            System.out.printf("%d %.4e %.4e\n", (int) numParticles,
                    100*(yieldRatios[0] - neutronReference)/neutronReference,
                    100*(yieldRatios[1] - protonReference )/protonReference);
        }
    }


    public static double[] runUniformModel(Plasma plasma, int numParticles) {

        // Create the particle distributions
        ParticleDistribution tritons = ParticleDistribution.ThermalDDtDistribution(
                plasma.getSpatialDDpBurnDistribution(), 1e-3);
        ParticleDistribution helium3s = ParticleDistribution.ThermalDD3HeDistribution(
                plasma.getSpatialDDnBurnDistribution(), 1e-3);


        // Create the models
        Model neutronModel = new Model(tritons, NuclearReaction.DTn, plasma);
        Model protonModel = new Model(helium3s, NuclearReaction.D3Hep, plasma);


        // Run the models
        double neutronYieldRatio = neutronModel.getYieldRatio(numParticles);
        double protonYieldRatio = protonModel.getYieldRatio(numParticles);


        // Return results
        double[] yieldRatios = new double[]{neutronYieldRatio, protonYieldRatio};
        return yieldRatios;
    }

    public static double[] runHotSpotModel(Plasma plasma, int numParticles){

        // Create the particle distributions
        ParticleDistribution tritons = ParticleDistribution.ThermalDDtDistribution(
                Distribution.deltaFunction(1e-32), 1e-3);
        ParticleDistribution helium3s = ParticleDistribution.ThermalDD3HeDistribution(
                Distribution.deltaFunction(1e-32), 1e-3);


        // Create the models
        Model neutronModel = new Model(tritons, NuclearReaction.DTn, plasma);
        Model protonModel = new Model(helium3s, NuclearReaction.D3Hep, plasma);


        // Run the models
        double neutronYieldRatio = neutronModel.getYieldRatio(numParticles);
        double protonYieldRatio = protonModel.getYieldRatio(numParticles);


        // Return results
        double[] yieldRatios = new double[]{neutronYieldRatio, protonYieldRatio};
        return yieldRatios;

    }



}
