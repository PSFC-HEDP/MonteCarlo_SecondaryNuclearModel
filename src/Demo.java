import MonteCarloParticleTracer.*;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.File;

public class Demo {

    public static void main(String... args) throws Exception {


	    //File cStopPowFile = new File("cStopPow/libcStopPow.so");        // Linux
        File cStopPowFile = new File("src/cStopPow/cStopPow.DLL");       // Windows
	    System.load(cStopPowFile.getAbsolutePath());

        testSpectraAsymmetry(0.2);
        //testConvergenceEffects(MonteCarloParticleTracer.Utils.linspace(1e-4*900.0, 1e-4*50.0, 2));

    }

    public static void testSpectraAsymmetry(double ... P2_Fractions){

        final double fillDensity_mgPerCc = 3.8;
        final double initialRadius_um = 900.0;

        final double P0 = 50.0 * 1e-4;
        final double Tion = 3.0;

        final int numParticles = (int) 1e5;
        final double[] energyNodes = Utils.linspace(11.0, 18.0, 100);

        final Vector3D specSP_LOS = new SphericalCoordinates(1800, Math.toRadians(56), Math.toRadians(161)).getCartesian();
        final Vector3D specE_LOS  = new SphericalCoordinates(1800, Math.toRadians(174), Math.toRadians(90)).getCartesian();

        final double totalVolume = (4.0/3.0) * Math.PI * Math.pow(1e-4* initialRadius_um, 3);
        final double totalMass = 1e-3 * fillDensity_mgPerCc * totalVolume;


        for (double fraction : P2_Fractions){

            // Calculate P2
            double P2 = P0 * fraction;


            // Name the model
            String name = String.format("Secondary_DTn_P2_%.2f", fraction);


            // Make the plasma
            PlasmaLayer plasma = PlasmaLayer.uniformPlasma(P0, totalMass, Tion);
            plasma.addDeuteriumSpecies(1.0);
            plasma.addLegendreMode(2, 0, P2);


            // Make the MonteCarloParticleTracer.Particle MonteCarloParticleTracer.Distribution
            ParticleDistribution tritons = ParticleDistribution.ThermalDDtDistribution(
                    plasma.getSpatialDDpBurnDistribution(), 1e-3);


            // Make the model
            Model model = new Model(name, tritons, NuclearReaction.DTn, plasma);


            // Run the model
            model.runSimulation(numParticles, specSP_LOS, energyNodes);
            model.runSimulation(numParticles, specE_LOS, energyNodes);

        }


    }

    public static void testConvergenceEffects(double ... P0s){

        final double fillDensity_mgPerCc = 3.8;
        final double initialRadius_um = 900.0;

        final double Tion = 1.0;

        final int numParticles = (int) 1e5;
        final double[] energyNodes = Utils.linspace(11.0, 18.0, 100);

        final Vector3D polarLineOfSight = new Vector3D(0.0, 0.0, 100.0);

        final double totalVolume = (4.0/3.0) * Math.PI * Math.pow(1e-4* initialRadius_um, 3);
        final double totalMass = 1e-3 * fillDensity_mgPerCc * totalVolume;


        for (double P0 : P0s){

            // Name the model
            String name = String.format("Secondary_DTn_P0_%.2f", P0);
            System.out.print(P0 + " ");


            // Make the plasma
            PlasmaLayer plasma = PlasmaLayer.uniformPlasma(P0, totalMass, Tion);
            plasma.addDeuteriumSpecies(1.0);


            // Make the MonteCarloParticleTracer.Particle MonteCarloParticleTracer.Distribution
            ParticleDistribution tritons = ParticleDistribution.ThermalDDtDistribution(
                    plasma.getSpatialDDpBurnDistribution(), 1e-3);


            // Make the model
            Model model = new Model(name, tritons, NuclearReaction.DTn, plasma);


            // Run the model
            model.runSimulation(numParticles, polarLineOfSight, energyNodes);

        }


    }

    public static void testConvergence(double Tion, double massDensity) {

        final double rhoR = 0.00001;
        //final double protonReference = 2.3501e-04;      // Determined at 1 keV, 1 g/cc, 1 mg/cm^2 (Hot Spot)
        //final double neutronReference = 1.2743e-04;     // Determined at 1 keV, 1 g/cc, 1 mg/cm^2 (Hot Spot)
        final double protonReference = 4.5779e-04;      // Determined at 1 keV, 1 g/cc, 10 mg/cm^2 (Uniform)
        final double neutronReference = 2.0326e-03;     // Determined at 1 keV, 1 g/cc, 10 mg/cm^2 (Uniform)

        /**
         * Make base plasma object
         */
        PlasmaLayer uniformPlasma = PlasmaLayer.uniformPlasma(1, 1, Tion);
        uniformPlasma.setCenterMassDensity(massDensity);
        uniformPlasma.addDeuteriumSpecies(1.0);

        double plasmaRadius = rhoR / massDensity;
        uniformPlasma.setP0(plasmaRadius);
        uniformPlasma.addLegendreMode(2, 0, 0.4*plasmaRadius);

        double[] particles = Utils.logspace(5, 5, 1000);
        for (double numParticles : particles) {
            runUniformModel(uniformPlasma, (int) numParticles);
        }
    }


    public static void runUniformModel(PlasmaLayer plasma, int numParticles) {

        final double[] energyNodes = Utils.linspace(11, 18, 100);

        // Create the particle distributions
        ParticleDistribution tritons = ParticleDistribution.ThermalDDtDistribution(
                plasma.getSpatialDDpBurnDistribution(), 1e-3);
        ParticleDistribution helium3s = ParticleDistribution.ThermalDD3HeDistribution(
                plasma.getSpatialDDnBurnDistribution(), 1e-3);


        // Create the models
        Model neutronModel = new Model("DT Neutron MonteCarloParticleTracer.Model", tritons, NuclearReaction.DTn, plasma);
        Model protonModel = new Model("D3He Proton MonteCarloParticleTracer.Model", helium3s, NuclearReaction.D3Hep, plasma);


        // Run the models
        neutronModel.runSimulation(numParticles, new Vector3D(100, 0, 0), energyNodes);
        protonModel.runSimulation(numParticles, new Vector3D(100, 0, 0), energyNodes);

    }



}
