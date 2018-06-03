import MonteCarloParticleTracer.*;

import java.io.File;

public class KO_Deuteron_Test {

    public static void main(String ... args) throws Exception{

        File cStopPowFile = new File("src/cStopPow/libcStopPow.so");        // Linux
        System.load(cStopPowFile.getAbsolutePath());

        runModel((int) 1e5, 2);

    }

    public static void runModel(int numParticles, int numCPUs) throws Exception{

        // Source is DT neutrons
        NuclearReaction sourceReaction = NuclearReaction.DT_n;
        Reactivity sourceReactivity = Reactivity.DTn_Reactivity;


        // Init the model
        String name = String.format("KOd_Test_%d.output", System.currentTimeMillis());
        Model model = new Model(name);


        // Create a "hot spot" layer
        Plasma hotspot = Plasma.UniformPlasmaOfKnownTotalMass(50.0, 1.0, 5.0);
        model.addPlasmaLayer(hotspot);


        // Create the "fuel layer"
        Plasma fuelLayer = Plasma.UniformPlasmaOfKnownMassDensity(70.0, 1.0, 5.0);
        model.addPlasmaLayer(fuelLayer);


        // Set the source info
        model.setSourceInformation(sourceReaction, sourceReactivity);


        // This is where we'll add d(n,n')d'


        // Run the simulation
        model.runSimulation(numParticles, numCPUs);
    }
}
