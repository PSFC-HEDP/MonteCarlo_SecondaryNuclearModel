import MonteCarloParticleTracer.*;
import SecondaryDTnAnalysisGUI.Capsule;

import java.io.FileWriter;
import java.util.ArrayList;

public class YieldRatioScript {

    public static void main(String ... args) throws Exception{

        // Grab a capsule
        Capsule capsule = getD2Capsule(1100, 100, 4.0);

        DoubleArray CR2s = DoubleArray.logspace(0, 3, 3);
        double gamma = 0;
        int numParticles = (int) 2e6;

        FileWriter w = new FileWriter("data.dat");
        for (int j = 0; j < CR2s.length(); j++) {

            // Build the fuel plasma
            Plasma fuelPlasma = capsule.getFuelPlasma(4.0, Math.sqrt(CR2s.get(j)), gamma);

            // Build the model

            String name = String.format("output/D2_CR_Scaling_CR2_%.2f_%.0e.output", CR2s.get(j), (double) numParticles);
            System.out.println(name);
            SecondaryDTnModel model = new SecondaryDTnModel(name);
            model.addPlasmaLayer(fuelPlasma);
            model.runSimulation(numParticles, 45);

            w.write(String.format("%.4e, %.4e, %.4e\n", fuelPlasma.getArealDensity(0.0, 0.0), CR2s.get(j), model.getYieldRatio()));

        }
        w.close();

    }

    public static Capsule getD2Capsule(double outerRadius, double thickness, double density){

        ArrayList<ParticleType> fillSpecies = new ArrayList<>();
        fillSpecies.add(ParticleType.deuteron);

        ArrayList<Double> proportions = new ArrayList<>();
        proportions.add(1.0);

        return new Capsule(outerRadius, thickness, density, fillSpecies, proportions);

    }
}
