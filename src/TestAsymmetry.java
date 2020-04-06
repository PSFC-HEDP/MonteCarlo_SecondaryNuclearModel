import MonteCarloParticleTracer.*;
import PlottingAPI.Figure;
import PlottingAPI.LineProperties;
import SecondaryDTnAnalysisGUI.Capsule;

import java.io.FileWriter;
import java.util.ArrayList;

public class TestAsymmetry {

    private final static PTOF_Detector SPEC_SP = PTOF_Detector.SPEC_SP01_DT();
    private final static PTOF_Detector SPEC_A  = PTOF_Detector.SPEC_A01_DT();
    private final static PTOF_Detector SPEC_E  = PTOF_Detector.SPEC_E01_DT();
    private final static PTOF_Detector SPEC_NP = PTOF_Detector.SPEC_NP01_DT();

    public static void main(String ... args) throws Exception{

        // Grab a capsule
        Capsule capsule = getD2Capsule(1100, 100, 4.0);

        // Detector LOS
        PTOF_Detector detector = SPEC_NP;

        /*
        Figure neutronFigure = new Figure("SPEC-SP DT Spectrum", "Time (s)", "Volts / Yield");
        neutronFigure.setLocation(100, 100);
        neutronFigure.setSize(600, 400);
        neutronFigure.setYLimits(0, 1e-12);
        neutronFigure.setXLimits(0, 1e-6);
         */


        DoubleArray CRs = new DoubleArray(new double[] {5.0});
        DoubleArray P2s = DoubleArray.linspace(-0.4, 0.4, 3);
        double gamma = 0;

        FileWriter w = new FileWriter(detector.getName() + ".dat");
        for (int i = 0; i < P2s.length(); i++) {
            for (int j = 0; j < CRs.length(); j++) {

                // Build the fuel plasma
                Plasma fuelPlasma = capsule.getFuelPlasma(4.0, CRs.get(j), gamma);
                fuelPlasma.addOuterLegendreMode(2, 0, P2s.get(i) * fuelPlasma.getOuterRadiusBound(0.0, 0.0));

                // Build the model
                String shape = "";
                if (P2s.get(i) > 0) shape = "prolate";
                else if (P2s.get(i) < 0) shape = "oblate";
                else shape = "round";


                String name = String.format("output/%s_%s_CR_%.1f_P2_%.1f_Gamma_%.4f_1e6.output", detector.getName(), shape, CRs.get(j), Math.abs(100 * P2s.get(i)), gamma);
                System.out.println(name);
                SecondaryDTnModel model = new SecondaryDTnModel(name);
                model.addPlasmaLayer(fuelPlasma);
                model.setProductDirection(detector);
                model.runSimulation((int) 2e6, 45);

                Tally data = detector.generateResponse(model.getSecondaryDTNeutronSpectrum());
                //System.out.println(data);
                w.write(name + "\n");
                w.write(data.toString() + "\n");
                //neutronFigure.plot(data.getBinCenters(), data.getWeights(), LineProperties.blackLine(2));

            }
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
