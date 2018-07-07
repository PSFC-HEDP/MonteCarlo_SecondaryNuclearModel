import MonteCarloParticleTracer.ParticleType;
import MonteCarloParticleTracer.Plasma;
import MonteCarloParticleTracer.SecondaryDTnModel;
import MonteCarloParticleTracer.Utils;
import SecondaryDTnAnalysisGUI.Capsule;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.ArrayList;

public class TestAsymmetry {

    private final static Vector3D SPEC_SP =
            Utils.getVectorFromSpherical(1798, 161, 56);

    private final static Vector3D SPEC_A =
            Utils.getVectorFromSpherical(2222, 116, 316);

    private final static Vector3D SPEC_E =
            Utils.getVectorFromSpherical(2009, 90, 174);

    private final static Vector3D SPEC_NP =
            Utils.getVectorFromSpherical(2161, 18, 303);



    public static void main(String ... args){

        // Grab a capsule
        Capsule capsule = getD2Capsule(1100, 100, 4.0);

        // Build the fuel plasma
        Plasma fuelPlasma = capsule.getFuelPlasma(4.0, 10.0, 0.0);
        fuelPlasma.addOuterLegendreMode(2, 0, 1.0*fuelPlasma.getOuterRadiusBound(0.0,0.0));

        // Build the model
        SecondaryDTnModel model = new SecondaryDTnModel("temp.OUTPUT");
        model.addPlasmaLayer(fuelPlasma);

        while (true) {
            model.setDetectorLineOfSight(SPEC_SP);
            model.runSimulation((int) 1e5, 1);
            System.out.printf("%.4e ", model.getYieldRatio());

            model.setDetectorLineOfSight(SPEC_A);
            model.runSimulation((int) 1e5, 24);
            System.out.printf("%.4e ", model.getYieldRatio());

            model.setDetectorLineOfSight(SPEC_E);
            model.runSimulation((int) 1e5, 24);
            System.out.printf("%.4e ", model.getYieldRatio());

            model.setDetectorLineOfSight(SPEC_NP);
            model.runSimulation((int) 1e5, 24);
            System.out.printf("%.4e ", model.getYieldRatio());

            System.out.println();
        }


    }

    public static Capsule getD2Capsule(double outerRadius, double thickness, double density){

        ArrayList<ParticleType> fillSpecies = new ArrayList<>();
        fillSpecies.add(ParticleType.deuteron);

        ArrayList<Double> proportions = new ArrayList<>();
        proportions.add(1.0);

        return new Capsule(1100.0, 100.0, 4.0, fillSpecies, proportions);

    }


}
