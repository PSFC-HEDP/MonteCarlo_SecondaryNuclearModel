import MonteCarloParticleTracer.*;
import PlottingAPI.Figure;
import PlottingAPI.LineProperties;
import SecondaryDTnAnalysisGUI.Capsule;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.ArrayList;

public class TestAsymmetry {

    private final static Vector3D SPEC_SP =
            Utils.getVectorFromSpherical(1798, Math.toRadians(161), Math.toRadians(56));

    private final static Vector3D SPEC_A =
            Utils.getVectorFromSpherical(2222, Math.toRadians(116), Math.toRadians(316));

    private final static Vector3D SPEC_E =
            Utils.getVectorFromSpherical(2009, Math.toRadians(90), Math.toRadians(174));

    private final static Vector3D SPEC_NP =
            Utils.getVectorFromSpherical(2161, Math.toRadians(18), Math.toRadians(303));



    public static void main(String ... args){

        // Grab a capsule
        Capsule capsule = getD2Capsule(1100, 100, 4.0);

        // Build the fuel plasma
        Plasma fuelPlasma = capsule.getFuelPlasma(4.0, 10.0, 0.0);
        fuelPlasma.addOuterLegendreMode(2, 0, -0.4*fuelPlasma.getOuterRadiusBound(0.0,0.0));
        //fuelPlasma.DEBUG_dumpOuterBounds();
        //fuelPlasma.DEBUG_printOuterModes();

        // Build the model
        SecondaryDTnModel model = new SecondaryDTnModel("temp.OUTPUT");
        model.addPlasmaLayer(fuelPlasma);

        Figure neutronFigure = new Figure("Secondary DT Spectrum", "Energy (MeV)", "Yield / MeV");

        neutronFigure.setLocation(100, 100 );
        neutronFigure.setSize(600, 400);


        model.setProductDirection(SPEC_SP);
        model.runSimulation((int) 1e5, 24);
        System.out.printf("%.4e ", model.getYieldRatio());
        Tally data = model.getSecondaryDTNeutronSpectrum();
        neutronFigure.plot(data.getBinCenters(), data.getWeights());


        model.setProductDirection(SPEC_A);
        model.runSimulation((int) 1e5, 24);
        System.out.printf("%.4e ", model.getYieldRatio());
        data = model.getSecondaryDTNeutronSpectrum();
        neutronFigure.plot(data.getBinCenters(), data.getWeights());


        model.setProductDirection(SPEC_E);
        model.runSimulation((int) 1e5, 24);
        System.out.printf("%.4e ", model.getYieldRatio());
        data = model.getSecondaryDTNeutronSpectrum();
        neutronFigure.plot(data.getBinCenters(), data.getWeights());


        model.setProductDirection(SPEC_NP);
        model.runSimulation((int) 1e5, 24);
        System.out.printf("%.4e ", model.getYieldRatio());data = model.getSecondaryDTNeutronSpectrum();
        neutronFigure.plot(data.getBinCenters(), data.getWeights());



    }

    public static Capsule getD2Capsule(double outerRadius, double thickness, double density){

        ArrayList<ParticleType> fillSpecies = new ArrayList<>();
        fillSpecies.add(ParticleType.deuteron);

        ArrayList<Double> proportions = new ArrayList<>();
        proportions.add(1.0);

        return new Capsule(outerRadius, thickness, density, fillSpecies, proportions);

    }


}
