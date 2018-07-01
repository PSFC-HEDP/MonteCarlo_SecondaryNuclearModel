package SecondaryDTnAnalysisGUI;

import MonteCarloParticleTracer.Plasma;
import MonteCarloParticleTracer.SecondaryDTnModel;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;

import java.util.ArrayList;
import java.util.Collections;


public class Analyzer {

    private final Double[] INTIAL_CR_GUESSES = new Double[] {10.0, 15.0};

    private GUI parentGUI;

    private Capsule capsule;
    private Double profileExponent;
    private Integer numParticles;
    private Integer numCPUs;

    private ArrayList<SimulationResult> results = new ArrayList<>();


    public Analyzer(GUI parentGUI, Capsule capsule, Double profileExponent, Integer numParticles, Integer numCPUs) {
        this.parentGUI = parentGUI;
        this.capsule = capsule;
        this.profileExponent = profileExponent;
        this.numParticles = numParticles;
        this.numCPUs = numCPUs;
    }


    public void doAnalysis(double measuredRatio) {

        parentGUI.setAnalysisEnabled(false);
        parentGUI.writeToLog("Beginning analysis ... ");


        // Our first guess at CR
        double CR = INTIAL_CR_GUESSES[0];
        double ratio = getYieldRatio(CR);
        double error = Math.abs(measuredRatio - ratio) / measuredRatio;
        results.add(new SimulationResult(CR, ratio));

        parentGUI.writeFormattedToLog("CR = %.2f, Ratio = %.2f (%.2f%% Off)",
                CR, ratio, 100*error);


        // Our second guess at CR
        CR = INTIAL_CR_GUESSES[1];
        ratio = getYieldRatio(CR);
        error = Math.abs(measuredRatio - ratio) / measuredRatio;
        results.add(new SimulationResult(CR, ratio));

        parentGUI.writeFormattedToLog("CR = %.2f, Ratio = %.2f (%.2f%% Off)",
                CR, ratio, 100*error);


        while(true){

            CR = getNextConvergenceGuess(measuredRatio);
            ratio = getYieldRatio(CR);
            error = Math.abs(measuredRatio - ratio) / measuredRatio;
            results.add(new SimulationResult(CR, ratio));

            parentGUI.writeFormattedToLog("CR = %.2f, Ratio = %.2f (%.2f%% Off)",
                    CR, ratio, 100*error);

        }

        //parentGUI.setAnalysisEnabled(true);




    }

    private double getYieldRatio(double CR){

        // Build the fuel plasma
        Plasma fuelPlasma = capsule.getFuelPlasma(CR, profileExponent);

        // Build the model
        SecondaryDTnModel model = new SecondaryDTnModel("temp.OUTPUT");
        model.addPlasmaLayer(fuelPlasma);

        // Run the simulation
        model.runSimulation(numParticles, numCPUs);

        // Return the yield ratio
        return model.getYieldRatio();
    }

    private double getNextConvergenceGuess(double desiredRatio){

        // Sort the results
        Collections.sort(results);

        // Convert the results to arrays
        double[] x = new double[results.size()];
        double[] y = new double[results.size()];
        for (int i = 0; i < results.size(); i++){
            x[i] = results.get(i).yieldRatio;
            y[i] = results.get(i).convergence;
        }

        // Return the interpolated value
        return new LinearInterpolator().interpolate(x, y).value(desiredRatio);

    }


    private class SimulationResult implements Comparable {

        private double convergence;
        private double yieldRatio;

        public SimulationResult(double convergence, double yieldRatio) {
            this.convergence = convergence;
            this.yieldRatio = yieldRatio;
        }

        @Override
        public int compareTo(Object o) {
            double otherRatio = ((SimulationResult) o).yieldRatio;

            if (otherRatio > this.yieldRatio) return -1;
            if (otherRatio < this.yieldRatio) return  1;
            return 0;
        }
    }
}
