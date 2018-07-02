package SecondaryDTnAnalysisGUI;

import MonteCarloParticleTracer.Plasma;
import MonteCarloParticleTracer.SecondaryDTnModel;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.ArrayList;
import java.util.Collections;


public class Analyzer {

    private Double INTIAL_CR_GUESS = 15.0;
    private Double INITAL_CR_DELTA = 5.0;

    private final Double   MAX_ERROR = 0.01;


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


    public void doAnalysis(double measuredRatio, double ratioUncertainty, double temperature, double temperatureUncertainty) {

        // Disable the GUI
        parentGUI.setAnalysisEnabled(false);


        // ***********************************
        // Run the simulations (5 loops total)
        // ***********************************

        // Nominal analysis
        results.clear();
        parentGUI.writeFormattedToLog("Beginning nominal analysis for Te = %.2f keV ... ",
                temperature);
        SimulationResult nominalResult = analysisLoop(temperature, measuredRatio);

        // Repeat analysis for plus uncertainty
        parentGUI.writeToLog("Calculating upper uncertainty ... ");
        SimulationResult plusResult = analysisLoop(temperature, measuredRatio + ratioUncertainty);

        // Repeat analysis for minus uncertainty
        parentGUI.writeToLog("Calculating lower uncertainty ... ");
        SimulationResult minusResult = analysisLoop(temperature, measuredRatio - ratioUncertainty);


        // Refine our guess
        INTIAL_CR_GUESS = nominalResult.fuelConvergence;
        INITAL_CR_DELTA = 1.0;


        // Hot analysis
        results.clear();
        parentGUI.writeFormattedToLog("Beginning hot analysis for Te = %.2f keV ... ",
                temperature + temperatureUncertainty);
        SimulationResult hotResult = analysisLoop(temperature + temperatureUncertainty, measuredRatio);


        // Cold analysis
        results.clear();
        parentGUI.writeFormattedToLog("Beginning cold analysis for Te = %.2f keV ... ",
                temperature - temperatureUncertainty);
        SimulationResult coldResult = analysisLoop(temperature - temperatureUncertainty, measuredRatio);


        // *************************************
        // Calculate the results / uncertainties
        // *************************************

        double[] rhoR_Yield_Uncertainties = nominalResult.calculateRhoRUncertainties(plusResult, minusResult);
        double[] CR_Yield_Uncertainties   = nominalResult.calculateCRUncertainties(plusResult, minusResult);
        double[] rho_Yield_Uncertainties  = nominalResult.calculateDensityUncertainties(plusResult, minusResult);
        double[] R_Yield_Uncertainties    = nominalResult.calculateRadiusUncertainties(plusResult, minusResult);

        double[] rhoR_Te_Uncertainties = nominalResult.calculateRhoRUncertainties(hotResult, coldResult);
        double[] CR_Te_Uncertainties   = nominalResult.calculateCRUncertainties(hotResult, coldResult);
        double[] rho_Te_Uncertainties  = nominalResult.calculateDensityUncertainties(hotResult, coldResult);
        double[] R_Te_Uncertainties    = nominalResult.calculateRadiusUncertainties(hotResult, coldResult);


        double[] rhoR_Uncertainties = new double[2];
        rhoR_Uncertainties[0] = Utils.quadSum(rhoR_Yield_Uncertainties[0], rhoR_Te_Uncertainties[0]);
        rhoR_Uncertainties[1] = Utils.quadSum(rhoR_Yield_Uncertainties[1], rhoR_Te_Uncertainties[1]);

        double[] CR_Uncertainties   = new double[2];
        CR_Uncertainties[0] = Utils.quadSum(CR_Yield_Uncertainties[0], CR_Te_Uncertainties[0]);
        CR_Uncertainties[1] = Utils.quadSum(CR_Yield_Uncertainties[1], CR_Te_Uncertainties[1]);

        double[] rho_Uncertainties  = new double[2];
        rho_Uncertainties[0] = Utils.quadSum(rho_Yield_Uncertainties[0], rho_Te_Uncertainties[0]);
        rho_Uncertainties[1] = Utils.quadSum(rho_Yield_Uncertainties[1], rho_Te_Uncertainties[1]);

        double[] R_Uncertainties    = new double[2];
        R_Uncertainties[0] = Utils.quadSum(R_Yield_Uncertainties[0], R_Te_Uncertainties[0]);
        R_Uncertainties[1] = Utils.quadSum(R_Yield_Uncertainties[1], R_Te_Uncertainties[1]);


        // **************************
        // Print final results to GUI
        // **************************

        parentGUI.writeToLog(" === Analysis Results === ");

        parentGUI.writeFormattedToLog("\u03c1R (mg/cm^2) = %.1f - %.2f + %.2f",
                1000*nominalResult.arealDensity, 1000*rhoR_Uncertainties[0], 1000*rhoR_Uncertainties[1]);

        parentGUI.writeFormattedToLog("Capsule CR = %.1f - %.2f + %.2f",
                nominalResult.capsuleConvergence, CR_Uncertainties[0], CR_Uncertainties[1]);

        parentGUI.writeFormattedToLog("\u03c1 (g/cc) = %.1f - %.2f + %.2f",
                nominalResult.density, rho_Uncertainties[0], rho_Uncertainties[1]);

        parentGUI.writeFormattedToLog("R (\u03bcm) = %.1f - %.2f + %.2f",
                10000*nominalResult.fuelRadius, 10000*R_Uncertainties[0], 10000*R_Uncertainties[1]);

        parentGUI.writeToLog(" === DONE! === ");


        // Re-enable the GUI
        parentGUI.setAnalysisEnabled(true);

    }

    private SimulationResult analysisLoop(double Te, double desiredRatio){

        SimulationResult result = null;

        double error = Double.MAX_VALUE;
        while(Math.abs(error) > MAX_ERROR){

            String direction = "";

            double CR = getNextConvergenceGuess(desiredRatio);
            result = doSimulation(Te, CR);
            results.add(result);

            error = (desiredRatio - result.yieldRatio) / desiredRatio;
            if (error > 0)  direction = "low";
            else            direction = "high";


            parentGUI.writeFormattedToLog("Fuel CR = %.2f, Ratio = %.4e (%.2f%% too %s)",
                    CR, result.yieldRatio, 100*Math.abs(error), direction);

        }

        return result;

    }

    private SimulationResult doSimulation(double Te, double CR){

        // Build the fuel plasma
        Plasma fuelPlasma = capsule.getFuelPlasma(Te, CR, profileExponent);

        // Build the model
        SecondaryDTnModel model = new SecondaryDTnModel("temp.OUTPUT");
        model.addPlasmaLayer(fuelPlasma);

        // Run the simulation
        model.runSimulation(numParticles, numCPUs);

        // Get the yield ratio
        double yieldRatio = model.getYieldRatio();

        // Return the simulation results
        return new SimulationResult(Te, CR, yieldRatio, fuelPlasma);

    }

    private double getNextConvergenceGuess(double desiredRatio){

        // If no results, use our first guess
        if (results.size() == 0) return INTIAL_CR_GUESS;

        // If only one result, we can only determine the direction
        if (results.size() == 1){
            if (desiredRatio < results.get(0).yieldRatio)
                return INTIAL_CR_GUESS - INITAL_CR_DELTA;
            else
                return INTIAL_CR_GUESS + INITAL_CR_DELTA;
        }

        // Sort the results
        Collections.sort(results);

        // Convert the results to arrays
        double[] x = new double[results.size()];
        double[] y = new double[results.size()];
        for (int i = 0; i < results.size(); i++){
            x[i] = results.get(i).yieldRatio;
            y[i] = results.get(i).fuelConvergence;
        }


        // Check if the answer is within our guesses before interpolating
        if (x[0]          > desiredRatio)   return y[0]          - INITAL_CR_DELTA;
        if (x[x.length-1] < desiredRatio)   return y[x.length-1] + INITAL_CR_DELTA;


        // Return the interpolated value
        return new LinearInterpolator().interpolate(x, y).value(desiredRatio);

    }


    private class SimulationResult implements Comparable {

        private double fuelConvergence;
        private double temperature;
        private double yieldRatio;

        private double arealDensity;
        private double capsuleConvergence;
        private double density;
        private double fuelRadius;


        SimulationResult(double temperature, double fuelConvergence, double yieldRatio, Plasma fuelPlasma) {

            this.fuelConvergence = fuelConvergence;
            this.temperature = temperature;
            this.yieldRatio = yieldRatio;

            this.arealDensity = fuelPlasma.getArealDensity(0.0, 0.0);
            this.capsuleConvergence = capsule.getCapsuleConvergence(fuelConvergence);
            this.density = fuelPlasma.getMassDensity(new Vector3D(0.0, 0.0, 0.0));
            this.fuelRadius = fuelPlasma.getOuterRadiusBound(0.0, 0.0);

        }

        private double[] calculateRhoRUncertainties(SimulationResult result1, SimulationResult result2){



            double minRhoR = Math.min(result1.arealDensity, result2.arealDensity);
            double maxRhoR = Math.max(result1.arealDensity, result2.arealDensity);

            double[] uncertainties = new double[2];
            uncertainties[0] = arealDensity - minRhoR;
            uncertainties[1] = maxRhoR - arealDensity;

            return uncertainties;
        }

        private double[] calculateCRUncertainties(SimulationResult result1, SimulationResult result2){

            double minCR = Math.min(result1.capsuleConvergence, result2.capsuleConvergence);
            double maxCR = Math.max(result1.capsuleConvergence, result2.capsuleConvergence);

            double[] uncertainties = new double[2];
            uncertainties[0] = capsuleConvergence - minCR;
            uncertainties[1] = maxCR - capsuleConvergence;

            return uncertainties;
        }

        private double[] calculateDensityUncertainties(SimulationResult result1, SimulationResult result2){

            double minRho = Math.min(result1.density, result2.density);
            double maxRho = Math.max(result1.density, result2.density);

            double[] uncertainties = new double[2];
            uncertainties[0] = density - minRho;
            uncertainties[1] = maxRho - density;

            return uncertainties;
        }

        private double[] calculateRadiusUncertainties(SimulationResult result1, SimulationResult result2){

            double minR = Math.min(result1.fuelRadius, result2.fuelRadius);
            double maxR = Math.max(result1.fuelRadius, result2.fuelRadius);

            double[] uncertainties = new double[2];
            uncertainties[0] = fuelRadius - minR;
            uncertainties[1] = maxR - fuelRadius;

            return uncertainties;
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
