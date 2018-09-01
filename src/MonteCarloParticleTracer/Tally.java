package MonteCarloParticleTracer;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

/**
 * Created by lahmann on 2017-04-05.
 */
public class Tally {

    double[] binEdges;
    double[] weights;               // Note that weights[i] is between binEdges[i] and binEdges[i+1]
    int   [] totalCounts;           // Number of histories that have contributed to this bin

    double totalWeight = 0.0;

    public Tally(double[] binEdges) {
        this.binEdges = binEdges;
        this.weights = new double[binEdges.length-1];
        this.totalCounts = new int[binEdges.length-1];
    }

    public void addValue(double value, double weight){
        addBiasedValue(value, weight, weight);
    }

    public void addBiasedValue(double value, double bias, double trueWeight){

        // Loop through all of the bin edges (ignoring the minimum)
        for (int i = 1; i < binEdges.length; i++) {

            // If this value is less than the current bin edge
            if (value < binEdges[i]) {

                // Verify that it's greater than the previous one (edge case / sanity check)
                if (value >= binEdges[i - 1]) {

                    // Add this weight to the corresponding tally
                    weights[i - 1] += bias;
                    totalCounts[i - 1]++;
                }

                totalWeight += trueWeight;
                return;
            }
        }
    }

    public void multiply(double scalar){
        for (int i = 0; i < weights.length; i++){
            weights[i] *= scalar;
        }
    }

    // TODO: More robust method for verifying the binEdges are the same
    public void addTally(Tally tally){
        for (int i = 0; i < weights.length; i++){
            weights[i]     += tally.weights[i];
            totalCounts[i] += tally.totalCounts[i];
        }

        totalWeight += tally.totalWeight;
    }

    public double[] getGaussFit(){
        WeightedObservedPoints points = new WeightedObservedPoints();
        double maxX = 0.0, maxY = 0.0;

        double[] binCenters   = getBinCenters();
        double[] weightPerBin = getWeightsPerBin();
        for (int i = 0; i < binCenters.length; i++) {
            double x = binCenters[i];
            double y = weightPerBin[i];

            points.add(x, y);

            if (y > maxY) {
                maxX = x;
                maxY = y;
            }
        }

        double[] startPoint = {maxY, maxX, 0.1};
        GaussianCurveFitter fitter = GaussianCurveFitter.create().withStartPoint(startPoint).withMaxIterations(1000);
        return fitter.fit(points.toList());
    }

    public double[] getGaussValues(){

        double[] binCenters = getBinCenters();
        double[] values = new double[binCenters.length];
        double[] fit = getGaussFit();

        NormalDistribution normalDistribution = new NormalDistribution(fit[1], fit[2]);
        for (int i = 0; i < binCenters.length; i++){
            values[i] = fit[0] * fit[2] * Math.sqrt(2*Math.PI) * normalDistribution.density(binCenters[i]);
        }

        return values;
    }


    public double getMean(){

        double mean = 0.0, totalWeight = 0.0;
        double[] binCenters = getBinCenters();
        for (int i = 0; i < binCenters.length; i++){
            totalWeight += weights[i];
            mean += weights[i] * binCenters[i];
        }

        return mean / totalWeight;

    }

    public double getTotalWeight(){
        return totalWeight;
    }

    public double getTotalOfWeights(){
        double total = 0.0;
        for (double weight : weights){
            total += weight;
        }
        return total;
    }

    public double getNormalization(){
        if (getTotalOfWeights() == 0) return 1.0;
        return getTotalWeight() / getTotalOfWeights();
    }

    public double[] getBinEdges() {
        return binEdges;
    }

    public double[] getBinCenters(){
        double[] binCenters = new double[binEdges.length - 1];
        for (int i = 1; i < binEdges.length; i++){
            binCenters[i-1] = 0.5*(binEdges[i] + binEdges[i-1]);
        }
        return binCenters;
    }

    public double[] getWeights() {
        return weights;
    }

    public double[] getWeightsPerBin(){

        double[] weightPerBin = new double[weights.length];
        for (int i = 0; i < weights.length; i++){
            weightPerBin[i] = weights[i] / (binEdges[i+1] - binEdges[i]);
        }
        return weightPerBin;

    }

    public double[] getUncertainties(){
        double[] uncertainties = new double[weights.length];
        for (int i = 0; i < uncertainties.length; i++){
            if (totalCounts[i] > 0) uncertainties[i] = weights[i] / Math.sqrt(totalCounts[i]);
        }
        return uncertainties;
    }

    public String toString(){
        String string = "Nodes, Weight, Uncertainty\n";
        double normalization = getNormalization();

        for (int i = 1; i < binEdges.length; i++){
            double binCenter   = 0.5*(binEdges[i] + binEdges[i-1]);
            double weight      = normalization * weights[i-1];
            double uncertainty = 0.0;

            if (totalCounts[i-1] > 0) {
                uncertainty = weights[i - 1] / Math.sqrt(totalCounts[i - 1]);
            }

            string += String.format("%.8e, %.8e, %.8e\n", binCenter, weight, uncertainty);
        }

        string += String.format("total, %.8e\n", totalWeight);
        return string;
    }
}
