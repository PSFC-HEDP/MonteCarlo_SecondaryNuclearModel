package MonteCarloParticleTracer;

/**
 * Created by lahmann on 2017-04-05.
 */
public class Tally {

    double[] nodes;
    double[] weights;

    public Tally(double[] nodes) {
        this.nodes = nodes;
        this.weights = new double[nodes.length];
    }

    public void addValue(double value, double weight){
        for (int i = 0; i < nodes.length; i++){
            if (value < nodes[i] && i > 0){
                weights[i-1] += weight;
                return;
            }
        }
    }

    // TODO: More robust method for verifying the nodes are the same
    public void addTally(Tally tally){
        for (int i = 0; i < weights.length; i++){
            weights[i] += tally.weights[i];
        }
    }

    public void setTotalWeight(double totalWeight){
        double currentTotal = 0.0;

        for (double weight : weights){
            currentTotal += weight;
        }

        for (int i = 0; i < weights.length; i++){
            weights[i] *= (totalWeight / currentTotal);
        }
    }

    public String toString(){
        String string = "   Nodes   |   Weight\n";

        for (int i = 0; i < nodes.length; i++){
            string += String.format("%.4e | %.4e\n", nodes[i], weights[i]);
        }

        return string;
    }
}
