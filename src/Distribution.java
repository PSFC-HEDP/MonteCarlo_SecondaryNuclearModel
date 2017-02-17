import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

/**
 * Class that represents an arbitrary distribution used for random sampling
 * Created by lahmann on 2016-06-15.
 */
public class Distribution {


    private static final int NUM_NODES_FOR_PRESETS = 100;      // TODO: Making this too large causes interpolation bugs? Need to understand this.

    private PolynomialSplineFunction probability;
    private PolynomialSplineFunction inverseCumulativeProbability;

    /**
     * Built in predefined distributions
     */

    static Distribution normDistribution(double mu, double sigma){
        final int NUM_SIGMAS_OUT = 5;
        double a = mu - NUM_SIGMAS_OUT*sigma;
        double b = mu + NUM_SIGMAS_OUT*sigma;

        double[] nodes = Utils.linspace(a, b, NUM_NODES_FOR_PRESETS);
        double[] probabilities = new double[NUM_NODES_FOR_PRESETS];

        for (int i = 0; i < nodes.length; i++){
            double power = -0.5*Math.pow((nodes[i] - mu)/sigma, 2);
            probabilities[i] = Math.exp(power);
        }

        return new Distribution(nodes, probabilities);
    }



    /**
     * Constructors
     */

    Distribution(double[] nodes, double[] frequencies) {

        // Normalize the frequencies
        frequencies = normalize(frequencies);


        // Calculate the cumulative probabilities
        double[] cumulativeProbabilities = calculateCumuProbability(frequencies);


        // Set the probability functions
        probability = new LinearInterpolator().interpolate(nodes, normalize(frequencies));
        inverseCumulativeProbability = new LinearInterpolator().interpolate(cumulativeProbabilities, nodes);

    }


    /**
     * Randomly sample from this distribution
     */
    double sample(){
        return inverseCumulativeProbability.value(Math.random());
    }

    double[] sample(int N){
        double[] samples = new double[N];
        for (int i = 0; i < samples.length; i++){
            samples[i] = sample();
        }
        return samples;
    }


    /**
     * Getters
     */

    public double[] getNodes() {
        return probability.getKnots();
    }

    /**
     * Private convenience functions
     */
    private double[] calculateCumuProbability(double[] frequencies){

        double[] cumuProbabilities = new double[frequencies.length];

        double sum = 0.0;
        for (int i = 0; i < frequencies.length; i++){
            sum += frequencies[i];
            cumuProbabilities[i] = sum;
        }

        return cumuProbabilities;
    }

    private double[] normalize(double[] probabilities){

        double[] norm = new double[probabilities.length];

        Sum sum = new Sum();
        double total = sum.evaluate(probabilities);

        for (int i = 0; i < probabilities.length; i++)
            norm[i] = probabilities[i] / total;

        return norm;
    }


}