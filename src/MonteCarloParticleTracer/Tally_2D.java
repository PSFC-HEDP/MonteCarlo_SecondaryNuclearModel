package MonteCarloParticleTracer;

public class Tally_2D {

    double[]   xBins;
    double[]   yBins;
    double[][] weights;               // Note that weights[i][j] is between xBins[i], xBins[i+1], yBins[j], and yBins[j+1]
    int   [][] totalCounts;           // Number of histories that have contributed to this bin

    double totalWeight = 0.0;

    public Tally_2D(double[] xBins, double[] yBins) {
        this.xBins = xBins;
        this.yBins = yBins;
        this.weights     = new double[xBins.length - 1][yBins.length - 1];
        this.totalCounts = new int   [xBins.length - 1][yBins.length - 1];
    }

    public void addValue(double xValue, double yValue, double weight){

        int xIndex = 0;
        for (int i = 1; i < xBins.length; i++){

            if (xValue < xBins[i] && xValue >= xBins[i-1]){
                xIndex = (i-1);
                break;
            }
        }

        int yIndex = 0;
        for (int j = 1; j < yBins.length; j++){

            if (yValue < yBins[j] && yValue >= yBins[j-1]){
                yIndex = (j-1);
                break;
            }
        }

        weights[xIndex][yIndex] += weight;
        totalCounts[xIndex][yIndex]++;

    }

    public String toString(){

        StringBuilder stringBuilder = new StringBuilder("0");

        // Write the y values
        String prefix = ",";
        for (int j = 0; j < weights[0].length; j++){
            stringBuilder.append(String.format("%s%.4e", prefix, 0.5*(yBins[j]+yBins[j+1])));
        }
        stringBuilder.append("\n");

        // Write the data
        for (int i = 0; i < weights.length; i++){

            stringBuilder.append(0.5*(xBins[i]+xBins[i+1]));
            for (int j = 0; j < weights[i].length; j++){
                stringBuilder.append(String.format("%s%.4e", prefix, weights[i][j]));
            }
            stringBuilder.append("\n");

        }
        return stringBuilder.toString();

    }
}
