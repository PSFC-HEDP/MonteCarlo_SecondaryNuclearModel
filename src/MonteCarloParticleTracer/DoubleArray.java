package MonteCarloParticleTracer;

public class DoubleArray {

    private double[] values;

    public DoubleArray(double[] values) {
        this.values = values;
    }

    public static DoubleArray linspace(double a, double b, int N){

        double[] values = new double[N];
        for (int i = 0; i < values.length; i++){
            values[i] = a + (b-a)*i/(values.length - 1);
        }

        return new DoubleArray(values);
    }

    public static DoubleArray logspace(double a, double b, int N){

        double[] values = new double[N];
        for (int i = 0; i < values.length; i++){
            double power =  a + (b-a)*i/(values.length - 1);
            values[i] = Math.pow(10, power);
        }

        return new DoubleArray(values);
    }

    public double[] getValues() {
        return values;
    }

    public double get(int i){
        return values[i];
    }

    public void set(int i, double value){
        this.values[i] = value;
    }

    public int length(){
        return values.length;
    }

    public double getMax(){
        double maxValue = Double.MIN_VALUE;
        for (double value : values) {
            maxValue = Math.max(maxValue, value);
        }
        return maxValue;
    }

    public double getMin(){
        double minValue = Double.MAX_VALUE;
        for (double value : values) {
            minValue = Math.min(minValue, value);
        }
        return minValue;
    }

    public void multiply(double scalar){
        for (int i = 0; i < values.length; i++){
            values[i] *= scalar;
        }
    }

    public void add(double scalar){
        for (int i = 0; i < values.length; i++){
            values[i] += scalar;
        }
    }

    public void subtract(double scalar){
        for (int i = 0; i < values.length; i++){
            values[i] -= scalar;
        }
    }





}
