package SecondaryDTnAnalysisGUI;

public class Utils {

    public static double quadSum(double ... values){

        double sum = 0.0;
        for (double value : values){
            sum += (value*value);
        }
        return Math.sqrt(sum);

    }
}
