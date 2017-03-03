
/**
 * Created by lahmann on 2016-06-15.
 */
public class Utils {



    public static final double MINIMUM_LI_PETRASSO_ENERGY_MeV = 0.01;


    public static final String D3He_ENDF_XS_FILE = "./data/d3He_BoschHale.dat";
    public static final String DT_ENDF_XS_FILE   = "./data/dT_BoschHale.dat";

    public static final String DDp_REACTIVITY_FILE = "./data/DDp_Reactivities.dat";
    public static final String DDn_REACTIVITY_FILE = "./data/DDn_Reactivities.dat";
    public static final String D3Hep_REACTIVITY_FILE = "./data/D3Hep_Reactivities.dat";


    public static double[] linspace(double a, double b, int N){
        double[] nodes = new double[N];

        for (int i = 0; i < nodes.length; i++){
            nodes[i] = a + (b-a)*i/(nodes.length - 1);
        }

        return nodes;
    }

    public static double[] logspace(double a, double b, int N){
        double[] nodes = new double[N];

        for (int i = 0; i < nodes.length; i++){
            double power =  a + (b-a)*i/(nodes.length - 1);
            nodes[i] = Math.pow(10, power);
        }

        return nodes;
    }

}
