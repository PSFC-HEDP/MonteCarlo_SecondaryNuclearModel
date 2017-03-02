
/**
 * Created by lahmann on 2016-06-15.
 */
public class Utils {

    public static final double protonMass_kg = 1.6726219e-27;

    public static final double DD_p_birthEnergy_MeV   = 3.02;
    public static final double DD_T_birthEnergy_MeV   = 1.01;

    public static final double DD_n_birthEnergy_MeV   = 2.45;
    public static final double DD_3He_birthEnergy_MeV = 0.82;

    public static final double MINIMUM_LI_PETRASSO_ENERGY_MeV = 0.01;


    public static final String D3He_ENDF_XS_FILE = "./data/d3He_BoschHale.dat";
    public static final String DT_ENDF_XS_FILE   = "./data/dT_BoschHale.dat";

    public static final String DDp_REACTIVITY_FILE = "./data/DDp_Reactivities.dat";
    public static final String DDn_REACTIVITY_FILE = "./data/DDn_Reactivities.dat";

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
