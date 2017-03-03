/**
 * Created by lahmann on 2017-03-02.
 */
public class Constants {

    /**
     * Physical constants
     */

    public static final double PROTON_MASS_KG = 1.6726219e-27;
    public static final double MEV_PER_AMU = 931.5;
    public static final double SPEED_OF_LIGHT_CM_PER_SEC = 2.998e+10;



    /**
     * Particle Types constants
     */

    public static final int NEUTRON_CHARGE = 0;
    public static final int PROTON_CHARGE = 1;
    public static final int DEUTERIUM_CHARGE = 1;
    public static final int TRITIUM_CHARGE = 1;
    public static final int HELIUM3_CHARGE = 2;
    public static final int ALPHA_CHARGE = 2;

    public static final double NEUTRON_MASS_AMU = 1.00866491588;
    public static final double PROTON_MASS_AMU = 1.00727647;
    public static final double DEUTERIUM_MASS_AMU = 2.014102 - 1*0.00054858;
    public static final double HELIUM3_MASS_AMU = 3.0160293 - 2*0.00054858;
    public static final double TRITIUM_MASS_AMU = 3.0160492 - 1*0.00054858;
    public static final double ALPHA_MASS_AMU = 4.001506179127;



    /**
     * Fusion reaction constants
     */

    // D + D -> p + T
    public static final double DD_P_BIRTH_ENERGY_MEV = 3.02;
    public static final double DD_T_BIRTH_ENERGY_MEV = 1.01;
    public static final double DD_P_ENERGY_RELEASE = DD_P_BIRTH_ENERGY_MEV + DD_T_BIRTH_ENERGY_MEV;

    // D + D -> n + 3He
    public static final double DD_N_BIRTH_ENERGY_MEV = 2.45;
    public static final double DD_3HE_BIRTH_ENERGY_MEV = 0.82;
    public static final double DD_N_ENERGY_RELEASE = DD_N_BIRTH_ENERGY_MEV + DD_3HE_BIRTH_ENERGY_MEV;

    // D + T -> n + a
    public static final double DT_N_BIRTH_ENERGY_MEV = 14.1;
    public static final double DT_A_BIRTH_ENERGY_MEV = 3.5;
    public static final double DT_N_ENERGY_RELEASE = DT_N_BIRTH_ENERGY_MEV + DT_A_BIRTH_ENERGY_MEV;




}
