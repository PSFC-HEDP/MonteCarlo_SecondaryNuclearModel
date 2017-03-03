/**
 * Created by lahmann on 2017-03-02.
 */
public class Constants {

    public static final double protonMass_kg = 1.6726219e-27;

    public static final double DD_p_birthEnergy_MeV   = 3.02;
    public static final double DD_T_birthEnergy_MeV   = 1.01;
    public static final double DDp_energyRelease = DD_p_birthEnergy_MeV + DD_T_birthEnergy_MeV;

    public static final double DD_n_birthEnergy_MeV   = 2.45;
    public static final double DD_3He_birthEnergy_MeV = 0.82;
    public static final double DDn_energyRelease = DD_n_birthEnergy_MeV + DD_3He_birthEnergy_MeV;

    public static final int neutronCharge = 0;
    public static final int protonCharge = 1;
    public static final int deuteriumCharge = 1;
    public static final int tritiumCharge = 1;
    public static final int helium3Charge = 2;
    public static final int alphaCharge = 2;

    public static final double neutronMass_amu = 1.00866491588;
    public static final double protonMass_amu = 1.00727647;
    public static final double deuteriumMass_amu = 2.014102 - 1*0.00054858;
    public static final double helium3Mass_amu = 3.0160293 - 2*0.00054858;
    public static final double tritiumMass_amu = 3.0160492 - 1*0.00054858;
    public static final double alphaMass_amu = 4.001506179127;

}
