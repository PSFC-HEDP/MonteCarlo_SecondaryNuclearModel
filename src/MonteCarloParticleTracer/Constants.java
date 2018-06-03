package MonteCarloParticleTracer;

/**
 * Created by lahmann on 2017-03-02.
 */
class Constants {
    
    // Physical constants
    // ******************

    static final double PROTON_MASS_KG = 1.6726219e-27;
    static final double MEV_PER_AMU = 931.5;
    static final double SPEED_OF_LIGHT_CM_PER_SEC = 2.998e+10;
    static final double BOLTZMANN_CONSTANT_J_PER_K = 1.38064852e-23;

    
    
    // Particle Type constants
    // ***********************

    static final int ELECTRON_CHARGE    = -1;
    static final int NEUTRON_CHARGE     = +0;
    static final int HYDROGEN_CHARGE    = +1;
    static final int HELIUM_CHARGE      = +2;
    static final int CARBON_CHARGE      = +6;

    static final double ELECTRON_MASS_AMU   = 0.00054858;
    static final double NEUTRON_MASS_AMU    = 1.00866491588;
    static final double PROTON_MASS_AMU     = 1.00727647;
    static final double DEUTERIUM_MASS_AMU  = 2.014102  - HYDROGEN_CHARGE*ELECTRON_MASS_AMU;
    static final double TRITIUM_MASS_AMU    = 3.0160492 - HYDROGEN_CHARGE*ELECTRON_MASS_AMU;
    static final double HELIUM3_MASS_AMU    = 3.0160293 -   HELIUM_CHARGE*ELECTRON_MASS_AMU;
    static final double ALPHA_MASS_AMU      = 4.002602  -   HELIUM_CHARGE*ELECTRON_MASS_AMU;
    static final double CARBON_MASS_AMU     = 12.0107   -   CARBON_CHARGE*ELECTRON_MASS_AMU;







}
